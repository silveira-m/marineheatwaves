import time
start_time_script = time.time()

import xarray as xr
import dask.array as da
import numpy as np

from netCDF4 import Dataset
from datetime import date, datetime

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import pandas as pd
import math

import geopandas as gpd
import regionmask
from shapely.ops import unary_union

# ============================================================
# Options
# ============================================================
print("OPTIONS")

pathdata = "DATA"
pathoutputfig = "PLOT"

# input data
model = 'oisst'  # oisst cmems
yearbeg = 1981
yearend = 2025  # 2024
region = 'Portugal_200m_bathymetry'  # Portugal 200m bathymetry from shapefile
choice_mask_region = 0  # -1=no mask 4=Med
choice_date = 0  # 0=no, 1=yes (choose start_date and end_date)
start_mth = 1
start_day = 1
end_mth = 12
end_day = 31
print("choice_mask_region", choice_mask_region)

# import_data
import_data = 'y'

# Calculation
mhw_calculation = 'y'
detrend = 0  # 0/1
clim_start_year = 1982
clim_end_year = 2024
percentile = 0.9
min_duration_mhw = 5
min_gap_mhw = 2
dt_clim = 'time.dayofyear'  # 'time.dayofyear'  # 'time.month'
dt = 'day'

print("OPTIONS/PLOT")
plot = 'y'
size_fonte = 14
format_fig = 'png'
dpi_value = 150

# ============================================================
# Import Data (common to all zones)
# ============================================================
if import_data == 'y':
    print("TIME IMPORT DATA - START", (time.time() - start_time_script) / 60)

    ds = xr.open_dataset(
        r"/Users/miguelsilveira/Documents/TeseFilesPython/marine-extreme-event-main/sst_merged_big_box_IbPen.nc",
        engine='netcdf4')
    sst_full = ds.sst
    print("type(var)", type(sst_full))
    print("sst.shape", sst_full.shape)
    print("sst.dims", sst_full.dims)

    ds_grid = xr.open_dataset(
        f'/Users/miguelsilveira/Documents/TeseFilesPython/marine-extreme-event-main/gridarea.nc',
        engine='netcdf4')
    print(ds_grid)
    cell_area = ds_grid.cell_area

    if choice_mask_region > 0:
        mask_region = np.load(f'{pathdata}/OceanMasks_oisst_{region}.npy')
        if model != 'oisst':
            mask_interp = xr.DataArray(mask_region,
                                        dims=('lat', 'lon'),
                                        coords={'lat': ds['lat'].values, 'lon': ds['lon'].values})
            mask_region = mask_interp.interp(lat=ds.lat, lon=ds.lon, method='nearest')
        print("mask_region", mask_region)
        print("np.shape(mask_region)", np.shape(mask_region))
        sst_full = sst_full.where(mask_region == choice_mask_region)

    # --------------------------------------------------------
    # ZONES: load shapefile and KEEP the 3 polygons separate
    # (previously this did unary_union to merge them into one)
    # --------------------------------------------------------
    shp_path = "200m"   # <- change this
    gdf = gpd.read_file(shp_path)

    if gdf.crs and gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(epsg=4326)

    # <<< If the automatic north->south sort below mislabels a zone, set this
    # to a 3-item list matching the RAW row order printed just below
    # (e.g. ["NW", "S", "SW"]), and it will be used instead of auto-sorting.
    ZONE_LABELS_OVERRIDE = None

    gdf["centroid_lat"] = gdf.geometry.centroid.y
    gdf["centroid_lon"] = gdf.geometry.centroid.x
    print("Raw shapefile zone centroids (original row order):")
    print(gdf[["centroid_lat", "centroid_lon"]])

    if ZONE_LABELS_OVERRIDE is not None and len(ZONE_LABELS_OVERRIDE) == len(gdf):
        gdf["zone_label"] = ZONE_LABELS_OVERRIDE
        gdf_sorted = gdf.reset_index(drop=True)
    else:
        # Sort north -> south by centroid latitude and assign NW / SW / S
        gdf_sorted = gdf.sort_values("centroid_lat", ascending=False).reset_index(drop=True)
        default_names = ["NW", "SW", "S"]
        if len(gdf_sorted) == len(default_names):
            gdf_sorted["zone_label"] = default_names
        else:
            gdf_sorted["zone_label"] = [f"zone_{i+1}" for i in range(len(gdf_sorted))]

    print("Assigned zone labels (verify these look right for your coastline!):")
    print(gdf_sorted[["zone_label", "centroid_lat", "centroid_lon"]])

    # regionmask needs a GeoDataFrame with one row per region
    regions = regionmask.from_geopandas(gdf_sorted, names="zone_label")

    # 2D mask: each cell gets the *index* (0,1,2,...) of the zone it falls
    # in, or NaN if it's outside every zone. Index order matches gdf_sorted.
    zone_mask_2d = regions.mask(sst_full.lon, sst_full.lat)

    zone_labels = gdf_sorted["zone_label"].tolist()  # ["NW", "SW", "S"]

    sst_full.isel(time=0).plot()
    plt.title("Full domain SST (t=0)")
    plt.show()

    zone_mask_2d.plot()
    plt.title("Zone index mask")
    plt.show()

mask_ocean = 1 * np.ones(sst_full.shape[1:]) * np.isfinite(sst_full.isel(time=0))
mask_land = 0 * np.ones(sst_full.shape[1:]) * np.isnan(sst_full.isel(time=0))
mask_full = mask_ocean + mask_land
mask_full.plot()
plt.title("Land/ocean mask")
plt.show()

if choice_date == 1:
    print("choice_date")
    if start_mth <= end_mth:
        mask_date = ((sst_full['time.month'] > start_mth) & (sst_full['time.month'] < end_mth)) | \
                    ((sst_full['time.month'] == start_mth) & (sst_full['time.day'] >= start_day)) | \
                    ((sst_full['time.month'] == end_mth) & (sst_full['time.day'] <= end_day))
    else:
        mask_date = ((sst_full['time.month'] > start_mth) | (sst_full['time.month'] < end_mth)) | \
                    ((sst_full['time.month'] == start_mth) & (sst_full['time.day'] >= start_day)) | \
                    ((sst_full['time.month'] == end_mth) & (sst_full['time.day'] <= end_day))
    sst_full = sst_full.sel(time=mask_date)

year_range = range(yearbeg, yearend + 1)

# ============================================================
# Iberian Peninsula basemap helper
# ============================================================
# lon_min, lon_max, lat_min, lat_max -- adjust if your zones extend
# outside this box (e.g. widen lon_min if you need more of Spain's
# Atlantic coast, or narrow it if you only need Portugal's west coast).
IBERIA_EXTENT = [-11.0, -5.5, 35.8, 44.0]


def make_iberia_map_grid(nrows, ncols, figsize):
    """
    Create a grid of cartopy GeoAxes, each pre-set to the Iberian
    Peninsula extent with coastlines, land shading, and country borders.
    Returns (fig, flattened_axes_array).
    """
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                              subplot_kw={'projection': ccrs.PlateCarree()})
    axes = np.atleast_1d(axes).flatten()
    for ax in axes:
        ax.set_extent(IBERIA_EXTENT, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, facecolor='0.85', zorder=2)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.6, zorder=3)
        ax.add_feature(cfeature.BORDERS, linewidth=0.4, linestyle=':', zorder=3)
    return fig, axes




# ============================================================
# Everything below runs ONCE PER ZONE
# ============================================================
def build_zone_month_summary(mhw_meanintensity_filtered, cell_area, event_df, zone_label):
    """
    Build a tidy (long-format) monthly summary table for one zone with columns:
    Year | Month | Zone | Parameter | Value

    Parameters included:
      - mean_intensity_degC     (field-based, from mhw_meanintensity_filtered)
      - activity_degC_days_m2   (field-based, from mhw_meanintensity_filtered * cell_area)
      - number_of_events        (event-based, from event_df, counted by event start month)
      - mean_duration_days      (event-based, from event_df)
      - mean_area_km2           (event-based, from event_df)
      - sum_area_km2            (event-based, from event_df)
    """
    rows = []

    # ---- Field-based: mean intensity ----
    intensity_da = mhw_meanintensity_filtered.where(mhw_meanintensity_filtered != 0)
    intensity_ts = intensity_da.mean(dim=('lat', 'lon'), skipna=True).to_series()
    if len(intensity_ts) > 0:
        intensity_monthly = intensity_ts.groupby([intensity_ts.index.year, intensity_ts.index.month]).mean()
        for (yr, mo), val in intensity_monthly.items():
            rows.append({"Year": yr, "Month": mo, "Zone": zone_label,
                         "Parameter": "mean_intensity_degC", "Value": val})

    # ---- Field-based: activity (intensity x cell area, summed over days in month) ----
    activity_da = mhw_meanintensity_filtered.fillna(0) * cell_area
    activity_ts = activity_da.sum(dim=('lat', 'lon'), skipna=True).to_series()
    if len(activity_ts) > 0:
        activity_monthly = activity_ts.groupby([activity_ts.index.year, activity_ts.index.month]).sum()
        for (yr, mo), val in activity_monthly.items():
            rows.append({"Year": yr, "Month": mo, "Zone": zone_label,
                         "Parameter": "activity_degC_days_m2", "Value": val})

    # ---- Event-based: number of events, mean duration, mean/sum area ----
    if not event_df.empty:
        ev = event_df.copy()
        ev["start_time"] = pd.to_datetime(ev["start_time"])
        ev["Year"] = ev["start_time"].dt.year
        ev["Month"] = ev["start_time"].dt.month

        n_events = ev.groupby(["Year", "Month"]).size()
        mean_duration = ev.groupby(["Year", "Month"])["duration"].mean()
        mean_area = ev.groupby(["Year", "Month"])["area"].mean() / 1e6   # m2 -> km2
        sum_area = ev.groupby(["Year", "Month"])["area"].sum() / 1e6     # m2 -> km2

        for (yr, mo), val in n_events.items():
            rows.append({"Year": yr, "Month": mo, "Zone": zone_label,
                         "Parameter": "number_of_events", "Value": val})
        for (yr, mo), val in mean_duration.items():
            rows.append({"Year": yr, "Month": mo, "Zone": zone_label,
                         "Parameter": "mean_duration_days", "Value": val})
        for (yr, mo), val in mean_area.items():
            rows.append({"Year": yr, "Month": mo, "Zone": zone_label,
                         "Parameter": "mean_area_km2", "Value": val})
        for (yr, mo), val in sum_area.items():
            rows.append({"Year": yr, "Month": mo, "Zone": zone_label,
                         "Parameter": "sum_area_km2", "Value": val})

    return pd.DataFrame(rows)


def run_mhw_analysis_for_zone(sst, mask, zone_id, zone_label):
    """
    Run the full MHW detection + plotting pipeline on `sst`, which has
    already been masked down to a single zone. `zone_label` is used to
    tag output filenames and plot titles.
    """
    print("=" * 70)
    print(f"ZONE {zone_id}: {zone_label}")
    print("=" * 70)

    region_tag = f"{region}.{zone_label}"

    print("TIME MHW_CALCULATION - START", (time.time() - start_time_script) / 60)
    print("np.shape(sst)", np.shape(sst))

    # --------------------------------------------------------
    # Detrending
    # --------------------------------------------------------
    if detrend == 1:
        time_axis = np.arange(len(sst['time']))
        sst_reshaped = sst.values.reshape(sst.sizes['time'], -1)
        trends = np.polyfit(time_axis, sst_reshaped, detrend)
        trend_lines = trends[0, :] * time_axis[:, np.newaxis] + trends[1, :]
        trend_lines_reshaped = trend_lines.reshape(sst.sizes['time'], sst.sizes['lat'], sst.sizes['lon'])
        sst_detrended = sst.values - trend_lines_reshaped
        sst = xr.DataArray(
            sst_detrended,
            dims=sst.dims,
            coords=sst.coords,
            attrs=sst.attrs,
            name=sst.name
        )

    # --------------------------------------------------------
    # Threshold
    # --------------------------------------------------------
    sst_clim_period = sst.sel(time=slice(f"{clim_start_year}-01-01", f"{clim_end_year}-12-31"))
    climatology = sst_clim_period.groupby(dt_clim).mean().rolling(dayofyear=11, center=True, min_periods=1).mean()
    sstanomalyclim = sst.groupby(dt_clim) - climatology
    threshold = sst.groupby(dt_clim).quantile(percentile, dim='time', keep_attrs=True, skipna=True)
    mhw_meanintensity_all = sstanomalyclim.groupby(dt_clim).where(sst.groupby(dt_clim) > threshold)
    print("TIME MHW_CALCULATION - AFTER THRESHOLD", (time.time() - start_time_script) / 60)
    del climatology, threshold

    # --------------------------------------------------------
    # min_duration & min_gap
    # --------------------------------------------------------
    def filter_short_sequences(da, min_duration=min_duration_mhw - 1, min_gap=min_gap_mhw):
        filtered_da = da.copy()
        for ii in range(da.sizes['lat']):
            for jj in range(da.sizes['lon']):
                ts = da.isel(lat=ii, lon=jj)
                if ts.sizes['time'] > 0:
                    ts_values = ts.notnull().values
                    changes = np.diff(ts_values.astype(int))
                    starts = np.where(changes == 1)[0] + 1
                    ends = np.where(changes == -1)[0]
                    if ts_values[0]:
                        starts = np.insert(starts, 0, 0)
                    if ts_values[-1]:
                        ends = np.append(ends, len(ts_values) - 1)

                    if starts.size > 0:
                        good_runs = []
                        for i_start, i_end in list(zip(starts, ends)):
                            if (i_end - i_start) >= min_duration:
                                good_runs.append((i_start, i_end))
                            else:
                                filtered_da.isel(lat=ii, lon=jj)[i_start:i_end + 1] = np.nan

                        if len(good_runs) == 0:
                            continue

                        merged = []
                        current_start, current_end = good_runs[0]
                        for i_start, i_end in good_runs[1:]:
                            gap = i_start - current_end - 1
                            if gap <= min_gap:
                                filtered_da.isel(lat=ii, lon=jj)[current_end + 1: i_start] = \
                                    sstanomalyclim.isel(lat=ii, lon=jj, time=slice(current_end + 1, i_start))
                                current_end = i_end
                            else:
                                merged.append((current_start, current_end))
                                current_start, current_end = i_start, i_end
                        merged.append((current_start, current_end))
        return filtered_da

    mhw_meanintensity_filtered = filter_short_sequences(mhw_meanintensity_all)
    out_nc = (f"{pathdata}/OUTPUT/MHWmeanintensity.threshold-percentile{percentile}-clim{clim_start_year}-"
              f"{clim_end_year}.minduration{min_duration_mhw}.min_gap{min_gap_mhw}.detrend{detrend}_"
              f"{model}.{yearbeg}{yearend}.{region_tag}.start{start_day}-{start_mth}to{end_day}-{end_mth}.nc")
    mhw_meanintensity_filtered.rename("mhw_meanintensity_filtered").to_netcdf(out_nc)
    del mhw_meanintensity_filtered

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------
    print("TIME PLOT - START", (time.time() - start_time_script) / 60)
    ds_out = xr.open_dataset(out_nc)
    mhw_meanintensity_filtered = ds_out["mhw_meanintensity_filtered"]
    print("mhw_meanintensity_filtered", mhw_meanintensity_filtered.shape)

    # ---- Intensity ----
    print("TIME PLOT - INTENSITY", (time.time() - start_time_script) / 60)
    mhw_meanintensity_per_year = (
        mhw_meanintensity_filtered
        .where(mhw_meanintensity_filtered != 0)
        .groupby('time.year')
        .mean(dim='time', skipna=True)
        .mean(dim=('lon', 'lat'), skipna=True)
    )

    plt.figure(figsize=(7, 4))
    plt.plot(year_range, mhw_meanintensity_per_year, c='red', marker='o')
    plt.ylabel("Mean intensity (°C) ", fontsize=size_fonte)
    plt.title(f"Zone {zone_label} - day {start_day} mth {start_mth} to day {end_day} mth {end_mth}",
              fontsize=size_fonte)
    plt.xticks(fontsize=size_fonte)
    plt.yticks(fontsize=size_fonte)
    plt.grid(True)
    plt.show()

    # ---- MHW features (event table) ----
    events_count = xr.DataArray(
        np.zeros((len(mhw_meanintensity_filtered.time), len(mhw_meanintensity_filtered.lat),
                   len(mhw_meanintensity_filtered.lon))),
        dims=('time', 'lat', 'lon'),
        coords={'time': mhw_meanintensity_filtered.time, 'lat': mhw_meanintensity_filtered.lat,
                'lon': mhw_meanintensity_filtered.lon}
    )

    event_records = []
    for ii in range(len(mhw_meanintensity_filtered.lat)):
        for jj in range(len(mhw_meanintensity_filtered.lon)):
            intensity_series = mhw_meanintensity_filtered.isel(lat=ii, lon=jj)
            mhw_events = (intensity_series > 0).astype(int)
            event_changes = np.diff(mhw_events, prepend=0)
            event_starts = np.where(event_changes == 1)[0]
            event_ends = np.where(event_changes == -1)[0]
            if len(event_starts) > len(event_ends):
                event_ends = np.append(event_ends, len(mhw_events) - 1)
            events_count[:, ii, jj] = event_changes
            for st, en in zip(event_starts, event_ends):
                start_time = mhw_meanintensity_filtered.time.isel(time=st).item()
                end_time = mhw_meanintensity_filtered.time.isel(time=en).item()
                duration = en - st
                area_cell = float(cell_area.isel(lat=ii, lon=jj).values)
                yr = pd.Timestamp(start_time).year
                event_records.append({
                    "year": yr,
                    "duration": duration,
                    "start_time": start_time,
                    "end_time": end_time,
                    "area": area_cell,
                    "lat": mhw_meanintensity_filtered.lat.isel(lat=ii).item(),
                    "lon": mhw_meanintensity_filtered.lon.isel(lon=jj).item()
                })

    event_df = pd.DataFrame(event_records)
    print(event_df.head())
    event_df.to_csv(
        f"{pathdata}/OUTPUT/MHWevents_{model}.{yearbeg}{yearend}.{region_tag}.csv", index=False)

    # ---- number of events ----
    print("TIME PLOT - NUMBER EVENTS", (time.time() - start_time_script) / 60)
    total_events = (events_count == 1).sum(dim=('lat', 'lon'))
    events_per_year = (
        total_events
        .where(total_events != 0)
        .groupby('time.year')
        .mean(dim='time', skipna=True)
    )

    plt.figure(figsize=(7, 4))
    plt.plot(year_range, events_per_year, c='red', marker='o')
    plt.ylabel("Mean number of events (-) ", fontsize=size_fonte)
    plt.title(f"Zone {zone_label} - day {start_day} mth {start_mth} to day {end_day} mth {end_mth}",
              fontsize=size_fonte)
    plt.xticks(fontsize=size_fonte)
    plt.yticks(fontsize=size_fonte)
    plt.grid(True)
    plt.show()

    # ---- map: number of events per year ----
    min_value_leg = 0
    max_value_leg = 2e11
    n_years = yearend - yearbeg + 1
    years = range(yearbeg, yearend + 1)
    ncols = 6
    nrows = math.ceil(n_years / ncols)

    fig, axes = make_iberia_map_grid(nrows, ncols, figsize=(ncols * 3, nrows * 2.5))
    im = None
    for i, yr in enumerate(years):
        startdate = f'{yr}-{start_mth}-{start_day}'
        enddate = f'{yr}-{end_mth}-{end_day}' if start_mth <= end_mth else f'{yr+1}-{end_mth}-{end_day}'
        nb_events_plot = (mhw_meanintensity_filtered > 0).sel(time=slice(startdate, enddate)).sum(dim='time')
        ax = axes[i]
        im = nb_events_plot.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap='Reds',
                                             add_colorbar=False, zorder=1)
        ax.set_title(f'{yr}', fontsize=12)
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    fig.subplots_adjust(bottom=0.15)
    if im is not None:
        cbar = fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.05, pad=0.05)
        cbar.set_label('Number of events', fontsize=12)
    fig.suptitle(f'MHW number of events, zone {zone_label} ({start_mth}-{start_day} - {end_mth}-{end_day}) '
                 f'from {yearbeg} to {yearend}', fontsize=12)
    plt.show()

    # ---- duration ----
    print("TIME PLOT - DURATION", (time.time() - start_time_script) / 60)
    durations_per_year_mean = event_df.groupby("year")["duration"].mean().reindex(year_range, fill_value=np.nan)

    plt.figure(figsize=(7, 4))
    plt.plot(durations_per_year_mean.index, durations_per_year_mean.values, c='red', marker='o')
    plt.ylabel("Mean Duration (days)", fontsize=size_fonte)
    plt.title(f"Zone {zone_label} - {start_day}/{start_mth} to {end_day}/{end_mth}", fontsize=size_fonte)
    plt.grid(True)
    plt.show()

    min_value_leg = 0
    max_value_leg = 60
    fig, axes = make_iberia_map_grid(nrows, ncols, figsize=(ncols * 3, nrows * 2.5))
    im = None
    for i, yr in enumerate(years):
        ax = axes[i]
        duration_2d = (
            event_df[event_df["year"] == yr]
            .groupby(["lat", "lon"])["duration"]
            .mean()
            .unstack(level="lon")
        )
        if duration_2d.empty or duration_2d.isna().all().all():
            ax.set_title(f'{yr} (no data)', fontsize=12)
            continue
        da_plot = xr.DataArray(duration_2d.values, dims=("lat", "lon"),
                                coords={"lat": duration_2d.index, "lon": duration_2d.columns})
        im = da_plot.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap='Reds',
                                      vmin=min_value_leg, vmax=max_value_leg, add_colorbar=False, zorder=1)
        ax.set_title(f'{yr}', fontsize=12)
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    fig.suptitle(f'MHW duration, zone {zone_label} ({start_mth}-{start_day} - {end_mth}-{end_day}) '
                 f'from {yearbeg} to {yearend}', fontsize=12)
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.02])
    if im is not None:
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar_ax.set_xlabel('days', fontsize=12)
    plt.show()

    # ---- area ----
    print("TIME PLOT - AREA", (time.time() - start_time_script) / 60)
    area_per_year_mean = event_df.groupby("year")["area"].mean().reindex(year_range, fill_value=0)
    area_per_year_sum = event_df.groupby("year")["area"].sum().reindex(year_range, fill_value=0)

    plt.figure(figsize=(7, 4))
    plt.plot(area_per_year_mean.index, area_per_year_mean.values / 1e6, c='red', marker='o')
    plt.ylabel("Mean area (km²)", fontsize=size_fonte)
    plt.title(f"Zone {zone_label} - {start_day}/{start_mth} to {end_day}/{end_mth}", fontsize=size_fonte)
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(7, 4))
    plt.plot(area_per_year_sum.index, area_per_year_sum.values / 1e6, c='red', marker='o')
    plt.ylabel("Sum area (km2)", fontsize=size_fonte)
    plt.title(f"Zone {zone_label} - {start_day}/{start_mth} to {end_day}/{end_mth}", fontsize=size_fonte)
    plt.grid(True)
    plt.show()

    fig, axes = make_iberia_map_grid(nrows, ncols, figsize=(ncols * 3, nrows * 2.5))
    da_plot = None
    for i, yr in enumerate(years):
        ax = axes[i]
        area_2d = (
            event_df[event_df["year"] == yr]
            .groupby(["lat", "lon"])["area"]
            .mean()
            .unstack(level="lon")
        )
        if area_2d.empty or area_2d.isna().all().all():
            ax.set_title(f'{yr} (no data)', fontsize=12)
            continue
        da_plot = xr.DataArray(area_2d.values, dims=("lat", "lon"),
                                coords={"lat": area_2d.index, "lon": area_2d.columns})
        da_plot.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap='Reds',
                                 add_colorbar=False, zorder=1)
        ax.set_title(f'{yr}', fontsize=12)
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    fig.subplots_adjust(bottom=0.15)
    if da_plot is not None:
        cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.02])
        da_plot.plot.pcolormesh(cmap='Reds', add_colorbar=True, cbar_ax=cbar_ax, cbar_kwargs={'orientation': 'horizontal'})
        cbar_ax.set_title('m²', fontsize=12)
        cbar_ax.set_xlabel('')
        cbar_ax.set_ylabel('')
    fig.suptitle(f'MHW area, zone {zone_label} ({start_mth}-{start_day} - {end_mth}-{end_day}) '
                 f'from {yearbeg} to {yearend}', fontsize=12)
    plt.show()

    # ---- activity ----
    print("TIME PLOT - ACTIVITY", (time.time() - start_time_script) / 60)
    mhw_meanintensity_weighted = mhw_meanintensity_filtered * cell_area
    activity_sum_per_year = mhw_meanintensity_weighted.groupby('time.year').sum(dim='time')
    total_activity = activity_sum_per_year.sum(dim=['lon', 'lat'])

    model_1 = np.polyfit(year_range, total_activity, 2)
    predict_1 = np.poly1d(model_1)
    y_lin_reg_1 = predict_1(year_range)

    plt.figure(figsize=(7, 4))
    plt.plot(year_range, y_lin_reg_1 / 10 ** 12, c='red', linestyle='--')
    plt.plot(year_range, total_activity / 10 ** 12, c='red', marker='o')
    plt.ylabel("MHWs Activity (°C.days.10⁶km²) ", fontsize=size_fonte)
    plt.title(f"Zone {zone_label} - day {start_day} mth {start_mth} to day {end_day} mth {end_mth}",
              fontsize=size_fonte)
    plt.xticks(fontsize=size_fonte)
    plt.yticks(fontsize=size_fonte)
    plt.grid(True)
    plt.show()

    min_value_leg = 0
    max_value_leg = 2e11
    fig, axes = make_iberia_map_grid(nrows, ncols, figsize=(ncols * 3, nrows * 2.5))
    activity_plot = None
    for i, yr in enumerate(years):
        startdate = f'{yr}-{start_mth}-{start_day}'
        enddate = f'{yr}-{end_mth}-{end_day}' if start_mth <= end_mth else f'{yr+1}-{end_mth}-{end_day}'
        activity_plot = mhw_meanintensity_filtered.sel(time=slice(startdate, enddate)).sum(dim='time') * cell_area
        ax = axes[i]
        activity_plot.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap='Reds', add_colorbar=False,
                                       vmin=min_value_leg, vmax=max_value_leg, zorder=1)
        ax.set_title(f'{yr}', fontsize=12)
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    fig.subplots_adjust(bottom=0.15)
    if activity_plot is not None:
        cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.02])
        activity_plot.plot.pcolormesh(cmap='Reds', add_colorbar=True, cbar_ax=cbar_ax,
                                       cbar_kwargs={'orientation': 'horizontal'}, vmin=min_value_leg, vmax=max_value_leg)
        cbar_ax.set_title('C.days.m²', fontsize=12)
        cbar_ax.set_xlabel('')
        cbar_ax.set_ylabel('')
    fig.suptitle(f'MHW Activity, zone {zone_label} ({start_mth}-{start_day} - {end_mth}-{end_day}) '
                 f'from {yearbeg} to {yearend}', fontsize=12)
    plt.show()

    summary_df = build_zone_month_summary(mhw_meanintensity_filtered, cell_area, event_df, zone_label)

    return event_df, summary_df


# ============================================================
# Loop over the 3 zones
# ============================================================
all_events = {}
all_summaries = []
for zone_id, zone_label in enumerate(zone_labels):
    sst_zone = sst_full.where(zone_mask_2d == zone_id)
    mask_zone = mask_full.where(zone_mask_2d == zone_id)
    event_df_zone, summary_df_zone = run_mhw_analysis_for_zone(sst_zone, mask_zone, zone_id, zone_label)
    all_events[zone_label] = event_df_zone
    all_summaries.append(summary_df_zone)

# ---- Combine all zones into one tidy summary table ----
summary_table = pd.concat(all_summaries, ignore_index=True)
summary_table = summary_table.sort_values(["Zone", "Year", "Month", "Parameter"]).reset_index(drop=True)

summary_csv_long = f"{pathdata}/OUTPUT/MHW_summary_table_long_{model}.{yearbeg}{yearend}.{region}.csv"
summary_table.to_csv(summary_csv_long, index=False)
print(f"Saved long-format summary table to: {summary_csv_long}")
print(summary_table.head(20))

# ---- Optional: wide format, one column per parameter ----
summary_table_wide = (
    summary_table
    .pivot_table(index=["Year", "Month", "Zone"], columns="Parameter", values="Value")
    .reset_index()
)
summary_csv_wide = f"{pathdata}/OUTPUT/MHW_summary_table_wide_{model}.{yearbeg}{yearend}.{region}.csv"
summary_table_wide.to_csv(summary_csv_wide, index=False)
print(f"Saved wide-format summary table to: {summary_csv_wide}")
print(summary_table_wide.head(20))

# End
duration_script = time.time() - start_time_script
print(f"Script executed in {int(duration_script / 60)} minutes or {int(duration_script / 3600)} hours.")