"""
Marine Heatwave (MHW) Descriptive Statistics & Plotting Script
================================================================
Loads MHW summary CSVs for multiple bathymetry depth cutoffs, combines them,
computes descriptive statistics, and generates comparison plots across
zones and bathymetry depths.

Input files expected (one per bathymetry depth), each with columns:
    Year, Month, Zone, Parameter, Value

Usage:
    1. Edit the INPUT_DIR and OUTPUT_DIR paths below.
    2. Run the script directly (no command-line arguments needed):
           python "Descriptive Statistics.py"

Outputs (each figure saved as both .png at 300 DPI and .pdf vector, for
publication use — set SAVE_PDF=False below to skip PDF export):
    - descriptive_statistics.csv        : full descriptive stats table
      (count, mean, std, min, 25%, 50%, 75%, max) per Parameter x Zone x Depth
    - cross_comparison_table.csv        : mean values, zone x depth, per parameter
    - MHW_all_parameters_comparison.(png|pdf) : mean +/- SD bar chart, all params
    - MHW_comparison_<PARAMETER>.(png|pdf)    : same comparison, one file per parameter
    - MHW_yearly_bar_zone_<ZONE>.(png|pdf)    : year-by-year bar chart per zone
    - MHW_yearly_<PARAMETER>_zone_<ZONE>.(png|pdf) : same, one file per zone x parameter
    - MHW_decadal_bar_zone_<ZONE>.(png|pdf)   : bar chart per zone, binned into decades
      (1981-1990, 1991-2000, 2001-2010, 2011-2020, 2021-2025)
    - MHW_decadal_<PARAMETER>_zone_<ZONE>.(png|pdf) : same, one file per zone x parameter
    - MHW_boxplots_distribution.(png|pdf)     : distribution boxplots, all params

Plotting is done entirely with matplotlib (no seaborn or other plotting libraries).
Uses the colorblind-safe Okabe-Ito palette for zones, and a sequential blue
scale for bathymetry depths.
"""

import glob
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# Publication-quality settings
# ----------------------------------------------------------------------
# - 300 DPI raster export (journal-standard minimum)
# - Vector PDF saved alongside every PNG (infinitely scalable, required by
#   most journals for line/bar charts)
# - Serif font (Times-like) matches typical journal body text; embed fonts
#   as actual text (not paths) so figures remain searchable/editable
# - SAVE_PDF = False to skip vector export and speed up iteration while
#   drafting; set True for the final publication run
SAVE_PDF = False
DPI = 300

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif", "Georgia"],
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 100,       # on-screen/preview resolution
    "savefig.dpi": DPI,      # exported raster resolution
    "pdf.fonttype": 42,      # embed fonts as text, not outlines (editable in Illustrator/Inkscape)
    "ps.fonttype": 42,
    "axes.linewidth": 0.8,
    "axes.grid": False,
})


def save_figure(fig, output_dir, filename_stem):
    """Save a figure as PNG (always) and PDF (if SAVE_PDF is True),
    at publication DPI with tight bounding boxes."""
    png_path = os.path.join(output_dir, f"{filename_stem}.png")
    fig.savefig(png_path, dpi=DPI, bbox_inches="tight")
    print(f"  Saved {png_path}")
    if SAVE_PDF:
        pdf_path = os.path.join(output_dir, f"{filename_stem}.pdf")
        fig.savefig(pdf_path, bbox_inches="tight")
        print(f"  Saved {pdf_path}")


# ----------------------------------------------------------------------
# >>> EDIT THESE TWO PATHS <<<
# ----------------------------------------------------------------------
INPUT_DIR = "/Users/miguelsilveira/Documents/GitHub/marineheatwaves/DATA/OUTPUT"
OUTPUT_DIR = "/Users/miguelsilveira/Documents/GitHub/marineheatwaves/outputs"


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
PARAMS = [
    "activity_degC_days_m2",
    "mean_intensity_degC",
    "mean_area_km2",
    "mean_duration_days",
    "number_of_events",
    "sum_area_km2",
]
PARAM_TITLES = {
    "activity_degC_days_m2": "Activity (°C·days·m²)",
    "mean_intensity_degC": "Mean Intensity (°C)",
    "mean_area_km2": "Mean Area per Event (km²)",
    "mean_duration_days": "Mean Duration (days)",
    "number_of_events": "Number of Events",
    "sum_area_km2": "Total (Sum) Area (km²)",
}
ZONE_COLORS = {"NW": "#471365", "SW": "#21918C", "S": "#BDDF26"}
DEPTH_COLORS = {
    "200m": "#440154", "400m": "#3B528B", "600m": "#21918C",
    "800m": "#5EC962", "1000m": "#FDE725",
}
# Colors above are sampled from matplotlib's 'viridis' colormap --
# perceptually uniform and colorblind-safe by design.


# ----------------------------------------------------------------------
# Data loading
# ----------------------------------------------------------------------
def load_combined_data(input_dir):
    """Find all MHW summary CSVs in input_dir, tag each with its bathymetry
    depth (parsed from the filename), and concatenate into one long DataFrame.

    Depth is detected by looking for a number immediately followed by 'm'
    in the filename (e.g. '200m', '_400m_', '600m.csv'). This works with
    any naming convention as long as each file's depth appears somewhere
    in its name, e.g.:
        MHW_..._200m_bathymetry.csv
        200m.csv
        depth_400m_summary.csv
    """
    pattern = os.path.join(input_dir, "*.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No .csv files found in {input_dir}")

    frames = []
    skipped = []
    for f in files:
        m = re.search(r"(\d{2,5})\s*m(?:[_\-.]|$)", os.path.basename(f), re.IGNORECASE)
        if not m:
            skipped.append(os.path.basename(f))
            continue
        depth_label = f"{m.group(1)}m"
        df = pd.read_csv(f)
        df["Depth"] = depth_label
        frames.append(df)
        print(f"  Loaded {os.path.basename(f)} -> Depth={depth_label}, rows={len(df)}")

    if skipped:
        print("  Skipped (no depth pattern like '200m' found in filename):")
        for s in skipped:
            print(f"    - {s}")

    if not frames:
        raise FileNotFoundError(
            f"Found {len(files)} CSV file(s) in {input_dir}, but none had a "
            "depth pattern (e.g. '200m') in the filename. Rename your files "
            "to include the depth (e.g. 'MHW_200m.csv'), or tell Claude the "
            "exact filenames so the matching pattern can be adjusted."
        )

    combined = pd.concat(frames, ignore_index=True)
    return combined


def natural_depth_order(depths):
    """Sort depth labels like ['200m','1000m','400m'] numerically."""
    return sorted(depths, key=lambda d: int(re.sub(r"\D", "", d)))


DECADE_BINS = [
    (1981, 1990, "1981-1990"),
    (1991, 2000, "1991-2000"),
    (2001, 2010, "2001-2010"),
    (2011, 2020, "2011-2020"),
    (2021, 2025, "2021-2025"),
]


def assign_decade(year):
    """Map a calendar year to one of the fixed decade labels above.
    Returns None if the year falls outside all defined ranges."""
    for start, end, label in DECADE_BINS:
        if start <= year <= end:
            return label
    return None


# ----------------------------------------------------------------------
# Descriptive statistics
# ----------------------------------------------------------------------
def compute_descriptive_stats(df):
    """Return a tidy descriptive-statistics table: one row per
    Parameter x Zone x Depth, with count/mean/std/min/quartiles/max."""
    stats = (
        df.groupby(["Parameter", "Zone", "Depth"])["Value"]
        .describe()
        .reset_index()
    )
    return stats


def compute_cross_comparison(df):
    """Mean value per Parameter, pivoted as Zone (rows) x Depth (columns)."""
    tables = {}
    mean_df = df.groupby(["Parameter", "Zone", "Depth"])["Value"].mean().reset_index()
    for param in PARAMS:
        sub = mean_df[mean_df.Parameter == param]
        depths = natural_depth_order(sub.Depth.unique())
        piv = sub.pivot(index="Zone", columns="Depth", values="Value")[depths]
        tables[param] = piv
    return tables


# ----------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------
def plot_overall_comparison(df, output_dir):
    """Bar chart (mean +/- SD) of all parameters, grouped by zone, x-axis = depth."""
    zones = [z for z in ["NW", "SW", "S"] if z in df.Zone.unique()]
    depths = natural_depth_order(df.Depth.unique())

    mean_df = df.groupby(["Parameter", "Zone", "Depth"])["Value"].mean().reset_index()
    std_df = df.groupby(["Parameter", "Zone", "Depth"])["Value"].std().reset_index()

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()
    x = np.arange(len(depths))
    width = 0.25

    for ax, param in zip(axes, PARAMS):
        m_piv = mean_df[mean_df.Parameter == param].pivot(
            index="Zone", columns="Depth", values="Value"
        ).reindex(index=zones, columns=depths)
        s_piv = std_df[std_df.Parameter == param].pivot(
            index="Zone", columns="Depth", values="Value"
        ).reindex(index=zones, columns=depths)

        for i, zone in enumerate(zones):
            vals = m_piv.loc[zone].values
            errs = s_piv.loc[zone].values
            ax.bar(
                x + (i - (len(zones) - 1) / 2) * width, vals, width,
                yerr=errs, capsize=3,
                error_kw={"elinewidth": 1, "ecolor": "black", "alpha": 0.6},
                label=zone, color=ZONE_COLORS.get(zone, None),
            )
        ax.set_title(PARAM_TITLES[param], fontsize=11, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(depths)
        ax.set_xlabel("Bathymetry depth cutoff")
        ax.grid(axis="y", alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(zones),
               bbox_to_anchor=(0.5, 1.03), fontsize=11, frameon=False)
    fig.suptitle(
        "Marine Heatwaves: Mean ± SD by Parameter, Zone, and Bathymetry Depth",
        fontsize=14, fontweight="bold", y=1.09,
    )
    plt.tight_layout()
    save_figure(fig, output_dir, "MHW_all_parameters_comparison")
    plt.close(fig)


def plot_individual_parameter_comparisons(df, output_dir):
    """Same comparison as plot_overall_comparison, but saved as one
    standalone PNG per parameter instead of a single combined grid."""
    zones = [z for z in ["NW", "SW", "S"] if z in df.Zone.unique()]
    depths = natural_depth_order(df.Depth.unique())

    mean_df = df.groupby(["Parameter", "Zone", "Depth"])["Value"].mean().reset_index()
    std_df = df.groupby(["Parameter", "Zone", "Depth"])["Value"].std().reset_index()

    x = np.arange(len(depths))
    width = 0.25

    for param in PARAMS:
        fig, ax = plt.subplots(figsize=(7, 5))

        m_piv = mean_df[mean_df.Parameter == param].pivot(
            index="Zone", columns="Depth", values="Value"
        ).reindex(index=zones, columns=depths)
        s_piv = std_df[std_df.Parameter == param].pivot(
            index="Zone", columns="Depth", values="Value"
        ).reindex(index=zones, columns=depths)

        for i, zone in enumerate(zones):
            vals = m_piv.loc[zone].values
            errs = s_piv.loc[zone].values
            ax.bar(
                x + (i - (len(zones) - 1) / 2) * width, vals, width,
                yerr=errs, capsize=3,
                error_kw={"elinewidth": 1, "ecolor": "black", "alpha": 0.6},
                label=zone, color=ZONE_COLORS.get(zone, None),
            )
        ax.set_title(PARAM_TITLES[param], fontsize=13, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(depths)
        ax.set_xlabel("Bathymetry depth cutoff")
        ax.grid(axis="y", alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(loc="best", frameon=False, fontsize=10)

        plt.tight_layout()
        save_figure(fig, output_dir, f"MHW_comparison_{param}")
        plt.close(fig)


def plot_yearly_bars_per_zone(df, output_dir):
    """Year-by-year bar chart (no grouping/binning) for each zone,
    all parameters, bars colored by bathymetry depth."""
    zones = [z for z in ["NW", "SW", "S"] if z in df.Zone.unique()]
    depths = natural_depth_order(df.Depth.unique())
    years = sorted(df.Year.unique())
    x = np.arange(len(years))
    width = 0.8 / len(depths)

    for zone in zones:
        zdf = df[df.Zone == zone]
        fig, axes = plt.subplots(2, 3, figsize=(26, 10))
        axes = axes.flatten()

        for ax, param in zip(axes, PARAMS):
            sub = zdf[zdf.Parameter == param]
            annual = sub.groupby(["Depth", "Year"])["Value"].mean().reset_index()
            piv = annual.pivot(index="Depth", columns="Year", values="Value").reindex(
                index=depths, columns=years
            )
            for i, depth in enumerate(depths):
                vals = piv.loc[depth].values
                ax.bar(
                    x + (i - (len(depths) - 1) / 2) * width, vals, width,
                    label=depth, color=DEPTH_COLORS.get(depth, None),
                )
            ax.set_title(PARAM_TITLES[param], fontsize=12, fontweight="bold")
            ax.set_xticks(x[::2])
            ax.set_xticklabels([years[i] for i in range(0, len(years), 2)],
                                rotation=90, fontsize=7)
            ax.grid(axis="y", alpha=0.3)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.set_xlim(-1, len(years))

        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=len(depths),
                   bbox_to_anchor=(0.5, 1.04), fontsize=11, frameon=False)
        fig.suptitle(
            f"Zone {zone}: All Parameters, Year by Year, by Bathymetry Depth",
            fontsize=15, fontweight="bold", y=1.1,
        )
        plt.tight_layout()
        save_figure(fig, output_dir, f"MHW_yearly_bar_zone_{zone}")
        plt.close(fig)


def plot_decadal_bars_per_zone(df, output_dir):
    """Bar chart (mean +/- SD) per zone, x-axis = fixed decade bins
    (1981-1990, 1991-2000, 2001-2010, 2011-2020, 2021-2025),
    bars colored by bathymetry depth."""
    zones = [z for z in ["NW", "SW", "S"] if z in df.Zone.unique()]
    depths = natural_depth_order(df.Depth.unique())
    decade_labels = [label for _, _, label in DECADE_BINS]

    ddf = df.copy()
    ddf["Decade"] = ddf["Year"].apply(assign_decade)
    ddf = ddf.dropna(subset=["Decade"])

    for zone in zones:
        zdf = ddf[ddf.Zone == zone]
        fig, axes = plt.subplots(2, 3, figsize=(18, 9))
        axes = axes.flatten()
        x = np.arange(len(decade_labels))
        width = 0.8 / len(depths)

        for ax, param in zip(axes, PARAMS):
            sub = zdf[zdf.Parameter == param]
            mean_p = sub.groupby(["Depth", "Decade"])["Value"].mean().reset_index()
            std_p = sub.groupby(["Depth", "Decade"])["Value"].std().reset_index()
            m_piv = mean_p.pivot(index="Depth", columns="Decade", values="Value").reindex(
                index=depths, columns=decade_labels
            )
            s_piv = std_p.pivot(index="Depth", columns="Decade", values="Value").reindex(
                index=depths, columns=decade_labels
            )
            for i, depth in enumerate(depths):
                vals = m_piv.loc[depth].values
                errs = s_piv.loc[depth].values
                ax.bar(
                    x + (i - (len(depths) - 1) / 2) * width, vals, width,
                    yerr=errs, capsize=2,
                    error_kw={"elinewidth": 0.7, "ecolor": "black", "alpha": 0.5},
                    label=depth, color=DEPTH_COLORS.get(depth, None),
                )
            ax.set_title(PARAM_TITLES[param], fontsize=11, fontweight="bold")
            ax.set_xticks(x)
            ax.set_xticklabels(decade_labels, rotation=20, ha="right", fontsize=8.5)
            ax.grid(axis="y", alpha=0.3)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=len(depths),
                   bbox_to_anchor=(0.5, 1.03), fontsize=10, frameon=False)
        fig.suptitle(
            f"Zone {zone}: All Parameters by Decade and Bathymetry Depth",
            fontsize=14, fontweight="bold", y=1.09,
        )
        plt.tight_layout()
        save_figure(fig, output_dir, f"MHW_decadal_bar_zone_{zone}")
        plt.close(fig)


def plot_individual_yearly_bars(df, output_dir):
    """Same as plot_yearly_bars_per_zone, but saved as one standalone PNG
    per (zone, parameter) combination instead of a combined grid."""
    zones = [z for z in ["NW", "SW", "S"] if z in df.Zone.unique()]
    depths = natural_depth_order(df.Depth.unique())
    years = sorted(df.Year.unique())
    x = np.arange(len(years))
    width = 0.8 / len(depths)

    for zone in zones:
        zdf = df[df.Zone == zone]
        for param in PARAMS:
            fig, ax = plt.subplots(figsize=(14, 5))
            sub = zdf[zdf.Parameter == param]
            annual = sub.groupby(["Depth", "Year"])["Value"].mean().reset_index()
            piv = annual.pivot(index="Depth", columns="Year", values="Value").reindex(
                index=depths, columns=years
            )
            for i, depth in enumerate(depths):
                vals = piv.loc[depth].values
                ax.bar(
                    x + (i - (len(depths) - 1) / 2) * width, vals, width,
                    label=depth, color=DEPTH_COLORS.get(depth, None),
                )
            ax.set_title(f"Zone {zone}: {PARAM_TITLES[param]}, Year by Year",
                         fontsize=13, fontweight="bold")
            ax.set_xticks(x[::2])
            ax.set_xticklabels([years[i] for i in range(0, len(years), 2)],
                                rotation=90, fontsize=8)
            ax.grid(axis="y", alpha=0.3)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.set_xlim(-1, len(years))
            ax.legend(loc="best", frameon=False, fontsize=9, ncol=len(depths))

            plt.tight_layout()
            save_figure(fig, output_dir, f"MHW_yearly_{param}_zone_{zone}")
            plt.close(fig)


def plot_individual_decadal_bars(df, output_dir):
    """Same as plot_decadal_bars_per_zone, but saved as one standalone PNG
    per (zone, parameter) combination instead of a combined grid."""
    zones = [z for z in ["NW", "SW", "S"] if z in df.Zone.unique()]
    depths = natural_depth_order(df.Depth.unique())
    decade_labels = [label for _, _, label in DECADE_BINS]

    ddf = df.copy()
    ddf["Decade"] = ddf["Year"].apply(assign_decade)
    ddf = ddf.dropna(subset=["Decade"])

    x = np.arange(len(decade_labels))
    width = 0.8 / len(depths)

    for zone in zones:
        zdf = ddf[ddf.Zone == zone]
        for param in PARAMS:
            fig, ax = plt.subplots(figsize=(8, 5))
            sub = zdf[zdf.Parameter == param]
            mean_p = sub.groupby(["Depth", "Decade"])["Value"].mean().reset_index()
            std_p = sub.groupby(["Depth", "Decade"])["Value"].std().reset_index()
            m_piv = mean_p.pivot(index="Depth", columns="Decade", values="Value").reindex(
                index=depths, columns=decade_labels
            )
            s_piv = std_p.pivot(index="Depth", columns="Decade", values="Value").reindex(
                index=depths, columns=decade_labels
            )
            for i, depth in enumerate(depths):
                vals = m_piv.loc[depth].values
                errs = s_piv.loc[depth].values
                ax.bar(
                    x + (i - (len(depths) - 1) / 2) * width, vals, width,
                    yerr=errs, capsize=2,
                    error_kw={"elinewidth": 0.7, "ecolor": "black", "alpha": 0.5},
                    label=depth, color=DEPTH_COLORS.get(depth, None),
                )
            ax.set_title(f"Zone {zone}: {PARAM_TITLES[param]}, by Decade",
                         fontsize=13, fontweight="bold")
            ax.set_xticks(x)
            ax.set_xticklabels(decade_labels, rotation=20, ha="right", fontsize=9)
            ax.grid(axis="y", alpha=0.3)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.legend(loc="best", frameon=False, fontsize=9, ncol=len(depths))

            plt.tight_layout()
            save_figure(fig, output_dir, f"MHW_decadal_{param}_zone_{zone}")
            plt.close(fig)


def plot_boxplots(df, output_dir):
    """Boxplots (pure matplotlib) showing the full distribution of each
    parameter, one subplot per parameter, boxes grouped by zone within
    each bathymetry depth."""
    zones = [z for z in ["NW", "SW", "S"] if z in df.Zone.unique()]
    depths = natural_depth_order(df.Depth.unique())

    fig, axes = plt.subplots(2, 3, figsize=(17, 9))
    axes = axes.flatten()
    n_zones = len(zones)
    group_width = 0.8
    box_width = group_width / n_zones

    for ax, param in zip(axes, PARAMS):
        sub = df[df.Parameter == param]
        positions, data, colors = [], [], []
        for di, depth in enumerate(depths):
            for zi, zone in enumerate(zones):
                vals = sub[(sub.Depth == depth) & (sub.Zone == zone)]["Value"].dropna().values
                data.append(vals)
                pos = di + (zi - (n_zones - 1) / 2) * box_width
                positions.append(pos)
                colors.append(ZONE_COLORS.get(zone, "#888888"))

        bp = ax.boxplot(
            data, positions=positions, widths=box_width * 0.9,
            patch_artist=True, showfliers=False,
            medianprops={"color": "black", "linewidth": 1.2},
        )
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.8)
            patch.set_edgecolor("black")
            patch.set_linewidth(0.5)

        ax.set_title(PARAM_TITLES[param], fontsize=11, fontweight="bold")
        ax.set_xticks(range(len(depths)))
        ax.set_xticklabels(depths)
        ax.set_xlabel("Bathymetry depth cutoff")
        ax.grid(axis="y", alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Manual legend (boxplot doesn't auto-generate one)
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=ZONE_COLORS.get(z, "#888888"), alpha=0.8, label=z,
              edgecolor="black", linewidth=0.5)
        for z in zones
    ]
    fig.legend(handles=legend_handles, loc="upper center", ncol=len(zones),
               bbox_to_anchor=(0.5, 1.03), fontsize=11, frameon=False)
    fig.suptitle(
        "Marine Heatwaves: Distribution (Boxplots) by Parameter, Zone, and Bathymetry Depth",
        fontsize=14, fontweight="bold", y=1.09,
    )
    plt.tight_layout()
    save_figure(fig, output_dir, "MHW_boxplots_distribution")
    plt.close(fig)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading data...")
    df = load_combined_data(INPUT_DIR)
    print(f"Combined dataset: {len(df)} rows, "
          f"{df.Zone.nunique()} zones, {df.Depth.nunique()} depths, "
          f"{df.Parameter.nunique()} parameters, years {df.Year.min()}-{df.Year.max()}")

    print("\nComputing descriptive statistics...")
    stats = compute_descriptive_stats(df)
    stats_path = os.path.join(OUTPUT_DIR, "descriptive_statistics.csv")
    stats.to_csv(stats_path, index=False)
    print(f"  Saved {stats_path}  ({len(stats)} rows)")

    print("\nComputing cross-comparison tables (Zone x Depth, per parameter)...")
    cross_tables = compute_cross_comparison(df)
    cross_rows = []
    for param, piv in cross_tables.items():
        piv_out = piv.copy()
        piv_out.insert(0, "Parameter", param)
        cross_rows.append(piv_out.reset_index())
    cross_df = pd.concat(cross_rows, ignore_index=True)
    cross_path = os.path.join(OUTPUT_DIR, "cross_comparison_table.csv")
    cross_df.to_csv(cross_path, index=False)
    print(f"  Saved {cross_path}")

    print("\nGenerating plots...")
    plot_overall_comparison(df, OUTPUT_DIR)
    plot_individual_parameter_comparisons(df, OUTPUT_DIR)
    plot_yearly_bars_per_zone(df, OUTPUT_DIR)
    plot_individual_yearly_bars(df, OUTPUT_DIR)
    plot_decadal_bars_per_zone(df, OUTPUT_DIR)
    plot_individual_decadal_bars(df, OUTPUT_DIR)
    plot_boxplots(df, OUTPUT_DIR)

    print("\nDone.")


if __name__ == "__main__":
    main()