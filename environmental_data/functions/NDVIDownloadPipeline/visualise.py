"""
Phase 5: Visualisation & validation for the NDVI data cube.

Usage:
    python visualise.py                     # run all plots
    python visualise.py --from netcdf       # load from NetCDF (default)
    python visualise.py --from zarr         # load from Zarr store
"""
import argparse
import os

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import config

import pipeline

OUTPUT_PLOTS_DIR = os.path.join(config.OUTPUT_DIR, "plots")


def load_ndvi(source: str) -> xr.DataArray:
    """Load the NDVI data cube from disk."""
    if source == "netcdf":
        path = os.path.join(config.OUTPUT_DIR, "ndvi.nc")
        return xr.open_dataarray(path, chunks=config.CHUNKS)
    elif source == "zarr":
        path = os.path.join(config.OUTPUT_DIR, "ndvi.zarr")
        return xr.open_zarr(path, chunks=config.CHUNKS)["NDVI"]
    else:
        raise ValueError(f"Unsupported source: {source!r}")


def load_aoi_projected():
    """Load AOI boundary reprojected to the output CRS for overlay."""
    gdf = gpd.read_file(config.SHAPEFILE_PATH)
    return gdf.to_crs(config.OUTPUT_CRS)


# ── Plot functions ───────────────────────────────────────────────────────────

def spot_check_map(ndvi: xr.DataArray, aoi: gpd.GeoDataFrame):
    """Plot a single time step to verify spatial extent, masking, and CRS."""
    t = ndvi.time.values[len(ndvi.time) // 2]  # pick a middle date
    date_str = str(np.datetime_as_string(t, unit="D"))

    fig, ax = plt.subplots(figsize=(10, 10))
    ndvi.sel(time=t).plot(
        ax=ax, cmap="RdYlGn", vmin=-0.2, vmax=0.9,
        cbar_kwargs={"label": "NDVI"},
    )
    aoi.boundary.plot(ax=ax, edgecolor="black", linewidth=1.5)
    ax.set_title(f"NDVI spot-check — {date_str}")
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")

    path = os.path.join(OUTPUT_PLOTS_DIR, "spot_check.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def time_series_plot(ndvi: xr.DataArray):
    """Plot spatial-mean NDVI over time with percentile bands."""
    # Load fully into memory — quantile requires single-chunk spatial dims
    ndvi_loaded = ndvi.load()
    mean = ndvi_loaded.mean(dim=["x", "y"])
    p10 = ndvi_loaded.quantile(0.1, dim=["x", "y"])
    p90 = ndvi_loaded.quantile(0.9, dim=["x", "y"])

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.fill_between(mean.time.values, p10.values, p90.values, alpha=0.25, color="green", label="10th–90th percentile")
    ax.plot(mean.time.values, mean.values, color="green", linewidth=1, label="Spatial mean")
    ax.set_ylabel("NDVI")
    ax.set_xlabel("Date")
    ax.set_title("NDVI Time Series — Doñana")
    ax.legend()
    ax.grid(True, alpha=0.3)

    path = os.path.join(OUTPUT_PLOTS_DIR, "time_series.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def multi_panel_map(ndvi: xr.DataArray):
    """Small-multiples grid showing NDVI at evenly spaced dates."""
    n_panels = min(12, len(ndvi.time))
    indices = np.linspace(0, len(ndvi.time) - 1, n_panels, dtype=int)
    cols = 4
    rows = int(np.ceil(n_panels / cols))

    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
    axes = axes.flatten()

    for i, idx in enumerate(indices):
        t = ndvi.time.values[idx]
        date_str = str(np.datetime_as_string(t, unit="D"))
        ndvi.sel(time=t).plot(
            ax=axes[i], cmap="RdYlGn", vmin=-0.2, vmax=0.9,
            add_colorbar=False,
        )
        axes[i].set_title(date_str, fontsize=9)
        axes[i].set_xlabel("")
        axes[i].set_ylabel("")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("NDVI Multi-Panel Overview", fontsize=13)
    fig.tight_layout()

    path = os.path.join(OUTPUT_PLOTS_DIR, "multi_panel.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def histogram_plot(ndvi: xr.DataArray):
    """Per-date NDVI histograms stacked as a 2-D heatmap (time vs NDVI bin)."""
    n_dates = len(ndvi.time)
    bins = np.linspace(-1, 1, 101)
    hist_matrix = np.empty((n_dates, len(bins) - 1))

    for i in range(n_dates):
        vals = ndvi.isel(time=i).values.ravel()
        vals = vals[~np.isnan(vals)]
        hist_matrix[i, :], _ = np.histogram(vals, bins=bins, density=True)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.imshow(
        hist_matrix.T, origin="lower", aspect="auto",
        extent=[0, n_dates, -1, 1], cmap="viridis",
    )
    # label x-axis with a few dates
    tick_idx = np.linspace(0, n_dates - 1, min(8, n_dates), dtype=int)
    ax.set_xticks(tick_idx)
    ax.set_xticklabels(
        [str(np.datetime_as_string(ndvi.time.values[j], unit="D")) for j in tick_idx],
        rotation=45, ha="right", fontsize=8,
    )
    ax.set_ylabel("NDVI")
    ax.set_title("NDVI Distribution Over Time")

    path = os.path.join(OUTPUT_PLOTS_DIR, "histogram.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def interactive_explorer(ndvi: xr.DataArray):
    """Launch an interactive hvplot map (requires hvplot)."""
    try:
        import hvplot.xarray  # noqa: F401
    except ImportError:
        print("  hvplot not installed — skipping interactive explorer.")
        return

    plot = ndvi.hvplot.image(
        x="x", y="y", groupby="time",
        cmap="RdYlGn", clim=(-0.2, 0.9),
        width=700, height=600,
        title="NDVI Interactive Explorer",
    )
    path = os.path.join(OUTPUT_PLOTS_DIR, "interactive.html")
    import holoviews as hv
    hv.save(plot, path)
    print(f"  Saved {path}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Visualise the NDVI data cube.")
    parser.add_argument(
        "--from", dest="source", default="netcdf",
        choices=["netcdf", "zarr"],
        help="Which output store to load (default: netcdf)",
    )
    parser.add_argument(
        "--output-dir", default=None, metavar="DIR",
        help="Output directory (default: from config.py)",
    )
    
    args = parser.parse_args()
    if args.output_dir is not None:
        config.OUTPUT_DIR = args.output_dir
    
    os.makedirs(OUTPUT_PLOTS_DIR, exist_ok=True)
    print(f"Loading NDVI from {args.source}...")
    ndvi = load_ndvi(args.source)
    aoi = load_aoi_projected()

    print("Generating plots:")
    spot_check_map(ndvi, aoi)
    time_series_plot(ndvi)
    multi_panel_map(ndvi)
    histogram_plot(ndvi)
    interactive_explorer(ndvi)
    print("Done.")


if __name__ == "__main__":
    main()
