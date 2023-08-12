import argparse
import re
import time
import shutil
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import ska_dbi

from .paths import FID_STATS_PATH

EXCLUDE_OBSIDS = "(2010, 2783, 1431, 1411)"


def get_fid_stats(db, det):
    query = (
        "select obsid, id_num, id_string, tstart, ang_y_med, ang_z_med, sim_z_offset"
        " FROM fid_stats"
        ' WHERE proc_status IS NULL AND id_string LIKE "{}%" '
        " and obsid not in {} ".format(det, EXCLUDE_OBSIDS)
    )
    vals = db.fetchall(query)
    return vals


def plotfids(detstats, det, data_dir):
    year0 = 0.0
    fids = {
        "ACIS-S": [1, 2, 3, 4, 5, 6],
        "ACIS-I": [1, 2, 3, 4, 5, 6],
        "HRC-I": [7, 8, 9, 10],
        "HRC-S": [11, 12, 13, 14],
    }
    fidcolor = " bgrcmkbgrcbgrc"
    sim_z_nom = np.median(detstats["sim_z_offset"])  # median for this detector
    plt.figure(1, figsize=(6, 4.5))
    plt.clf()
    plt.subplot(2, 1, 1)
    plt.ylabel("Y offset (arcsec)")
    maxyear = int(time.strftime("%Y")) + 1
    DEFAULT_MAX = 10  # arcsec
    DEFAULT_MIN = -25  # arcsec
    plt.axis([1999, maxyear, DEFAULT_MIN, DEFAULT_MAX])
    plt.title(det + " Fid Drift")
    plt.grid()
    plt.subplot(2, 1, 2)
    plt.xlabel("Year")  # MET (years)')
    plt.ylabel("Z offset (arcsec)")
    plt.axis([1999, maxyear, DEFAULT_MIN, 10])
    plt.grid()
    min_plot_y = DEFAULT_MIN
    max_plot_y = DEFAULT_MAX
    for fid in fids[det]:
        fidstats = detstats[detstats["id_num"] == fid]
        year = fidstats["tstart"] / 86400.0 / 365.25 + 1998.0
        normmask = np.logical_and(year > 2010.0, year < 2020.0)
        fidstatsnorm = fidstats[normmask]
        if len(fidstatsnorm) < 3:
            print("Fid %s %d" % (det, fid))
            print("Not enough readouts for the median")
            continue
        y0 = np.median(fidstatsnorm["ang_y_med"])
        z0 = np.median(
            fidstatsnorm["ang_z_med"] + fidstatsnorm.sim_z_offset - sim_z_nom
        )
        fid_min_y = np.min(fidstats["ang_y_med"] - y0)
        fid_min_z = np.min(
            fidstats["ang_z_med"] - z0 + fidstats["sim_z_offset"] - sim_z_nom
        )
        fid_max_y = np.max(fidstats["ang_y_med"] - y0)
        fid_max_z = np.max(
            fidstats["ang_z_med"] - z0 + fidstats["sim_z_offset"] - sim_z_nom
        )
        fid_min = np.min([fid_min_y, fid_min_z])
        fid_max = np.max([fid_max_y, fid_max_z])
        # if the plot min value is less than we've seen for this detector
        # or less than the default, update the minimum (rounded down to the nearest 5)
        if fid_min < min_plot_y:
            min_plot_y = fid_min - (fid_min % 5)
        if fid_max > max_plot_y:
            max_plot_y = fid_max + (5 - fid_max % 5)
        plt.subplot(2, 1, 1)
        plt.plot(
            year - year0,
            fidstats["ang_y_med"] - y0,
            ".",
            markersize=4,
            markerfacecolor=fidcolor[fid],
            mew=0,
            scaley=False,
            scalex=False,
        )
        plt.ylim(ymin=min_plot_y, ymax=max_plot_y)
        plt.subplot(2, 1, 2)
        plt.plot(
            year - year0,
            fidstats["ang_z_med"] - z0 + fidstats["sim_z_offset"] - sim_z_nom,
            ".",
            markersize=4,
            markerfacecolor=fidcolor[fid],
            mew=0,
            scaley=False,
            scalex=False,
        )
        plt.ylim(ymin=min_plot_y, ymax=max_plot_y)

    det = re.sub(r"-", "_", det.lower())
    plt.savefig(Path(data_dir) / f"drift_{det}.png")


def parse_args():
    parser = argparse.ArgumentParser(description="Plot drift data")
    parser.add_argument(
        "--data-dir", type=str, default=".", help="Fid drift data directory"
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    matplotlib.use("Agg")

    dets = ("ACIS-S", "ACIS-I", "HRC-S", "HRC-I")

    with ska_dbi.DBI(dbi="sqlite", server=FID_STATS_PATH(args.data_dir)) as dbh:
        for det in dets:
            # Some filtering here?
            detstats = get_fid_stats(dbh, det)
            plotfids(detstats, det, args.data_dir)

    # Copy the index template into the output directory
    shutil.copyfile(
        Path(__file__).parent / "data" / "index_template.html",
        Path(args.data_dir) / "index.html")


if __name__ == "__main__":
    main()
