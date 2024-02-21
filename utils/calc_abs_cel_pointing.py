"""
Compute the 99th percentile radius for fid light offsets over the last year.

  % ipython --pylab

  >>> run calc_abs_cel_pointing.py
  Fid ACIS-S 3
  Not enough readouts for the median
  99th percentile radius for ACIS-S fid light offset from median is 12.1
  Fid ACIS-I 3
  Not enough readouts for the median
  99th percentile radius for ACIS-I fid light offset from median is 11.4
  99th percentile radius for HRC-S fid light offset from median is 11.5
  99th percentile radius for HRC-I fid light offset from median is 13.3

  99th percentile radius for fid light offset from median is 11.8 arcsec
"""
import argparse

import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import ska_dbi
from cxotime import CxoTime

from fid_drift_mon.paths import FID_STATS_PATH


def parse_args(args):
    parser = argparse.ArgumentParser(description="Plot drift model")
    parser.add_argument(
        "--data-dir", type=str, default=".", help="Output data directory"
    )
    args = parser.parse_args(args)
    return args


def get_fid_stats(db, det, start=None):
    if start is None:
        start = CxoTime.now() - 365.0 * u.day
    tstart = CxoTime(start).secs
    query = (
        "select id_num, id_string, tstart, ang_y_med, ang_z_med, sim_z_offset"
        " FROM fid_stats"
        " WHERE proc_status IS NULL"
        '   AND id_string LIKE "{}%"'
        "   AND tstart > {}".format(det, tstart)
    )
    vals = db.fetchall(query)
    db.conn.close()

    return vals


def calc_abs_cel_pointing(detstats, det):
    year0 = 0.0
    fids = {
        "ACIS-S": [1, 2, 3, 4, 5, 6],
        "ACIS-I": [1, 2, 3, 4, 5, 6],
        "HRC-I": [7, 8, 9, 10],
        "HRC-S": [11, 12, 13, 14],
    }
    fidcolor = " bgrcmkbgrcbgrc"
    sim_z_nom = np.median(detstats["sim_z_offset"])  # median for this detector
    plt.figure(figsize=(6, 4.5))
    plt.clf()
    plt.subplot(2, 1, 1)
    plt.ylabel("Y offset (arcsec)")
    plt.title(det + " Fid Drift")
    plt.grid()
    plt.subplot(2, 1, 2)
    plt.xlabel("Year")  # MET (years)')
    plt.ylabel("Z offset (arcsec)")
    plt.grid()
    drs = []
    for fid in fids[det]:
        fidstats = detstats[detstats["id_num"] == fid]
        year = fidstats["tstart"] / 86400.0 / 365.25 + 1998.0
        fidstatsnorm = fidstats
        if len(fidstatsnorm) < 3:
            print("Fid %s %d" % (det, fid))
            print("Not enough readouts for the median")
            continue
        y0 = np.median(fidstatsnorm["ang_y_med"])
        z0 = np.median(
            fidstatsnorm["ang_z_med"] + fidstatsnorm.sim_z_offset - sim_z_nom
        )
        plt.subplot(2, 1, 1)
        dy = fidstats["ang_y_med"] - y0
        dz = fidstats["ang_z_med"] - z0 + fidstats["sim_z_offset"] - sim_z_nom
        dr = np.sqrt(dy**2 + dz**2)
        plt.plot(year - year0, dy, ",", markerfacecolor=fidcolor[fid], mew=0)
        plt.subplot(2, 1, 2)
        plt.plot(year - year0, dz, ",", markerfacecolor=fidcolor[fid], mew=0)
        drs.append(dr)

    dr = np.concatenate(drs)
    dr99 = np.percentile(dr, 99.0)
    print(
        "99th percentile radius for {} fid light offset from median is {:.1f}".format(
            det, dr99
        )
    )

    return dr


def main(args=None):
    args = parse_args(args)

    matplotlib.use("Agg")

    dets = ("ACIS-S", "ACIS-I", "HRC-S", "HRC-I")

    drs = []
    for det in dets:
        with ska_dbi.DBI(dbi="sqlite", server=FID_STATS_PATH(args.data_dir)) as dbh:
            detstats = get_fid_stats(dbh, det)
        dr = calc_abs_cel_pointing(detstats, det)
        drs.append(dr)

    dr = np.concatenate(drs)
    dr99 = np.percentile(dr, 99.0)
    print()
    print(
        "99th percentile radius for fid light offset from median is {:.1f} arcsec".format(
            dr99
        )
    )


if __name__ == "__main__":
    main()
