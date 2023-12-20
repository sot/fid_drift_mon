#!/usr/bin/env python
"""
Plot positions of fid lights from commands versus the actual seen in telemetry.
::

  Usage: plot_commands_vs_telem.py [-h] [--start START] [--stop STOP]
                                    [--out OUT]

  Commanded vs. telemetry fid positions

  optional arguments:
    -h, --help     show this help message and exit
    --start START  Start date (default=NOW - 90 days)
    --stop STOP    Stop date (default=NOW)
    --out OUT      Output plot file
"""

import argparse
from collections import OrderedDict

import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import ska_dbi
from astropy import table
from astropy.table import Table
from cxotime import CxoTime
from kadi import events
from kadi.commands import observations
from pyyaks.logger import get_logger
from ska_matplotlib import plot_cxctime

from .paths import FID_STATS_PATH

logger = get_logger()


def get_opt():
    parser = argparse.ArgumentParser(
        description="Commanded vs. telemetry fid positions"
    )
    parser.add_argument("--start", type=str, help="Start date (default=NOW - 90 days)")
    parser.add_argument("--stop", type=str, help="Stop date (default=NOW)")
    parser.add_argument(
        "--data-dir", type=str, default=".", help="Fid drift data directory"
    )
    parser.add_argument("--out", type=str, help="Output plot file")
    return parser


def get_fids_commands(dwells):
    """
    Get fid command values corresponding to ``dwells`` (keyed by obsid)

    :returns: dict of tables keyed by obsid
    """

    logger.info("Getting fids from commands")
    fids_commands = {}

    for obsid, dwell in dwells.items():
        # Only accept the first dwell of science observations
        if obsid in fids_commands:
            logger.info("Skipping obsid {} already in fids_commands".format(obsid))
            continue

        try:
            starcats = observations.get_starcats(obsid=obsid)
        except ValueError:
            logger.info("Skipping obsid {} with no starcat".format(obsid))
            continue

        if len(starcats) == 0:
            logger.info("Skipping obsid {} with no starcat".format(obsid))
            continue

        starcat = starcats[0]
        fids = starcat[starcat["type"] == "FID"]
        fids["obsid"] = obsid
        if len(fids) > 0:
            fids_commands[obsid] = Table(fids)
        else:
            logger.info("No fids found for obsid {} in kadi starcats".format(obsid))

    return fids_commands


def get_dwells_with_fids(start, stop):
    """
    Get telemetry yag, zag values for each fid in ``fids_commands`` at the start
    of the corresponding ``dwells``.

    :returns: dict of tables keyed by obsid (like commands fids)
    """
    # Only get dwells that have an obsid
    stop = min(CxoTime(stop).date, events.obsids.all().reverse()[0].stop)

    logger.info("Getting dwells betweeen {} and {}".format(start, stop))
    dwells = OrderedDict()
    for dwell in events.dwells.filter(start, stop):
        obsid = dwell.get_obsid()
        if obsid > 38000:
            continue
        if obsid in dwells:
            logger.info("Skipping duplicate obsid {} for dwell {}".format(obsid, dwell))
            continue
        dwells[obsid] = dwell

    return dwells


def get_fids_telem(dwells, fids_commands, dbh):
    """
    Get telemetry yag, zag values for each fid in ``fids_commands`` at the start
    of the corresponding ``dwells``.

    :returns: dict of tables keyed by obsid (like commands fids)
    """

    fids_telem = {}

    for obsid, dwell in dwells.items():
        logger.debug("Get_fids_telem for obsid {}".format(obsid))
        if obsid not in fids_commands:
            logger.info("Skipping obsid {} not in fids_commands".format(obsid))
            continue

        rows = dbh.fetchall(
            f"""SELECT tstart, obsid, slot, ang_y_med, ang_z_med from fid_stats
                             where obsid = {obsid}"""
        )

        if len(rows) > 0:
            fids_telem[obsid] = Table(
                rows, names=["tstart", "obsid", "slot", "aoacyan", "aoaczan"]
            )
        else:
            logger.info("No fids found for obsid {}".format(obsid))

    return fids_telem


def join_commands_telem(fids_commands, fids_telem):
    """
    Remake dict of tables into a single table for each structure
    """
    # Stack the dict of tables into a single table
    t_fids_commands = table.vstack(
        [fids_commands[obsid] for obsid in sorted(fids_commands)],
        metadata_conflicts="silent",
    )
    t_fids_telem = table.vstack(
        [fids_telem[obsid] for obsid in sorted(fids_telem)], metadata_conflicts="silent"
    )

    # Join on obsid and slot columns into a single table
    commands_telem = table.join(
        t_fids_commands,
        t_fids_telem,
        keys=["obsid", "slot"],
        metadata_conflicts="silent",
    )

    # Reject unacquired fids
    ok = commands_telem["aoacyan"] > -3276
    return commands_telem[ok]


def plot_commands_telem(commands_telem, savefig=None):
    plt.close(1)
    plt.figure(1, figsize=(6, 4))
    tstart = commands_telem["tstart"]
    dyag = commands_telem["aoacyan"] - commands_telem["yang"]
    dzag = commands_telem["aoaczan"] - commands_telem["zang"]
    plot_cxctime(tstart, dyag, ".r", label="Yag")
    plot_cxctime(tstart, dzag, ".b", label="Zag")
    plt.ylim(-45, 45)
    plt.grid()
    plt.legend(fontsize="small", numpoints=1, loc="lower left")
    plt.ylabel("Offset (arcsec)")
    plt.title("Fid light commanded and observed angles")
    x0, x1 = plt.xlim()
    dx = (x1 - x0) * 0.05
    x0, x1 = x0 - dx, x1 + dx
    plt.xlim(x0, x1)
    plt.hlines([-35, 35], x0, x1, colors="g", linestyles="--")
    plt.hlines([-40, 40], x0, x1, colors="r", linestyles="--")
    if savefig is not None:
        logger.info("Writing {}".format(savefig))
        plt.savefig(savefig)


def main(sys_args=None):
    opt = get_opt().parse_args(sys_args)
    matplotlib.use("Agg")

    start = CxoTime(opt.start) - 90 * u.day if opt.start is None else CxoTime(opt.start)
    stop = CxoTime(opt.stop)

    dwells = get_dwells_with_fids(start.date, stop.date)
    fids_commands = get_fids_commands(dwells)
    with ska_dbi.DBI(dbi="sqlite", server=FID_STATS_PATH(opt.data_dir)) as dbh:
        fids_telem = get_fids_telem(dwells, fids_commands, dbh)
    commands_telem = join_commands_telem(fids_commands, fids_telem)
    plot_commands_telem(commands_telem, savefig=opt.out)


if __name__ == "__main__":
    main()
