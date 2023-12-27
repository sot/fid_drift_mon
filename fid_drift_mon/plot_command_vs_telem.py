# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Plot positions of fid lights from commands and the actual values seen in telemetry.
::

  Usage: plot_commands_vs_telem.py [-h] [--start START] [--stop STOP]
                                    [--out OUT]

  Commanded vs. telemetry fid positions

  optional arguments:
    -h, --help     show this help message and exit
    --start START  Start date (default=NOW - 90 days)
    --stop STOP    Stop date (default=NOW)
    --data-dir DATA_DIR (default ".")
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
from ska_matplotlib import plot_cxctime
from ska_helpers.logging import basic_logger

from .paths import FID_STATS_PATH

LOGGER = basic_logger(__name__)
LOGGER.setLevel("DEBUG")


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
    Get fid commanded values from star catalogs corresponding to ``dwells`` (keyed by obsid).

    Parameters:
    -----------
    dwells : dict
        Dictionary of dwells keyed by obsid.

    Returns:
    --------
    dict
        Dictionary of astropy Tables keyed by obsid including slot, yang, zang etc for each fid.

    """
    LOGGER.info("Getting fids from commands")
    fids_commands = {}

    for obsid, dwell in dwells.items():
        # Only accept the first dwell of science observations
        if obsid in fids_commands:
            LOGGER.info(f"Skipping obsid {obsid} already in fids_commands")
            continue

        try:
            starcats = observations.get_starcats(obsid=obsid)
        except ValueError:
            LOGGER.info(f"Skipping obsid {obsid} with no starcat")
            continue

        if len(starcats) == 0:
            LOGGER.info(f"Skipping obsid {obsid} with no starcat")
            continue

        starcat = starcats[0]
        fids = starcat[starcat["type"] == "FID"]
        fids["obsid"] = obsid
        if len(fids) > 0:
            fids_commands[obsid] = Table(fids)
        else:
            LOGGER.info(f"No fids found for obsid {obsid} in kadi starcats")

    return fids_commands


def get_dwells_with_fids(start, stop):
    """
    Get a dictionary of kadi OR dwell events with fid lights from start to stop, keyed by obsid.

    Parameters
    ----------
    start : str
        The start time of the dwells.
    stop : str
        The stop time of the dwells.

    Returns
    -------
    dict
        A dictionary of kadi dwell events keyed by obsid.

    Notes:
    ------
    - Only the first science observation dwell with the commanded obsid is included.
      In the case of scs107 and multiple dwells labeled with the first obsid, only the first
      can have fid lights turned on anyway.
    """

    # Limit to the range that has processed obsids
    stop = min(CxoTime(stop).date, events.obsids.all().reverse()[0].stop)

    LOGGER.info(
        f"Getting dwells between {CxoTime(start).date} and {CxoTime(stop).date}"
    )
    dwells = OrderedDict()
    for dwell in events.dwells.filter(start, stop):
        obsid = dwell.get_obsid()
        if obsid > 38000:
            continue
        if obsid in dwells:
            LOGGER.info(f"Skipping duplicate obsid {obsid} for dwell {dwell}")
            continue
        dwells[obsid] = dwell

    return dwells


def get_fids_telem(dwells, fids_commands, dbh):
    """
    Retrieve the fid position yag, zag values from the fid stats database for each commanded
    fid in ``fids_commands`` at the start of the corresponding dwell in ``dwells``.

    Parameters
    ----------
    dwells : dict
        A dictionary of dwells keyed by obsid.
    fids_commands : dict
        A dictionary of the fid positions from star catalog commands fids keyed by obsid.
    dbh : object
        The database handler object.

    Returns
    -------
    dict
        A dictionary of astropy Tables keyed by obsid, containing telemetry aoacyan aoaczan for each fid.
    """
    fids_telem = {}

    for obsid, dwell in dwells.items():
        LOGGER.debug(f"Get_fids_telem for obsid {obsid}")
        if obsid not in fids_commands:
            LOGGER.info(f"Skipping obsid {obsid} not in fids_commands")
            continue

        rows = dbh.fetchall(
            f"""SELECT tstart, obsid, slot, ang_y_start_med, ang_z_start_med from fid_stats
                             where obsid = {obsid}"""
        )

        if len(rows) > 0:
            fids_telem[obsid] = Table(
                rows, names=["tstart", "obsid", "slot", "aoacyan", "aoaczan"]
            )
        else:
            LOGGER.info(f"No fids found for obsid {obsid}")

    return fids_telem


def join_commands_telem(fids_commands, fids_telem):
    """
    Join data from ``fids_commands`` and ``fids_telem`` into a single table.

    The ``fids_commands`` and ``fids_telem`` dictionaries are keyed by obsid and contain astropy Tables.
    The returned table contains the columns from both tables, joined on obsid and slot.

    Parameters
    ----------
    fids_commands : dict
        A dictionary of astropy Tables containing star catalog fid positions for each obsid.
    fids_telem : dict
        A dictionary of astropy Tables containing observed fid positions for each obsid.

    Returns
    -------
    commands_telem : table
        A single astropy Table containing joined command and telemetry data.
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
    """
    Plot the commanded and observed angles for the fid lights.

    Parameters
    ----------
    commands_telem : astropy Table
        A dTable containing fid light commanded and observed angles.
    savefig : str, optional
        The filepath to save the plot as an image file. Defaults to None.

    Returns
    -------
    None
        This function does not return anything.

    """
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
        LOGGER.info(f"Writing {savefig}")
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
