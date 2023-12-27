# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Gather fid statistics and store in the fid_stats table in a sqlite3 database.
"""

from pathlib import Path

import argparse
import astropy.units as u
import numpy as np
import ska_dbi
from astropy.table import Table, vstack
from chandra_aca.transform import count_rate_to_mag
from cxotime import CxoTime
from kadi.commands import get_observations
from mica.archive import asp_l1
from ska_helpers.logging import basic_logger

from . import __version__

LOGGER = basic_logger(__name__)
LOGGER.setLevel("DEBUG")

TABLE_SQL = """
create table fid_stats (
    obsid int,
    slot int,
    id_num int,
    id_string varchar(10),
    tstart float,
    date_obs varchar(21),
    ang_y_med float,
    ang_z_med float,
    ang_y_start_med float,
    ang_z_start_med float,
    ang_y_5th float,
    ang_y_95th float,
    ang_z_5th float,
    ang_z_95th float,
    mag_med  float,
    mag_i_avg float,
    exp_time  float,
    sim_z_offset float,
) ;
create index fid_stats_obsid_index on fid_stats ( obsid ) ;
"""

# Define reasonable precision for fid stats numerical values
FORMAT_ROUND = dict(
    ang_y_med=2,
    ang_z_med=2,
    ang_y_start_med=2,
    ang_z_start_med=2,
    ang_y_95th=2,
    ang_z_95th=2,
    ang_y_5th=2,
    ang_z_5th=2,
    mag_med=3,
    mag_i_avg=3,
    exp_time=2,
    sim_z_offset=1,
    tstart=1,
)


def process_obs(dbh, obs):
    """
    Process the observation data and insert the fid statistics into the database.

    Parameters
    ----------
    dbh : database handle
        The database handle object.
    obs : dict
        The observation data.

    Returns
    -------
    None
    """
    obsid = obs["obsid"]
    if obsid > 38000:
        return

    # Check if this is in the database already. Unlike the perl version this
    # does not typically update existing entries.
    fids = dbh.fetchall(f"select * from fid_stats where obsid={obsid}")
    if len(fids) > 0:
        LOGGER.debug(f"Obsid {obsid} already in database, skipping")
        return

    LOGGER.info(f"***** PROCESSING OBSID {obsid} at {obs['obs_start']} *****")

    # Read FIDPROPS and ACACENT data tables. This vstacks the files as needed.
    # If anything went wrong an error message is logged and we carry on.
    fidprops = get_archive_file_data(obs, "FIDPROPS", max_files=1)
    if fidprops is None:
        return
    acen = get_archive_file_data(obs, "ACACENT")
    if acen is None:
        return

    LOGGER.info(f"Got {len(fidprops)} rows fidprops and {len(acen)} rows acacen")

    # Calculate fid statistics, returning a list of dict for each fid
    stats = calc_fid_stats(obs, fidprops, acen)

    # Write the fid statistics to the database
    for stat in stats:
        LOGGER.debug(f"Inserting row {stat}")
        dbh.insert(stat, "fid_stats")


def get_archive_file_data(obs, content, max_files=None):
    """
    Get archive file data for a given observation.

    This is used for FIDPROPS and ACACENT content type files.

    Parameters
    ----------
    obs : dict
        The observation details.
    content : str
        The mica asp_l1 CONTENT of the files to retrieve.
    max_files : int, optional
        The maximum number of files to retrieve.

    Returns
    -------
    out : Table
        The combined data from the retrieved files.
    """
    files = asp_l1.get_files(
        start=obs["obs_start"], stop=obs["obs_stop"], content=[content]
    )

    if len(files) == 0:
        LOGGER.info(
            f"No {content} files found for obsid {obs['obsid']} between "
            f"{obs['obs_start']} and {obs['obs_stop']}"
        )
        return None

    if max_files is not None and len(files) > max_files:
        LOGGER.error(
            f"Too many {content} files found for obsid {obs['obsid']} between "
            f"{obs['obs_start']} and {obs['obs_stop']}"
        )
        return None

    dats = [Table.read(file) for file in files]
    LOGGER.debug(f"Got {len(dats)} {content} files")
    if len(dats) == 1:
        out = dats[0]
    else:
        out = vstack(dats, join_type="exact", metadata_conflicts="silent")

    return out


def delete_obs(dbh, obsid):
    """
    Delete the fid_stats entries for a given obsid.

    Parameters
    ----------
    dbh : object
        The database connection object.
    obsid : int
        The obsid for which fid_stats entries need to be deleted.

    Returns
    -------
    None

    """
    LOGGER.info(f"Deleting fid_stats entries for obsid {obsid}")
    dbh.execute(f"delete from fid_stats where obsid={obsid}")


def calc_fid_stats(obs, fidprops, acen):
    """
    Calculate fid statistics.

    Parameters
    ----------
    obs : dict
        The observation dictionary.
    fidprops : astropy Table
        The FIDPROPS data.
    acen : astropy Table
        aspect L1 centroids

    Returns
    -------
    stats: list
        The list of dicts of fid statistics.

    """
    stats = []
    for fidpr in fidprops:
        if fidpr["id_status"].strip() != "GOOD":
            LOGGER.info(
                f"WARNING: obsid {obs['obsid']} has bad fid status "
                f"{fidpr['id_status']!r} in slot {fidpr['slot']}"
            )
            continue
        stat = calc_stats_for_fidpr(obs, acen, fidpr)
        if stat is not None:
            stats.append(stat)

    return stats


def calc_stats_for_fidpr(obs, acen, fidpr):
    """
    Calculate statistics for an individual fid light

    Parameters
    ----------
    obs : dict
        Dictionary containing observation information.
    acen : astropy Table
        astropy Table of aspect L1 centroid data
    fidpr : astropy Row
        astropy Row of fidprops data

    Returns
    -------
    stat: dict
        Dictionary containing calculated statistics for single fid light in an observation

    """
    t_int = 1.696  # Assume 1.696 sec integration time
    slot = fidpr["slot"]
    ok = (acen["slot"] == slot) & (acen["alg"] == 8) & (acen["status"] == 0)
    cen = acen[ok]
    n_cen = len(cen)

    if n_cen < 200:
        LOGGER.info(
            f"WARNING: obsid {obs['obsid']} has {n_cen} < 200 "
            f"fid counts in slot {slot}, skipping"
        )
        return

    ang_y_sm = cen["ang_y_sm"] * 3600
    ang_z_sm = cen["ang_z_sm"] * 3600
    # SIM_Z offset in arcsec (using 20.493 arcsec/mm)
    sim_z_offset = (fidpr.meta["LSI0STT3"] + fidpr.meta["STT0STF3"]) * 20.493
    stat = {
        "slot": slot,
        "id_num": fidpr["id_num"],
        "id_string": fidpr["id_string"],
        "tstart": CxoTime(obs["obs_start"]).secs,
        "obsid": obs["obsid"],
        "date_obs": obs["obs_start"],
        "ang_y_med": np.median(ang_y_sm[-1000:]),
        "ang_z_med": np.median(ang_z_sm[-1000:]),
        "ang_y_start_med": np.median(ang_y_sm[:1000]),
        "ang_z_start_med": np.median(ang_z_sm[:1000]),
        "ang_y_5th": np.percentile(ang_y_sm, 5),
        "ang_y_95th": np.percentile(ang_y_sm, 95),
        "ang_z_5th": np.percentile(ang_z_sm, 5),
        "ang_z_95th": np.percentile(ang_z_sm, 95),
        "mag_med": count_rate_to_mag(np.median(cen["counts"]) / t_int),
        "mag_i_avg": fidpr["mag_i_avg"],
        "exp_time": acen.meta["TSTOP"] - acen.meta["TSTART"],
        "sim_z_offset": sim_z_offset,
    }

    # Make the numbers prettier
    for key, precision in FORMAT_ROUND.items():
        stat[key] = round(stat[key], precision)

    return stat


def get_options(sys_args=None):
    parser = argparse.ArgumentParser(description=f"get_fid_data {__version__}")
    parser.add_argument(
        "--data-dir", type=str, default=".", help="Data directory (default='.')"
    )
    parser.add_argument(
        "--stop", type=str, help="Stop date for processing (default=NOW)"
    )
    parser.add_argument("--obsid", type=int, help="Obsid to process manually")
    parser.add_argument(
        "--lookback",
        type=float,
        default=30,
        help="Number of days to look back before --stop (default=30)",
    )
    parser.add_argument(
        "--delete",
        type=int,
        action="append",
        default=[],
        help="Obsid to delete (can be supplied multiple times)",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        help="logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL, default=INFO)",
    )
    return parser


def main(sys_args=None):
    opt = get_options().parse_args(sys_args)
    LOGGER.setLevel(opt.log_level)

    data_dir = Path(opt.data_dir)
    db3_path = data_dir / "fid_stats.db3"

    if not db3_path.exists():
        with ska_dbi.DBI(dbi="sqlite", server=db3_path) as dbh:
            dbh.execute(TABLE_SQL)

    if opt.obsid:
        obss = get_observations(obsid=opt.obsid)
    else:
        stop = CxoTime(opt.stop)
        start = stop - opt.lookback * u.day
        obss = get_observations(start, stop)

    with ska_dbi.DBI(dbi="sqlite", server=db3_path) as dbh:
        # Delete the specified (comma-sep) list of obsids first.  Mostly for testing.
        for obsid in opt.delete:
            delete_obs(dbh, obsid)

        for obs in obss:
            try:
                process_obs(dbh, obs)
            except Exception:
                # Some unhandled exception occurred. Log traceback and continue.
                import traceback

                msg = traceback.format_exc()
                LOGGER.error(f"ERROR\n{msg}")
                # Delete any stats for this obsid that have been inserted
                delete_obs(dbh, obs["obsid"])


if __name__ == "__main__":
    main()
