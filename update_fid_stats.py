# ##########################################################################################
# #
# # Gather fid statistics and store in the fid_stats table in the sybase aca_ops database
# #
# ##########################################################################################

from pathlib import Path

import numpy as np
import Ska.DBI
from astropy.table import Table, vstack
from chandra_aca.transform import count_rate_to_mag
from cxotime import CxoTime
from kadi.commands import get_observations
from mica.archive import asp_l1
from ska_helpers.utils import basic_logger

# use File::Basename;
# use Getopt::Long;
# use PDL;
# use PDL::NiceSlice;
# # use CXC::Envs;
# # use Chandra::Tools::Common;
# use Ska::Process qw(get_archive_files);
# use Ska::IO qw(read_param_file cfitsio_read_header);
# use Ska::HashTable;
# use CFITSIO::Simple qw(fits_read_bintbl fits_read_hdr);
# use Data::Dumper;
# use Ska::DatabaseUtil qw(:all);
# use File::Temp qw/tempdir tempfile/;

# use App::Env;
# my $ftools = App::Env->new('HEADAS');        # Set up HEADAS (aka FTOOLS) environment

# # Set up some constants

# my $TASK = "fid_drift_mon";
# my $SKA = $ENV{SKA} || die "Must set SKA env var\n";
# my $TASKDATA = "$SKA/data/$TASK";
# our @row;
# our @loglines;                  # All logged messages (per obsid) get stored here

# my $tmp_data_dir = tempdir(CLEANUP => 1); #
# chdir $tmp_data_dir or die "Could not change to tmp dir '$tmp_data_dir'\n";
# # Use END block to always change out of that temp directory at the end of the program
# # This should help with cleanup of the directory
# END {
#     chdir("/");
# }

# # Mostly for reference, here is the table creation command
# my $create_table_sql = <<CREATE_TABLE_SQL;
#   create table fid_stats (obsid int,
#                           id_num int,
#                           id_string varchar(10),
#                           tstart float,
#                           date_obs varchar(20),
#                           ang_y_med float,
#                           ang_z_med float,
#                           ang_y_5th float,
#                           ang_y_95th float,
#                           ang_z_5th float,
#                           ang_z_95th float,
#                           mag_med  float,
#                           mag_i_avg float,
#                           exp_time  float,
#                           sim_z_offset float,
#                           ap_date datetime,
#                           proc_status varchar(40)
#                           ) ;
#   create index fid_stats_obsid_index on fid_stats ( obsid ) ;
# CREATE_TABLE_SQL

__version__ = "0.1"

LOGGER = basic_logger(__name__)
LOGGER.setLevel("DEBUG")

TABLE_SQL = """
create table fid_stats (
    obsid int,
    id_num int,
    id_string varchar(10),
    tstart float,
    date_obs varchar(21),
    ang_y_med float,
    ang_z_med float,
    ang_y_5th float,
    ang_y_95th float,
    ang_z_5th float,
    ang_z_95th float,
    mag_med  float,
    mag_i_avg float,
    exp_time  float,
    sim_z_offset float,
    ap_date datetime,
    proc_status varchar(40)
) ;
create index fid_stats_obsid_index on fid_stats ( obsid ) ;
"""

# @fidprops_cols = qw(slot mag_i_avg id_string id_num id_status);
# %formats = (detector => '%s',
# 	    ang_y_med => '%.2f',
# 	    ang_z_med => '%.2f',
# 	    ang_y_95th => '%.2f',
# 	    ang_z_95th => '%.2f',
# 	    ang_y_5th => '%.2f',
# 	    ang_z_5th => '%.2f',
# 	    mag_med => '%.3f',
# 	    mag_i_avg => '%.3f',
# 	    exp_time => '%.2f',
# 	    sim_z_offset => '%.1f',
# 	    tstart => '%.1f',
#             obsid => '%d',
# 	    );

FIDPROPS_COLS = ("slot", "mag_i_avg", "id_string", "id_num", "id_status")
FORMAT_ROUND = dict(
    ang_y_med=2,
    ang_z_med=2,
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


class NoArchiveFilesFound(FileNotFoundError):
    pass


def main():
    # our $fid_table = 'fid_stats';
    # get_options();			# Parse command line options
    # my $dbh = sql_connect('sybase-aca-aca_ops');
    opt = get_options()

    data_dir = Path(opt.data_dir)
    db3_path = data_dir / "fid_stats.db3"

    if not db3_path.exists():
        with Ska.DBI.DBI(dbi="sqlite", server=db3_path) as dbh:
            dbh.execute(TABLE_SQL)

    # my @obs = get_obs($dbh);
    # my @obsids = map { $_->{obsid} } @obs;
    # # Extract the list of observations
    # log_message("Processing obsids @obsids");
    obss = get_observations(opt.start, opt.stop)

    with Ska.DBI.DBI(dbi="sqlite", server=db3_path) as dbh:
        # # Delete the specified (comma-sep) list of obsids first.  Mostly for testing.
        # delete_obs($dbh) if $par{delete};
        if opt.delete:
            delete_obs(dbh, opt.delete)

        for obs in obss:
            process_obs(dbh, obs)


def process_obs(dbh, obs):
    # my $obsid = $obs->{obsid};
    obsid = obs["obsid"]
    if obsid > 38000:
        return

    # foreach $obs (@obs) {
    #     clean_up($obsid);           # Make completely sure no old files are hanging around
    #     @loglines = ();             # Clear logged messages for this obsid
    #     my $error;
    #     $obsid = $obs->{obsid};
    #     $date = localtime;
    #     log_message("********* PROCESSING OBSID $obsid at $date ***********\n");

    #     my @rows = sql_fetchall_array_of_hashref($dbh,
    #                                              qq{select obsid from $fid_table
    #                                                 where obsid=$obsid});
    #     # Delete existing values if needed
    #     if (@rows) {
    #         my $cmd = qq{delete from $fid_table where obsid=$obsid};
    #         $par{readonly} ? print("$cmd\n") : sql_do($dbh, $cmd);
    #     }

    # Check if this is in the database already. Unlike the perl version this
    # does not typically update existing entries.
    fids = dbh.fetchall(f"select * from fid_stats where obsid={obsid}")
    if len(fids) > 0:
        LOGGER.debug(f"Obsid {obsid} already in database, skipping")
        return

    LOGGER.info(f"***** PROCESSING OBSID {obsid} at {obs['obs_start']} *****")

    #     # Do some processing for each obsid.  If any step along the way fails,
    #     # bail out to the clean-up step
    fidprops = get_archive_file_data(obs, "FIDPROPS", max_files=1)
    if fidprops is None:
        return

    acen = get_archive_file_data(obs, "ACACENT")
    if acen is None:
        return

    LOGGER.info(f"Got {len(fidprops)} rows fidprops and {len(acen)} rows acacen")

    stats = calc_fid_stats(obs, fidprops, acen)

    for stat in stats:
        LOGGER.debug(f"Inserting row {stat}")
        dbh.insert(stat, "fid_stats")


def get_archive_file_data(obs, content, max_files=None):
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

    #     unless (merge_files()) {
    #         $error = "Could not merge acacen files";
    # 	log_message($error);
    # 	goto CLEAN_UP;
    #     }

    #     unless (@fidprops = read_fidprops()) {
    #         $error = "Could not read fidprops";
    # 	log_message($error);
    # 	goto CLEAN_UP;
    #     }

    #     unless (@fid_stats = calc_stats()) {
    #         $error = "Could not calculate statistics files";
    # 	log_message($error);
    # 	goto CLEAN_UP;
    #     }

    #     log_message("Adding fid stats to database");
    #     if (not $par{readonly}) {
    #         for my $fid_stat (@fid_stats) {
    #             $fid_stat->{ap_date} = $obs->{last_ap_date};

    #             eval { sql_insert_hash($dbh, $fid_table, $fid_stat) };
    #             if ($@) {
    #                 log_message($@);
    #                 $error = substr($@, 0, 39);
    #                 goto CLEAN_UP;
    #             }
    #         }
    #     }

    #   CLEAN_UP:
    #     # In case of any error, remove database rows (which may be corrupted)
    #     # and put in an entry with the error message.
    #     if (defined $error) {
    #         my $cmd = qq{delete from $fid_table where obsid=$obsid};
    #         if ($par{readonly}) {
    #             print "$cmd\n";
    #             print "Insert process status error message\n";
    #         } else {
    #             sql_do($dbh, $cmd);
    #             sql_insert_hash($dbh, $fid_table, {obsid => $obsid,
    #                                                proc_status => $error} );
    #         }
    #     }

    #     clean_up($obsid, $error);
    # }

    # # Close down arc5gl properly (somehow normal garbage collection doesn't work here)
    # if (defined $Ska::Process::arc5gl) {
    #     undef $Ska::Process::arc5gl;
    #     sleep 2;
    # }


# ##****************************************************************************
# sub get_options {
# ##****************************************************************************
#     $par{loud}     = 1;		# Run loudly

#     GetOptions( \%par,
# 	       'help!',
# 	       'loud!',
# 	       'tstart=s',
# 	       'tstop=s',
#                'delete=s',
#                'readonly!',
# 	       ) ||
# 		   exit( 1 );

#     if ($par{loud}) {
# 	log_message("COMMAND LINE PARAMETERS");
# 	foreach (sort keys %par) {
# 	    log_message(sprintf("  %-16s = %s", $_, $par{$_}));
# 	}
# 	log_message("");
#     }

# }


def get_options(sys_args=None):
    import argparse

    # RESTORE THIS LINE
    # from . import __version__

    parser = argparse.ArgumentParser(description=f"get_fid_data {__version__}")
    parser.add_argument("--start", type=str, help="Start date")
    parser.add_argument("--stop", type=str, help="Stop date")
    parser.add_argument(
        "--data-dir", type=str, default=".", help="Data directory (default='.')"
    )
    parser.add_argument(
        "--delete", type=str, help="Comma-separated list of obsids to delete"
    )
    parser.add_argument(
        "--log-level",
        type=str,
        help="logging level (debug, info, warning, error, critical)",
    )
    args = parser.parse_args(sys_args)
    return args


# ##****************************************************************************
# sub delete_obs {
# ##****************************************************************************
#     my $dbh = shift;
#     for my $obsid (split(',', $par{delete})) {
#         next unless $obsid =~ /^ \s* \d+ \s* $/x;
#         log_message("Deleting $fid_table entries for obsid $obsid");
#         if (not $par{readonly}) {
#             sql_do($dbh, "delete from $fid_table where obsid=$obsid");
#         }
#     }
# }


def delete_obs(dbh, delete_obsids):
    pass


# ##****************************************************************************
# sub get_obs {
# ##****************************************************************************
#     my $dbh = shift;


# ##****************************************************************************
# sub get_fidprops_acacen_files {
# ##****************************************************************************

# # Get event file and source, either from current directory or from archive
#     log_message("Getting acacen file(s)");
#     @acen = get_archive_files(obsid 	 => $obsid,
# 			      prod  	 => "asp1{acacent}",
# 			      file_glob  => "pcad*acen1.fits*",
# 			      version    => [4,3,2,1],
#                               dir       => $tmp_data_dir,
#                               loud      => 0,
#                               guestuser => 1,
#         );
#     return 0 unless (@acen);

#     log_message("Getting fidprops file(s)");
#     @fidpr = get_archive_files(obsid 	 => $obsid,
# 				prod  	 => "asp1{fidprops}",
# 			      file_glob => "pcad*fidpr1.fits*",
# 			      version    => [4,3,2,1],
#                               dir       => $tmp_data_dir,
#                               loud      => 0,
#                               guestuser => 1,
#         );
#     return 0 unless (@fidpr);
# }

# ##***************************************************************************
# sub merge_files {
# ##***************************************************************************
# # If necessary, merge centroid files
#     $" = "\n";			# Set the list separator to "\n"
#     foreach (qw(acen fidpr)) {
# 	if (@{$_} == 1) {
# 	    $cat = ($ {$_}[0] =~ /\.gz$/) ? "gunzip --stdout" : "cat";
# 	    $cmd = "$cat $ {$_}[0] > obs_${_}1.fits";
# 	    log_message($cmd);
# 	    system($cmd);
#         } else {
# 	    open (MERGE, "> merge.lis") || die "ERROR - Couldn't open merge.lis for writing\n";
# 	    print MERGE "@{$_}\n";
# 	    close MERGE;
# 	    $cmd = "fmerge \@merge.lis obs_${_}1.fits - clobber=yes";
# 	    log_message($cmd);
# 	    $ftools->system($cmd);
# 	    unlink "merge.lis";
# 	}
#     }
#     $" = ' ';			# Reset the list separator to " "
#     return (-e "obs_acen1.fits" && -e "obs_fidpr1.fits");
# }

# ##***************************************************************************
# sub read_fidprops {
# ##***************************************************************************
#     my @entries;
#     my $i;
#     my $j;
#     my @vals;
#     my %tmp_vals;

#     @vals = CFITSIO::Simple::read_bintbl_cols("obs_fidpr1.fits", @fidprops_cols, { extname => 'FIDPROPS' });

#     # Convert PDL elements into an array of hashes, for easier manipulation

#     @entries = ();
#     for $j (0 .. nelem($vals[0])-1) {
# 	%tmp_vals = ();
# 	for $i (0 .. $#fidprops_cols) {
# 	    $tmp_vals{$fidprops_cols[$i]} = (UNIVERSAL::isa($vals[$i],'PDL')) ?
# 		$vals[$i]->at($j) : $vals[$i]->[$j];
# 	}
# 	push @entries, { %tmp_vals };
#     }
#     return @entries;
# }

# ##***************************************************************************
# sub clean_up {
# ##***************************************************************************
#     my $obsid = shift;
#     my $error = shift;
#     my @files;

#     @files = glob("pcad*fidpr*.fits* pcad*acen*fits* obs_fidpr1.fits* obs_acen1.fits*");
#     if ($error) {
#         my $dir = "$TASKDATA/failed/obs$obsid";
#         log_message("Copying files for obsid to $dir");
#         if (not -e $dir) {
#             mkdir($dir) or die "Could not make dir $dir";
#         }
#         open(LOG, "> $dir/log") or die "Could not open $dir/log";
#         print LOG join('', @loglines);
#         close(LOG);
#         if (@files) {
#             system("cp", "-p", @files, "$dir/") and die "Could not copy files to $dir";
#             system("gzip -f $dir/*.fits");
#         }
#     }

#     unlink @files if (@files);
# }

# ##***************************************************************************
# sub calc_stats {
# ##***************************************************************************
#     my $M0 = 10.32;                 # ACA magnitude corresponding to C0 electrons/second
#     my $C0 = 5263.0;                     # Electrons/second
#     my $t_int = 1.696;                   # Assume fixed readout time of 1.696 seconds

#     my $acen = Ska::HashTable->new('obs_acen1.fits');
#     my $acen_hdr = fits_read_hdr('obs_acen1.fits');
#     my $fidpr_hdr = fits_read_hdr("obs_fidpr1.fits[FIDPROPS]");
#     my @fid_stats;
#     log_message("Calculating fid statistics");

#   FID: foreach $fidpr (@fidprops) {
#         my $data;
#         if ($fidpr->{id_status} ne 'GOOD') {
#             log_message("WARNING: Obsid $obsid has bad fid in slot $fidpr->{slot}") ;
#             next FID;
#         }
#         my $slot = $fidpr->{slot};
#         my $cen = $acen->row("slot == $slot && alg == 8 && status == 0");
#         my $ncen = $cen->rows();
#         if ($ncen > 200) {
#             my $ind_y = qsorti($cen->ang_y_sm);
#             my $ind_z = qsorti($cen->ang_z_sm);
#             my $cts_med = median($cen->counts);
#             my $mag_med = $M0 - 2.5 * log10($cts_med/$t_int / $C0);
#             my $ok0 = ($ncen > 1500) ? $ncen - 1000 : 0;
#             my $ang_y = $cen->ang_y_sm->($ind_y) * 3600;
#             my $ang_z = $cen->ang_z_sm->($ind_z) * 3600;

#             $data = { id_num => $fidpr->{id_num},
#                       id_string => $fidpr->{id_string},
#                       tstart => $cen->time->at(0),
#                       obsid => $acen_hdr->{OBS_ID},
#                       date_obs => $acen_hdr->{'DATE-OBS'},
#                       ang_y_med => median($cen->ang_y_sm->($ok0:)) * 3600,
#                       ang_z_med => median($cen->ang_z_sm->($ok0:)) * 3600,
#                       ang_y_5th => $ang_y->at(int($ncen*0.05)),
#                       ang_y_95th => $ang_y->at(int($ncen*0.95)),
#                       ang_z_5th => $ang_z->at(int($ncen*0.05)),
#                       ang_z_95th => $ang_z->at(int($ncen*0.95)),
#                       mag_med  => $mag_med,
#                       mag_i_avg => $fidpr->{mag_i_avg},
#                       exp_time => ($acen_hdr->{TSTOP} - $acen_hdr->{TSTART}) / 1000.,
#                       sim_z_offset => ($fidpr_hdr->{LSI0STT3} + $fidpr_hdr->{STT0STF3}) * 20.493, # arcsec/mm
#                     };
#             foreach (keys %{$data}) {
# 		$data->{$_} = sprintf($formats{$_} || "%s", $data->{$_});
# 	    }
#             push @fid_stats, $data;
#         }
#     }

#     return @fid_stats;
# }


def calc_fid_stats(obs, fidprops, acen):
    """
    Calculate fid statistics
    """
    stats = []
    for fidpr in fidprops:
        if fidpr["id_status"].strip() != "GOOD":
            LOGGER.info(
                f"WARNING: obsid {obs['obsid']} has bad fid status {fidpr['id_status']!r} in slot {fidpr['slot']}"
            )
            continue
        stat = calc_stats_for_fidpr(obs, acen, fidpr)
        if stat is not None:
            stats.append(stat)

    return stats


def calc_stats_for_fidpr(obs, acen, fidpr):
    t_int = 1.696
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
        "id_num": fidpr["id_num"],
        "id_string": fidpr["id_string"],
        "tstart": CxoTime(obs["obs_start"]).secs,
        "obsid": obs["obsid"],
        "date_obs": obs["obs_start"],
        # CHANGE ME TO :1000 (use first 4 ksec not last 4 ksec)
        "ang_y_med": np.median(ang_y_sm[-1000:]),
        "ang_z_med": np.median(ang_z_sm[-1000:]),
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


#            $data = { id_num => $fidpr->{id_num},
#                       id_string => $fidpr->{id_string},
#                       tstart => $cen->time->at(0),
#                       obsid => $acen_hdr->{OBS_ID},
#                       date_obs => $acen_hdr->{'DATE-OBS'},
#                       ang_y_med => median($cen->ang_y_sm->($ok0:)) * 3600,
#                       ang_z_med => median($cen->ang_z_sm->($ok0:)) * 3600,
#                       ang_y_5th => $ang_y->at(int($ncen*0.05)),
#                       ang_y_95th => $ang_y->at(int($ncen*0.95)),
#                       ang_z_5th => $ang_z->at(int($ncen*0.05)),
#                       ang_z_95th => $ang_z->at(int($ncen*0.95)),
#                       mag_med  => $mag_med,
#                       mag_i_avg => $fidpr->{mag_i_avg},
#                       exp_time => ($acen_hdr->{TSTOP} - $acen_hdr->{TSTART}) / 1000.,
#                       sim_z_offset => ($fidpr_hdr->{LSI0STT3} + $fidpr_hdr->{STT0STF3}) * 20.493, # arcsec/mm


if __name__ == "__main__":
    main()
