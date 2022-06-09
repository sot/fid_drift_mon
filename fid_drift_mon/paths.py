# Licensed under a 3-clause BSD style license - see LICENSE.rst
from pathlib import Path


def FID_STATS_PATH(data_dir):
    return Path(data_dir) / 'fid_stats.db3'
