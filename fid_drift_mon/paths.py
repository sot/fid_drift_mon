# Licensed under a 3-clause BSD style license - see LICENSE.rst
from pathlib import Path


def FID_STATS_PATH(data_dir):
    return Path(data_dir) / "fid_stats.db3"


def INDEX_TEMPLATE_PATH():
    return Path(__file__).parent / "data" / "index_template.html"
