import os
from pathlib import Path
import shutil
import argparse

FILEDIR = os.path.dirname(__file__)


def parse_args(args):
    parser = argparse.ArgumentParser(description="Copy index.html")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory")
    args = parser.parse_args(args)
    return args


def main(args=None):
    args = parse_args(args)
    if not Path(args.outdir).exists():
        Path(args.outdir).mkdir(parents=True)
    shutil.copy(
        Path(FILEDIR) / "data" / "index_template.html", Path(args.outdir) / "index.html"
    )


if __name__ == "__main__":
    main()
