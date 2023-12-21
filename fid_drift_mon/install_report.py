#!/usr/bin/env python
import argparse
from pathlib import Path
from ska_helpers.logging import basic_logger
import shutil


# Constants and file path definitions
FILE_DIR = Path(__file__).parent

logger = basic_logger(__name__, level="INFO")


def INDEX_TEMPLATE_PATH():
    return FILE_DIR / "data" / "index_template.html"


def get_opt():
    parser = argparse.ArgumentParser(description="Install report and plots")
    parser.add_argument(
        "--data-dir", type=str, default=".", help="Fid drift data directory"
    )
    parser.add_argument("--web-dir", type=str, default=".", help="Output web directory")
    return parser


def main(sys_args=None):
    opt = get_opt().parse_args(sys_args)
    data_dir = Path(opt.data_dir)
    web_dir = Path(opt.web_dir)

    if not web_dir.exists():
        logger.info(f"Creating web directory {web_dir}")
        web_dir.mkdir(parents=True)

    index_template_html = INDEX_TEMPLATE_PATH().read_text()
    out_html = web_dir / "index.html"
    logger.info(f"Writing HTML to {out_html}")
    out_html.write_text(index_template_html)

    # Copy the pngs from data to web
    for png in data_dir.glob("*.png"):
        shutil.copy(png, web_dir / png.name)


if __name__ == "__main__":
    main()
