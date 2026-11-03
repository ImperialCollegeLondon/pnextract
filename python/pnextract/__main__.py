import os
import sys
from pathlib import Path

PACKAGE_DIR = Path(__file__).parent
BINARY_SUFFIX = ".exe" if os.name == "nt" else ""

BINARY_DIR = PACKAGE_DIR
if not (BINARY_DIR / f"pnextract{BINARY_SUFFIX}").exists():
    BINARY_DIR = PACKAGE_DIR.parent.parent / "bin"


def pnextract() -> None:
    os.execv(BINARY_DIR / f"pnextract{BINARY_SUFFIX}", sys.argv)


def voxel_image_process() -> None:
    os.execv(BINARY_DIR / f"voxelImageProcess{BINARY_SUFFIX}", sys.argv)


if __name__ == "__main__":
    pnextract()
