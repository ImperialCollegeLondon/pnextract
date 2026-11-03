import os
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import Distribution, setup
from setuptools.command.bdist_wheel import bdist_wheel
from setuptools.command.build_py import build_py

BINARIES = (
    "pnextract",
    "voxelImageProcess",
)
WINDOWS_WHEEL = sys.platform == "win32"
BINARY_SUFFIX = ".exe" if WINDOWS_WHEEL else ""

if sys.platform == "darwin":
    os.environ.setdefault("MACOSX_DEPLOYMENT_TARGET", "11.0")


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


class BuildPy(build_py):
    def run(self):
        binaries = [Path("bin") / f"{name}{BINARY_SUFFIX}" for name in BINARIES]
        command = (
            ["cmd.exe", "/d", "/c", "msys2", "-c", "make mgw"]
            if WINDOWS_WHEEL
            else ["make"]
        )
        subprocess.run(command, check=True)

        super().run()

        package_dir = Path(self.build_lib) / "pnextract"
        package_dir.mkdir(parents=True, exist_ok=True)

        for binary in binaries:
            shutil.copy2(binary, package_dir / binary.name)


class BdistWheel(bdist_wheel):
    def get_tag(self):
        _, _, platform = super().get_tag()
        # Override Python version and ABI compatibility since this project does not link against Python
        return self.python_tag, "none", platform


setup(
    distclass=BinaryDistribution,
    cmdclass={
        "build_py": BuildPy,
        "bdist_wheel": BdistWheel,
    },
)
