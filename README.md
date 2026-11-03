# pnextract

[![Package](https://github.com/lmmp-puc-rio/pnextract/actions/workflows/package.yml/badge.svg)](https://github.com/lmmp-puc-rio/pnextract/actions/workflows/package.yml)
[![PyPI](https://img.shields.io/pypi/v/pnextract)](https://pypi.org/project/pnextract/)


**pnextract** is a C++ application for extracting pore networks from images and statistical reconstructions of porous materials. The extracted networks preserve the topology of the pore space and provide geometric information used in pore-network modelling and related analyses.

The project originates from the **Imperial College London pore-scale modelling group**. This repository is a fork of [`ImperialCollegeLondon/pnextract`](https://github.com/ImperialCollegeLondon/pnextract).

## About this fork

This fork is maintained by [**LMMP/PUC-Rio**](https://github.com/LMMP-PUC-Rio) to provide packaging, automated builds, and binary distribution of **pnextract**. Distribution through this fork is done with permission.

Where appropriate, changes that are generally useful to **pnextract** rather than specific to the **LMMP/PUC-Rio** distribution should be proposed back to the [upstream project](https://github.com/ImperialCollegeLondon/pnextract).

## Installation

### PyPI

Despite not strictly being a Python application, **pnextract** is distributed as a package on [PyPI](https://pypi.org/project/pnextract/)

Thus, it can be installed with `pip`:

```bash
pip install pnextract
```

or with [**uv**](https://github.com/astral-sh/uv):

```bash
uv tool install pnextract
```

This installs the command-line tools:

```text
pnextract
voxelImageProcess
```

### Standalone binaries

Pre-built binaries produced by **LMMP/PUC-Rio** are attached to releases of this repository. Currently these target Linux only.
