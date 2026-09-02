# Figure Files

This directory contains the code, input data, and final artwork for the PRX manuscript.

## Directory layout

- `code/`: MATLAB scripts and other figure-generation code.
- `data/`: input data required by the corresponding plotting scripts.
- `Figs/`: exported manuscript figures in publication formats (`.png`, `.eps`, and/or `.pdf`).

Within `code/`, `data/`, and `Figs/`, folders are organized by manuscript figure:

- `Fig1`--`Fig5`: main-text figures and their associated resources.
- `Fig6`: final exported main-text figure artwork. Its code and data have not yet been reorganized into this directory structure.
- `AppendixA`--`AppendixD`: figures and resources for the corresponding appendices.

## Usage

When moving or archiving figures, keep each plotting script together with the data in
the same named subfolder under `code/` and `data/`. Exported files should be written to
the matching subfolder under `Figs/` so that the LaTeX manuscript can use stable relative
paths.

The manuscript references the finalized assets from `Figs/Figs/`.
