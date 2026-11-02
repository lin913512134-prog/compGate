# -*- coding: utf-8 -*-
"""Signature scoring conventions for composition-aware audits.

Two mean-Z conventions are provided because their choice can flip audit
conclusions (R12 methodological finding, OSF 10.17605/OSF.IO/CQDP8):

- ``frozen``: z-score each gene *within a sample* across the panel, then
  average (the R10/R11 registry convention). Removes per-sample panel level;
  damps composition-driven level shifts; suited to pattern/shape signatures.
- ``level``: z-score each gene *across samples*, then average. Preserves
  level effects; at signature training margins it recovers published
  discrimination magnitudes (d 0.6-1.7 for 18 TB signatures) and yields
  materially higher composition-attribution fractions.

The convention must be stated with every reported CAF.
"""
import pandas as pd


def mean_z(X, genes, convention="level"):
    """Mean-Z score over the measurable subset of `genes`.

    Parameters
    ----------
    X : DataFrame, genes x samples
    genes : iterable of gene symbols
    convention : {"frozen", "level"}

    Returns
    -------
    pd.Series over samples, or None if no gene is measurable.
    """
    g = [x for x in genes if x in X.index]
    if not g:
        return None
    sub = X.loc[g].astype(float)
    if convention == "frozen":
        z = (sub - sub.mean(axis=0)) / (sub.std(axis=0) + 1e-12)
        return z.mean(axis=0, skipna=True)
    if convention == "level":
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-12, axis=0)
        return z.mean(axis=0, skipna=True)
    raise ValueError(f"unknown convention: {convention!r}")
