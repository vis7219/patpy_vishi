import os
from pathlib import Path

import pytest
import scanpy as sc
from anndata import AnnData

import patpy

CACHE_DIR = Path(os.environ["PATPY_TEST_DATASETDIR"]) if os.environ.get("PATPY_TEST_DATASETDIR") else None


def _datasetdir(tmp_path: Path) -> Path:
    """Return a persistent cache dir when available, otherwise fall back to tmp_path."""
    if CACHE_DIR is not None:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        return CACHE_DIR
    return tmp_path


@pytest.mark.dataset
def test_combat_preprocessed_shape(tmp_path):
    original = sc.settings.datasetdir
    sc.settings.datasetdir = _datasetdir(tmp_path)
    try:
        adata = patpy.dt.combat_preprocessed()

        assert isinstance(adata, AnnData)
        assert adata.n_obs == 783677
        assert adata.n_vars == 3000
        assert "PCs" in adata.varm
    finally:
        sc.settings.datasetdir = original


@pytest.mark.dataset
def test_hlca_preprocessed_shape(tmp_path):
    original = sc.settings.datasetdir
    sc.settings.datasetdir = _datasetdir(tmp_path)
    try:
        adata = patpy.dt.hlca_preprocessed()

        assert isinstance(adata, AnnData)
        assert adata.n_obs == 1687127
        assert adata.n_vars == 3000
        assert "PCs" in adata.varm
    finally:
        sc.settings.datasetdir = original


@pytest.mark.dataset
def test_onek1k_preprocessed_shape(tmp_path):
    original = sc.settings.datasetdir
    sc.settings.datasetdir = _datasetdir(tmp_path)
    try:
        adata = patpy.dt.onek1k_preprocessed()

        assert isinstance(adata, AnnData)
        assert adata.n_obs == 1248980
        assert adata.n_vars == 3000
        assert "PCs" in adata.varm
    finally:
        sc.settings.datasetdir = original


@pytest.mark.dataset
def test_stephenson_preprocessed_shape(tmp_path):
    original = sc.settings.datasetdir
    sc.settings.datasetdir = _datasetdir(tmp_path)
    try:
        adata = patpy.dt.stephenson_preprocessed()

        assert isinstance(adata, AnnData)
        assert adata.n_obs == 639482
        assert adata.n_vars == 3000
        assert "PCs" in adata.varm
    finally:
        sc.settings.datasetdir = original


@pytest.mark.dataset
def test_ticatlas_preprocessed_shape(tmp_path):
    original = sc.settings.datasetdir
    sc.settings.datasetdir = _datasetdir(tmp_path)
    try:
        adata = patpy.dt.ticatlas_preprocessed()

        assert isinstance(adata, AnnData)
        assert adata.n_obs == 267547
        assert adata.n_vars == 3000
        assert "PCs" in adata.varm
    finally:
        sc.settings.datasetdir = original
