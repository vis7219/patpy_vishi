import os
from pathlib import Path

import pytest
from anndata import AnnData

import patpy
from patpy.datasets import DatasetInfo

CACHE_DIR = Path(os.environ["PATPY_TEST_DATASETDIR"]) if os.environ.get("PATPY_TEST_DATASETDIR") else None


def _datasetdir(tmp_path: Path) -> Path:
    """Return a persistent cache dir when available, otherwise fall back to tmp_path."""
    if CACHE_DIR is not None:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        return CACHE_DIR
    return tmp_path


def _check_dataset_info(adata: AnnData, info: DatasetInfo) -> None:
    """Common assertions tying a returned DatasetInfo back to the AnnData."""
    assert isinstance(info, DatasetInfo)
    assert info.n_cells == adata.n_obs
    assert info.n_features == adata.n_vars
    assert info.sample_key in adata.obs.columns
    assert info.cell_type_key in adata.obs.columns
    for col in info.sample_metadata_columns:
        assert col in adata.obs.columns


@pytest.mark.dataset
def test_combat_preprocessed_shape(tmp_path):
    original = patpy.settings.datasetdir
    patpy.settings.datasetdir = _datasetdir(tmp_path)
    try:
        adata = patpy.datasets.combat_preprocessed()

        assert isinstance(adata, AnnData)
        assert adata.n_obs == 783677
        assert adata.n_vars == 3000
        assert "PCs" in adata.varm

        adata, info = patpy.datasets.combat_preprocessed(return_dataset_info=True)
        _check_dataset_info(adata, info)
    finally:
        patpy.settings.datasetdir = original


@pytest.mark.dataset
def test_hlca_preprocessed_shape(tmp_path):
    original = patpy.settings.datasetdir
    patpy.settings.datasetdir = _datasetdir(tmp_path)
    try:
        adata = patpy.datasets.hlca_preprocessed()

        assert isinstance(adata, AnnData)
        assert adata.n_obs == 1687127
        assert adata.n_vars == 3000
        assert "PCs" in adata.varm

        adata, info = patpy.datasets.hlca_preprocessed(return_dataset_info=True)
        _check_dataset_info(adata, info)
    finally:
        patpy.settings.datasetdir = original


@pytest.mark.dataset
def test_onek1k_preprocessed_shape(tmp_path):
    original = patpy.settings.datasetdir
    patpy.settings.datasetdir = _datasetdir(tmp_path)
    try:
        adata = patpy.datasets.onek1k_preprocessed()

        assert isinstance(adata, AnnData)
        assert adata.n_obs == 1248980
        assert adata.n_vars == 3000
        assert "PCs" in adata.varm

        adata, info = patpy.datasets.onek1k_preprocessed(return_dataset_info=True)
        _check_dataset_info(adata, info)
    finally:
        patpy.settings.datasetdir = original


@pytest.mark.dataset
def test_stephenson_preprocessed_shape(tmp_path):
    original = patpy.settings.datasetdir
    patpy.settings.datasetdir = _datasetdir(tmp_path)
    try:
        adata = patpy.datasets.stephenson_preprocessed()

        assert isinstance(adata, AnnData)
        assert adata.n_obs == 639482
        assert adata.n_vars == 3000
        assert "PCs" in adata.varm

        adata, info = patpy.datasets.stephenson_preprocessed(return_dataset_info=True)
        _check_dataset_info(adata, info)
    finally:
        patpy.settings.datasetdir = original


@pytest.mark.dataset
def test_ticatlas_preprocessed_shape(tmp_path):
    original = patpy.settings.datasetdir
    patpy.settings.datasetdir = _datasetdir(tmp_path)
    try:
        adata = patpy.datasets.ticatlas_preprocessed()

        assert isinstance(adata, AnnData)
        assert adata.n_obs == 267547
        assert adata.n_vars == 3000
        assert "PCs" in adata.varm

        adata, info = patpy.datasets.ticatlas_preprocessed(return_dataset_info=True)
        _check_dataset_info(adata, info)
    finally:
        patpy.settings.datasetdir = original
