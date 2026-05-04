from pathlib import Path

import scanpy as sc
from scanpy import settings

from patpy.datasets._dataloader import _download


def combat_preprocessed(
    overwrite: bool = False,
):
    """Processed COvid-19 Multi-omics Blood ATlas (COMBAT) dataset.

    Here, data was preprocessed according to the standard scanpy pipeline, cells annotated as "nan" removed, PCA, scVI, scANVI and scPoli dimensionality reduction and reduction were applied.

    The dataset contains 783'677 cells and 3,000 features.
    The function downloads a zip compressed file of approximately 1.5 GB, which is unzipped to an :class:`~anndata.AnnData` file of approximately 5 GB. Download with good internet connection takes approximately 2 min.

    References
    ----------
        Ahern, D. J., Ai, Z., Ainsworth, M., Allan, C., Allcock, A., Angus, B., ... & Salio, M. (2022). A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. Cell, 185(5), 916-938. https://doi.org/10.1016/j.cell.2022.01.012.
        COvid-19 Multi-omics Blood ATlas (COMBAT) Consortium. (2021). A blood atlas of COVID-19 defines hallmarks of disease severity and specificity: Associated data (1.0.1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6120249

    Returns
    -------
        :class:`~anndata.AnnData` object of scRNA-seq profiles.

    Examples
    --------
        >>> import patpy
        >>> adata = patpy.dt.combat_preprocessed()
        >>> adata
        AnnData object with n_obs × n_vars = 783677 × 3000
            obs: 'Annotation_cluster_id', 'Annotation_cluster_name', 'Annotation_minor_subset', 'Annotation_major_subset', 'Annotation_cell_type', 'GEX_region', 'QC_ngenes', 'QC_total_UMI', 'QC_pct_mitochondrial', 'QC_scrub_doublet_scores', 'TCR_chain_composition', 'TCR_clone_ID', 'TCR_clone_count', 'TCR_clone_proportion', 'TCR_contains_unproductive', 'TCR_doublet', 'TCR_chain_TRA', 'TCR_v_gene_TRA', 'TCR_d_gene_TRA', 'TCR_j_gene_TRA', 'TCR_c_gene_TRA', 'TCR_productive_TRA', 'TCR_cdr3_TRA', 'TCR_umis_TRA', 'TCR_chain_TRA2', 'TCR_v_gene_TRA2', 'TCR_d_gene_TRA2', 'TCR_j_gene_TRA2', 'TCR_c_gene_TRA2', 'TCR_productive_TRA2', 'TCR_cdr3_TRA2', 'TCR_umis_TRA2', 'TCR_chain_TRB', 'TCR_v_gene_TRB', 'TCR_d_gene_TRB', 'TCR_j_gene_TRB', 'TCR_c_gene_TRB', 'TCR_productive_TRB', 'TCR_chain_TRB2', 'TCR_v_gene_TRB2', 'TCR_d_gene_TRB2', 'TCR_j_gene_TRB2', 'TCR_c_gene_TRB2', 'TCR_productive_TRB2', 'TCR_cdr3_TRB2', 'TCR_umis_TRB2', 'BCR_umis_HC', 'BCR_contig_qc_HC', 'BCR_functionality_HC', 'BCR_v_call_HC', 'BCR_v_score_HC', 'BCR_j_call_HC', 'BCR_j_score_HC', 'BCR_junction_aa_HC', 'BCR_total_mut_HC', 'BCR_s_mut_HC', 'BCR_r_mut_HC', 'BCR_c_gene_HC', 'BCR_clone_per_replicate_HC', 'BCR_clone_global_HC', 'BCR_clonal_abundance_HC', 'BCR_locus_LC', 'BCR_umis_LC', 'BCR_contig_qc_LC', 'BCR_functionality_LC', 'BCR_v_call_LC', 'BCR_v_score_LC', 'BCR_j_call_LC', 'BCR_j_score_LC', 'BCR_junction_aa_LC', 'BCR_total_mut_LC', 'BCR_s_mut_LC', 'BCR_r_mut_LC', 'BCR_c_gene_LC', 'COMBAT_ID', 'scRNASeq_sample_ID', 'COMBAT_participant_timepoint_ID', 'Source', 'Age', 'Sex', 'Race', 'BMI', 'Hospitalstay', 'Death28', 'Institute', 'PreExistingHeartDisease', 'PreExistingLungDisease', 'PreExistingKidneyDisease', 'PreExistingDiabetes', 'PreExistingHypertension', 'PreExistingImmunocompromised', 'Smoking', 'Symptomatic', 'Requiredvasoactive', 'Respiratorysupport', 'SARSCoV2PCR', 'Outcome', 'TimeSinceOnset', 'Ethnicity', 'Tissue', 'DiseaseClassification', 'Pool_ID', 'Channel_ID', 'ifn_1_score', '_scvi_batch', '_scvi_labels', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'vertex', 'eigenvector_centrality'
            var: 'gene_ids', 'feature_types', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'
            uns: 'Institute', 'ObjectCreateDate', 'Source_colors', 'Technology', 'X_gloscope_cuml_distances', 'X_gloscope_pynndescent_distances', 'X_scpoli', '_scvi_manager_uuid', '_scvi_uuid', 'genome_annotation_version', 'gloscope_representation', 'gloscope_scpoli_distances', 'hvg', 'log1p', 'neighbors', 'pca', 'scpoli_distances', 'scpoli_parameters', 'scpoli_samples'
            obsm: 'X_pca', 'X_scANVI_batch', 'X_scANVI_sample', 'X_scVI_batch', 'X_scVI_sample', 'X_scpoli', 'X_umap', 'X_umap_source'
            varm: 'PCs'
            layers: 'X_raw_counts'
            obsp: 'connectivities', 'distances'
    """
    output_file_name = "combat_processed.h5ad"
    output_file_path = settings.datasetdir / output_file_name
    if not Path(output_file_path).exists() or overwrite:
        _download(
            url="https://ndownloader.figshare.com/files/64217586",
            output_file_name="combat_processed.h5ad.zip",
            output_path=settings.datasetdir,
            is_zip=True,
        )
        (settings.datasetdir / "combat_processed.h5ad.zip").unlink(missing_ok=True)
    adata = sc.read_h5ad(output_file_path)

    return adata

def hlca_preprocessed(
    overwrite: bool = False,
):
    """Processed Human Lung Cell Atlas (HLCA) dataset.

    Here, data was preprocessed according to the standard scanpy pipeline, cells annotated as "nan" removed, PCA, scVI, scANVI and scPoli dimensionality reduction and reduction were applied.

    The dataset contains 1,687,127 cells and 3,000 features.
    The function downloads a zip compressed file of approximately 3 GB, which is unzipped to an :class:`~anndata.AnnData` file of approximately 6.5 GB. Download with good internet connection takes approximately 3 min.

    References:
        Sikkema, L., Ramírez-Suástegui, C., Strobl, D. C., Gillett, T. E., Zappia, L., Madissoon, E., ... & Theis, F. J. (2023). An integrated cell atlas of the lung in health and disease. Nature medicine, 29(6), 1563-1577. https://doi.org/10.1038/s41591-023-02327-2

    Returns:
        :class:`~anndata.AnnData` object of scRNA-seq profiles.

    Examples:
        >>> import patpy
        >>> adata = patpy.dt.hlca_preprocessed()
        >>> adata
        AnnData object with n_obs × n_vars = 1687127 × 3000
            obs: 'suspension_type', 'donor_id', 'is_primary_data', 'assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 'disease_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'tissue_ontology_term_id', 'organism_ontology_term_id', 'sex_ontology_term_id', "3'_or_5'", 'BMI', 'age_or_mean_of_age_range', 'age_range', 'anatomical_region_ccf_score', 'ann_coarse_for_GWAS_and_modeling', 'ann_finest_level', 'ann_level_1', 'ann_level_2', 'ann_level_3', 'ann_level_4', 'ann_level_5', 'cause_of_death', 'core_or_extension', 'dataset', 'fresh_or_frozen', 'log10_total_counts', 'lung_condition', 'mixed_ancestry', 'original_ann_level_1', 'original_ann_level_2', 'original_ann_level_3', 'original_ann_level_4', 'original_ann_level_5', 'original_ann_nonharmonized', 'reannotation_type', 'sample', 'scanvi_label', 'sequencing_platform', 'smoking_status', 'study', 'subject_type', 'tissue_coarse_unharmonized', 'tissue_detailed_unharmonized', 'tissue_dissociation_protocol', 'tissue_level_2', 'tissue_level_3', 'tissue_sampling_method', 'total_counts', 'transf_ann_level_1_label', 'transf_ann_level_1_uncert', 'transf_ann_level_2_label', 'transf_ann_level_2_uncert', 'transf_ann_level_3_label', 'transf_ann_level_3_uncert', 'transf_ann_level_4_label', 'transf_ann_level_4_uncert', 'transf_ann_level_5_label', 'transf_ann_level_5_uncert', 'tissue_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid', '_scvi_batch', '_scvi_labels', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'vertex', 'eigenvector_centrality'
            var: 'feature_is_filtered', 'original_gene_symbols', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'
            uns: 'X_gloscope_cuml_distances', '_scvi_manager_uuid', '_scvi_uuid', 'batch_condition', 'citation', 'default_embedding', 'gloscope_representation', 'hvg', 'log1p', 'neighbors', 'pca', 'schema_reference', 'schema_version', 'title'
            obsm: 'X_pca', 'X_scANVI_batch', 'X_scANVI_sample', 'X_scVI_batch', 'X_scVI_sample', 'X_scanvi_emb', 'X_umap'
            varm: 'PCs'
            layers: 'X_raw_counts', 'soupX'
            obsp: 'connectivities', 'distances'
    """

    output_file_name = "hlca_processed.h5ad"
    output_file_path = settings.datasetdir / output_file_name
    if not Path(output_file_path).exists() or overwrite:
        _download(
            url="https://ndownloader.figshare.com/files/64225983",
            output_file_name="hlca_processed.h5ad.zip",
            output_path=settings.datasetdir,
            is_zip=True,
        )
        (settings.datasetdir / "hlca_processed.h5ad.zip").unlink(missing_ok=True)
    adata = sc.read_h5ad(output_file_path)

    return adata

def onek1k_preprocessed(
    overwrite: bool = False,
):
    """Processed onek1k dataset.

    Here, data was preprocessed according to the standard scanpy pipeline, cells annotated as "nan" removed, PCA, scVI, scANVI and scPoli dimensionality reduction and reduction were applied.

    The dataset contains 1,248,980 cells and 3,000 features.
    The function downloads a zip compressed file of approximately 2.5 GB, which is unzipped to an :class:`~anndata.AnnData` file of approximately 4 GB. Download with good internet connection takes approximately 3 min.

    References
    ----------
        Yazar, S., Alquicira-Hernandez, J., Wing, K., Senabouth, A., Gordon, M. G., Andersen, S., ... & Powell, J. E. (2022). Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease. Science, 376(6589), eabf3041. https://doi.org/10.1126/science.abf3041
        https://onek1k.org/

    Returns
    -------
        :class:`~anndata.AnnData` object of scRNA-seq profiles.

    Examples
    --------
        >>> import patpy
        >>> adata = patpy.dt.onek1k_preprocessed()
        >>> adata
        AnnData object with n_obs × n_vars = 1248980 × 3000
            obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'donor_id', 'pool_number', 'predicted.celltype.l2', 'predicted.celltype.l2.score', 'age', 'organism_ontology_term_id', 'tissue_ontology_term_id', 'assay_ontology_term_id', 'disease_ontology_term_id', 'cell_type_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'is_primary_data', 'suspension_type', 'tissue_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid', '_scvi_batch', '_scvi_labels', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'vertex', 'eigenvector_centrality'
            var: 'vst.mean', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized', 'vst.variable', 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'
            uns: 'X_gloscope_cuml_distances', 'X_scpoli', '_scvi_manager_uuid', '_scvi_uuid', 'cell_type_ontology_term_id_colors', 'citation', 'default_embedding', 'gloscope_representation', 'hvg', 'log1p', 'neighbors', 'pca', 'schema_reference', 'schema_version', 'scpoli_distances', 'scpoli_parameters', 'scpoli_samples', 'title'
            obsm: 'X_azimuth_spca', 'X_azimuth_umap', 'X_harmony', 'X_pca', 'X_scANVI_batch', 'X_scANVI_sample', 'X_scVI_batch', 'X_scVI_sample', 'X_scpoli', 'X_umap'
            varm: 'PCs'
            layers: 'X_raw_counts'
            obsp: 'connectivities', 'distances'
    """
    output_file_name = "onek1k_processed.h5ad"
    output_file_path = settings.datasetdir / output_file_name
    if not Path(output_file_path).exists() or overwrite:
        _download(
            url="https://ndownloader.figshare.com/files/64225884",
            output_file_name="onek1k_processed.h5ad.zip",
            output_path=settings.datasetdir,
            is_zip=True,
        )
        (settings.datasetdir / "onek1k_processed.h5ad.zip").unlink(missing_ok=True)
    adata = sc.read_h5ad(output_file_path)

    return adata

def stephenson_preprocessed(
    overwrite: bool = False,
):
    """Processed multi-omics immune response in COVID-19 (Stephenson) dataset.

    Here, data was preprocessed according to the standard scanpy pipeline, cells annotated as "nan" removed, PCA, scVI, scANVI and scPoli dimensionality reduction and reduction were applied.

    The dataset contains 639,482 cells and 3,000 features.
    The function downloads a zip compressed file of approximately 1.5 GB, which is unzipped to an :class:`~anndata.AnnData` file of approximately 4.5 GB. Download with good internet connection takes approximately 2 min.

    References:
        Stephenson, E., Reynolds, G., Botting, R. A., Calero-Nieto, F. J., Morgan, M. D., Tuong, Z. K., ... & Haniffa, M. (2021). Single-cell multi-omics analysis of the immune response in COVID-19. Nature medicine, 27(5), 904-916. https://doi.org/10.1038/s41591-021-01329-2

    Returns:
        :class:`~anndata.AnnData` object of scRNA-seq profiles.

    Examples:
        >>> import patpy
        >>> adata = patpy.dt.stephenson_preprocessed()
        >>> adata
        AnnData object with n_obs × n_vars = 639482 × 3000
            obs: 'sample_id', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'initial_clustering', 'Resample', 'Collection_Day', 'Swab_result', 'Status', 'Smoker', 'Status_on_day_collection', 'Status_on_day_collection_summary', 'Days_from_onset', 'Site', 'time_after_LPS', 'Worst_Clinical_Status', 'Outcome', 'donor_id', 'assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 'disease_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'is_primary_data', 'organism_ontology_term_id', 'sex_ontology_term_id', 'tissue_ontology_term_id', 'author_cell_type', 'suspension_type', 'tissue_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid', '_scvi_batch', '_scvi_labels', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'vertex', 'eigenvector_centrality'
            var: 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'
            uns: 'X_gloscope_cuml_distances', 'X_gloscope_pynndescent_distances', 'X_scpoli', '_scvi_manager_uuid', '_scvi_uuid', 'antibody_X', 'antibody_features', 'antibody_raw.X', 'citation', 'default_embedding', 'gloscope_representation', 'gloscope_scpoli_distances', 'hvg', 'leiden', 'log1p', 'neighbors', 'pca', 'schema_reference', 'schema_version', 'scpoli_distances', 'scpoli_parameters', 'scpoli_samples', 'title', 'umap'
            obsm: 'X_pca', 'X_pca_harmony', 'X_scANVI_batch', 'X_scANVI_sample', 'X_scVI_batch', 'X_scVI_sample', 'X_scpoli', 'X_umap'
            varm: 'PCs'
            layers: 'X_raw_counts'
            obsp: 'connectivities', 'distances'
    """

    output_file_name = "stephenson_processed.h5ad"
    output_file_path = settings.datasetdir / output_file_name
    if not Path(output_file_path).exists() or overwrite:
        _download(
            url="https://ndownloader.figshare.com/files/64226109",
            output_file_name="stephenson_processed.h5ad.zip",
            output_path=settings.datasetdir,
            is_zip=True,
        )
        (settings.datasetdir / "stephenson_processed.h5ad.zip").unlink(missing_ok=True)
    adata = sc.read_h5ad(output_file_path)

    return adata

def ticatlas_preprocessed(
    overwrite: bool = False,
):
    """Processed Tumor Immune Cell Atlas (TICAtlas) dataset.

    Here, data was preprocessed according to the standard scanpy pipeline, cells annotated as "nan" removed, PCA, scVI, scANVI and scPoli dimensionality reduction and reduction were applied.

    The dataset contains 267,547 cells and 3,000 features.
    The function downloads a zip compressed file of approximately 0.5 GB, which is unzipped to an :class:`~anndata.AnnData` file of approximately 1.8 GB. Download with good internet connection takes approximately 2 min.

    References:
        Nieto, P., Elosua-Bayes, M., Trincado, J. L., Marchese, D., Massoni-Badosa, R., Salvany, M., ... & Heyn, H. (2021). A single-cell tumor immune atlas for precision oncology. Genome research, 31(10), 1913-1926. https://doi.org/10.1101/gr.273300.120

    Returns:
        :class:`~anndata.AnnData` object of scRNA-seq profiles.

    Examples:
        >>> import patpy
        >>> adata = patpy.dt.ticatlas_preprocessed()
        >>> adata
        AnnData object with n_obs × n_vars = 267547 × 3000
            obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'patient', 'percent.mt', 'gender', 'subtype', 'source', 'lv1_annot', 'lv2_annot', 'kmeans_cluster', 'n_genes', '_scvi_batch', '_scvi_labels', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo'
            var: 'features', 'n_cells', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'
            uns: 'X_scpoli', '_scvi_manager_uuid', '_scvi_uuid', 'hvg', 'log1p', 'pca', 'scpoli_distances', 'scpoli_parameters', 'scpoli_samples'
            obsm: 'X_pca', 'X_scANVI_batch', 'X_scANVI_sample', 'X_scVI_batch', 'X_scVI_sample', 'X_scpoli', 'X_umap'
            varm: 'PCs'
            layers: 'X_raw_counts', 'shifted_log_counts'
    """

    output_file_name = "ticatlas_processed.h5ad"
    output_file_path = settings.datasetdir / output_file_name
    if not Path(output_file_path).exists() or overwrite:
        _download(
            url="https://ndownloader.figshare.com/files/64226097",
            output_file_name="ticatlas_processed.h5ad.zip",
            output_path=settings.datasetdir,
            is_zip=True,
        )
        (settings.datasetdir / "ticatlas_processed.h5ad.zip").unlink(missing_ok=True)
    adata = sc.read_h5ad(output_file_path)

    return adata

