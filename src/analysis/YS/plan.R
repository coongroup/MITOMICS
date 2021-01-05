# Drake plan for running everything.
# This version loads everything into the global environment.
#
# Author: Yuriy Sverchkov

# Input files:
sample_matching_file <- here::here(
    "data-raw",
    "DataForYuriy",
    "h3k_Samples_RespectiveControls.txt")

proteins_file <- here::here(
    "data-raw",
    "DataForYuriy",
    "2020-05-12_Combined_LFQ_imputed_filtered_ComBat.tsv")

proteins_file_w_metadata <- here::here(
    "data-raw",
    "H3KDataForSpotChecking",
    "20200514_H3K_proteomics_avg_log2_FC_w_metadata.csv")

metabolites_file <- here::here(
    "data-raw",
    "DataForYuriy",
    "20191106_IJM_metabolomics_data_forCombat.csv")

lipids_file <- here::here(
    "data-raw",
    "DataForYuriy",
    "2020-06-04_combined_lipidomics_data_collapsed.tsv")

# Deliverable:
excel_file <- here::here(
    "results",
    "molecule_data.xlsx"
)

## JWR's excel file and related constants
xls_jwr_file <- here::here(
    "data-raw",
    "molecule_data_201001_JWR.xlsx"
)

# Protein metadata column types selector
xls_jwr_pmct <- c(rep("text", 14), rep("skip", 203))

# Column renaming
xls_jwr_rename <- c("CI assmebly factor" = "CI assembly factor")

# Columns to smartly convert to logical (after rename)
xls_jwr_lgl_cols <- c(
    "KO target",
    "MitoCarta2.0",
    "mtDNA encoded",
    "CI assembly factor",
    "Complex Q")

# Protein metadata file
protein_metadata_file <- here::here(
    "data-clean",
    "protein_metadata.csv"
)

## tSNE-related output files
embedding_html_template <- here::here("templates", "embedding_plot.html")
mrdm_tsne_html2 <- here::here("figures", "tsne_20201106.html")

tsne_neighbors <- here::here("results", "tsne_neighbors.csv")
mxp_genes_file <- here::here("data_manual", "mxp_genes.csv")


############
# The plan #
############
plan <- drake_plan(

    # Reshape into standardized long table
    protein_intensities =
        make_protein_intensity_table(file_in(!!proteins_file)),
    metabolite_intensities =
        make_metabolite_intensity_table(file_in(!!metabolites_file)),
    lipid_intensities =
        make_lipid_intensity_table(file_in(!!lipids_file)),

    # Build batch metadata
    batch_metadata = extract_batch_metadata(
        get_protein_condition_matches(file_in(!!sample_matching_file)),
        metabolite_intensities,
        lipid_intensities,
        extract_condition_mapping(file_in(!!metabolites_file))),

    # And compute statistics
    long_protein_df = compute_protein_l2fc(
        protein_intensities,
        get_protein_condition_matches(file_in(!!sample_matching_file)),
        test_strategy = "technical avg"),
    long_metabolite_df = compute_l2fc_per_batch_wt(metabolite_intensities),
    long_lipid_df = compute_l2fc_per_batch_wt(lipid_intensities),

    # Compose one big table, normalize condition names & run fdr corrections
    long_multiomics = compose_long_multiomic(
        list(
            Protein = long_protein_df,
            Metabolite = long_metabolite_df,
            Lipid = long_lipid_df),
        extract_condition_mapping(file_in(!!metabolites_file))),

    # Build metadata
    molecule_metadata = build_metadata(
        distinct(long_multiomics, ID, `molecule type`),
        file_in(!!proteins_file_w_metadata),
        uniprot_genes),

    # Write excel file
    wrote_excel = write_excel_tables(
        long_multiomics,
        molecule_metadata,
        file_out(!!excel_file)),

    # Protein-specific metadata
    wrote_protein_metadata = protein_metadata_from_jwr(
        file_in(!!xls_jwr_file),
        xls_jwr_pmct,
        xls_jwr_rename,
        xls_jwr_lgl_cols,
        file_out(!!protein_metadata_file)
    ),
    
    mrdm = derive_mrdm_tuple(long_multiomics),

    mrdm_tsne_perplexity = sqrt(nrow(mrdm$matrix)),
    mrdm_tsne = tsne::tsne(
        X = mrdm$matrix,
        k = 2,
        perplexity = mrdm_tsne_perplexity),

    mrdm_tsne_plot2 = plot_embedding2(
        x = mrdm_tsne[, 1],
        y = mrdm_tsne[, 2],
        id = rownames(mrdm$matrix),
        metadata = metadata_20201007(molecule_metadata, read_protein_metadata(file_in(!!protein_metadata_file))),
        colors = c(
            "Lipid" = "#2cb781",
            "Metabolite" = "#e3b74a",
            "Protein" = "#99b0b5",
            "Protein (Mitochondria)" = "#003947",
            "Protein (Mitoribosome)" = "#9d2063",
            "Protein (OxPhos)" = "#ed2c51"),
        html_template = file_in(!!embedding_html_template),
        html_out = file_out(!!mrdm_tsne_html2)
    ),

    tsne_neighbors_df = nearest_embedding_neighbors(
        radius = 1.0,
        embedding = mrdm_tsne[, 1:2],
        id = rownames(mrdm$matrix),
        protein_metadata = read_protein_metadata(file_in(!!protein_metadata_file)),
        molecule_types = molecule_metadata %>% select(ID, `molecule type`),
        mxp_genes = read_csv(file_in(!!mxp_genes_file)),
        out_file = file_out(!!(tsne_neighbors))
    )
)
