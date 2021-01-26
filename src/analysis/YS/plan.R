# Drake plan for running everything.
# This version loads everything into the global environment.
#
# Author: Yuriy Sverchkov

# Check that required variables are set
if (!exists('src_folder')) stop("Set the src_folder variable before running this.")
if (!exists('data_folder')) stop("Set the data_folder variable before running this.")
if (!exists('results_folder')) stop("Set the results_folder variable before running this.")

metadata_folder <- file.path(data_folder, "metadata")

# Input files:
sample_matching_file <- file.path(metadata_folder, "sample_control_matches.txt")

proteins_file <- file.path(data_folder,
    "combined_LFQ_filtered_imputed_combat.tsv")

metabolites_file <- file.path(data_folder,
    "metabolomics_data_w_batch.csv")

lipids_file <- file.path(data_folder,
    "combined_lipidomics_data_filtered.tsv")

# Deliverable:
excel_file <- file.path(results_folder, "molecule_data.xlsx")

# Protein metadata file
protein_metadata_file <- file.path(metadata_folder, "protein_metadata.csv")

## tSNE-related output files
embedding_html_template <- file.path(src_folder, "templates", "embedding_plot.html")
mrdm_tsne_html2 <- file.path(results_folder, "tsne_figure.html")

tsne_neighbors <- file.path(results_folder, "tsne_neighbors.csv")
mxp_genes_file <- file.path(metadata_folder, "mxp_genes.csv")


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
        long_multiomics,
        read_protein_metadata(file_in(!!protein_metadata_file))),

    # Write excel file
    wrote_excel = write_excel_tables(
        long_multiomics,
        molecule_metadata,
        file_out(!!excel_file)),

    # Make a matrix of q-adjusted relative differences    
    mrdm = derive_mrdm_tuple(long_multiomics),

    # Compute perplexity for tSNE
    mrdm_tsne_perplexity = sqrt(nrow(mrdm$matrix)),

    # Compute tSNE embedding
    mrdm_tsne = tsne::tsne(
        X = mrdm$matrix,
        k = 2,
        perplexity = mrdm_tsne_perplexity),

    # Generate tSNE plot
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

    # Generate tSNE neighbors table
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
