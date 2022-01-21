# Main target script for running everything.
#
# Author: Yuriy Sverchkov

src_folder <- here::here("src", "analysis", "YS")
data_folder <- here::here("data")
results_path <- here::here("results")

metadata_folder <- file.path(data_folder, "metadata")

library(targets)

source(file.path(src_folder, "R/functions.R"))

options(tidyverse.quiet = TRUE)
tar_option_set(packages= c("magrittr", "tidyverse", "ggplot2", "patchwork", "qvalue"))

# Targets list
list(

  # Templates

  tar_target(
    embedding_html_template,
    file.path(src_folder, "templates", "embedding_plot.html"),
    format = "file"
  ),

  tar_target(
    embedding_with_clusters_html_template,
    file.path(src_folder, "templates", "embedding_with_clusters_plot.html"),
    format = "file"
  ),

  # Metadata files

  tar_target(
    mxp_genes_file,
    file.path(metadata_folder, "mxp_genes.csv"),
    format = "file"
  ),

  tar_target(
    metadata_xlsx,
    file.path(metadata_folder, 'H3K Project Master Lists.xlsx'),
    format = "file"
  ),

  # Data files

  tar_target(
    data_xlsx,
    file.path(data_folder, 'Supplementary Table 3.xlsx'),
    format = "file"
  ),

  # QC

  tar_target(
    excluded_lines,
    c(
      'FUNDC2-KO2', 'HDHD3-KO1', 'HDHD3-KO2', 'MIPEP-KO', 'MPC1-KO2', 'MTRES1-KO', 'PISD-KO1', 'PISD-KO2',
      'SFXN3-KO1', 'SFXN3-KO2', 'SPATA20-KO1', 'SPATA20-KO2', 'TOMM7-KO2'
    )
  ),

  # Output folders
  tar_target(
    results_folder,
    {
      dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
      results_path
    },
    format = "file"  
  ),

  tar_target(
    tsne_folder,
    {
      path <- file.path(results_folder, "tsne")
      dir.create(path, showWarnings = FALSE, recursive = TRUE)
      path
    },
    format = "file"  
  ),

  tar_target(
    sister_kos_folder,
    {
      path <- file.path(results_folder, "sister-ko-analysis")
      dir.create(path, showWarnings = FALSE, recursive = TRUE)
      path
    },
    format = "file"  
  ),

  # Load data
  tar_target(lipids, read_lipids(data_xlsx, "Lipidomics log2FCsPvaluesSDs")),
  tar_target(metabolites, read_metabolites(data_xlsx, "Metabolomics log2FCsPvaluesSDs")),
  tar_target(proteins, read_proteins(data_xlsx, "Proteomics log2FCsPvaluesSDs")),

  # Make a combined molecule type dataframe
  tar_target(
    molecule_types,
    bind_rows(
      lipids$metadata %>% select(ID) %>% mutate(`molecule type` = 'Lipid'),
      metabolites$metadata %>% select(ID) %>% mutate(`molecule type` = 'Metabolite'),
      proteins$metadata %>% select(ID) %>% mutate(`molecule type` = 'Protein')
    )
  ),

  # Make combined molecules frame, and compute q-values
  tar_target(
    molecules_df,
    bind_rows(lipids$data, metabolites$data, proteins$data) %>%
      filter(
        !str_detect(`cell line`, '^P?WT'),
        !(`cell line` %in% excluded_lines)) %>%
      mutate(`Q Value` = qvalue(`P Value`)$qvalues)
  ),

  # Make a matrix of q-adjusted relative differences    
  tar_target(
    qardm,
    derive_QARD_matrix(molecules_df)
  ),

  # Compute tSNE embedding with Rtsne
  tar_target(
    qardm_tsne,
    compute_tsne_layout(
      qardm,
      n_tsne_runs = 10,
      n_pca_perm = 200,
      theta = 0.1,
      exaggeration_factor = 8.0,
      seed = 2109021610,
      scale = 'tsne_pk'
    )
  ),

  # Generate rSNE plot
  tar_target(
    qardm_tsne_plot,
    plot_embedding2(
      x = qardm_tsne$best$Y[, 1],
      y = qardm_tsne$best$Y[, 2],
      id = rownames(qardm),
      metadata = metadata_for_plot(lipids$metadata, metabolites$metadata, proteins$metadata),
      colors = c(
          "Lipid" = "#2cb781",
          "Metabolite" = "#e3b74a",
          "Protein" = "#99b0b5",
          "Protein (Mitochondria)" = "#003947",
          "Protein (Mitoribosome)" = "#9d2063",
          "Protein (OxPhos)" = "#ed2c51"),
      html_template = embedding_html_template,
      html_out = file.path(tsne_folder, "qardm-tsne.html")
    ),
    format = "file"
  ),

  # Generate tSNE neighbors table
  tar_target(
    tsne_neighbors_df,
    nearest_embeddings_neighbors(
      radius = 1.0,
      embeddings = qardm_tsne$tsnes,
      id = rownames(qardm),
      protein_metadata = proteins$metadata,
      molecule_types = molecule_types,
      mxp_genes = read_csv(mxp_genes_file),
      out_file = file.path(tsne_folder, "tsne-neighbors.csv")
    )
  ),

  # Clustering
  tar_target(
    tsne_hdbscan_selected,
    dbscan::hdbscan(qardm_tsne$best$Y, minPts = 6),
  ),

  # tSNE plot with core and noiseless HDBSCAN clusters
  tar_target(
    tsne_both_dbscan_selected_html,
    plot_embedding_clustered2(
      x = qardm_tsne$best$Y[, 1],
      y = qardm_tsne$best$Y[, 2],
      id = rownames(qardm),
      metadata = metadata_for_plot(lipids$metadata, metabolites$metadata, proteins$metadata),
      colors = c(
          "Lipid" = "#2cb781",
          "Metabolite" = "#e3b74a",
          "Protein" = "#99b0b5",
          "Protein (Mitochondria)" = "#003947",
          "Protein (Mitoribosome)" = "#9d2063",
          "Protein (OxPhos)" = "#ed2c51"),
      core_clusters = tsne_hdbscan_selected$cluster,
      extended_clusters = cluster_noise_points(qardm_tsne$best$Y, tsne_hdbscan_selected$cluster),
      html_template = embedding_with_clusters_html_template,
      html_out = file.path(tsne_folder, "dbscan-selected.html")
    ),
    format = "file"
  ),

  # KO-specific phenotype detection
  tar_target(
    outliers_csv,
    {
      csv_file <- file.path(results_folder, 'ksp_outliers.csv')
      molecules_df %>%
        outlier_distances(
            `Control Normalized Log2 Fold Change`,
            `P Value`,
            k = 3
        ) %>%
        write_csv(csv_file)
      csv_file
    },
    format = "file"
  ),

  # KO-KO similarity charactarization
  ###################################

  tar_target(
    cell_line_metadata,
    read_cell_line_metadata(metadata_xlsx),
  ),

  # Make calls
  tar_target(
    long_calls,
    molecules_df %>% make_calls()
  ),

  # Make call profile
  tar_target(
    call_profile,
    long_calls %>% long_to_profile_matrix(call, remove_zeros = 'both'),
  ),

  # Make filtered l2fc profile
  tar_target(
    filtered_profile,
    long_calls %>% long_to_profile_matrix(sigbigL2FC, remove_zeros = 'both')
  ),

  # Make correlation matrix
  tar_target(
    spearman_ko_cor,
    filtered_profile %>% t() %>% cor(method = "spearman")
  ),

  # Save correlation matrix
  tar_target(
    saved_spearman_ko_cor,
    {
    csv_file <- file.path(sister_kos_folder, "spearman_cor.csv")
    spearman_ko_cor %>%
      as_tibble(rownames = 'KO') %>%
      write_csv(file_out(csv_file))
    csv_file
    },
    format = "file"
  ),

  # Make plots
  tar_target(
    plot_spearman_ranks,
    plot_sister_ko_cor_ranks(
      correlation_matrix = spearman_ko_cor,
      plot_file = file.path(sister_kos_folder, "spearman_ranks.pdf"),
      #full_match_table_file = file.path(sister_kos_folder, "full_match_table.csv"),
      neat_match_table_file = file.path(sister_kos_folder, "neat_match_table.csv"),
      profile_matrix = call_profile,
      sort_by = 'best rank',
      title = 'Spearman correlations between sister KO profiles',
      pair_properties = get_pair_properties(cell_line_metadata)
    ),
    format = "file"
  )
)
