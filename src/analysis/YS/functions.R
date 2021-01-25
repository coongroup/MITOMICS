# Function definitions
#
# Licensed under the BSD 3-Clause License
# Copyright (c) 2019, Yuriy Sverchkov

library(magrittr)
library(drake)
library(tidyverse)
library(tidygraph)
library(visNetwork)
library(plotly)
library(qvalue)


make_protein_intensity_table <- function(proteins_file) {
    # Assumption: Proteins are rows, replicates are columns.
    # Assumption: Protein ID is first column
    # Assumption: Replicates are encoded as LFQ.intensity.[condition]_rep[X]
    # Assumption: Values in table are log2-transformed

    protein_table <- read_tsv(proteins_file) %>% rename(ID = 1)

    intensity_table <- pivot_longer(
        protein_table,
        cols = -ID,
        names_to = "replicate",
        values_to = "log2 intensity")

    # Split out replicate IDS from cell lines
    rep_matrix <- str_match(
        pull(intensity_table, replicate),
        "^LFQ[.]intensity[.](.+)_rep([A-Z])$")

    intensity_table <- mutate(intensity_table,
        condition = rep_matrix[, 2],
        rep = rep_matrix[, 3])

    return(intensity_table)
}


make_metabolite_intensity_table <- function(filename) {
    # Assumptions:
    # * Rows are samples and columns are molecules
    # * Sample names are Bath*_Day*_[condition]M for KO strains and
    #   Batch*_Day*_PWT for pooled wild type
    # * First column is sample id
    # * Also contains metadata columns "batch" and "H3K_cell_line"

    metabolites_table <- read_csv(filename) %>%
        rename(replicate = 1) %>%
        select(-batch, -H3K_cell_line)

    intensity_table <- pivot_longer(
        metabolites_table,
        cols = -replicate,
        names_to = "ID",
        values_to = "log2 intensity")

    # Split out replicate IDS from cell lines
    rep_matrix <- str_match(
        pull(intensity_table, replicate),
        "^(Batch.+)_Day.+_([0-9]+M|PWT)$")

    c_vector <- rep_matrix[, 3] %>% str_replace("M$", "")

    # Annotate conditions
    return(mutate(intensity_table,
        condition = c_vector,
        batch = rep_matrix[, 2]))
}


make_lipid_intensity_table <- function(filename) {
    # Assumptions:
    # * Rows are samples and columns are molecules
    # * Sample names are "A[condition]L{-#} Quant Values" for KO strains and
    #   "PWT{-#-#} Quant Values" for pooled wild type
    # * First column is sample id
    # * Also contains metadata columns "batch" and "H3K_cell_line"

    metabolites_table <- read_tsv(filename)

    intensity_table <- pivot_longer(
        metabolites_table,
        cols = c(-replicate, -batch),
        names_to = "ID",
        values_to = "log2 intensity")

    # Split out replicate IDS from cell lines
    rep_matrix <- str_match(
        pull(intensity_table, replicate),
        "^(PWT|[A-Z]([0-9]+)L).* Quant Values$")

    c_vector <- if_else(
        is.na(rep_matrix[, 3]),
        rep_matrix[, 2],
        rep_matrix[, 3])

    # Annotate conditions
    return(intensity_table %>% mutate(condition = c_vector))
}


#' Compute intensity summary statistics
#' 
summarize_intensities <- function(df)
    summarize(df,
        samples = n(),
        `mean of log2 intensity` = mean(`log2 intensity`),
        `SD of log2 intensity` = sd(`log2 intensity`),
        `log2 intensities` = list(`log2 intensity`))


#' Function to extract metadata about batches
#' 
#' @param protein_batch_matches Dataframe with columns "condition" and "control"
#' @param lipid_intensities Lipid intensities table with columns "condition"
#'      and "batch"
#' @param metabolite_intensities Metabolite intensities table with columns
#'      "condition" and "batch"
#' @param condition_mapping Dataframe with columns "condition" and "cell line"
extract_batch_metadata <- function(
    protein_batch_matches,
    metabolite_intensities,
    lipid_intensities,
    condition_mapping
) {
    # Proteins
    protein_batch_matches %<>%
        mutate(batch = str_replace(control, "_pWT$", "")) %>%
        distinct(condition, batch) %>%
        mutate(
            batch = batch %>%
                str_replace_all("(\\D)(\\d)", "\\1 \\2") %>%
                str_replace_all("(\\d)(\\D)", "\\1 \\2") %>%
                str_replace_all("( )(\\d)( |$)", "\\1 \\2\\3"),
            condition = as.integer(condition),
            `molecule type` = "Protein")
    # Metabolites
    metabolite_intensities %<>%
        distinct(condition, batch) %>%
        mutate(
            batch = batch %>%
                str_replace_all("(\\D)(\\d)", "\\1 \\2") %>%
                str_replace_all("(\\d)(\\D)", "\\1 \\2") %>%
                str_replace_all("( )(\\d)( |$)", "\\1 \\2\\3"),
            condition = as.integer(condition),
            `molecule type` = "Metabolite")
    # Lipids
    lipid_intensities %<>%
        distinct(condition, batch) %>%
        mutate(batch = paste("Batch", batch),
        condition = as.integer(condition),
        `molecule type` = "Lipid")
    # Join cell line names
    bind_rows(
        protein_batch_matches,
        metabolite_intensities,
        lipid_intensities) %>%
        left_join(condition_mapping, by = "condition")
}


#' Compose long multiomics table
#' 
#' Composes a combined table with adjusted p-values
#' 
#' @param list_of_dfs A named list of long-shape dataframes that are expected
#' to contain the LFC values and p-values.
#' List element names will become molecule types in the output.
#' Required columns: ID, condition, `t-test p-value` (other columns will be
#' preserved as-is in the output)
#' @param condition_map A dataframe containing condition-to-cell-line-name
#' mappings. Required columns: condition (expected to be integers at least for
#' KO conditions), `cell line`, and `WT condition` (logical).
#' @param method string specifying method for producing adjusted p-values or q-values.
#' @return A dataframe consisting of a row-bind of the dataframes from
#' list_of_dfs with added columns
#' `molecule type` (taken from list names),
#' `cell line` (taken from a join with condition_map),
#' `q-value (multi-omic)` or `BH-adjusted multi-omic p-value`, and
#' `q-value (single-omic)` or `BH-adjusted single-omic p-value`
compose_long_multiomic <- function(list_of_dfs, condition_map, method = "q-value") {
    full_df <- imap_dfr(list_of_dfs, function(df, m_type) {
        df %>%
            mutate(
                `molecule type` = m_type,
                condition = as.integer(condition)) %>%
            left_join(condition_map, by = "condition")
    })

    adjusted <- switch(method,
        "q-value" = full_df %>%
            filter(!`WT condition`) %>%
            mutate(`q-value (multi-omic)` =
                qvalue(`t-test p-value`)$qvalues) %>%
            group_by(`molecule type`) %>%
            mutate(`q-value (single-omic)` =
                qvalue(`t-test p-value`)$qvalues) %>%
            ungroup() %>%
            select(
                ID,
                `cell line`,
                `q-value (multi-omic)`,
                `q-value (single-omic)`),
        "BH" = full_df %>%
            filter(!`WT condition`) %>%
            mutate(`BH-adjusted multi-omic p-value` =
                p.adjust(`t-test p-value`, method = "BH")) %>%
            group_by(`molecule type`) %>%
            mutate(`BH-adjusted single-omic p-value` =
                p.adjust(`t-test p-value`, method = "BH")) %>%
            ungroup() %>%
            select(
                ID,
                `cell line`,
                `BH-adjusted multi-omic p-value`,
                `BH-adjusted single-omic p-value`)
    )

    left_join(full_df, adjusted, by = c("ID", "cell line"))
}


get_protein_condition_matches <- function(sample_matching_file) {
    # Assumption: Sample column is called `Sample` and has the form
    #   [Condition]_rep[Replicate]
    # Assumption: Control column is called `Respective pooled control`
    read_tsv(sample_matching_file) %>%
        rename(control = `Respective pooled control`) %>%
        mutate(condition = str_match(Sample, "^(.+)_rep[A-Z]$")[, 2]) %>%
        distinct(condition, control)
}


extract_condition_mapping <- function(file) {
    # Assumptions:
    # * First column is sample ID
    # * Sample ID format is Batch*_Day*_[condition #]M OR ends in PWT
    # * Cell line column is H2K_cell_line
    # * Cell line condition mapping is one-to-one except for WT-# cell lines
    #   which have multiple conditions per cell line.

    read_csv(file) %>%
        rename(sample = 1, `cell line` = H3K_cell_line) %>%
        select(sample, `cell line`) %>%
        mutate(condition =
            str_match(sample, "^Batch.+_Day.+_([0-9]+)M$")[, 2] %>%
            as.integer()) %>%
        filter(!is.na(condition)) %>%
        distinct(condition, `cell line`) %>%
        mutate(
            `WT condition` = str_detect(`cell line`, "^WT-"),
            `cell line` = if_else(
                `WT condition`,
                paste0(`cell line`, " (", condition, ")"),
                `cell line`))
}


#' Function for computing protein log-fold-changes and statistics.
#'
#' @param intensities Table of intensities containing columns:
#' condition, ID, rep, and `log2 intensity`
#' @param pwt_matches Dataframe with columns condition and control that match
#' conditions to their in-batch pwt.
#' @param test_strategy String specifying the strategy to use for the
#' statistical test to produce p-values.
#' @return Data table with L2FC computed as a difference in mean
#' `log2 intensity` in each non-pWT condition and all pWT-conditions
#' And `t-statistic` and `t-test p-value` computed from a 2-sample
#' independent (Welch) two-sample t-test where the non-pWT conditions'
#' mean, SD, and DoF are estimated the usual way, and the pWT conditions'
#' mean and SD are estimated from the full pool of all pWT considions and
#' the DoF is set to the number of distinct rep values in the pool of all
#' pWT conditions.
compute_protein_l2fc <- function (
    intensities, pwt_matches, test_strategy
) {

    wt_intensities <- switch(test_strategy,
        "reduced dof" = intensities %>%
            filter(
                condition %in%
                (pwt_matches %>% distinct(control) %>% pull())
            ) %>%
            group_by(ID) %>%
            summarize(
                `PWT effective DoF` = n_distinct(rep),
                `PWT mean of log2 intensity` = mean(`log2 intensity`),
                `PWT SD of log2 intensity` = sd(`log2 intensity`)) %>%
            ungroup()
        ,
        "technical avg" = intensities %>%
                filter(
                    condition %in%
                    (pwt_matches %>% distinct(control) %>% pull())
                ) %>%
                group_by(ID, rep) %>%
                summarize(`log2 intensity` = mean(`log2 intensity`)) %>%
                group_by(ID) %>%
                summarize(
                    `PWT effective DoF` = n() - 1,
                    `PWT mean of log2 intensity` = mean(`log2 intensity`),
                    `PWT SD of log2 intensity` = sd(`log2 intensity`)) %>%
                ungroup()
        ,
        "pooled sd" = intensities %>%
            filter(
                condition %in%
                (pwt_matches %>% distinct(control) %>% pull())
            ) %>%
            group_by(ID, condition) %>%
            summarize(
                mean = mean(`log2 intensity`),
                sd = sd(`log2 intensity`),
                df = n() - 1) %>%
            group_by(ID) %>%
            summarize(
                `PWT mean of log2 intensity` = mean(mean),
                `PWT SD of log2 intensity` = sqrt(sum(sd^2 * df) / sum(df)),
                `PWT effective DoF` = sum(df)
            ) %>%
            ungroup()
    )

    condition_intensities <- intensities %>%
        filter(
            condition %in%
            (pwt_matches %>% distinct(condition) %>% pull())
        ) %>%
        group_by(ID, condition) %>%
        summarize(
            `effective DoF` = n() - 1,
            `mean of log2 intensity` = mean(`log2 intensity`),
            `SD of log2 intensity` = sd(`log2 intensity`)) %>%
        ungroup()

    matched_table <- condition_intensities %>% left_join(wt_intensities, by = "ID")

    # Return value
    matched_table %>%
        mutate(
            L2FC = `mean of log2 intensity` - `PWT mean of log2 intensity`,
            `sd/n` = `SD of log2 intensity`^2 / (`effective DoF` + 1),
            `PWT sd/n` = `PWT SD of log2 intensity`^2 / (`PWT effective DoF` + 1),
            `t statistic` = L2FC / sqrt(`sd/n` + `PWT sd/n`),
            `t-test DoF` = (`sd/n` + `PWT sd/n`)^2 /
                (`sd/n`^2 / `effective DoF` + `PWT sd/n`^2 / `PWT effective DoF`),
            `t-test p-value` = 2 * pt(-abs(`t statistic`), `t-test DoF`),
            `t-test type` = "Independent 2-sample t-test with adjusted degrees of freedom"
        ) %>%
        select(-`sd/n`, -`PWT sd/n`, -`t-test DoF`) %>%
        group_by(condition) %>%
        mutate(`q-value (within condition)` = qvalue(`t-test p-value`)$qvalues) %>%
        ungroup()
}


compute_l2fc_per_batch_wt <- function(intensities) {
    # Assumptions:
    # * Molecule ID column is "ID"
    # * PWT conditions are called PWT
    # * batch column exists

    condition_intensities <- intensities %>%
        group_by(ID, condition, batch) %>%
        summarize_intensities() %>%
        ungroup()

    matched_table <- condition_intensities %>%
        filter(condition != "PWT") %>%
        mutate(control = "PWT") %>%
        left_join(
            condition_intensities %>% rename_all(~ paste("PWT", .)),
            by = c(
                "ID" = "PWT ID",
                "control" = "PWT condition",
                "batch" = "PWT batch")) %>%
        mutate(control = paste(control, "batch", batch)) %>%
        select(-batch)

    return(compute_l2fc_from_(matched_table))
}


graceful_t_test <- function(x, y) {
    if (length(x) < 2 || length(y) < 2) {
       return(t.test(x - y))
    } else {
       return(t.test(x, y))
    }
}

graceful_q_values <- function(pvalues) {
    tryCatch(
        qvalue(pvalues)$qvalues,
        error = function(e) qvalue(pvalues, pi0 = 1)$qvalues
    )
}

compute_l2fc_from_ <- function(df)
    df %>%
        mutate(
            `t test` = map2(
                `log2 intensities`,
                `PWT log2 intensities`,
                graceful_t_test),
            `t-test p-value` = map_dbl(`t test`, ~ .$p.value),
            `t-statistic` = map_dbl(`t test`, ~ .$statistic),
            `t-test type` = map_chr(`t test`, ~ .$method),
            L2FC = `mean of log2 intensity` - `PWT mean of log2 intensity`) %>%
        select(
            -`log2 intensities`,
            -`PWT log2 intensities`,
            -`t test`) %>%
        group_by(condition) %>%
        mutate(`q-value (within condition)` =
            graceful_q_values(`t-test p-value`)) %>%
        ungroup()


#' Reads the clean protein metadata file into a variable
#'
#' @param path Path to the protein metadata file
#' @return Data frame
read_protein_metadata <- function(path) {
    readr::read_csv(path, col_types = readr::cols(`MICOS/MIB` = readr::col_character()))
}

#' Build metadata table with rows for all molecules
#' 
#' @param molecule_df Dataframe with all molecule IDs and types
#' @param protein_metadata Dataframe with additional protein metadata
#' @return Dataframe with metadata linked to every molecule
build_metadata <- function(molecule_df, protein_metadata) {
    molecule_df %>%
        distinct(ID, `molecule type`) %>%
        left_join(protein_metadata, by = c("ID" = "Molecule ID"))
}


write_excel_tables <- function(long_multiomics, molecule_metadata, out_file) {

    molecule_types <- molecule_metadata %>% distinct(`molecule type`) %>% pull()

    table_names <- c(
        "L2FC",
        "t-test p-value",
        "q-value in-cond",
        "q-value multi-o",
        "q-value single-o")
    table_values <- c(
        "L2FC",
        "t-test p-value",
        "q-value (within condition)",
        "q-value (multi-omic)",
        "q-value (single-omic)")
    #    "BH-adjusted within condition p-value",
    #    "BH-adjusted multi-omic p-value",
    #    "BH-adjusted single-omic p-value")

    nested_list <- map(molecule_types, function(m_type) {
        meta_df <- molecule_metadata %>%
            filter(`molecule type` == m_type) %>%
            select_if(~ length(unique(.)) > 1)
        df <- long_multiomics %>%
            filter(`molecule type` == m_type) %>%
            arrange(`WT condition`)

        tvs <- setNames(table_values, paste(m_type, table_names))

        map(tvs, function(val_col) {
            if(str_detect(val_col, "adjusted")){
                my_df <- df %>% filter(!`WT condition`)
            } else {
                my_df <- df
            }
            right_join(meta_df,
                pivot_wider(
                    my_df,
                    id_cols = ID,
                    names_from = `cell line`,
                    values_from = one_of(val_col)
                ),
                by = "ID")
        })
    })

    table_list <- flatten(nested_list)

    writexl::write_xlsx(table_list, path = out_file)
    return(TRUE)
}


#' Derive the q-adjusted relative difference matrix tuple
#' 
#' The q-adjusted relative difference is:
#' (1 - q-value) * (KO - WT) / (KO + WT)
#' This value has the property of being bounded between -1 and 1, tending towards 0 when
#' wither the q-value is large or the KO is close to the WT
#' It is computed from the L2FC as:
#' (1 - q-value) * (2 / (1 + 2^-L2FC) - 1)
#' 
#' @return A named list with member matrix (the matrix of values),
#' and one member for each molecule type (indices of those molecules)
derive_mrdm_tuple <- function(long_multiomics) {

    df <- long_multiomics %>%
        filter(!`WT condition`) %>%
        mutate(MRD = (1 - `q-value (multi-omic)`) * (2 / (1 + 2^L2FC) - 1)) %>%
        pivot_wider(
            id_cols = c(ID, `molecule type`),
            names_from = `cell line`,
            values_from = MRD)

    m <- df %>%
        select(-`molecule type`) %>%
        column_to_rownames("ID") %>%
        as.matrix()

    mtypes <- c("Protein", "Metabolite", "Lipid")
    names(mtypes) <- mtypes

    df %<>% mutate(i = 1:n())

    result <- map(mtypes,
        function(mtype) {
            df %>% filter(`molecule type` == mtype) %>% pull(i)
        }
    )

    result$matrix <- m

    return(result)
}


#' Prepare molecule metadata for plots
#' 
#' @param molecule_metadata contains molecule type information
#' @param protein_metadata contains additional protein metadata
#' @return Metadata table for plot
metadata_20201007 <- function(molecule_metadata, protein_metadata) {
    molecule_metadata %>%
        select(ID, `molecule type`) %>%
        left_join(protein_metadata %>% rename(ID = `Molecule ID`), by="ID") %>%
        mutate(
            mt = `molecule type`,
            `molecule type` = if_else(`molecule type` == "Protein",
                if_else(!is.na(`OxPhos complex`), "Protein (OxPhos)",
                    if_else(!is.na(`Mitoribosome subunit`), "Protein (Mitoribosome)",
                        if_else(`MitoCarta2.0`, "Protein (Mitochondria)", "Protein")
                    )
                ),
                `molecule type`
            ),
            searchstr = paste(`ID`, `Protein groups`, `HGNC symbols`),
            tt1 = if_else(
                    `MitoCarta2.0`,
                    paste0("In MitoCarta2.0",
                        if_else(is.na(`OxPhos complex`),
                            if_else(!is.na(`Mitoribosome subunit`),
                                paste("<br />Mitoribosome subunit", `Mitoribosome subunit`),
                                ""
                            ),
                            paste("<br />OxPhos", `OxPhos complex`)
                        )
                    ),
                as.character(NA)),
            tooltip = if_else(
                mt != "Protein",
                paste(ID, mt, sep = "<br />"),
                paste(
                    `Protein groups`,
                    if_else(is.na(tt1),
                        `HGNC symbols`,
                        paste(`HGNC symbols`, tt1, sep = "<br />")),
                    sep = "<br />"
                )
            )
        ) %>%
        select(ID, `molecule type`, searchstr, tooltip)
}


#' Plot molecule embedding
#' 
#' @param x x-coordinate
#' @param y y-coordinate
#' @param id ID array
#' @param metadata contains metadata to join to table by the "ID" column
#' @param colors array of colors to use
#' @param html_template path to the html template file
#' @param html_out html output file
plot_embedding2 <- function(x, y, id, metadata, colors, html_template, html_out){

    data_json <- embedding_to_json(x, y, id, metadata, colors)

    read_file(html_template) %>%
        str_replace(fixed("/*DATA*/"), data_json) %>%
        write_file(html_out)
}


#' Plot molecule embedding
#' 
#' @param x x-coordinate
#' @param y y-coordinate
#' @param id ID array
#' @param metadata contains metadata to join to table by the "ID" column
#' @param colors a named array of colors to use
embedding_to_json <- function(x, y, id, metadata, colors){

    color_df <- if(is.null(names(colors))){
        metadata %>%
            distinct(`molecule type`) %>%
            arrange() %>%
            mutate(color = colors)
    } else {
        tibble(`molecule type` = names(colors), color = colors)
    }

    df <- tibble(
        ID = id,
        x = x,
        y = y) %>%
        left_join(metadata, by = "ID")

    result <- df %>%
        group_by(`molecule type`) %>%
        group_map(function(g_df, g_key){
            g_key %<>% left_join(color_df, by = "molecule type")
            return (list(
                molecule_type = g_key %>% pull(`molecule type`) %>% jsonlite::unbox(),
                color = g_key %>% pull(`color`) %>% jsonlite::unbox(),
                data = g_df
            ))
        })

    jsonlite::toJSON(result, dataframe = "columns")
}


#' Find embedding neighbors of MXP molecules
#' 
nearest_embedding_neighbors <- function(
    radius,
    embedding,
    id,
    protein_metadata,
    molecule_types,
    mxp_genes,
    out_file
){

    distances <- as.matrix(dist(embedding))

    mxp_hgnc <- mxp_genes %>% filter(MXP == "MXP") %>% pull(`HGNC Symbol`)

    molecules <- tibble(ID = id) %>%
        mutate(i = 1:n())

    mxp_molecules <- molecules %>%
        left_join(protein_metadata %>% select(ID = `Molecule ID`, `HGNC symbols`), by="ID") %>%
        filter(map_lgl(`HGNC symbols`, function(sym){
            any(mxp_hgnc %in% str_split(sym, ";")[[1]])
        }))

    result_df <- map_dfr(mxp_molecules %>% pull(i),
    function(index){
        mxp_molecule_row = mxp_molecules %>% filter(i==index)
        molecules %>%
            mutate(
                `MXP Molecule` = mxp_molecule_row$ID,
                `MXP Molecule HGNC Symbols` = mxp_molecule_row$`HGNC symbols`,
                `Embedding Distance` = distances[index,]) %>%
            filter(`Embedding Distance` <= radius, i != index) %>%
            select(
                `MXP Molecule`,
                `MXP Molecule HGNC Symbols`,
                `Neighbor Molecule` = ID,
                `Embedding Distance`)
    })

    result_df %<>%
        left_join(molecule_types, by = c("Neighbor Molecule" = "ID")) %>%
        rename(`Neighbor Molecule Type` = `molecule type`)

    result_df %>% write_csv(out_file)

    return(result_df)
}