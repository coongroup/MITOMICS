# Function definitions
#
# Licensed under the BSD 3-Clause License
# Copyright (c) 2019, Yuriy Sverchkov

library(magrittr)
library(drake)
library(tidyverse)
library(ggplot2)
library(patchwork)
#library(tidygraph)
#library(visNetwork)
#library(plotly)
library(qvalue)

# Data loading for 2021-08-19 data

read_lipids <- function(file_path, sheet) {

    readxl::read_excel(file_path, sheet = sheet) %>%
        rename(ID = `Molecule.ID`) %>%
        make_data_object()
}

read_metabolites <- function(file_path, sheet) {

    readxl::read_excel(file_path, sheet = sheet) %>%
        make_data_object()
}

read_proteins <- function(file_path, sheet) {

    readxl::read_excel(file_path, sheet = sheet, guess_max=1E7) %>%
        rename(ID = `Majority protein ID`) %>%
        make_data_object()
}

make_data_object <- function(wide_df) {
    list(
        metadata = wide_df %>% select(!where(is.numeric)),
        data = wide_df %>%
            select(ID, where(is.numeric)) %>%
            pivot_longer(
                where(is.numeric),
                names_to = c("cell line", ".value"),
                names_pattern = "(\\S+) (.+)")
    )
}

# If we're computing the p-values
process_protein_intensities <- function(file_path) {
    proteins_csv <- read_csv(
        file_path,
        col_types =
            readr::cols(
                `MICOS/MIB` = readr::col_character(),
                `mtDNA encoded` = readr::col_character()
            )
        ) %>%
        rename(ID = `Majority protein ID`)
    
    metadata <- proteins_csv %>% select(!where(is.numeric))

    protein_intensities <- proteins_csv %>%
        select(ID, where(is.numeric)) %>%
        pivot_longer(
            where(is.numeric),
            names_to = c("cell line", "replicate"),
            names_pattern = "^(.+)-([A-Z])$"
        )

    wt_intensities <- protein_intensities %>%
        filter(`cell line` == "WT") %>%
        group_by(ID) %>%
        summarize(wt_values = list(value), .groups = "drop")

    fold_changes <- protein_intensities %>%
        group_by(ID, `cell line`) %>%
        summarize(values = list(value), .groups = "drop") %>%
        left_join(wt_intensities, by="ID") %>%
        mutate(
            ttest = purrr::map2(values, wt_values, t.test),
            `P Value` = purrr::map_dbl(ttest, ~ .$p.value),
            `Control Normalized Log2 Fold Change` = purrr::map_dbl(ttest, ~ (.$estimate[1] - .$estimate[2]))
        ) %>%
        select(ID, `cell line`, `Control Normalized Log2 Fold Change`, `P Value`)
    
    list(
        metadata = metadata,
        data = fold_changes
    )
}


#' Derive the q-adjusted relative difference matrix
#' 
#' The q-adjusted relative difference is:
#' (1 - q-value) * (KO - WT) / (KO + WT)
#' This value has the property of being bounded between -1 and 1, tending towards 0 when
#' wither the q-value is large or the KO is close to the WT
#' It is computed from the L2FC as:
#' (1 - q-value) * (2 / (1 + 2^-L2FC) - 1)
#' 
#' @return A matrix of the QARD values where rows are molecules and columns are cell lines.
derive_QARD_matrix <- function(long_multiomics) {

    long_multiomics %>%
        mutate(QARD = (1 - `Q Value`) * (2 / (1 + 2^(-`Control Normalized Log2 Fold Change`)) - 1)) %>%
        pivot_wider(
            id_cols = c(ID),
            names_from = `cell line`,
            values_from = QARD) %>%
        column_to_rownames("ID") %>%
        as.matrix()
}


##################
# Embedding plot #
##################

#' Create metadata for plot
#'
#' @param lipid_metadata
#' @param metabolite_metadata
#' @param protein_metadata
#' @return Data frame with columns ID, molecule type, searchstr, tooltip
metadata_for_plot <- function(lipid_metadata, metabolite_metadata, protein_metadata) {

    br = "<br />"

    bind_rows(
        lipid_metadata %>%
            mutate(
                `molecule type` = "Lipid",
                searchstr = ID,
                tooltip = paste(ID, "Lipid", sep = br)
            ),
        metabolite_metadata %>%
            mutate(
                `molecule type` = "Metabolite",
                searchstr = paste(ID, HMDB),
                tooltip = paste(ID, HMDB, "Metabolite", sep = br)
            ),
        protein_metadata %>%
            mutate(
                `molecule type` = if_else(!is.na(`OxPhos complex`), "Protein (OxPhos)",
                    if_else(!is.na(`Mitoribosome subunit`), "Protein (Mitoribosome)",
                        if_else(!is.na(`MitoCarta2.0`), "Protein (Mitochondria)", "Protein")
                    )
                ),
                searchstr = paste(`HGNC symbols`, `Protein groups`),
                tt1 = if_else(
                    !is.na(`MitoCarta2.0`),
                    paste0("In MitoCarta2.0",
                        if_else(is.na(`OxPhos complex`),
                            if_else(!is.na(`Mitoribosome subunit`),
                                paste0(br, "Mitoribosome subunit ", `Mitoribosome subunit`),
                                ""
                            ),
                            paste0(br, "OxPhos ", `OxPhos complex`)
                        )
                    ),
                    as.character(NA)
                ),
                tooltip = paste(
                    if_else(is.na(tt1),
                        `HGNC symbols`,
                        paste(`HGNC symbols`, tt1, sep = br)),
                    `Protein groups`,
                    `molecule type`,
                    sep = br
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

    html_out
}

#' Plot molecule embedding
#' 
#' @param x x-coordinate
#' @param y y-coordinate
#' @param id ID array
#' @param metadata contains metadata to join to table by the "ID" column
#' @param colors a named array of colors to use
embedding_to_json <- function(x, y, id, metadata, colors, cluster = NULL){

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
    
    if (!is.null(cluster)) {
        df$size <- if_else(cluster == 0, 4, 8)
    }
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


#' Find embedding neighbors of MXP molecules consistent across layouts
#' 
nearest_embeddings_neighbors <- function(
    radius,
    embeddings,
    id,
    protein_metadata,
    molecule_types,
    mxp_genes,
    out_file
){

    distances <- do.call(pmax, lapply(embeddings, function(embedding){
        as.matrix(dist(embedding$Y))
    }))

    mxp_hgnc <- mxp_genes %>% filter(MXP == "MXP") %>% pull(`HGNC Symbol`)

    molecules <- tibble(ID = id) %>%
        mutate(i = 1:n())

    mxp_molecules <- molecules %>%
        left_join(protein_metadata %>% select(ID, `HGNC symbols`), by="ID") %>%
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

    return(out_file)
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
        left_join(protein_metadata %>% select(ID, `HGNC symbols`), by="ID") %>%
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


#####################################
# KO-KO similarity characterization #
#####################################

#' Add calls to long-form table
make_calls <- function (long_multiomics, q_cutoff = 0.05, l2fc_cutoff = 1) {
    long_multiomics %>%
        mutate(
            sigL2FC = if_else(`Q Value` < q_cutoff, `Control Normalized Log2 Fold Change`, 0),
            sigbigL2FC = if_else(abs(sigL2FC) > l2fc_cutoff, sigL2FC, 0),
            call = sign(sigbigL2FC))
}

#' Convert long_multiomics (or a subset) to a KOs-by-molecules feature matrix
#'
#' @param long_multiomics is the original long-form data frame 
#' @param values_from is the variable from which to take the values (e.g. L2FC)
#' @param remove_zeros (character) is the policy on which zeros to remove:
#'  'none' - don't remove zeros (default)
#'  'cols' - columns where all values are zero
#'  'rows' - rows where all values are zero
#'  'both' - both rows and columns
long_to_profile_matrix <- function (long_multiomics, values_from, remove_zeros = 'none'){

    x <- long_multiomics %>%
        pivot_wider(
            id_cols = `cell line`,
            names_from = `ID`,
            values_from = {{ values_from }} # "embrace"d variable promise
        ) %>%
        column_to_rownames("cell line") %>%
        as.matrix()
    
    if (remove_zeros %in% c('cols', 'both')) {
        x <- x[,colSums(x != 0) > 0]
    }

    if (remove_zeros %in% c('rows', 'both')) {
        x <- x[rowSums(x != 0) > 0,]
    }

    return (x)
}

read_cell_line_metadata <- function(excel_file) {
    readxl::read_excel(excel_file, sheet = 'H3K Cell Lines', na = c('', 'N/A'))
}

get_pair_properties <- function(cell_line_metadata) {
    df <- cell_line_metadata %>%
        dplyr::select(
            `HGNC Symbol`,
            `H3K Cell Line Name`,
            `Target Exon`,
            `Genomic Location`,
            `Guide RNA Sequence`,
            `Indel Type`,
            `Indel Size`) %>%
        dplyr::rename(target = `HGNC Symbol`, line = `H3K Cell Line Name`) %>%
        filter(!is.na(target)) %>%
        pivot_longer(!c(target, line), names_to = "Property", values_to = "value")
    
    inner_join(df, df, by = c('target', 'Property')) %>%
        filter(line.x < line.y, value.x == value.y) %>%
        dplyr::rename(`Property Match` = Property)
}

#' Make 'sister KO' correlation/similarity rank plot
#' 
#' Sort-by values:
#'  * rank - default - rank of worse match
#'  * distance - rank of correlation between sisters
#'  * best rank - rank of better match
plot_sister_ko_cor_ranks <- function (
    correlation_matrix,
    plot_file,
    full_match_table_file = NULL,
    neat_match_table_file = NULL,
    profile_matrix=NULL,
    sort_by='rank',
    title=NULL,
    pair_properties=NULL){

    correlations_df <- correlation_matrix %>%
        as_tibble(rownames = 'a') %>%
        pivot_longer(!a, names_to = 'b', values_to = 'correlation') %>%
        filter(a != b) %>%
        mutate(
            match = str_remove(a, "-.*") == str_remove(b, "-.*"),
            target = str_remove(a, "-.*"),
            `correlation profile` = str_replace(a, ".*-", ''))
    
    df <- correlations_df %>%
        mutate(distance = abs(correlation - 1)) %>%
        group_by(a) %>%
        arrange(distance, .by_group = TRUE) %>%
        mutate(`rank of sister KO` = 1:n()) %>%
        ungroup() %>%
        filter(match) %>%
        mutate(`rank of` = paste(str_replace(b, '.*-', ''), 'in', str_replace(a, '.*-', '')))

    if(!is.null(full_match_table_file)) df %>% write_csv(full_match_table_file)

    # Impose order over target names
    if (sort_by == 'rank') {
        targets_order <- df %>%
            group_by(target) %>%
            summarize(rank = max(`rank of sister KO`)) %>%
            ungroup() %>%
            arrange(desc(rank)) %>%
            pull(target)
    } else if (sort_by == 'best rank') {
        targets_order <- df %>%
            group_by(target) %>%
            summarize(rank = min(`rank of sister KO`)) %>%
            ungroup() %>%
            arrange(desc(rank)) %>%
            pull(target)
    } else if (sort_by == 'distance') {
        targets_order <- df %>%
            group_by(target) %>%
            summarize(rank = dplyr::first(distance)) %>%
            ungroup() %>%
            arrange(desc(rank)) %>%
            pull(target)
    }

    correlations_df %<>%
        filter(target %in% targets_order) %>%
        mutate(target = factor(target, levels = targets_order))

    df$target <- factor(df$target, levels = targets_order)

    if(!is.null(neat_match_table_file)) df %>%
        arrange(desc(target), `rank of sister KO`) %>%
        select(target, correlation, `rank of sister KO`, `rank in` = a, `rank of` = b) %>%
        write_csv(neat_match_table_file)

    match_distance_df <- df %>% group_by(target) %>% summarize(correlation = dplyr::first(correlation))

    # ggplot
    g1 <- ggplot(data = df, aes(x = `rank of sister KO`, y = target, shape = `rank of`)) +
        geom_point() +
        scale_x_continuous(sec.axis = dup_axis())+
        theme_light()
    
    g2 <- ggplot(data = correlations_df, aes(x = correlation, y = target, fill = `correlation profile`)) +
        geom_col(data = match_distance_df, fill = 'grey', alpha = 0.7) +
        geom_violin(size=0.1) +
        scale_x_continuous(sec.axis = dup_axis())+
        theme_light()

    g <- g1 | g2
    # Relative widths of subplots
    plot_widths = c(1,1)

    # If we have a profile matrix also plot raw values
    if (!is.null(profile_matrix)) {
        profile_df <- profile_matrix %>%
            as_tibble(rownames = 'ko') %>%
            pivot_longer(!ko, names_to = 'molecule', values_to = 'value') %>%
            mutate(
                target = str_remove(ko, "-.*"),
                `call profile` = str_replace(ko, ".*-", '')) %>%
            filter(target %in% targets_order) %>%
            mutate(target = factor(target, levels = targets_order))
        
        g3 <- ggplot(data = profile_df, aes(y=target, fill=`call profile`)) +
            geom_bar(position = 'dodge2') +
            facet_grid(
                cols = vars(sign(value)),
                scales = 'free_x',
                switch = 'x') +
            scale_x_continuous(
                guide = guide_axis(check.overlap = TRUE),
                sec.axis = dup_axis())+
            theme_light()
        
        g <- g | g3
        plot_widths <- c(plot_widths, 1)
    }

    # If we have pair properties also plot those
    if (!is.null(pair_properties)) {
        pair_properties %<>%
            filter(target %in% targets_order) %>%
            mutate(target = factor(target, levels = targets_order))

        g4 <- ggplot(data = pair_properties, aes(y=target, x=`Property Match`)) +
            geom_point() +
            scale_y_discrete(drop = FALSE) +
            #scale_x_discrete(position = 'top') +
            theme_light() +
            theme(axis.text.x = element_text(angle=90))
        
        g <- g | g4
        plot_widths <- c(plot_widths, 0.6)
    }
    
    g <- g +
        plot_layout(widths = plot_widths, guides = 'collect') &
        theme(axis.title.y = element_blank())

    if (!is.null(title)) g <- g + plot_annotation(title = title)

    ggplot2::ggsave(plot_file, g, width=3 + 3 * sum(plot_widths), height=16)

    plot_file
}


#' Plot a 1-d distribution with annotations
#'
plot_distribution <- function(values, annotate_by_prop, annotate_by_val, plot_file){

    g <- ggplot(mapping = aes(x=values)) + geom_density(adjust = 0.25)

    sorted_values <- sort(values)
    n <- length(sorted_values)

    text_df <- bind_rows(
        tibble(prop = annotate_by_prop) %>%
            mutate(
                index = ceiling(n * prop),
                x = sorted_values[index]),
        tibble(x = annotate_by_val, prop = sapply(x, function(v) mean(sorted_values < v)))
        ) %>%
        mutate(
            label = paste0(
                format(100 * prop, digits = 3),
                "% of values\n< ",
                format(x, digits = 4))
        ) %>%
        arrange(desc(x)) %>%
        mutate(y = 0.5 + 0.4 * 1:n())
    
    g <- g +
        geom_vline(aes(xintercept=x), data=text_df) +
        geom_label(aes(x=x, y=y, label=label), data=text_df, hjust= "right")

    ggplot2::ggsave(plot_file, g)
}

################
# Cluster plot #
################

#' Plot molecule embedding with hdbscan overlay
#' 
#' @param x x-coordinate
#' @param y y-coordinate
#' @param id ID array
#' @param metadata contains metadata to join to table by the "ID" column
#' @param colors array of colors to use
#' @param html_template path to the html template file
#' @param html_out html output file
plot_embedding_clustered <- function(x, y, id, metadata, colors, clusters, html_template, html_out){

    contour_json <- mk_cluster_shapes_json(x, y, clusters)
    data_json <- embedding_to_json(x, y, id, metadata, colors, cluster = clusters)

    read_file(html_template) %>%
        str_replace(fixed("/*SHAPES*/"), contour_json) %>%
        str_replace(fixed("/*DATA*/"), data_json) %>%
        write_file(html_out)
}


#' Plot molecule embedding with hdbscan overlay
#' 
#' @param x x-coordinate
#' @param y y-coordinate
#' @param id ID array
#' @param metadata contains metadata to join to table by the "ID" column
#' @param colors array of colors to use
#' @param html_template path to the html template file
#' @param html_out html output file
plot_embedding_clustered2 <- function(x, y, id, metadata, colors, core_clusters, extended_clusters, html_template, html_out){

    clusters_df <- tibble(ID = id, core_cluster = core_clusters, extended_cluster = extended_clusters) %>%
        mutate(cluster_text = if_else(
            core_cluster != 0,
            paste("Core cluster", core_cluster),
            paste("Extended cluster", extended_clusters)))

    metadata <- metadata %>%
        left_join(clusters_df, by="ID") %>%
        mutate(tooltip = paste(tooltip, cluster_text, sep="<br />"))

    contour_json <- mk_cluster_ashapes_json(x, y, core_clusters, extended_clusters)

    data_json <- embedding_to_json(x, y, id, metadata, colors)#, cluster = core_clusters)

    read_file(html_template) %>%
        str_replace(fixed("/*SHAPES*/"), contour_json) %>%
        str_replace(fixed("/*DATA*/"), data_json) %>%
        write_file(html_out)

    html_out
}


#' Make Plotly shapes to indicate clusters
#'
#' Returns a specification of the convex hull
#'
#' @param x x-coordinate vector
#' @param y y-coordinate vector
#' @param clusters cluster numbers (0 is noise)
#' @param alpha alpha value for alpha shapes calculation
mk_cluster_shapes_json <- function(x, y, clusters){
    shapes <- list()
    for (c in unique(clusters)) if (c != 0) {
        cx <- x[clusters == c]
        cy <- y[clusters == c]
        inds <- chull(cx, cy)
        path_str <- paste(
            'M',
            paste(cx[inds], cy[inds], collapse = ' L '),
            'Z'
        )
        shapes[[length(shapes)+1]] <- list(
            type = 'path',
            path = path_str,
            #fillcolor = 'rgba(128, 128, 128, 0.2)',
            line = list(color = 'rgba(210, 128, 0, 0.6)')
        )
    }
    
    jsonlite::toJSON(shapes, auto_unbox = TRUE)
}


#' Make Plotly shapes to indicate clusters
#'
#' Returns a specification of the convex hull
#'
#' @param x x-coordinate vector
#' @param y y-coordinate vector
#' @param core_clusters cluster numbers (0 is noise)
#' @param extended_clusters cluster numbers (0 is noise)
#' @param alpha alpha value for alpha shapes calculation
mk_cluster_shapes2_json <- function(x, y, core_clusters, extended_clusters){
    
    makeshapes <- function(clusters, linecolor, fillcolor = NULL){
        shapes <- list()

        for (c in unique(clusters)) if (c != 0) {
            cx <- x[clusters == c]
            cy <- y[clusters == c]
            inds <- chull(cx, cy)
            path_str <- paste(
                'M',
                paste(cx[inds], cy[inds], collapse = ' L '),
                'Z'
            )
            shapes[[length(shapes)+1]] <- list(
                type = 'path',
                path = path_str,
                line = list(color = linecolor),
                layer = 'below'
            )
            if(!is.null(fillcolor)) shapes[[length(shapes)]]$fillcolor <- fillcolor
        }

        shapes
    }
    
    jsonlite::toJSON(
        c(
            makeshapes(core_clusters, 'rgba(97, 159, 208, 0.9)', 'rgba(97, 159, 208, 0.1)'),
            makeshapes(extended_clusters, 'rgba(0, 0, 0, 0.3)')
        ),
        auto_unbox = TRUE
    )
}


#' Make Plotly shapes to indicate clusters
#'
#' Returns a specification of the convex alpha-hulls
#'
#' @param x x-coordinate vector
#' @param y y-coordinate vector
#' @param core_clusters cluster numbers (0 is noise)
#' @param extended_clusters cluster numbers (0 is noise)
#' @param alpha alpha value for alpha shapes calculation
mk_cluster_ashapes_json <- function(x, y, core_clusters, extended_clusters){

    #message("Computing distance matrix")
    #dist_m <- as.matrix(dist(cbind(x,y)))
    
    makeshapes <- function(clusters, linecolor, fillcolor = NULL){
        shapes <- list()

        for (c in unique(clusters)) if (c != 0) {

            cx <- x[clusters == c]
            cy <- y[clusters == c]

            # Alpha = distance to nearest out-of-cluster point guarantees no overlap between clusters
            #alpha <- min(dist_m[clusters == c, clusters != c])

            # Use maximum 6th nearest neighbor?
            ldm <- as.matrix(dist(cbind(cx, cy)))
            #ldm <- dist_m[clusters == c, clusters == c]
            #diag(ldm) <- max(ldm)
            alpha <- max(apply(ldm, 1, function(v) sort(v)[min(6, length(v))]))

            message("Computing alpha=", alpha, " hull for cluster ", c)

            ash <- alphahull::ashape(cx, cy, alpha=alpha)

            # Get outline, assuming edges are always in sequence although sometimes they have reversed
            # point order?
            edges <- ash$edges[, c(1,2)]
            inds <- edges[1,]
            edges <- edges[-1, ]
            while(!is.null(dim(edges)) && nrow(edges) > 1){
                last <- inds[length(inds)]
                # Find
                if (last %in% edges[,1]){
                    j <- which(edges[,1] == last)
                    inds <- c(inds, edges[j,2])
                    edges <- edges[-j,]
                }else{
                    j <- which(edges[,2] == last)
                    inds <- c(inds, edges[j,1])
                    edges <- edges[-j,]
                }
            }

            #for(i in 2:(nrow(ash$edges)-1)){
            #    if (ash$edges$ind1 == inds[length(inds))
            #        inds <- c(inds, ash$edges$ind2)
            #    else
            #        inds <- c(inds, ash$edges$ind1)
            #}

            #inds <- ash$alpha.extremes
            path_str <- paste(
                'M',
                paste(cx[inds], cy[inds], collapse = ' L '),
                'Z'
            )
            shapes[[length(shapes)+1]] <- list(
                type = 'path',
                path = path_str,
                line = list(color = linecolor),
                layer = 'below'
            )
            if(!is.null(fillcolor)) shapes[[length(shapes)]]$fillcolor <- fillcolor
        }

        shapes
    }
    
    jsonlite::toJSON(
        c(
            makeshapes(core_clusters, 'rgba(97, 159, 208, 0.9)', 'rgba(97, 159, 208, 0.1)'),
            makeshapes(extended_clusters, 'rgba(0, 0, 0, 0.3)')
        ),
        auto_unbox = TRUE
    )
}


#####################
# Outlier detection #
#####################


#' Function for identifying outliers using a modified Y3K uGPS algorithm
outlier_distances <- function(long_multiomics, L2FC, p_var, p_cutoff = 0.05, k = 3) {

    # For each molecule
    # For +ve and -ve direction
    # Consider log-2fc vs log-10 p-value
    # Normalize so that the largest is 1 in either dimension
    # Consider points with 2*k largest fold-changes at p < 0.05
    # For each point consider distances to all other points for which
    # lfc is smaller and the ko condition is targeting a different gene
    # Consider the smallest distance for each (of top 2k) points
    # Return top k

    # Return
    long_multiomics %>%
        # Split up data into upregulated and downregulated
        mutate(regulation = if_else({{ L2FC }} > 0, "up", "down")) %>%
        # For each molecule and regulation mode
        group_by(ID, regulation) %>%
        group_modify(function(df, key) {
            # Create table of normalized volcano plot coordinates
            df %<>%
                transmute(
                    `cell line`,
                    fc = if(key$regulation == "up")
                        {{ L2FC }} / max({{ L2FC }}) else
                        {{ L2FC }} / min({{ L2FC }}),
                    significant = {{ p_var }} < p_cutoff,
                    l10pv = -log10({{ p_var }}),
                    pv = l10pv / max(l10pv),
                    target = str_remove(`cell line`, "-.*")
                )

            # Pull out top significant conditions, up to 3 distinct targets
            candidates <- df %>%
                filter(significant) %>%
                arrange(desc(fc))

            if (nrow(candidates) < 1) {
                return(tibble(`cell line` = character(), distance = numeric()))
            }

            if (nrow(candidates) > 2*k) {
                n <- k*2;
                while (n_distinct(candidates$target[1:n]) < k*2) n <- n + 1
                candidates <- candidates[1:n, ]
            }

            # For each candidate
            candidates %>%
                group_by(`cell line`) %>%
                group_modify(function(candidate, key){
                    df %>%
                        filter(
                            fc < candidate$fc, # Less extreme change
                            target != candidate$target # Different targets
                        ) %>%
                        summarize(distance = min(sqrt(
                            (candidate$fc - fc)^2
                            + (candidate$pv - pv)^2)))
                }) %>%
                ungroup() %>%
                filter(!is.na(`cell line`)) %>% # NAs introduced if there were empty groups above
                slice_max(distance, n = k)
        }) %>%
        ungroup() %>%
        #filter(!is.infinite(distance)) %>%
        arrange(desc(distance))
}


##############
# Tuned tSNE #
##############

compute_tsne_layout <- function(
    x_mat,
    n_tsne_runs = 10,
    n_pca_perm = 50,
    theta = 0.5,
    pc_size = NULL,
    exaggeration_factor = 4.0,
    seed = NULL,
    scale = 'none'
){

    if(!is.null(seed)){
        set.seed(seed)
    }

    # Scale
    if (scale == 'tsne_pk'){
        # Scale as in tsne package
        x_mat <- (x_mat - min(x_mat)) / (max(x_mat) - min(x_mat))
        x_mat <- scale(x_mat, center = TRUE, scale = FALSE)
    } else if (scale == 'znorm'){
        # Use z-normalization
        x_mat <- scale(x_mat, center = TRUE, scale = TRUE)
    }

    pval = NA
    if(is.null(pc_size)) {
        # Find optimal PCA dimension
        message("Finding optimal PCA dimension...")
        # Get full PCA
        pcs <- prcomp(x_mat, center = FALSE, scale = FALSE)
        # Get explained variation per component
        expl_var <- pcs$sdev^2 / sum(pcs$sdev^2)

        message("Running ", n_pca_perm, " permutations...")
        # Get permuted PCA explained variations
        perm_vars <- sapply(1:n_pca_perm, function(i){
            message("Permutation ", i)
            x_perm <- t(apply(t(x_mat), 2, sample))
            perm_pcs <- prcomp(x_perm, center = FALSE, scale = FALSE)
            perm_pcs$sdev^2 / sum(perm_pcs$sdev^2)
        })

        # Compute p-value per component
        pval <- apply(perm_vars >= expl_var, 1, sum) / n_pca_perm

        # Last significant p-val gives optimum pca dimension
        pc_size <- tail(which(pval < 0.05), 1)
        message("Optimal PCA dimension is ", pc_size)
    }

    message('Whitening to ', pc_size, ' dimensions.')
    # Do PCA as in tsne package
    n <- nrow(x_mat)
    p <- ncol(x_mat)
    n.comp = pc_size

    X <- t(x_mat)
    V <- X %*% t(X)/n
    s <- La.svd(V)
    D <- diag(c(1/sqrt(s$d)))
    K <- D %*% t(s$u)
    K <- matrix(K[1:n.comp, ], n.comp, p)
    X = t(K %*% X)

    # set perplexity = sqrt(nrows)
    perplexity <- sqrt(nrow(x_mat))

    print(paste("Perplexity is", perplexity))

    # run n_tsne_runs tsnes
    tsnes <- lapply(1:n_tsne_runs, function(i){
        # run optimized tsne
        tsne_1 <- Rtsne::Rtsne(
            X,
            dims = 2,
            initial_dims = pc_size,
            perplexity = perplexity,
            theta = theta,
            check_duplicates = FALSE,
            pca = FALSE,
            max_iter = 2000,
            verbose = TRUE,
            #pca_center = FALSE,
            #pca_scale = FALSE,
            #normalize = FALSE,
            stop_lying_iter = 250,
            momentum = 0.5,
            final_momentum = 0.8,
            eta = 500,
            exaggeration_factor = exaggeration_factor
        )
    })

    # pick best scoring layout
    best_tsne <- tsnes[[which.min(sapply(tsnes, function(tsne_obj){
        sum(tsne_obj$costs)
    }))]]

    # Return
    list(
        best = best_tsne,
        tsnes = tsnes,
        perplexity = perplexity,
        pc_size = pc_size,
        pc_pval = pval,
        seed = seed
    )
}

compute_l1_tsne_layout <- function(
    x_mat,
    n_tsne_runs = 10,
    seed = NULL
){

    if(!is.null(seed)){
        set.seed(seed)
    }

    # set perplexity = sqrt(nrows)
    perplexity <- sqrt(nrow(x_mat))

    message("Perplexity is ", perplexity)

    # Get distance matrix
    x_dist <- dist(x_mat, method="manhattan")

    # run n_tsne_runs tsnes
    tsnes <- lapply(1:n_tsne_runs, function(i){
        # run optimized tsne
        Rtsne::Rtsne(
            x_dist,
            dims = 2,
            #initial_dims = pc_size,
            perplexity = perplexity,
            theta = 0.5,
            check_duplicates = FALSE,
            pca = FALSE,
            max_iter = 2000,
            verbose = TRUE,
            #pca_center = TRUE,
            is_distance = TRUE,
            stop_lying_iter = 250,
            momentum = 0.5,
            final_momentum = 0.8,
            eta = 200, #500,
            exaggeration_factor = 12.0
        )
    })

    # pick best scoring layout
    best_tsne <- tsnes[[which.min(sapply(tsnes, function(tsne_obj){
        sum(tsne_obj$costs)
    }))]]

    # Return
    list(
        best = best_tsne,
        tsnes = tsnes,
        perplexity = perplexity,
        x_dist = x_dist,
        seed = seed
    )
}

compute_dist_tsne_layout <- function(
    x_dist,
    n_tsne_runs = 10,
    seed = NULL
){

    if(!is.null(seed)){
        set.seed(seed)
    }

    n <- if(class(x_dist) == 'dist') attr(x_dist, 'Size') else nrow(x_dist)

    # set perplexity = sqrt(n)
    perplexity <- sqrt(n)

    message("Perplexity is ", perplexity)

    # run n_tsne_runs tsnes
    tsnes <- lapply(1:n_tsne_runs, function(i){
        # run optimized tsne
        Rtsne::Rtsne(
            x_dist,
            dims = 2,
            #initial_dims = pc_size,
            perplexity = perplexity,
            theta = 0.5,
            check_duplicates = FALSE,
            pca = FALSE,
            max_iter = 2000,
            verbose = TRUE,
            #pca_center = TRUE,
            is_distance = TRUE,
            stop_lying_iter = 250,
            momentum = 0.5,
            final_momentum = 0.8,
            eta = 200, #500,
            exaggeration_factor = 12.0
        )
    })

    # pick best scoring layout
    best_tsne <- tsnes[[which.min(sapply(tsnes, function(tsne_obj){
        sum(tsne_obj$costs)
    }))]]

    # Return
    list(
        best = best_tsne,
        tsnes = tsnes,
        perplexity = perplexity,
        x_dist = x_dist,
        seed = seed
    )
}

####################
# Distance metrics #
####################

cosine_distance <- function(x, transform = 'acos', return_type = 'matrix'){

    # Compute cosines
    dot_product <- x %*% t(x)
    squares <- diag(dot_product)
    cosims <- dot_product / sqrt(squares %*% t(squares))

    # Clamp values
    cosims[cosims < -1] <- -1
    cosims[cosims > 1] <- 1

    if (transform == 'acos')
        dist_mat <- acos(cosims)
    else if (transform == 'logarithmic')
        dist_mat <- 1 - log2(cosims + 1)

    if(return_type == 'dist'){
        return(as.dist(dist_mat))
    }

    return(dist_mat)
}

#################
# Tuned HDBSCAN #
#################

compute_hdbscan <- function(x_mat, minPts_start = 3, minPts_end = 20, step = 1) {

    print("Searching for best minPts parameter for HDBSCAN")

    searching <- TRUE
    minPts <- 3
    clusterings <- NULL

    for (minPts in seq(minPts_start, minPts_end, step)) {

        print(paste("Trying minPts =", minPts))

        clustered <- dbscan::hdbscan(x_mat, minPts)
        clusterings[[length(clusterings) + 1]] <- clustered

        n_clusters <- length(unique(clustered$cluster))
        n_noise <- sum(clustered$cluster == 0)

        print(paste(
            n_clusters,
            "clusters,",
            n_noise,
            "noise points"
        ))

        if (length(unique(clustered$cluster)) < 3)
            break
    }

    # Return
    list(
        best_clustering =
            clusterings[[head(which.min(sapply(
                clusterings[1:(length(clusterings)-1)],
                function(obj){
                    sum(obj$cluster == 0)
                }
            )), 1)]],
        clusterings = clusterings
    )
}

#' Noiseless cluster assignment
#'
#' Gven a cluster assignment with 0's (noise points) assign those points to their nearest clusters
cluster_noise_points <- function(in_x, in_clust) {
    dist_mat <-
        if (class(in_x) == "dist") as.matrix(in_x)
        else as.matrix(dist(in_x))
    non_noise <- which(in_clust != 0)
    in_clust[non_noise[max.col(-dist_mat[, non_noise])]]
}

#' Evaluate cluster stability across tsnes
#' @param tsnes list of Rtsne outputs
#' @param minPts minPts parameter for hdbscan
#' @param n_perms number of permutations for estimating p-value
eval_tsne_cluster_stability <- function(tsnes, minPts, n_perms) {

    message('Computing clusters')
    clusterings <- lapply(tsnes, function(tsne_obj) noiseless_hdbscan(tsne_obj$Y, minPts = minPts))

    message('Computing MI of clusterings')
    mi <- compute_clustering_mi(clusterings)
    message('MI is ', mi)

    perm_mi <- sapply(1:n_perms, function(i) {
        message('Permutation ', i, ' of ', n_perms)
        # Create permuted clustering
        perm_clustering <- lapply(clusterings, sample)
        compute_clustering_mi(perm_clustering)
    })

    pval <- mean(perm_mi >= mi)
    message('p value is ', pval)

    list(
        mi = mi,
        perm_mi = perm_mi,
        pval = pval
    )
}

#' Compute normalized mutual information of a list of clusterings
#'
compute_clustering_mi <- function(clusterings) {

    # Combine clusterings into one df
    c_df <- data.frame(do.call(cbind, clusterings))

    # Compute mutual information

    # Count
    n_pts = nrow(c_df)

    # Get intersections
    intersections <- c_df %>% count(across())

    # Get marginals
    marginals <- lapply(1:ncol(c_df), function(i) {c_df %>% count(across({{ i }}))})

    work_df <- intersections %>% mutate(lg_denominator = 0)
    for (marginal in marginals)
        work_df <- work_df %>%
            left_join(marginal %>% rename(m = n)) %>%
            mutate(lg_denominator = lg_denominator + log2(n_pts / m)) %>%
            select(-m)
    
    # MI:
    mi <- work_df %>%
        mutate(term = (n / n_pts) * (log2(n / n_pts) + lg_denominator)) %>%
        summarize(mi = sum(term)) %>%
        pull(mi)

    # Entropies:
    entropies <- sapply(marginals, function(m_df) {
        m_df %>%
            mutate(term = - n / n_pts * log2(n / n_pts)) %>%
            summarize(entropy = sum(term)) %>%
            pull(entropy)
    })

    # Normalized MI:
    mi / max(entropies)
}
