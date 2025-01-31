library(multiMiR)
library(dplyr)
library(igraph)
library(visNetwork)
library(ReactomePA)
library(DT)  # For interactive tables

# Define miRNAs and fetch targets
miRNAs <- c("hsa-miR-223-3p", "hsa-let-7f-5p", "hsa-miR-1307-3p", "hsa-miR-223-5p", "hsa-miR-221-3p", "hsa-miR-199a-3p")
results_list <- lapply(miRNAs, function(mirna) tryCatch(get_multimir(mirna = mirna, summary = TRUE), error = function(e) NULL))
results_list <- Filter(Negate(is.null), results_list)
all_targets <- do.call(rbind, lapply(results_list, function(x) x@data))  # Assume @data holds the relevant data

# Assume that you have some previous steps correctly defined.

# Combining and deduplicating (be careful with NA values and ensure all columns exist)
if (exists("all_targets")) {
  unique_targets <- distinct(all_targets, mature_mirna_id, target_symbol, .keep_all = TRUE)
} else {
  cat("all_targets does not exist. Ensure that your data fetching and previous processing steps are correct.\n")
}

# Create network data only if unique_targets is correctly formed
if (exists("unique_targets") && "mature_mirna_id" %in% names(unique_targets) && "target_symbol" %in% names(unique_targets)) {
  nodes <- data.frame(
    id = unique(c(unique_targets$mature_mirna_id, unique_targets$target_symbol)),
    label = unique(c(unique_targets$mature_mirna_id, unique_targets$target_symbol)),
    stringsAsFactors = FALSE
  )
  nodes$group <- ifelse(nodes$id %in% unique_targets$mature_mirna_id, "miRNA", "Gene")

  edges <- data.frame(
    from = unique_targets$mature_mirna_id,
    to = unique_targets$target_symbol,
    stringsAsFactors = FALSE
  )

  # Initialize visNetwork
  network_plot <- visNetwork(nodes, edges, width = "100%") %>%
    visNodes(color = list(background = ifelse(nodes$group == "miRNA", "green", "blue"), border = "black")) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visEdges(arrows = "to")

  print(network_plot)
} else {
  cat("Check your unique_targets data frame for correct columns and non-NA entries.\n")
}


# Step 2: Reactome Pathway Analysis
if ("target_symbol" %in% names(unique_targets)) {
    entrez_ids <- bitr(unique_targets$target_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    valid_entrez_ids <- entrez_ids[!is.na(entrez_ids$ENTREZID), "ENTREZID"]

    # Perform enrichment analysis
    reactome_results <- enrichPathway(as.character(valid_entrez_ids), organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "bonferroni")

    # Interactive Table with Enrichment Results
    reactome_table <- datatable(as.data.frame(reactome_results), options = list(pageLength = 5, searchHighlight = TRUE))
}

# Step 3: Interactive Highlighting in Network
observeEvent(input$reactome_table_rows_selected, {
  selected_pathway <- reactome_results[input$reactome_table_rows_selected, ]
  selected_genes <- entrez_ids[entrez_ids$ENTREZID %in% selected_pathway$geneID, "SYMBOL"]

  # Highlight selected genes in the network
  network_plot <- visNetwork(nodes, edges, width = "100%") %>%
    visNodes(color = list(background = ifelse(nodes$id %in% selected_genes, "yellow", ifelse(nodes$group == "miRNA", "green", "blue")), border = "black")) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visEdges(arrows = "to")
  
  output$networkPlot <- renderVisNetwork(network_plot)
})
