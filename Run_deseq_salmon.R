suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

# =========================================================
# EDIT THESE TWO LINES EACH RUN
# =========================================================
dataset_name <- "Carboplatin"   # "Cisplatin" or "Carboplatin"
quant_dir    <- "/home/a1321/Desktop/Bladder/carp/Salmon"  # folder containing *.quant.sf
# =========================================================

# -----------------------------
# 1) Collect Salmon quant.sf files
# -----------------------------
files <- list.files(quant_dir, pattern = "\\.quant\\.sf$", full.names = TRUE)
if (length(files) == 0) stop("No *.quant.sf found in: ", quant_dir)

sample_names <- sub("\\.quant\\.sf$", "", basename(files))
names(files) <- sample_names

condition <- dplyr::case_when(
  stringr::str_detect(sample_names, "^Untreated_") ~ "untreated",
  stringr::str_detect(sample_names, "^Treated_")   ~ "treated",
  TRUE ~ NA_character_
)

if (any(is.na(condition))) {
  bad <- sample_names[is.na(condition)]
  stop("These files do not start with Treated_ or Untreated_:\n",
       paste(bad, collapse = "\n"))
}

sample_table <- data.frame(
  sample = sample_names,
  condition = factor(condition, levels = c("untreated", "treated")),
  stringsAsFactors = FALSE
)
rownames(sample_table) <- sample_table$sample

cat("Dataset:", dataset_name, "\n")
cat("Quant dir:", quant_dir, "\n")
cat("Found", length(files), "quant.sf files:",
    sum(condition == "untreated"), "untreated +",
    sum(condition == "treated"), "treated\n")

# output prefix
out_prefix <- file.path(quant_dir, paste0("DESeq2_gene_", dataset_name, "_Treated_vs_Untreated"))

# -----------------------------
# 2) Build tx2gene (GENCODE v44 headers) — remove transcript/gene versions
# -----------------------------
ex <- readr::read_tsv(files[1], show_col_types = FALSE)

tx2gene <- tibble::tibble(TXFULL = ex$Name) %>%
  dplyr::mutate(
    TXNAME = sub("\\|.*", "", TXFULL),      # ENST....version
    TXNAME = sub("\\..*$", "", TXNAME),     # ENST.... (remove version)
    GENEID = sub("^.*\\|(ENSG[^|]+)\\|.*$", "\\1", TXFULL),  # ENSG....version
    GENEID = sub("\\..*$", "", GENEID)      # ENSG.... (remove version)
  ) %>%
  dplyr::select(TXNAME, GENEID) %>%
  dplyr::distinct()

cat("tx2gene rows:", nrow(tx2gene), "\n")
print(head(tx2gene, 3))

# sanity check
ex_ids <- sub("\\..*$", "", sub("\\|.*", "", ex$Name))
match_n <- sum(ex_ids %in% tx2gene$TXNAME)
cat("Sanity check matches between quant and tx2gene TXNAME =", match_n, "\n")
if (match_n == 0) stop("0 matches. tx2gene parsing failed.")

# -----------------------------
# 3) tximport
# -----------------------------
txi <- tximport::tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreAfterBar = TRUE,
  ignoreTxVersion = TRUE
)

# -----------------------------
# 4) DESeq2
# -----------------------------
dds <- DESeq2::DESeqDataSetFromTximport(
  txi,
  colData = sample_table,
  design = ~ condition
)

dds <- dds[rowSums(DESeq2::counts(dds)) >= 10, ]
dds <- DESeq2::DESeq(dds)

res <- DESeq2::results(dds, contrast = c("condition", "treated", "untreated"))

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::mutate(ENSEMBL = sub("\\..*$", "", gene_id))

# -----------------------------
# 5) Map Ensembl -> SYMBOL / ENTREZID (clean 1-to-1, no select warnings)
# -----------------------------
res_df$SYMBOL <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = res_df$ENSEMBL,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

res_df$ENTREZID <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = res_df$ENSEMBL,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

# -----------------------------
# 6) Save full results
# -----------------------------
readr::write_csv(res_df, paste0(out_prefix, "_full.csv"))
cat("Saved full results to:\n", paste0(out_prefix, "_full.csv"), "\n")

# -----------------------------
# 7) Significant genes + Top10 up/down
# -----------------------------
sig <- res_df %>%
  dplyr::filter(!is.na(padj), padj < 0.05)

cat("Significant genes (padj < 0.05):", nrow(sig), "\n")

top10_up <- sig %>%
  dplyr::filter(log2FoldChange > 0) %>%
  dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::slice_head(n = 10)

top10_down <- sig %>%
  dplyr::filter(log2FoldChange < 0) %>%
  dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::slice_head(n = 10)

readr::write_csv(top10_up,   paste0(out_prefix, "_Top10_Up.csv"))
readr::write_csv(top10_down, paste0(out_prefix, "_Top10_Down.csv"))

cat("\nTop10 UP (padj<0.05):\n")
print(top10_up %>% dplyr::select(gene_id, ENSEMBL, SYMBOL, baseMean, log2FoldChange, padj))

cat("\nTop10 DOWN (padj<0.05):\n")
print(top10_down %>% dplyr::select(gene_id, ENSEMBL, SYMBOL, baseMean, log2FoldChange, padj))

cat("\nSaved Top10 files:\n",
    paste0(out_prefix, "_Top10_Up.csv\n"),
    paste0(out_prefix, "_Top10_Down.csv\n"))

# -----------------------------
# 8) Volcano plot
# -----------------------------
plot_df <- res_df %>%
  dplyr::mutate(
    status = dplyr::case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      !is.na(padj) & padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    neglog10padj = dplyr::if_else(is.na(padj), NA_real_, -log10(padj))
  )

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FoldChange, y = neglog10padj, color = status)) +
  ggplot2::geom_point(alpha = 0.75, size = 1.6) +
  ggplot2::scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not significant" = "grey70"
  )) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggplot2::labs(
    title = paste0("Volcano plot: ", dataset_name, " Treated vs Untreated — Gene expression"),
    x = "log2 Fold Change (treated / untreated)",
    y = "-log10(adjusted p-value)",
    color = "Legend"
  ) +
  ggplot2::theme_minimal(base_size = 13)

volcano_path <- paste0(out_prefix, "_volcano.png")
cat("Saving volcano plot to:\n", volcano_path, "\n")

ggplot2::ggsave(volcano_path, p, width = 7.5, height = 6.2, dpi = 300)

cat("DONE \nAll outputs saved with prefix:\n", out_prefix, "\n")
