cat > skmel28_deseq2_10_vs_01.R <<'EOF'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(readr)
})

# -----------------------------
# SETTINGS (EDIT ONLY THIS PATH)
# -----------------------------
quant_dir <- "~/Desktop/New_Skml28/gene expreission/skmel28"
out_dir   <- file.path(quant_dir, "DESeq2_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# INPUT FILES (Skmel28 only)
# -----------------------------
files <- list.files(quant_dir, pattern="^Skmel28-.*\\.quant\\.sf$", full.names=TRUE)
if (length(files) < 2) stop("Found <2 Skmel28 quant.sf files. Check quant_dir path/pattern.")

sample_id <- sub("\\.quant\\.sf$", "", basename(files))

# Condition extracted from filename like: Skmel28-01_... or Skmel28-10_...
condition <- sub("^Skmel28-([0-9]+)_.*$", "\\1", sample_id)
condition <- factor(condition, levels = c("01", "10"))  # reference = 01

sampleTable <- data.frame(
  sample = sample_id,
  condition = condition,
  file = files,
  stringsAsFactors = FALSE
)

# Check both groups exist
tab <- table(sampleTable$condition)
print(tab)
if (!all(c("01","10") %in% names(tab))) {
  stop("Missing a condition group. Need BOTH Skmel28-01_*.quant.sf and Skmel28-10_*.quant.sf in quant_dir.")
}

write.csv(sampleTable, file.path(out_dir, "sampleTable.csv"), row.names = FALSE)

# tximport expects named vector
names(files) <- sample_id

# -----------------------------
# BUILD tx2gene + gene_symbol map from Name column
# Name format:
# ENST | ENSG | ... | transcript_name | gene_symbol | ...
# So gene symbol = field 6 (based on your files)
# -----------------------------
ex <- read_tsv(files[1], show_col_types = FALSE)
parts <- strsplit(ex$Name, "\\|")

tx_id <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
gn_id <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
gsym  <- vapply(parts, function(x) if (length(x) >= 6) x[6] else NA_character_, character(1))

tx2gene <- unique(data.frame(TXNAME = tx_id, GENEID = gn_id, stringsAsFactors = FALSE))
tx2gene <- tx2gene[!is.na(tx2gene$TXNAME) & !is.na(tx2gene$GENEID), ]

# one gene -> one symbol
gene_map <- unique(data.frame(GENEID = gn_id, gene_symbol = gsym, stringsAsFactors = FALSE))
gene_map <- gene_map[!is.na(gene_map$GENEID), ]
gene_map <- gene_map[!duplicated(gene_map$GENEID), ]

# -----------------------------
# TXIMPORT (Salmon -> gene-level)
# IMPORTANT: Name contains "|" so use ignoreAfterBar=TRUE
# -----------------------------
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM",
  ignoreAfterBar = TRUE
)

# -----------------------------
# DESEQ2
# -----------------------------
dds <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ condition)

# filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

# Compare 10 vs 01
res <- results(dds, contrast = c("condition", "10", "01"))
res_df <- as.data.frame(res)
res_df$ENSG <- rownames(res_df)

# attach gene symbols safely
res_df <- merge(res_df, gene_map, by.x = "ENSG", by.y = "GENEID", all.x = TRUE)
res_df <- res_df[!duplicated(res_df$ENSG), ]
res_df$gene_symbol[is.na(res_df$gene_symbol) | res_df$gene_symbol == ""] <- res_df$ENSG

write.csv(res_df, file.path(out_dir, "DESeq2_results_10_vs_01.csv"), row.names = FALSE)

# -----------------------------
# Volcano prep (handle NA for plotting)
# -----------------------------
res_df$padj_plot <- res_df$padj
res_df$padj_plot[is.na(res_df$padj_plot)] <- 1
res_df$neglog10padj <- -log10(res_df$padj_plot)

res_df$group <- "NS"
res_df$group[is.na(res_df$padj)] <- "NA"
res_df$group[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "Up"
res_df$group[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

# top 10 up/down (tested only)
res_ok <- res_df[!is.na(res_df$padj), ]
top_up   <- head(res_ok[order(res_ok$padj, -res_ok$log2FoldChange), ], 10)
top_down <- head(res_ok[order(res_ok$padj,  res_ok$log2FoldChange), ], 10)

write.csv(top_up, file.path(out_dir, "top10_up.csv"), row.names = FALSE)
write.csv(top_down, file.path(out_dir, "top10_down.csv"), row.names = FALSE)

lab <- unique(rbind(
  top_up[,   c("log2FoldChange", "neglog10padj", "gene_symbol")],
  top_down[, c("log2FoldChange", "neglog10padj", "gene_symbol")]
))

# -----------------------------
# COLORFUL VOLCANO
# -----------------------------
p <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10padj, color = group)) +
  geom_point(alpha = 0.7, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_minimal() +
  labs(
    title = "Skmel28 DESeq2: condition 10 vs 01",
    x = "log2 fold change (10 vs 01)",
    y = "-log10(adj p-value)"
  ) +
  scale_color_manual(values = c(
    "Up"   = "red3",
    "Down" = "royalblue3",
    "NS"   = "grey60",
    "NA"   = "grey85"
  ))

p <- p + geom_text_repel(
  data = lab,
  inherit.aes = FALSE,
  aes(x = log2FoldChange, y = neglog10padj, label = gene_symbol),
  size = 3,
  max.overlaps = 30
)

ggsave(file.path(out_dir, "volcano_10_vs_01_color.png"), plot = p, width = 8, height = 6, dpi = 300)

cat("Done! Results in:", out_dir, "\n")
EOF
