suppressWarnings(suppressMessages({
  library(data.table)
  library(ape)
}))

# ----- Arguments -----
args <- commandArgs(trailingOnly = TRUE)

mdist_file <- args[1]
ids_file   <- args[2]
out_prefix <- args[3]     # <- ВАЖНО: здесь переменная определяется!

# ----- Load data -----
M   <- as.matrix(fread(mdist_file, header = FALSE))
ids <- fread(ids_file, header = FALSE)
labs <- ids$V2

rownames(M) <- colnames(M) <- labs
D <- as.dist(M)
tree <- nj(D)

# Ladderize improves visualization
tree <- ladderize(tree)

n_tips <- length(tree$tip.label)

# ----- Plot -----
pdf(paste0(out_prefix), width = 10, height = 15)

plot(
  tree,
  type = "phylogram",
  cex = 0.4,
  no.margin = TRUE,
  x.lim = max(tree$edge.length) * 1.5,
  y.lim = c(0, n_tips * 1.2),
  label.offset = 0.01,
  edge.width = 0.6
)

title("Neighbor-Joining tree (1 - IBS)", cex.main = 1.5)

dev.off()
