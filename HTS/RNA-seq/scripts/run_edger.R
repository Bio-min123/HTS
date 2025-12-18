args <- commandArgs(trailingOnly=TRUE)
counts_file <- args[1]
meta_file <- args[2]
out_file <- args[3]

library(edgeR)

# Load counts and metadata
counts <- read.table(counts_file, header=TRUE, row.names=1, check.names=FALSE)
meta <- read.csv(meta_file, row.names=1)

group <- factor(meta$condition)
dge <- DGEList(counts=counts, group=group)

# Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Normalize
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~group)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

res <- topTags(lrt, n=Inf)
write.csv(res$table, file=out_file)
