# Libraries
library(data.table)
library(BiocFileCache)
library(GenomicRanges)

# Parameters
chromos <- paste0("chr", c(1:22, "X"))

# Cache file
bfc <- BiocFileCache()
resource <- "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"
path <- bfcrpath(bfc, resource)

# Process file
ideo <- fread(path)
ideo <- ideo[V1 %in% chromos & V5 == "acen"]
ideo$V1 <- factor(ideo$V1, chromos)
ideo <- ideo[, .(V2 = min(V2), V3 = max(V3)), by = "V1"]
gr   <- ideo[, GRanges(V1, IRanges(V2 + 1, V3))]
gr   <- sort(gr)

# Saving
saveRDS(gr, here("references", "centromeres.rds"))
