
library(magick)
library(data.table)
library(here)

meta <- readRDS(here("rds", "translocation_images.rds"))
meta <- rbindlist(meta, idcol = "sample")

# Shuffle
set.seed(42)
meta <- meta[sample(nrow(meta)),]

answer <- character(nrow(meta))

for (i in 1:length(answer)) {
  img <- image_read(meta$file[i])
  x <- print(img)
  Sys.sleep(1)
  answer[i] <- readline(prompt = "Translocation? T/F: ")
}

meta$label <- answer == "T"

saveRDS(meta, here("rds", "translocation_confirmed.rds"))
