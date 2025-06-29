
library(scales)
library(here)
library(data.table)
library(GENOVA)
library(grid)

calls <- list.files(here("processed_data", "translocations"), "_computed.tsv", full.names = TRUE)
names(calls) <- do.call(paste, c(tstrsplit(basename(calls), "_")[1:2], list(sep = "_")))

files <- list.files(
  here("data", "hic", "pools"),
  full.names = TRUE
)
names(files) <- tstrsplit(basename(files), "\\.")[[1]]

stopifnot(names(calls) == names(files))


pal <- colour_ramp(myplotdefaults:::fire_colours)

window <- CJ(-300:300, -300:300)

all <- list()

for (i in seq_along(files)) {
  
  trans <- fread(calls[i])
  exp   <- load_contacts(files[i], balancing = FALSE, resolution = 200e3)
  
  trans[, x_i := bed2idx(exp$IDX, data.frame(chr_x, start_x, end_x))]
  trans[, y_i := bed2idx(exp$IDX, data.frame(chr_y, start_y, end_y))]
  trans$file <- ""
  
  template <- matrix("#FFFFFF", nrow = 601, ncol = 601)
  
  bins <- exp$IDX
  
  rect <- rectGrob(width = unit(0.9, "npc"), height = unit(0.9, "npc"),
                   gp = gpar(col = 1, fill = NA))
  circ <- circleGrob(r = 0.05, gp = gpar(col = "limegreen", fill = NA))
  pxl  <- round(601 * (1 / 0.9))
  fname <- here("processed_data", "translocations", "images", names(files)[i])
  
  for (j in seq_len(nrow(trans))) {
    
    idx <- copy(window)
    idx[, V1 := V1 + trans$x_i[j]]
    idx[, V2 := V2 + trans$y_i[j]]
    idx[, c("V1", "V2") := list(pmin(V1, V2), pmax(V1, V2))]
    
    mins <- idx[, .(V1 = min(V1), V2 = min(V2))]
    
    bin_left  <- bins[V4 %in% idx$V1]
    bin_right <- bins[V4 %in% idx$V2]
    
    mat <- exp$MAT[idx]
    mat[is.na(V3), V3 := 0]
    mat[, V3 := pmin(V3, 300) / 300]
    mat[, V3 := pal(V3)]
    img <- template
    img[mat[, cbind(V1 - mins$V1 + 1, V2 - mins$V2 + 1)]] <- mat$V3
    
    x_breaks <- bin_right[which(tail(V1, -1) != head(V1, -1)), V4] - 0.5
    y_breaks <- bin_left[ which(tail(V1, -1) != head(V1, -1)), V4] - 0.5
    if (length(x_breaks) > 0) {
      x_breaks <- rescale((x_breaks - mins$V2 + 1) / 601, to = c(0.05, 0.95), from = c(0, 1))
      x_breaks <- polylineGrob(
        y = unit(rep(c(0.05, 0.95), length(x_breaks)), "npc"), 
        x = unit(rep(x_breaks, each = 2), "npc"),
        gp = gpar(lty = "dotted", lwd = 0.5),
        id.lengths = rep(2, length(x_breaks))
      )
    } else {
      x_breaks <- ggplot2::zeroGrob()
    }
    if (length(y_breaks) > 0) {
      y_breaks <- rescale((y_breaks - mins$V1 + 1) / 601, to = c(0.05, 0.95), from = c(0, 1))
      y_breaks <- polylineGrob(
        x = unit(rep(c(0.05, 0.95), length(y_breaks)), "npc"), 
        y = unit(rep(y_breaks, each = 2), "npc"),
        gp = gpar(lty = "dotted", lwd = 0.5),
        id.lengths = rep(2, length(y_breaks))
      )
    } else {
      y_breaks <- ggplot2::zeroGrob()
    }
    
    lab <- trans[j, paste(chr_x, "&", chr_y)]
    lab <- textGrob(x = unit(0.05, "npc"), y = unit(0.975, "npc"), hjust = 0, label = lab)
    
    grb <- rasterGrob(img, width = unit(0.9, "npc"), height = unit(0.9, "npc"),
                      interpolate = FALSE)
    
    tree <- gTree(
      children = gList(grb, y_breaks, x_breaks, rect, circ, lab)
    )
    filename <- paste0(fname, "_", j, ".png")
    ragg::agg_png(
      filename,
      width = pxl, height = pxl
    )
    grid.newpage(); grid.draw(tree)
    dev.off()
    
    trans$file[j] <- filename
  }
  all[[names(files)[i]]] <- trans
}

saveRDS(all, here("rds", "translocation_images.rds"))

