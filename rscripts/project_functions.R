
coolers <- list.files(here::here("data", "hic", "pools"), full.names = TRUE)
names(coolers) <- data.table::tstrsplit(basename(coolers), "\\.")[[1]]

sample_lut <- c(
  "P14_HB" = "P01_HB",
  "P14_PB" = "P01_PB",
  "P14_LM" = "P01_LM",
  "P12_LM" = "P01_??",
  "P12_HB" = "P02_HB",
  "P16_PB" = "P02_PB",
  "P12_PB" = "P02_BR",
  "P11_PB" = "P03_PB",
  "P11_LM" = "P03_LM",
  "P15_PB" = "P04_PB",
  "P15_LM" = "P04_LM",
  "P13_PB" = "P05_PB",
  "P04_PE" = "P06_PE",
  "P05_PE" = "P07_PE",
  "P06_PE" = "P08_PE",
  "P07_PE" = "P09_PE",
  "P08_PE" = "P10_PE",
  "P09_PE" = "P11_PE",
  "P10_PE" = "P12_PE",
  "P01_PE" = "P13_PE", 
  "P02_PE" = "P14_PE",
  "P03_PE" = "P15_PE"
  # Previous numbering system
  # "P01_PE" = "P01_PE", 
  # "P02_PE" = "P02_PE",
  # "P03_PE" = "P03_PE",
  # "P04_PE" = "P04_PE",
  # "P05_PE" = "P05_PE",
  # "P06_PE" = "P06_PE",
  # "P07_PE" = "P07_PE",
  # "P08_PE" = "P08_PE",
  # "P09_PE" = "P09_PE",
  # "P10_PE" = "P10_PE",
  # "P11_LM" = "P11_LM",
  # "P11_PB" = "P11_PB",
  # "P12_HB" = "P12_HB",
  # "P12_LM" = "P14_??",
  # "P12_PB" = "P12_BR",
  # "P13_PB" = "P13_PB",
  # "P14_HB" = "P14_HB",
  # "P14_PB" = "P14_PB",
  # "P14_LM" = "P14_LM",
  # "P15_PB" = "P15_PB",
  # "P15_LM" = "P15_LM",
  # "P16_PB" = "P12_PB"
)

names(coolers) <- sample_lut[names(coolers)]

load_experiments <- function(
  resolution = 1e6,
  sample_subset = NULL, raw = FALSE,
  ...
) {
  if (is.null(sample_subset)) {
    files <- coolers
  } else {
    files <- coolers[sample_subset]
  }
  if (length(files) == 0 || any(!file.exists(files))) {
    stop("Cannot find samples.")
  }
  if (!rlang::is_named(files)) {
    stop("Invalidly named files")
  }
  
  exps <- mapply(
    GENOVA::load_contacts,
    signal_path = files,
    sample_name = names(files),
    MoreArgs  = c(list(resolution = resolution), list(...)),
    SIMPLIFY  = FALSE,
    USE.NAMES = FALSE
  )
  setNames(exps, names(files))
}

# Plotting ----------------------------------------------------------------

set_figure_theme <- function() {
  mycolour <- "#000000FF" # Opaque Black
  
  # Geom defaults
  
  text_defaults <- list(family = "Arial", size = 6 / .pt, colour = mycolour)
  
  ggplot2::update_geom_defaults("text",  text_defaults)
  ggplot2::update_geom_defaults("label", text_defaults)
  
  # Theme defaults 
  
  ggplot2::theme_set(ggplot2::theme_gray())
  ggplot2::theme_update(
    text              = ggplot2::element_text(colour = mycolour, size = 8,
                                              family = "Arial"),
    line              = ggplot2::element_line(colour = mycolour),
    axis.line         = ggplot2::element_line(colour = mycolour, lineend = "square"),
    axis.ticks        = ggplot2::element_line(colour = mycolour),
    axis.text         = ggplot2::element_text(colour = mycolour, size = 6),
    legend.title      = ggplot2::element_text(colour = mycolour, size = 6),
    legend.text       = ggplot2::element_text(colour = mycolour, size = 6),
    legend.key        = ggplot2::element_blank(),
    legend.background = ggplot2::element_rect(colour = NA, fill = NA),
    legend.box.spacing = unit(5.5, "pt"),
    legend.spacing    = unit(5.5, "pt"),
    panel.background  = ggplot2::element_blank(),
    panel.grid.major  = ggplot2::element_blank(),
    panel.grid.minor  = ggplot2::element_blank(),
    plot.background   = ggplot2::element_blank(),
    strip.background  = ggplot2::element_blank(),
    strip.text        = ggplot2::element_text(colour = mycolour, size = 8),
    strip.placement   = "outside"
  )
}

thin_line <- ggplot2::element_line(colour = "grey95", linewidth = 0.25)
dark_outline <- ggplot2::aes(colour = after_scale(colorspace::darken(fill, 0.2)))

number_signed <- function(
  x, ..., big.mark = ",",
  style_positive = "plus",
  style_negative = "minus"
) {
  scales::number(
    x, ..., big.mark = big.mark,
    style_positive = style_positive,
    style_negative = style_negative
  )
}

number_integer <- function(x, ..., accuracy = 1, big.mark = ",") {
  scales::number(x, accuracy = accuracy, big.mark = big.mark, ...)
}

save_svg <- function(filename, plot = ggplot2::last_plot(),
                     width = NA, height = NA, units = "mm") {
  if (!grepl("\\.svg$", filename)) {
    filename <- paste0(filename, ".svg")
  }
  file <- here("figures", "manuscript", filename)
  ggplot2::ggsave(
    filename = file,
    plot = plot, width = width, height = height, units = units,
    device = svglite::svglite, fix_text_size = FALSE
  )
  knitr::include_graphics(file)
}

vcolourbar <- function(..., barwidth = unit(5.5, "pt"), 
                       direction = "vertical",
                       title.position = "top") {
  ggplot2::guide_colourbar(
    ..., barwidth = barwidth, direction = direction, 
    title.position = title.position
  )
}

hcolourbar <- function(..., barheight = unit(5.5, "pt"),
                       direction = "horizontal",
                       title.position = "top") {
  ggplot2::guide_colourbar(
    ..., barheight = barheight, direction = direction,
    title.position = title.position
  )
}

sample_type <- c("HB" = "Healthy Breast", "PB" = "Primary Breast",
                 "BR" = "Breast Recurrent", "LM" = "Liver Metastasis",
                 "PE" = "Pleural Effusion")
sample_colours <- c("HB" = "#00E6FF", "PB" = "#8F9FFF", "BR" = "#88FF99",
                    "LM" = "#F546B4", "PE" = "red")

div_colours <- div_colours <- 
  c('#0065b2', '#2399de', '#67cdfe', '#f5f5f5', '#fcb09d', '#ed6855', '#be2a21')
