#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Transcriptomic timing.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(svglite)
  library(ragg)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else "hitplot_input.tsv"
output_dir <- if (length(args) >= 2) args[[2]] else "."
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Figure contract
# Conclusion: the DET-weighted timing of transcriptomic ageing differs among tissues.
# Evidence: preserve the original lollipop plot, interval shading, reference lines,
#           tissue palette and total-DET size encoding; only resize and typeset it.

width_mm <- 69.5
height_mm <- 58.5
pt_to_mm <- 25.4 / 72

text_col <- "#333333"
axis_col <- "#555555"
guide_col <- "#D1D6DA"
reference_col <- "#9AA1A8"

tissue_palette <- c(
  "Heart" = "#CD0000", "Lung" = "#7FFFD4", "Trachea" = "#00EEEE",
  "Pectoral" = "#EED5B7", "Thigh_muscle" = "#CDAA7D", "Spleen" = "#CDC9C9",
  "Hypothalamus" = "#6A0C8A", "Pituitary" = "#9B59B6", "Cerebrum" = "#A74B9A",
  "Cerebellum" = "#DDA0DD", "Medulla_oblongata" = "#DA70D6",
  "Olfactory_bulb" = "#BA55D3", "Kidney" = "#8B0000", "Ovary" = "#EE9A49",
  "Infundibulum" = "#FF4500", "Magnum" = "#FF8C00", "Isthmus" = "#FFA500",
  "Uterine" = "#FFD700", "Vignia" = "#FFA07A", "Esophagus" = "#0000FF",
  "Crop" = "#1E90FF", "Proventriculus" = "#4682B4", "Gizzard" = "#87CEEB",
  "Duodenum" = "#008000", "Jejunum" = "#32CD32", "Ileum" = "#00FF7F",
  "Cecum" = "#8FBC8F", "Rectum" = "#00FF00", "Liver" = "#FFC0CB",
  "Gallbladder" = "#FF69B4", "Pancreas" = "#CD6090", "Adipose" = "#EEEE00",
  "Retina" = "#000000"
)

raw <- as_tibble(read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE))
required <- c("Tissue", "Group", "DEGs_number")
if (!all(required %in% names(raw))) {
  stop("Input must contain: ", paste(required, collapse = ", "))
}

timing <- raw %>%
  transmute(
    Tissue = as.character(Tissue),
    Comparison = as.character(Group),
    DET = as.numeric(DEGs_number)
  ) %>%
  filter(Comparison %in% c("20vs50", "50vs80", "80vs100")) %>%
  group_by(Tissue, Comparison) %>%
  summarise(DET = sum(DET, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Comparison, values_from = DET, values_fill = 0) %>%
  transmute(
    Tissue,
    D_20_50 = `20vs50`,
    D_50_80 = `50vs80`,
    D_80_100 = `80vs100`,
    Total_DET = D_20_50 + D_50_80 + D_80_100,
    p_20_50 = D_20_50 / Total_DET,
    p_50_80 = D_50_80 / Total_DET,
    p_80_100 = D_80_100 / Total_DET,
    Timing_score = 35 * p_20_50 + 65 * p_50_80 + 90 * p_80_100
  ) %>%
  arrange(Timing_score, Tissue)

missing_colours <- setdiff(timing$Tissue, names(tissue_palette))
if (length(missing_colours) > 0) {
  stop("Missing tissue colours for: ", paste(missing_colours, collapse = ", "))
}

display_labels <- setNames(gsub("_", " ", timing$Tissue), timing$Tissue)
display_labels["Vignia"] <- "Vagina"

timing <- timing %>%
  mutate(
    Tissue_factor = factor(Tissue, levels = rev(Tissue)),
    Tissue_colour = Tissue
  )

write_tsv(
  timing %>%
    select(Tissue, D_20_50, D_50_80, D_80_100, Total_DET,
           p_20_50, p_50_80, p_80_100, Timing_score),
  file.path(output_dir, "transcriptomic_aging_timing_score_NG_source_data.tsv")
)

p <- ggplot(timing, aes(x = Timing_score, y = Tissue_factor)) +
  annotate("rect", xmin = 50, xmax = 60, ymin = -Inf, ymax = Inf,
           fill = "#D9E7F2", alpha = 0.36) +
  annotate("rect", xmin = 60, xmax = 70, ymin = -Inf, ymax = Inf,
           fill = "#E8E8E8", alpha = 0.42) +
  annotate("rect", xmin = 70, xmax = 80, ymin = -Inf, ymax = Inf,
           fill = "#F2D8D6", alpha = 0.34) +
  geom_vline(
    xintercept = c(60, 70),
    linetype = "22",
    linewidth = 0.5 * pt_to_mm,
    colour = "#6F7782"
  ) +
  geom_segment(
    aes(x = 50, xend = Timing_score, yend = Tissue_factor),
    linewidth = 0.5 * pt_to_mm,
    colour = "#C8CDD3"
  ) +
  geom_point(
    aes(size = Total_DET, fill = Tissue_colour),
    shape = 21,
    stroke = 0.18,
    colour = "#202124",
    alpha = 0.96
  ) +
  scale_fill_manual(values = tissue_palette, guide = "none") +
  scale_size_continuous(
    range = c(1.05, 2.60),
    breaks = c(15000, 25000, 35000),
    labels = comma,
    name = "Total DETs"
  ) +
  scale_x_continuous(
    limits = c(50, 80),
    breaks = c(50, 55, 60, 65, 70, 75, 80),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_discrete(labels = display_labels, expand = expansion(add = c(0.35, 0.35))) +
  labs(
    x = "Weighted estimate of aging onset (weeks)",
    y = NULL
  ) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      title.hjust = 0,
      override.aes = list(
        size = c(1.15, 1.85, 2.55),
        fill = "white", colour = "#202124", alpha = 1, stroke = 0.18
      )
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "Arial", base_size = 5.2) +
  theme(
    text = element_text(family = "Arial", colour = text_col),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_line(colour = "#E5E8EC", linewidth = 0.25 * pt_to_mm),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(colour = axis_col, linewidth = 0.5 * pt_to_mm),
    axis.line.y = element_line(colour = axis_col, linewidth = 0.5 * pt_to_mm),
    axis.ticks = element_line(colour = axis_col, linewidth = 0.5 * pt_to_mm),
    axis.ticks.length = unit(1.3, "pt"),
    axis.text.x = element_text(size = 5.2, colour = axis_col, margin = margin(t = 1.5, unit = "pt")),
    axis.text.y = element_text(size = 5.0, colour = text_col, lineheight = 0.86,
                               margin = margin(r = 1.8, unit = "pt")),
    axis.title.x = element_text(size = 6.0, colour = text_col, margin = margin(t = 3.0, unit = "pt")),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.92), colour = "#777777",
                                     linewidth = 0.4 * pt_to_mm),
    legend.margin = margin(1.5, 2.0, 1.5, 2.0, unit = "pt"),
    legend.title = element_text(size = 5.2, face = "plain", colour = text_col),
    legend.text = element_text(size = 5.0, colour = text_col),
    legend.key.height = unit(5.5, "pt"),
    legend.key.width = unit(7.0, "pt"),
    legend.spacing.y = unit(0.1, "pt"),
    plot.margin = margin(t = 2.5, r = 2.5, b = 2.0, l = 2.0, unit = "pt")
  )

pdf_file <- file.path(output_dir, "transcriptomic_aging_timing_score_NG.pdf")
svg_file <- file.path(output_dir, "transcriptomic_aging_timing_score_NG.svg")
png_file <- file.path(output_dir, "transcriptomic_aging_timing_score_NG.png")

ggsave(
  pdf_file, p,
  device = cairo_pdf,
  # Windows cairo truncates the PDF media box to whole points; +0.10 mm
  # yields a 166-pt page (58.56 mm), the closest representable height to 58.5 mm.
  width = width_mm, height = height_mm + 0.10, units = "mm",
  family = "Arial", bg = "white", limitsize = FALSE
)

ggsave(
  svg_file, p,
  device = function(filename, width, height, ...) {
    svglite(filename, width = width, height = height, bg = "white",
            system_fonts = list(sans = "Arial"))
  },
  width = width_mm, height = height_mm, units = "mm",
  bg = "white", limitsize = FALSE
)

ggsave(
  png_file, p,
  device = agg_png,
  width = width_mm, height = height_mm, units = "mm",
  dpi = 600, background = "white", limitsize = FALSE
)

message("Saved: ", pdf_file)
message("Saved: ", svg_file)
message("Saved: ", png_file)
