#!/usr/bin/env Rscript
# Uso:
# Rscript abundance_species_plot.R input=merged_abundance_table_species.tsv output=test2pct.png title=test min_pct=2
# filtrando solo especies con ≥2% en alguna muestra.

args <- commandArgs(trailingOnly = TRUE) #turn on arguments
kv <- function(a){o<-list();for(x in a) if(grepl("=",x)){k<-sub("=.*","",x);v<-sub("^[^=]*=","",x);o[[k]]<-v};o}
opt <- kv(args); get <- function(n,d) if(is.null(opt[[n]])) d else opt[[n]]

##VARIABLES##
input <- get("input") #Load input table
output <- get("output") #output name
plot_title <- get("title", "Species abundance") #graph title
dpi    <- as.integer(get("dpi", "200")) #png output dpi
width  <- as.numeric(get("width","18")) #png output width
height <- as.numeric(get("height","14")) #png output height
samples_order <- get("samples_order", NA)  # optional: txt file with sample order
min_pct <- as.numeric(get("min_pct","0"))  # minimum percentage of species abundance to show in the graph. Pon 2 para filtrar ≥2% si quieres.

###Installing necessary libraries if the case
need <- c("readr","dplyr","tidyr","ggplot2","forcats","stringr","digest", "colorspace", "ggsci", "RColorBrewer")
miss <- setdiff(need, rownames(installed.packages()))
if(length(miss)>0){ install.packages(miss, repos="https://cloud.r-project.org") }
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(forcats); library(stringr); library(digest); library(colorspace); library(ggsci); library(RColorBrewer)
})

#Set working director
#setwd("~/Documentos/candida_auris_project")

# --- Lectura TSV (primera col = clade_name) ---
df <- readr::read_tsv(input, show_col_types = FALSE)
#df <- readr::read_tsv(file = "merged_abundance_table_species_contaminated_samples.tsv", show_col_types = FALSE)
stopifnot(ncol(df) >= 2)
colnames(df)[1] <- "Species"
num_cols <- setdiff(colnames(df), "Species")

# Numeric + NAs -> 0
df[num_cols] <- lapply(df[num_cols], function(x) suppressWarnings(as.numeric(x)))
df[is.na(df)] <- 0

# (Opcional) filtra especies con ≥ min_pct en alguna muestra
if (min_pct > 0) {
  df <- df %>% filter(apply(across(all_of(num_cols)), 1, max) >= min_pct)
}

# Elimina especies todo 0
df <- df %>% filter(rowSums(across(all_of(num_cols))) > 0)

# orden de muestras opcional
if (!is.na(samples_order) && file.exists(samples_order)) {
  desired <- readr::read_lines(samples_order)
  desired <- desired[desired != ""]
  exist <- intersect(desired, num_cols)
  if (length(exist) > 0) { df <- df %>% select(Species, all_of(exist)); num_cols <- exist }
}

# Ordena species por abundancia total (asc; las más grandes quedan arriba en el apilado)
df <- df %>%
  mutate(Total = rowSums(across(all_of(num_cols)))) %>%
  arrange(Total) %>% select(-Total)

# Long table fpr plotting
long <- df %>%
  pivot_longer(all_of(num_cols), names_to="Sample", values_to="Abundance") %>%
  mutate(Sample = factor(Sample, levels = num_cols),
         Species = fct_inorder(Species))

#################
# Opción 1: Set3
#palette_vec <- setNames(colorRampPalette(brewer.pal(12, "Spectral"))(length(species_all)), species_all)
# Paleta Set3 (expandida si hay >12 especies)
#n_sp <- length(levels(long$Species))
#base_set3 <- brewer.pal(12, "Spectral")
#set3_extended <- colorRampPalette(base_set3)(max(n_sp, length(base_set3)))
#palette_vec <- setNames(set3_extended[seq_len(n_sp)], levels(long$Species))
###################

## Paleta de 25 colores (Basada en Paired + Extended)
microbiome_palette <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928", "#00FFFF", "#1B9E77", "#D95F02",
  "#D9D9D9", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
  "#666666", "#F781BF", "#8DD3C7", "#FFFFB3", "#BEBADA"
)

## Función para obtener n colores
microbiome_pal <- function(n) {
  if (n > length(microbiome_palette)) {
    warning("Se han pedido más de 25 colores, se reciclarán colores.")
    return(rep(microbiome_palette, length.out = n))
  }
  microbiome_palette[seq_len(n)]
}

# Construir paleta para las Species
n_sp <- length(levels(long$Species))
palette_vec <- setNames(microbiome_pal(n_sp), levels(long$Species))


# Plot
p <- ggplot(long, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = palette_vec) +
  labs(title = plot_title, x = "Sample", y = "Relative Abundance (%)", fill = "Species") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# Plot con ggsci solo funciona con hasta 8 colores!!!!!
#p <- ggplot(long, aes(x = Sample, y = Abundance, fill = Species)) +
#  geom_col() +
#  coord_flip() +
#  labs(title = plot_title, x = "Sample", y = "Relative abundance (%)", fill = "Species") +
#  theme_bw(base_size = 11) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1),
#        legend.position = "right") +
#  scale_fill_jco()    # ← puedes cambiar por: _aaas(), _lancet(), _nejm(), etc.

# Exportar PNG (dispositivo base)
png(filename = output, width = width * 100, height = height * 100, res = dpi)
print(p)
dev.off()

# Exportar SVG (dispositivo base, no requiere svglite)
svg_out <- sub("\\.[A-Za-z0-9]+$", ".svg", output)
svg(filename = svg_out, width = width, height = height)
print(p)
dev.off()

message("[OK] Guardado: ", output, " y ", svg_out)