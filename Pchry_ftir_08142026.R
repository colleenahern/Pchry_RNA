# C. B. Ahern
# Plotting GFCA and SACA ftir data

library(tidyverse)
library(patchwork)

# ---- File paths ----
gf_path <- "/Users/colleenahern/Documents/Pchry/ftir/GF_CA_film.CSV"
sa_path <- "/Users/colleenahern/Documents/Pchry/ftir/SA_CA_powder_nocg.CSV"

# ---- Read data ----
read_ftir <- function(path) {
  read_csv(path, col_names = c("wavenumber", "transmittance"), show_col_types = FALSE) %>%
    filter(transmittance > 0)  # trim zero-padded baseline artifacts outside true scan range
}

gf <- read_ftir(gf_path) %>% mutate(sample = "GF CA (lowly-acetylated film)")
sa <- read_ftir(sa_path) %>% mutate(sample = "SA CA (highly-acetylated powder)")

# Shared x-axis limits across both panels
xlims <- range(c(gf$wavenumber, sa$wavenumber))

# ---- Build each panel ----
p_gf <- ggplot(gf, aes(x = wavenumber, y = transmittance)) +
  geom_line(color = "#2166ac", linewidth = 0.4) +
  scale_x_reverse(limits = c(4000, 600)) +
  labs(title = "GF CA (lowly-acetylated cellulose acetate)", x = NULL, y = "Transmittance (%)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

p_sa <- ggplot(sa, aes(x = wavenumber, y = transmittance)) +
  geom_line(color = "#b2182b", linewidth = 0.4) +
  scale_x_reverse(limits = c(4000, 600)) +
  labs(title = "SA CA (highly-acetylated cellulose acetate)",
       x = expression("Wavenumber (cm"^-1*")"), y = "Transmittance (%)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )
# ---- Stack panels vertically ----
combined <- p_gf / p_sa 

combined

n# ---- Save ----
ggsave("/Users/colleenahern/Documents/Pchry/Figures/GF_SA_CA_FTIR_stacked.pdf", combined, width = 8, height = 8, dpi = 300)
