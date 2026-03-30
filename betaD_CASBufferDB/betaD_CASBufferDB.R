# =============================================================================
# Beta diversity analysis — Dung beetles
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(vegan)
library(betapart)
library(ggplot2)
library(ggforce)

# --- Data --------------------------------------------------------------------
dataDB <- read.csv("dataCASBufferDB.csv", header = TRUE, sep = ";")
dataDB$vegcover <- factor(dataDB$vegcover,
                          levels = c("SF", "CAS", "PA"),
                          labels = c("Secondary forest",
                                     "Cacao agroforestry system",
                                     "Pasture"))

# =============================================================================
# 1. COMMUNITY MATRIX — mean abundances aggregated at site level
# =============================================================================

matriz_abund_site <- dataDB %>%
  group_by(vegcover, site, species) %>%
  summarise(abund_mean = mean(abund), .groups = "drop") %>%
  unite("site_id", vegcover:site, remove = FALSE) %>%
  pivot_wider(names_from = species, values_from = abund_mean, values_fill = 0)

# Remove descriptor columns; keep species matrix only
matriz_site <- as.data.frame(matriz_abund_site[, -c(1:2)])
rownames(matriz_site) <- matriz_abund_site$site_id

# =============================================================================
# 2. SITE METADATA
# =============================================================================

metadata_site <- matriz_abund_site %>%
  distinct(site_id, vegcover, site)

# =============================================================================
# 3. SQUARE-ROOT TRANSFORMATION AND BRAY-CURTIS DISTANCE
# =============================================================================

matriz_site_sqrt <- sqrt(matriz_site)
dist_bray_site   <- vegdist(matriz_site_sqrt, method = "bray")

# =============================================================================
# 4. PERMANOVA — vegetation cover as grouping factor
# =============================================================================

permanova_site <- adonis2(dist_bray_site ~ vegcover,
                          data = metadata_site, permutations = 999)
print(permanova_site)

# =============================================================================
# 5. PERMDISP — within-group multivariate dispersion
# =============================================================================

mod_disp_site <- betadisper(dist_bray_site, metadata_site$vegcover)
print(anova(mod_disp_site))
print(permutest(mod_disp_site))
plot(mod_disp_site)

# =============================================================================
# 6. NMDS — site scores coloured by vegetation cover
# =============================================================================

nmdsDB <- metaMDS(matriz_site_sqrt, distance = "bray", k = 2, trymax = 100)

nmds_scores <- as.data.frame(scores(nmdsDB)$sites)
nmds_data   <- nmds_scores

nmds_data$sites <- c(rep("Secondary forest", 4),
                     rep("Cacao agroforestry system", 4),
                     rep("Pasture", 4))

nmds_data$sites <- factor(nmds_data$sites,
                          levels = c("Secondary forest",
                                     "Cacao agroforestry system",
                                     "Pasture"))

# Group centroids and spider segments
centroids <- nmds_data %>%
  group_by(sites) %>%
  summarise(across(c(NMDS1, NMDS2), mean), .groups = "drop")

segments <- left_join(nmds_data, centroids, by = "sites",
                      suffix = c(".x", ".y"))

# Custom theme
theme_custom <- theme_minimal(base_family = "serif") +
  theme(
    axis.title        = element_text(size = 18, face = "bold"),
    axis.text         = element_text(size = 14),
    legend.text       = element_text(size = 14),
    legend.title      = element_text(size = 18, face = "bold"),
    legend.position   = "top",
    axis.line         = element_line(colour = "black"),
    panel.grid        = element_blank(),
    axis.title.x      = element_text(margin = margin(t = 10)),
    axis.title.y      = element_text(margin = margin(t = 10)),
    panel.grid.major  = element_line(color = "gray90", linetype = "dashed")
  )

# NMDS plot
stress_label <- paste("Stress:", round(nmdsDB$stress, 3))

ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
  geom_segment(data = segments,
               aes(x = NMDS1.x, y = NMDS2.x,
                   xend = NMDS1.y, yend = NMDS2.y,
                   color = sites),
               linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(shape = sites, fill = sites),
             size = 5, color = "black", stroke = 1) +
  geom_point(data = centroids,
             aes(x = NMDS1, y = NMDS2, shape = sites),
             size = 5, color = "black", stroke = 1.5, fill = "white") +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  geom_mark_ellipse(expand = 0, aes(color = sites, fill = sites), alpha = 0.3) +
  scale_shape_manual(values = c("Secondary forest"          = 21,
                                "Cacao agroforestry system" = 22,
                                "Pasture"                   = 24)) +
  labs(x = "NMDS1", y = "NMDS2",
       color = "Vegetation cover",
       shape = "Vegetation cover",
       fill  = "Vegetation cover") +
  theme_custom +
  annotate("text",
           x = max(nmds_scores$NMDS1), y = max(nmds_scores$NMDS2),
           label = stress_label, fontface = "italic", size = 5, hjust = 1)

# =============================================================================
# 7. BETA DIVERSITY DECOMPOSITION — presence/absence, per vegetation cover
# =============================================================================

matriz_pa_site <- ifelse(matriz_site > 0, 1, 0)
rownames(matriz_pa_site) <- rownames(matriz_site)

beta_results <- metadata_site %>%
  group_split(vegcover) %>%
  lapply(function(df) {
    sub_mat <- matriz_pa_site[df$site_id, ]
    beta_m  <- beta.multi(sub_mat, index.family = "jaccard")
    data.frame(
      vegcover        = unique(df$vegcover),
      beta_total      = beta_m$beta.JAC,
      beta_turnover   = beta_m$beta.JTU,
      beta_nestedness = beta_m$beta.JNE
    )
  }) %>%
  do.call(rbind, .)

# --- Butterfly plot (turnover vs. nestedness) ---------------------------------
beta_long <- beta_results[, c("vegcover", "beta_turnover", "beta_nestedness")] %>%
  tidyr::pivot_longer(
    cols      = c(beta_turnover, beta_nestedness),
    names_to  = "component",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    value     = ifelse(component == "beta_nestedness", -value, value),
    component = dplyr::recode(component,
                              beta_turnover   = "β-JTU",
                              beta_nestedness = "β-JNE")
  )

ggplot(beta_long, aes(y = vegcover, x = value, fill = component)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("β-JTU" = "cyan4",
                               "β-JNE" = "cyan3")) +
  labs(y    = "Vegetation cover",
       x    = "Beta diversity (Jaccard)",
       fill = "Component") +
  theme(
    axis.title.x      = element_text(margin = margin(t = 10),
                                     family = "serif", size = 20, face = "bold"),
    axis.title.y      = element_text(margin = margin(r = 10),
                                     family = "serif", size = 20, face = "bold"),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.border      = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "gray95", linetype = "dashed"),
    axis.text         = element_text(family = "serif", size = 18),
    axis.line         = element_line(colour = "black"),
    legend.position   = "top",
    legend.text       = element_text(family = "serif", size = 18),
    legend.title      = element_text(family = "serif", size = 20, face = "bold"),
    strip.background  = element_blank()
  )

ggsave("Beta/Beta.svg", width = 12, height = 8, dpi = 600)
