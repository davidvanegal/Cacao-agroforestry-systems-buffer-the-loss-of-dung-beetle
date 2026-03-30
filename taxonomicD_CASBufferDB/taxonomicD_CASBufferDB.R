# =============================================================================
# Hill numbers (diversity orders q = 0, 1, 2) — Dung beetles
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(iNEXT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(multcompView)
library(multcomp)
library(patchwork)

# --- Data --------------------------------------------------------------------
dataDB <- read.csv("dataCASBufferDB.csv", header = TRUE, sep = ";")
dataDB$vegcover <- factor(dataDB$vegcover,
                          levels = c("SF", "CAS", "PA"),
                          labels = c("Secondary forest",
                                     "Cacao agroforestry system",
                                     "Pasture"))

# =============================================================================
# 1. COMMUNITY MATRIX — aggregated at site level
# One row per locality × vegcover × site
# =============================================================================

comm_site <- dataDB %>%
  group_by(locality, vegcover, site, species) %>%
  summarise(abund = sum(abund), .groups = "drop") %>%
  mutate(siteID = paste(locality, vegcover, site, sep = "_")) %>%
  pivot_wider(names_from = species, values_from = abund, values_fill = 0)

# Site metadata
metadata <- comm_site %>%
  distinct(siteID, locality, vegcover)

# Species-only community matrix
comm_mat2 <- as.data.frame(comm_site[, setdiff(colnames(comm_site),
                                               c("locality", "vegcover",
                                                 "site", "siteID"))])
rownames(comm_mat2) <- comm_site$siteID

# =============================================================================
# 2. HILL NUMBERS PER SITE
# q = 0: species richness; q = 1: Shannon diversity; q = 2: Simpson diversity
# =============================================================================

hill_results <- do.call(rbind, lapply(rownames(comm_mat2), function(sid) {
  x  <- as.numeric(comm_mat2[sid, ])
  N  <- sum(x)
  p  <- x[x > 0] / N

  q0 <- sum(x > 0)                   # species richness
  q1 <- exp(-sum(p * log(p)))        # Shannon entropy exponentiated
  q2 <- 1 / sum(p^2)                 # inverse Simpson

  data.frame(siteID = sid, q0 = q0, q1 = q1, q2 = q2)
})) %>%
  left_join(metadata, by = "siteID")

print(hill_results)

# =============================================================================
# 3. GLMMs — locality as random intercept
# =============================================================================

# q0: Poisson (count data)
model_q0 <- glmmTMB(q0 ~ vegcover + (1 | locality),
                    data = hill_results, family = poisson)

# q1: Gaussian (continuous, approximately log-normal)
model_q1 <- glmmTMB(q1 ~ vegcover + (1 | locality),
                    data = hill_results, family = gaussian)

# q2: Gaussian
model_q2 <- glmmTMB(q2 ~ vegcover + (1 | locality),
                    data = hill_results, family = gaussian)

# =============================================================================
# 4. POST-HOC — CLD letters (Bonferroni)
# =============================================================================

cld_q0 <- cld(emmeans(model_q0, ~ vegcover), Letters = letters, 
              adjust = "bonferroni")
cld_q1 <- cld(emmeans(model_q1, ~ vegcover), Letters = letters, 
              adjust = "bonferroni")
cld_q2 <- cld(emmeans(model_q2, ~ vegcover), Letters = letters, 
              adjust = "bonferroni")

# =============================================================================
# 5. FIGURE — three stacked panels (q0, q1, q2)
# =============================================================================

cols <- c("Secondary forest"          = "forestgreen",
          "Cacao agroforestry system" = "darkorange4",
          "Pasture"                   = "springgreen")

upper_whisker <- function(x) quantile(x, 0.75, na.rm = TRUE) + 
  1.5 * IQR(x, na.rm = TRUE)

theme_hill <- theme_classic() +
  theme(
    axis.title.x     = element_text(margin = margin(t = 10),
                                    family = "serif", size = 20, face = "bold"),
    axis.title.y     = element_text(margin = margin(r = 10),
                                    family = "serif", size = 20, face = "bold"),
    axis.text        = element_text(family = "serif", size = 18),
    axis.line        = element_line(colour = "black"),
    panel.grid.major = element_line(color = "gray95", linetype = "dashed"),
    legend.position  = "none"
  )

make_hill_panel <- function(data, var, cld_obj, ylab, show_x = FALSE) {
  pos <- data %>%
    group_by(vegcover) %>%
    summarise(ymax = upper_whisker(.data[[var]]), .groups = "drop")

  cld_plot <- as.data.frame(cld_obj) %>%
    left_join(pos, by = "vegcover") %>%
    mutate(.group = trimws(.group))

  p <- ggplot(data, aes(x = vegcover, y = .data[[var]], fill = vegcover)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_text(data = cld_plot,
              aes(x = vegcover, y = ymax * 1.08, label = .group),
              inherit.aes = FALSE, family = "serif", size = 6) +
    scale_fill_manual(values = cols) +
    labs(x = NULL, y = ylab) +
    theme_hill

  if (!show_x) p <- p + theme(axis.text.x = element_blank())
  p
}

p_q0 <- make_hill_panel(hill_results, "q0", cld_q0,
                        expression(bold(italic("q") * " = 0  (Species richness)")))
p_q1 <- make_hill_panel(hill_results, "q1", cld_q1,
                        expression(bold(italic("q") * " = 1  (Shannon diversity)")))
p_q2 <- make_hill_panel(hill_results, "q2", cld_q2,
                        expression(bold(italic("q") * " = 2  (Simpson diversity)")),
                        show_x = TRUE)

fig_hill <- (p_q0 / p_q1 / p_q2) +
  plot_annotation(
    tag_levels = "a",
    caption    = "Vegetation cover",
    theme      = theme(plot.caption = element_text(hjust   = 0.5,
                                                   family  = "serif",
                                                   size    = 20,
                                                   face    = "bold"))
  )

ggsave("Diversity/Fig_Hill.png", fig_hill, width = 10, height = 12, dpi = 600)
ggsave("Diversity/Fig_Hill.svg", fig_hill, width =  9, height = 14, dpi = 600)

# =============================================================================
# 6. DIVERSITY PROFILE — Hill numbers q = 0 to 2
# =============================================================================

q_vals <- seq(0, 2, by = 0.1)

hill_profile <- do.call(rbind, lapply(rownames(comm_mat2), function(sid) {
  x <- as.numeric(comm_mat2[sid, ])
  N <- sum(x)
  p <- x[x > 0] / N

  do.call(rbind, lapply(q_vals, function(q) {
    D <- if (q == 1) exp(-sum(p * log(p))) else sum(p^q)^(1 / (1 - q))
    data.frame(siteID = sid, q = q, D = D)
  }))
})) %>%
  left_join(metadata, by = "siteID")

# Mean ± SE per vegcover × q (q = 0, 1, 2 only)
hill_profile_summary <- hill_profile %>%
  filter(q %in% c(0, 1, 2)) %>%
  mutate(q_label = factor(q,
                          levels = c(0, 1, 2),
                          labels = c("q0", "q1", "q2"))) %>%
  group_by(vegcover, q_label) %>%
  summarise(
    mean_D = mean(D, na.rm = TRUE),
    se_D   = sd(D,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_profile <- ggplot(hill_profile_summary,
                    aes(x = vegcover, y = mean_D,
                        color = q_label, fill = q_label,
                        group = q_label, shape = q_label)) +
  geom_errorbar(aes(ymin = mean_D - se_D, ymax = mean_D + se_D),
                width = 0.15, linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4, stroke = 1.2) +
  scale_color_manual(values = c("q0" = "#3B4DA8",
                                "q1" = "#2A9D8F",
                                "q2" = "#8AB733"),
                     name = "Diversity order (q)") +
  scale_fill_manual(values  = c("q0" = "#3B4DA8",
                                "q1" = "#2A9D8F",
                                "q2" = "#8AB733"),
                    name = "Diversity order (q)") +
  scale_shape_manual(values = c("q0" = 19, "q1" = 15, "q2" = 17),
                     name = "Diversity order (q)") +
  labs(x = "Vegetation cover",
       y = "Effective number of species") +
  theme_classic() +
  theme(
    axis.title.x     = element_text(margin = margin(t = 10),
                                    family = "serif", size = 20, face = "bold"),
    axis.title.y     = element_text(margin = margin(r = 10),
                                    family = "serif", size = 20, face = "bold"),
    axis.text        = element_text(family = "serif", size = 16),
    axis.line        = element_line(colour = "black"),
    panel.grid.major = element_line(color = "gray95", linetype = "dashed"),
    legend.position  = "top",
    legend.text      = element_text(family = "serif", size = 14),
    legend.title     = element_text(family = "serif", size = 16, face = "bold")
  )

ggsave("Diversity/Fig_DiversityProfile.svg", p_profile,
       width = 10, height = 7, dpi = 600)
