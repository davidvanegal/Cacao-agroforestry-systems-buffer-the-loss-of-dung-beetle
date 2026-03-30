# =============================================================================
# Dung beetle analysis
# GLMMs by vegetation cover — general and by foraging strategy
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(dplyr)
library(tibble)
library(tidyr)
library(emmeans)
library(multcomp)
library(multcompView)
library(glmmTMB)
library(ggplot2)
library(FD)
library(viridis)
library(patchwork)
library(broom.mixed)
library(knitr)
library(openxlsx)
library(FSA)

# --- Data --------------------------------------------------------------------
dataDB <- read.csv("dataCASBufferDB.csv", header = TRUE, sep = ";")
dataDB$vegcover <- factor(dataDB$vegcover, levels = c("SF", "CAS", "PA"))

# =============================================================================
# 1. DATA PREPARATION
# =============================================================================

# General aggregation: one row per locality × vegcover × site
general_data <- dataDB %>%
  group_by(locality, vegcover, site) %>%
  summarise(
    richness        = n_distinct(species),
    abund_total     = sum(abund, na.rm = TRUE),
    biom_total      = sum(biom,  na.rm = TRUE),
    length_weighted = sum(length * biom, na.rm = TRUE) / sum(biom, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(site = as.integer(site))

# Full sampling design (locality × vegcover combinations present in the field)
design <- tibble::tribble(
  ~locality, ~vegcover,
  "EH", "CAS",
  "ET", "PA",
  "LP", "SF",
  "LP", "CAS",
  "LV", "SF",
  "LV", "PA",
  "LV", "CAS",
  "SE", "SF",
  "SE", "CAS",
  "SR", "SF",
  "SR", "PA",
  "SR", "CAS",
  "VN", "PA"
)

# Expand to 4 sites per combination; fill zeros for sites with no captures
general_data_complete <- design %>%
  crossing(site = 1:4) %>%
  mutate(site_id = paste(locality, site, sep = "_")) %>%
  left_join(general_data, by = c("locality", "vegcover", "site")) %>%
  mutate(
    across(c(richness, abund_total, biom_total), ~ replace_na(., 0)),
    length_weighted = replace_na(length_weighted, NA)
  )

# Strategy-level aggregation: one row per locality × vegcover × site × strategy
strategy_data <- dataDB %>%
  group_by(locality, vegcover, site, strategy) %>%
  summarise(
    richness        = n_distinct(species),
    abund_total     = sum(abund, na.rm = TRUE),
    biom_total      = sum(biom,  na.rm = TRUE),
    length_weighted = sum(length * biom, na.rm = TRUE) / sum(biom, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(site = as.integer(site))

# Expand with zeros per strategy
strategy_final <- design %>%
  crossing(site = 1:4) %>%
  crossing(strategy = unique(dataDB$strategy)) %>%
  mutate(site_id = paste(locality, site, sep = "_")) %>%
  left_join(strategy_data, by = c("locality", "vegcover", "site", "strategy")) %>%
  mutate(
    across(c(richness, abund_total, biom_total), ~ replace_na(., 0)),
    length_weighted = replace_na(length_weighted, NA)
  )

# Exclude Roller: detected in only 11.5% of sites
strategy_final_sub <- strategy_final %>% filter(strategy != "Roller")

# Body length: exclude NAs and zeros (no biological meaning without captures)
data_length          <- general_data_complete %>%
  filter(!is.na(length_weighted) & length_weighted > 0)

data_length_strategy <- strategy_final_sub %>%
  filter(!is.na(length_weighted) & length_weighted > 0)

# Vegetation cover factor: ordered from least to most disturbed
veg_levels <- c("SF", "CAS", "PA")

general_data_complete <- general_data_complete %>%
  mutate(vegcover = factor(vegcover, levels = veg_levels))

strategy_final_sub <- strategy_final_sub %>%
  mutate(vegcover = factor(vegcover, levels = veg_levels))

data_length <- data_length %>%
  mutate(vegcover = factor(vegcover, levels = veg_levels))

data_length_strategy <- data_length_strategy %>%
  mutate(vegcover = factor(vegcover, levels = veg_levels))

# =============================================================================
# 2. GENERAL MODELS
# =============================================================================

# Species richness — Poisson + (1 | locality)
model_richG <- glmmTMB(richness ~ vegcover + (1 | locality),
                       data = general_data_complete, family = poisson)
summary(model_richG)
drop1(model_richG, test = "Chisq")

# Total abundance — nbinom1 + (1 | locality)
model_abundG <- glmmTMB(abund_total ~ vegcover + (1 | locality),
                        data = general_data_complete, family = nbinom1)
summary(model_abundG)
drop1(model_abundG, test = "Chisq")

# Total biomass — Tweedie + (1 | locality)
model_biomG <- glmmTMB(biom_total ~ vegcover + (1 | locality),
                       data = general_data_complete,
                       family = tweedie(link = "log"))
summary(model_biomG)
drop1(model_biomG, test = "Chisq")

# Body length — Kruskal-Wallis (residuals not normally distributed)
kruskal_lengthG <- kruskal.test(length_weighted ~ vegcover, data = data_length)
kruskal_lengthG

dunnTest(length_weighted ~ vegcover, data = data_length, method = "bonferroni")
dunnTest(length_weighted ~ vegcover, data = data_length, method = "holm")
dunnTest(length_weighted ~ vegcover, data = data_length, method = "bh")

# =============================================================================
# 3. STRATEGY-LEVEL MODELS
# =============================================================================

# Species richness — additive model (interaction not significant)
model_rich_strategy <- glmmTMB(richness ~ vegcover + strategy + (1 | locality),
                               data = strategy_final_sub, family = poisson)

# Total abundance — additive nbinom1 model
model_abund_strategy <- glmmTMB(abund_total ~ vegcover + strategy + (1 | locality),
                                data = strategy_final_sub, family = nbinom1)

# Total biomass — additive Tweedie model
model_biom_strategy <- glmmTMB(biom_total ~ vegcover + strategy + (1 | locality),
                               data = strategy_final_sub,
                               family = tweedie(link = "log"))

# Body length — additive Gamma model (neither term significant)
model_length_strategy <- glmmTMB(length_weighted ~ vegcover + strategy + (1 | locality),
                                 data = data_length_strategy,
                                 family = Gamma(link = "log"))

# =============================================================================
# 4. POST-HOC — pairwise contrasts (Bonferroni)
# =============================================================================

# General models
pairs(emmeans(model_richG,  ~ vegcover), adjust = "bonferroni")
pairs(emmeans(model_abundG, ~ vegcover), adjust = "bonferroni")
pairs(emmeans(model_biomG,  ~ vegcover), adjust = "bonferroni")

# Strategy models — vegcover and strategy separately
pairs(emmeans(model_rich_strategy,  ~ vegcover),  adjust = "bonferroni")
pairs(emmeans(model_rich_strategy,  ~ strategy),  adjust = "bonferroni")
pairs(emmeans(model_biom_strategy,  ~ vegcover),  adjust = "bonferroni")
pairs(emmeans(model_biom_strategy,  ~ strategy),  adjust = "bonferroni")
pairs(emmeans(model_abund_strategy, ~ vegcover),  adjust = "bonferroni")
pairs(emmeans(model_abund_strategy, ~ strategy),  adjust = "bonferroni")

# Marginal means on the response scale
emmeans(model_richG,  ~ vegcover, type = "response")
emmeans(model_abundG, ~ vegcover, type = "response")
emmeans(model_biomG,  ~ vegcover, type = "response")

# =============================================================================
# 5. POST-HOC — CLD letters (for figures)
# =============================================================================

# General models
cld_richG  <- cld(emmeans(model_richG,  ~ vegcover),
                  Letters = letters, adjust = "bonferroni")
cld_abundG <- cld(emmeans(model_abundG, ~ vegcover),
                  Letters = letters, adjust = "bonferroni")
cld_biomG  <- cld(emmeans(model_biomG,  ~ vegcover),
                  Letters = letters, adjust = "bonferroni")

# Strategy models
cld_rich_veg      <- cld(emmeans(model_rich_strategy,  ~ vegcover),
                         Letters = letters, adjust = "bonferroni")
cld_abund_veg     <- cld(emmeans(model_abund_strategy, ~ vegcover),
                         Letters = letters, adjust = "bonferroni")
cld_abund_strategy <- cld(emmeans(model_abund_strategy, ~ strategy),
                          Letters = letters, adjust = "bonferroni")
cld_biom_veg      <- cld(emmeans(model_biom_strategy,  ~ vegcover),
                         Letters = letters, adjust = "bonferroni")

# =============================================================================
# 6. FIGURES
# =============================================================================

# Colour palette (consistent with other scripts)
cols <- setNames(viridis(3, option = "D"),
                 c("SF", "CAS", "PA"))

# Upper whisker: robust CLD label position (avoids outlier distortion)
upper_whisker <- function(x) quantile(x, 0.75) + 0.5 * IQR(x)

# Full vegetation cover names for x-axis labels
veg_labels <- c("SF"  = "Secondary forest",
                "CAS" = "Cacao agroforestry system",
                "PA"  = "Pasture")

# Base theme
theme_base <- theme_classic() +
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
    legend.position   = "none",
    legend.text       = element_text(family = "serif", size = 18),
    legend.title      = element_text(family = "serif", size = 20, face = "bold"),
    strip.background  = element_blank()
  )

# Strategy theme: adds facet strip styling
theme_strategy <- theme_base +
  theme(strip.background = element_rect(fill = "grey90", color = NA),
        strip.text       = element_text(size = 10))

# --- General figures ---------------------------------------------------------

pos_richG  <- general_data_complete %>%
  group_by(vegcover) %>% summarise(ymax = upper_whisker(richness),    
                                   .groups = "drop")
pos_abundG <- general_data_complete %>%
  group_by(vegcover) %>% summarise(ymax = upper_whisker(abund_total), 
                                   .groups = "drop")
pos_biomG  <- general_data_complete %>%
  group_by(vegcover) %>% summarise(ymax = upper_whisker(biom_total),  
                                   .groups = "drop")

cld_richG_plot  <- as.data.frame(cld_richG)  %>%
  left_join(pos_richG,  by = "vegcover") %>% mutate(.group = trimws(.group))
cld_abundG_plot <- as.data.frame(cld_abundG) %>%
  left_join(pos_abundG, by = "vegcover") %>% mutate(.group = trimws(.group))
cld_biomG_plot  <- as.data.frame(cld_biomG)  %>%
  left_join(pos_biomG,  by = "vegcover") %>% mutate(.group = trimws(.group))

p_richG <- ggplot(general_data_complete,
                  aes(x = vegcover, y = richness, fill = vegcover)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_text(data = cld_richG_plot,
            aes(x = vegcover, y = ymax * 1.08, label = .group),
            inherit.aes = FALSE, family = "serif", size = 8) +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = "Species richness") +
  theme_base +
  theme(axis.text.x = element_blank())

p_abundG <- ggplot(general_data_complete,
                   aes(x = vegcover, y = abund_total, fill = vegcover)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_text(data = cld_abundG_plot,
            aes(x = vegcover, y = ymax * 1.08, label = .group),
            inherit.aes = FALSE, family = "serif", size = 8) +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = "Total abundance") +
  theme_base +
  theme(axis.text.x = element_blank())

p_biomG <- ggplot(general_data_complete,
                  aes(x = vegcover, y = biom_total, fill = vegcover)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_text(data = cld_biomG_plot,
            aes(x = vegcover, y = ymax * 1.08, label = .group),
            inherit.aes = FALSE, family = "serif", size = 8) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(labels = veg_labels) +
  labs(x = NULL, y = "Total biomass (g)") +
  theme_base

fig_general <- (p_richG / p_abundG / p_biomG) +
  plot_annotation(
    caption = "Vegetation cover",
    theme   = theme(plot.caption = element_text(hjust  = 0.5,
                                                family = "serif",
                                                size   = 20, face = "bold"))
  )

ggsave("GLM/Fig_General.svg", fig_general, width = 12, height = 12, dpi = 600)

# --- Strategy figures --------------------------------------------------------

pos_rich_strategy  <- strategy_final_sub %>%
  group_by(vegcover, strategy) %>%
  summarise(ymax = upper_whisker(richness),    .groups = "drop")
pos_abund_strategy <- strategy_final_sub %>%
  group_by(vegcover, strategy) %>%
  summarise(ymax = upper_whisker(abund_total), .groups = "drop")
pos_biom_strategy  <- strategy_final_sub %>%
  group_by(vegcover, strategy) %>%
  summarise(ymax = upper_whisker(biom_total),  .groups = "drop")

cld_rich_veg_plot <- as.data.frame(cld_rich_veg) %>%
  left_join(pos_rich_strategy %>%
              group_by(vegcover) %>%
              summarise(ymax = max(ymax), .groups = "drop"),
            by = "vegcover") %>%
  mutate(.group = trimws(.group))

cld_abund_veg_plot <- as.data.frame(cld_abund_veg) %>%
  left_join(pos_abund_strategy %>%
              group_by(vegcover) %>%
              summarise(ymax = max(ymax), .groups = "drop"),
            by = "vegcover") %>%
  mutate(.group = trimws(.group))

cld_biom_veg_plot <- as.data.frame(cld_biom_veg) %>%
  left_join(pos_biom_strategy %>%
              group_by(vegcover) %>%
              summarise(ymax = max(ymax), .groups = "drop"),
            by = "vegcover") %>%
  mutate(.group = trimws(.group))

p_rich_strategy <- ggplot(strategy_final_sub,
                          aes(x = vegcover, y = richness, fill = vegcover)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_text(data = cld_rich_veg_plot,
            aes(x = vegcover, y = ymax * 0.3, label = .group),
            inherit.aes = FALSE, family = "serif", size = 8) +
  scale_fill_manual(values = cols) +
  facet_wrap(~ strategy, scales = "free_y") +
  labs(x = NULL, y = "Species richness") +
  theme_strategy +
  theme(strip.text = element_blank(), axis.text.x = element_blank())

p_abund_strategy <- ggplot(strategy_final_sub,
                           aes(x = vegcover, y = abund_total, fill = vegcover)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_text(data = cld_abund_veg_plot,
            aes(x = vegcover, y = ymax * 0.2, label = .group),
            inherit.aes = FALSE, family = "serif", size = 8) +
  scale_fill_manual(values = cols) +
  facet_wrap(~ strategy, scales = "free_y") +
  labs(x = NULL, y = "Total abundance") +
  theme_strategy +
  theme(strip.text = element_blank(), axis.text.x = element_blank())

p_biom_strategy <- ggplot(strategy_final_sub,
                          aes(x = vegcover, y = biom_total, fill = vegcover)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_text(data = cld_biom_veg_plot,
            aes(x = vegcover, y = ymax * 0.2, label = .group),
            inherit.aes = FALSE, family = "serif", size = 8) +
  scale_fill_manual(values = cols) +
  facet_wrap(~ strategy, scales = "free_y") +
  scale_x_discrete(labels = veg_labels) +
  labs(x = NULL, y = "Total biomass (g)") +
  theme_strategy +
  theme(strip.text = element_blank())

fig_strategy <- (p_rich_strategy / p_abund_strategy / p_biom_strategy) +
  plot_annotation(
    caption = "Vegetation cover",
    theme   = theme(plot.caption = element_text(hjust  = 0.5,
                                                family = "serif",
                                                size   = 20, face = "bold"))
  )

ggsave("GLM/Fig_strategy.svg", fig_strategy, width = 17, height = 12, dpi = 600)

# =============================================================================
# 7. EXPORT RESULTS TO EXCEL
# =============================================================================

wb <- createWorkbook()

# Helper: write a titled block (title row + data) starting at a given row
write_block <- function(wb, sheet_name, title, df, start_row) {
  title_style <- createStyle(
    textDecoration = "bold", fontSize = 12,
    fontColour = "#FFFFFF", fgFill = "#4472C4",
    halign = "left"
  )
  header_style <- createStyle(textDecoration = "bold", border = "Bottom")

  writeData(wb, sheet_name, title, startRow = start_row, startCol = 1)
  addStyle(wb, sheet_name, title_style,
           rows = start_row, cols = 1:ncol(df), stack = TRUE)
  mergeCells(wb, sheet_name, cols = 1:ncol(df), rows = start_row)

  data_start <- start_row + 1
  writeData(wb, sheet_name, df,
            startRow = data_start, startCol = 1,
            headerStyle = header_style)
  setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")

  data_start + nrow(df) + 3
}

# --- Sheet 1: Summary — tidy model coefficients ------------------------------
addWorksheet(wb, "Summary")
row <- 1
summaries <- list(
  "Abundance General"     = model_abundG,
  "Richness General"      = model_richG,
  "q1"                    = model_q1,
  "q2"                    = model_q2,
  "Biomass General"       = model_biomG,
  "Richness × Strategy"   = model_rich_strategy,
  "Abundance × Strategy"  = model_abund_strategy,
  "Biomass × Strategy"    = model_biom_strategy,
  "Length × Strategy"     = model_length_strategy,
  "FDis"                  = model_FDis,
  "RaoQ"                  = model_RaoQ
)
for (nm in names(summaries)) {
  row <- write_block(wb, "Summary",
                     paste("Summary —", nm),
                     as.data.frame(tidy(summaries[[nm]])),
                     row)
}

# --- Sheet 2: Drop1 — global Chi-sq tests ------------------------------------
addWorksheet(wb, "Drop1")
row <- 1
for (nm in names(summaries)) {
  row <- write_block(wb, "Drop1",
                     paste("Drop1 —", nm),
                     as.data.frame(broom::tidy(
                       drop1(summaries[[nm]], test = "Chisq"))),
                     row)
}

# --- Sheet 3: Post-hoc — pairwise Bonferroni + Kruskal-Wallis ----------------
addWorksheet(wb, "Post-hoc")
row <- 1

posthoc_veg <- list(
  "Abundance General" = model_abundG,
  "Richness General"  = model_richG,
  "q1"                = model_q1,
  "q2"                = model_q2,
  "Biomass General"   = model_biomG,
  "FDis"              = model_FDis,
  "RaoQ"              = model_RaoQ
)
for (nm in names(posthoc_veg)) {
  row <- write_block(wb, "Post-hoc",
                     paste("Post-hoc vegcover —", nm),
                     as.data.frame(pairs(
                       emmeans(posthoc_veg[[nm]], ~ vegcover),
                       adjust = "bonferroni")),
                     row)
}

posthoc_both <- list(
  "Richness × Strategy"  = model_rich_strategy,
  "Abundance × Strategy" = model_abund_strategy,
  "Biomass × Strategy"   = model_biom_strategy,
  "Length × Strategy"    = model_length_strategy
)
for (nm in names(posthoc_both)) {
  row <- write_block(wb, "Post-hoc",
                     paste("Post-hoc vegcover —", nm),
                     as.data.frame(pairs(
                       emmeans(posthoc_both[[nm]], ~ vegcover),
                       adjust = "bonferroni")),
                     row)
  row <- write_block(wb, "Post-hoc",
                     paste("Post-hoc strategy —", nm),
                     as.data.frame(pairs(
                       emmeans(posthoc_both[[nm]], ~ strategy),
                       adjust = "bonferroni")),
                     row)
}

# Kruskal-Wallis for body length (no GLMM)
row <- write_block(wb, "Post-hoc",
                   "Kruskal-Wallis — Length General",
                   data.frame(
                     test      = "Kruskal-Wallis",
                     variable  = "length_weighted",
                     statistic = kruskal_lengthG$statistic,
                     df        = kruskal_lengthG$parameter,
                     p.value   = kruskal_lengthG$p.value),
                   row)

# --- Sheet 4: CLD — group letters for figures --------------------------------
addWorksheet(wb, "CLD")
row <- 1
clds <- list(
  "Abundance General — vegcover"    = cld_abundG,
  "q1 — vegcover"                   = cld_q1,
  "q2 — vegcover"                   = cld_q2,
  "Biomass General — vegcover"      = cld_biomG,
  "Richness × Strategy — vegcover"  = cld_rich_veg,
  "Abundance × Strategy — vegcover" = cld_abund_veg,
  "Abundance × Strategy — strategy" = cld_abund_strategy,
  "Biomass × Strategy — vegcover"   = cld_biom_veg,
  "FDis — vegcover"                 = cld_FDis,
  "RaoQ — vegcover"                 = cld_RaoQ
)
for (nm in names(clds)) {
  row <- write_block(wb, "CLD",
                     paste("CLD —", nm),
                     as.data.frame(clds[[nm]]),
                     row)
}

saveWorkbook(wb, "GLMM/DungBeetlesGLMM.xlsx", overwrite = TRUE)
