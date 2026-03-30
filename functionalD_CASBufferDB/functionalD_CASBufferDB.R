# =============================================================================
# Functional diversity analysis — Dung beetles
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(FD)
library(viridis)
library(glmmTMB)
library(emmeans)
library(multcomp)
library(multcompView)
library(ggplot2)
library(patchwork)
library(scales)

# --- Data --------------------------------------------------------------------
dataDB <- read.csv("dataCASBufferDB.csv", header = TRUE, sep = ";")
dataDB$vegcover <- factor(dataDB$vegcover,
                          levels = c("SF", "CAS", "PA"),
                          labels = c("Secondary forest",
                                     "Cacao agroforestry system",
                                     "Pasture"))

# =============================================================================
# 1. COMMUNITY MATRIX — aggregated at site level
# One row per locality × vegcover × site (traps are pseudoreplicates)
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
comm_mat <- as.data.frame(comm_site[, setdiff(colnames(comm_site),
                                              c("locality", "vegcover",
                                                "site", "siteID"))])
rownames(comm_mat) <- comm_site$siteID

# =============================================================================
# 2. TRAIT MATRIX — mean values per species
# =============================================================================

# Collapse individual-level variation to species means
traits <- dataDB %>%
  group_by(species) %>%
  summarise(
    length = mean(length, na.rm = TRUE),
    biom   = mean(biom,   na.rm = TRUE),
    guild  = first(guild),
    .groups = "drop"
  ) %>%
  column_to_rownames("species")

# Retain only species present in both matrices
common_sp <- intersect(colnames(comm_mat), rownames(traits))
comm_mat2 <- as.matrix(comm_mat[, common_sp])
mode(comm_mat2) <- "numeric"
traits2   <- traits[common_sp, ]
traits2$guild <- as.factor(traits2$guild)

# =============================================================================
# 3. FUNCTIONAL DIVERSITY INDICES
# =============================================================================

fd_indices <- dbFD(traits2, comm_mat2, calc.FRic = TRUE)

# Note: FRic and FEve excluded — 47% NA in Pasture (insufficient species
# to compute functional volume at low-capture sites)
fd_results <- data.frame(
  siteID = rownames(comm_mat2),
  FRic   = fd_indices$FRic,
  FEve   = fd_indices$FEve,
  FDiv   = fd_indices$FDiv,
  FDis   = fd_indices$FDis,
  FRaoQ  = fd_indices$RaoQ
) %>%
  left_join(metadata, by = "siteID")

print(fd_results)

# =============================================================================
# 4. GLMMs — locality as random intercept
# =============================================================================

model_FDis <- glmmTMB(FDis  ~ vegcover + (1 | locality),
                      data = fd_results, family = gaussian)
summary(model_FDis)
drop1(model_FDis, test = "Chisq")

model_RaoQ <- glmmTMB(FRaoQ ~ vegcover + (1 | locality),
                      data = fd_results, family = gaussian)
summary(model_RaoQ)
drop1(model_RaoQ, test = "Chisq")

# =============================================================================
# 5. POST-HOC — pairwise contrasts and CLD letters (Bonferroni)
# =============================================================================

pairs(emmeans(model_FDis, ~ vegcover), adjust = "bonferroni")
cld_FDis <- cld(emmeans(model_FDis, ~ vegcover),
                Letters = letters, adjust = "bonferroni")

cld_RaoQ <- cld(emmeans(model_RaoQ, ~ vegcover),
                Letters = letters, adjust = "bonferroni")

# =============================================================================
# 6. FIGURE — FDis and FRaoQ boxplots
# =============================================================================

cols <- setNames(viridis(3, option = "D"),
                 c("Secondary forest", "Cacao agroforestry system", "Pasture"))

upper_whisker <- function(x) quantile(x, 0.75, na.rm = TRUE) + 0.5 * IQR(x, na.rm = TRUE)

theme_fd <- theme_classic() +
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

# FDis panel
pos_FDis <- fd_results %>%
  group_by(vegcover) %>%
  summarise(ymax = upper_whisker(FDis), .groups = "drop")

cld_FDis_plot <- as.data.frame(cld_FDis) %>%
  left_join(pos_FDis, by = "vegcover") %>%
  mutate(.group = trimws(.group))

p_FDis <- ggplot(fd_results, aes(x = vegcover, y = FDis, fill = vegcover)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_text(data = cld_FDis_plot,
            aes(x = vegcover, y = ymax * 1.08, label = .group),
            inherit.aes = FALSE, family = "serif", size = 8) +
  scale_fill_manual(values = cols) +
  labs(x = "Vegetation cover", y = "FDis") +
  theme_fd

# FRaoQ panel
pos_RaoQ <- fd_results %>%
  group_by(vegcover) %>%
  summarise(ymax = upper_whisker(FRaoQ), .groups = "drop")

cld_RaoQ_plot <- as.data.frame(cld_RaoQ) %>%
  left_join(pos_RaoQ, by = "vegcover") %>%
  mutate(.group = trimws(.group))

p_RaoQ <- ggplot(fd_results, aes(x = vegcover, y = FRaoQ, fill = vegcover)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_text(data = cld_RaoQ_plot,
            aes(x = vegcover, y = ymax * 1.08, label = .group),
            inherit.aes = FALSE, family = "serif", size = 8) +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = "FRaoQ") +
  theme_fd

# Combined figure (FDis on top, x labels hidden)
p_FDis_top <- p_FDis +
  theme(axis.text.x  = element_blank(),
        axis.title.x = element_blank())

fig_fd <- (p_FDis_top / p_RaoQ) +
  plot_annotation(
    caption = "Vegetation cover",
    theme   = theme(plot.caption = element_text(hjust  = 0.5,
                                                family = "serif",
                                                size   = 20,
                                                face   = "bold"))
  )

ggsave("FD/Fig_FD.svg",   fig_fd,  width = 10, height = 12, dpi = 600)
ggsave("FD/Fig_FDis.svg", p_FDis,  width = 10, height =  7, dpi = 600)

# =============================================================================
# 7. COMMUNITY WEIGHTED MEAN (CWM)
# =============================================================================

cwm_df <- as.data.frame(fd_indices$CWM) %>%
  rownames_to_column("siteID") %>%
  left_join(metadata, by = "siteID")

# --- 7a. Guild proportions (categorical trait) — stacked bar -----------------

# Build dummy matrix for guild levels
guild_dummies <- model.matrix(~ guild - 1, data = traits2) %>%
  as.data.frame()

colnames(guild_dummies)   # verify column names

# Relative abundance per site; CWM = rel. abund. × guild dummy
comm_rel       <- comm_mat2 / rowSums(comm_mat2)
cwm_guild_manual <- as.data.frame(
  as.matrix(comm_rel) %*% as.matrix(guild_dummies)
) %>%
  rownames_to_column("siteID") %>%
  left_join(metadata, by = "siteID")

guild_cols <- grep("^guild", colnames(cwm_guild_manual), value = TRUE)

cwm_guild_long <- cwm_guild_manual %>%
  dplyr::select(siteID, vegcover, all_of(guild_cols)) %>%
  pivot_longer(cols      = all_of(guild_cols),
               names_to  = "guild",
               values_to = "cwm_prop") %>%
  mutate(guild = sub("^guild", "", guild))   # remove "guild" prefix

guild_pal <- c("Dweller"  = "#4e9a8c",
               "Roller"   = "#e8a838",
               "Tunneler" = "#7b4fa6")

cwm_guild_summary <- cwm_guild_long %>%
  dplyr::group_by(vegcover, guild) %>%
  dplyr::summarise(
    mean_prop = mean(cwm_prop, na.rm = TRUE),
    se_prop   = sd(cwm_prop,   na.rm = TRUE) / sqrt(n()),
    .groups   = "drop"
  )

p_cwm_guild <- ggplot(cwm_guild_summary,
                      aes(x = vegcover, y = mean_prop, fill = guild)) +
  geom_bar(stat = "identity", position = "stack",
           colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = guild_pal, name = "Functional guild") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Vegetation cover",
       y = "CWM — guild proportion") +
  theme_classic() +
  theme(
    axis.title        = element_text(family = "serif", size = 16, face = "bold"),
    axis.text         = element_text(family = "serif", size = 14),
    axis.text.x       = element_text(angle = 15, hjust = 1),
    legend.title      = element_text(family = "serif", size = 14, face = "bold"),
    legend.text       = element_text(family = "serif", size = 13),
    panel.grid.major.y = element_line(color = "gray92", linetype = "dashed")
  )

ggsave("FD/Fig_CWM_guild.svg", p_cwm_guild, width = 8, height = 6, dpi = 600)

# --- 7b. Summary table — CWM for continuous traits (length, biom) ------------
cwm_cont <- cwm_df %>%
  dplyr::select(siteID, vegcover, length, biom) %>%
  pivot_longer(cols      = c(length, biom),
               names_to  = "trait",
               values_to = "cwm_val") %>%
  group_by(vegcover, trait) %>%
  summarise(mean = mean(cwm_val, na.rm = TRUE),
            se   = sd(cwm_val,   na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

print(cwm_cont)
