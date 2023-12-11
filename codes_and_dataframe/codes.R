## install libraries
library(tidyverse)
library(janitor)
library(viridis)
library(here)

# import .tsv file
tsv <- read_tsv(here("df.tsv")) 

# make a data frame which only contain genus
tsv_g <- tsv %>%
  clean_names() %>%
  select(-phylum, -class, -order) 

# eliminate "%" from value
tsv_g = tsv_g %>%
  mutate_all(~ str_replace_all(., "%", ""))

# make a longer data = combine 0.22um, 5um, 20um, 2 seines into "filter_size" column
tsv_long <- tsv_g %>%
  clean_names() %>%
  select(genus, seine_river, seine_river_2, wk20um, wk5um, wk0_2um) %>%
  rename("wk_0.2μm" = wk0_2um,
         "wk_5μm" = wk5um,
         "wk_20μm" =  wk20um) %>%
  pivot_longer(cols = -c(genus), # don't include "genus" into "filter_size" column
               names_to = "filter_size",
               values_to = "ratio")

# make stacked bar plots
sb <- tsv_long %>%
  clean_names() %>%
  mutate(tsv_long, ratio = as.numeric(ratio)) %>%
  ggplot() +
  geom_col(aes(x = filter_size, y = ratio, fill = genus), position = "fill") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), 
                     labels = scales::percent) + 
  scale_x_discrete(expand = expansion(mult = c(0, 0)),
                   limits = c("seine_river", "seine_river_2", "wk_20μm", "wk_5μm", "wk_0.2μm")) +
  scale_fill_viridis_d(option = "D") +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.x  = element_text(angle = 90, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 0))

print(sb)

# calculate Shannon, Simpson and InvSimpson diversity index (alpha diversity)
library(vegan)

# make wide data for calculating diversity indices
tsv_t <- t(tsv_g) 
tsv_t <- as.data.frame(tsv_t) 
colnames(tsv_t) <- tsv_g$genus
tsv_t <- tsv_t %>%
  filter(row_number() > 1)
fil_vec <- unique(tsv_long$filter_size)
rownames(tsv_t) <- fil_vec
tsv_num <- tsv_t %>%
  mutate_all(as.numeric)

# calculate Shannon index
shannon <- round(diversity(tsv_num, index = "shannon"), 2)
print(shannon)

# draw rarefaction curves
rarecurve <- rarecurve(round(tsv_num*200), 
                       label = FALSE,
                       col = viridis(5), 
                       lwd = 3) 

axis(1, at = seq(0, 20000, by = 2500))

legend("bottomright",
       legend = rownames(tsv_num), 
       title = "Samples",
       lty = 1, 
       lwd = 3, 
       cex = 1,
       col = viridis(5),
       text.width = 600, 
       inset = c(0.01, 0.01),
       bg = "transparent")

# calculate beta diversity with distance matrix
beta_dis <- vegdist(tsv_num, method = "bray")

# calculate dissimilarities and conduct clustering analysis (beta diversity)
library(gplots) 

tsv_heat <- tsv_g %>%
  select(-genus) 

tsv_heat <- as.data.frame(tsv_heat) 
fil_vec2 <- unique(tsv_long$genus) 
rownames(tsv_heat) <- fil_vec2
colnames(tsv_heat) <- c("seine_river", "seine_river_2", "wk_20μm", "wk_5μm", "wk_0.2μm")
tsv_heat <- tsv_heat %>%
  mutate_all(as.numeric) 

heat2 <- heatmap.2(as.matrix(tsv_heat), 
                   distfun=function(x) vegdist(x, method="bray"),
                   hclustfun=function(x) hclust(x, method="ward.D2"),
                   ColSideColors = viridis(5), 
                   scale="none",
                   trace="none",
                   margins = c(10,9), 
                   Rowv = FALSE, 
                   
                   col = viridis_pal(), 
                   breaks = c(0:15, seq(25, max(tsv_num)+10, 10)),
                   cexRow = 1
                   )

# calculate species richness in each sample
richness <- specnumber(tsv_num)
print(richness)

# calculate Pielou’s Evenness by shannon / log(richness)
evenness <- round(shannon/log(richness), 2) 
print(evenness)

# make scatter plots for species richness, species evenness and Shannon index
alpha_mat <- cbind(shannon = shannon, richness = richness, pielou = evenness,
                   filter_size = c("seine_river", "seine_river_2", 
                                   "wk_20μm", "wk_5μm", "wk_0.2μm"), 
                   )
rownames(alpha_mat) <- c(1,2,3,4,5)

alpha_df <- alpha_mat %>%
  as.data.frame(alpha_mat) %>%
  select(numbering, filter_size, richness, pielou, shannon) %>%
  mutate(shannon = as.numeric(shannon)) %>% 
  mutate(richness = as.numeric(richness)) %>%
  mutate(pielou = as.numeric(pielou))

# make plots for alpha diversity
plot_shannon <- alpha_df %>%
  ggplot() +
  geom_point(aes(x = filter_size, y = shannon, colour = filter_size), size = 10) +
  scale_color_viridis_d() + 
  labs(x = "", y = "Shannon index (H')") +
  scale_y_continuous(breaks = seq(1.79, 2.1, 0.05)) + 
  theme(legend.position = "none", 
        axis.text.x  = element_text(angle = 90, size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 30)) 
print(plot_shannon)

plot_richness <- alpha_df %>%
  ggplot() +
  geom_point(aes(x = filter_size, y = richness, colour = filter_size), size = 10) +
  scale_color_viridis_d() + 
  labs(x = "", y = "Species richness") +
  theme(legend.position = "none",
        axis.text.x  = element_text(angle = 90, size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 30))
print(plot_richness)

plot_evenness <- alpha_df %>%
  ggplot() +
  geom_point(aes(x = filter_size, y = pielou, colour = filter_size), size = 10) +
  scale_color_viridis_d() + # filter_sizeはdiscreteなので、_dにする
  labs(x = "", y = "Pielou's evenness") +
  scale_y_continuous(breaks = seq(0.43, 0.63, 0.05)) + #y軸いじり
  theme(legend.position = "none",
        axis.text.x  = element_text(angle = 90, size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 30))
print(plot_evenness)



