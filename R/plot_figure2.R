#################
### written by LS
### Plotting repeat content 

#Clear the workspace
rm(list=ls())

# Setup Load packages
require(tidyverse)
library(viridis)
library(ggh4x)
plotdir="plots/"
d<-Sys.Date()
format(d, format="%Y-%m-%d")

################################################################################
# Input files
del_file="data/repeats/DEL.repeat.summary.txt"
dup_file="data/repeats/DUP.repeat.summary.txt"
inv_file="data/repeats/INV.repeat.summary.txt"
del_peak_file<-"data/repeats/DEL.peak190.repeat.summary.txt"
dup_peak_file<-"data/repeats/DUP.peak160.repeat.summary.txt"

# Read in the data
del_tib<-del_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(Family=V1, Length=V2) %>% mutate(Set="All deletions") %>%
              mutate(Frac=Length/sum(Length))
dup_tib<-dup_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(Family=V1, Length=V2) %>% mutate(Set="All duplications") %>%
              mutate(Frac=Length/sum(Length))
inv_tib<-inv_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(Family=V1, Length=V2) %>% mutate(Set="All inversions") %>%
              mutate(Frac=Length/sum(Length))
del_peak_tib<-del_peak_file %>% read.table(header=FALSE) %>% as_tibble() %>%
            rename(Family=V1, Length=V2) %>% mutate(Set="Deletions 180-200bp") %>%
            mutate(Frac=Length/sum(Length))
dup_peak_tib<-dup_peak_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(Family=V1, Length=V2) %>% mutate(Set="Duplications 140-180bp") %>%
              mutate(Frac=Length/sum(Length))

# Combine
comb_tib <-bind_rows(del_tib, dup_tib,inv_tib, del_peak_tib, dup_peak_tib) %>%
          mutate(Family=case_when(Family=="Low_complexity" ~ "Low complexity",
                                  Family=="Simple_repeat" ~ "Simple repeat",
                                  TRUE ~ Family))
comb_tib$Set<-factor(comb_tib$Set, levels=c("All deletions","All duplications","All inversions", "Deletions 180-200bp", "Duplications 140-180bp"))

# Plot
p<-ggplot(comb_tib, aes(x="", y=Frac, fill=Family)) +
  geom_bar(stat="identity", width=1, alpha=0.75) +
  coord_polar("y", start=0) +
  facet_wrap2(vars(Set), nrow = 2, ncol = 3, trim_blank = FALSE) +
  scale_color_viridis(discrete=TRUE, na.value="lightgrey") +
	scale_fill_viridis(discrete=TRUE, na.value="lightgrey",labels = function(breaks) {breaks[is.na(breaks)] <- "Not repetetive"; breaks}) +
  labs(fill= "Repeat class") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        strip.background = element_rect(fill = 'lightgrey', colour = 'lightgrey'),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.key.size = unit(0.35, 'cm'))
outfile=paste(plotdir,"Figure2.",d,".png", sep="")
ggsave(outfile, plot = p, scale = 1, dpi = 600, limitsize = TRUE, width=5.5, height=4)
