#################
### written by LS
### Plotting Length distributions

#Clear the workspace
rm(list=ls())

# Setup Load packages
require(tidyverse)
library(viridis)
library(ggh4x)
library(ggplot2)
library(tibble)
library(ggpmisc)
plotdir="plots/"
d<-Sys.Date()
format(d, format="%Y-%m-%d")


################################################################################
# Input files and read in the data
len_file="data/lengths/all_lengths.txt"

comb_tib<-len_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(Type=V1, Set=V2, Length=V3) %>%
              mutate(Length=as.double(Length)) %>%
              filter(Set=="Curated")

comb_tib$Type<-factor(comb_tib$Type, levels=c("Deletions","Duplications","Inversions"))

################################################################################
# Plot with zoom in insets
acc_tib<-comb_tib %>% filter(Length<7500)

# Make main plot
mainplot<- ggplot() +
  geom_histogram(data=acc_tib, aes(x=Length), alpha=0.5, bins=100, colour="#440154", fill="#440154") +
  facet_wrap2(vars(Type), scales="free", ncol = 3, trim_blank = FALSE) +
  theme(panel.grid.major = element_line(colour = '#f5f4e6'),
      panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_blank(),
      axis.text=element_text(size=8))+
      labs(x="SV length (bp)", y="Count")

# Inset plot 1, Deletions
p2<-ggplot(subset(acc_tib, Type=="Deletions" & Length<300), aes(x=Length)) +
  geom_histogram(alpha=0.5, bins=30, colour="#440154", fill="#440154") +
  theme(panel.grid.major = element_line(colour = '#FFFDF8'),
        panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(margin=margin(t=0), size=6),
        plot.margin=unit(c(0.5,0.5,0.5,0.5), "mm")) +
  scale_x_continuous(breaks=c(0,100,200,300), labels=c(0,100,200,300), limits=c(0,300))

# Inset plot 2, Duplications
p3<-ggplot(subset(acc_tib, Type=="Duplications" & Length<300), aes(x=Length)) +
  geom_histogram(alpha=0.5, bins=30, colour="#440154", fill="#440154") +
  theme(panel.grid.major = element_line(colour = '#FFFDF8'),
        panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(margin=margin(t=0), size=6),
        plot.margin=unit(c(0.5,0.5,0.5,0.5), "mm")) +
  scale_x_continuous(breaks=c(0,100,200,300), labels=c(0,100,200,300), limits=c(0,300))

# Inset plot 3, Inversions
p4<-ggplot(subset(acc_tib, Type=="Inversions" & Length<300), aes(x=Length)) +
  geom_histogram(alpha=0.5, bins=30, colour="#440154", fill="#440154") +
  theme(panel.grid.major = element_line(colour = '#FFFDF8'),
      panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(margin=margin(t=0), size=6),
        plot.margin=unit(c(0.5,0.5,0.5,0.5), "mm")) +
  scale_x_continuous(breaks=c(0,100,200,300), labels=c(0,100,200,300), limits=c(0,300))


# Combine inset plots into a data tibble
data.tb <- tibble(x = c(8000, 8000, 8000),
                  y = c(14000, 120, 25),
                  Type = factor(c("Deletions", "Duplications", "Inversions")),
                  plot = list(p2, p3, p4))

# Combine main plot with inset plots
combplot<-mainplot +
  geom_plot(data=data.tb, aes(x, y, label = plot), vp.width = 0.5, vp.height = 0.6)
# Create output file
outfile=paste(plotdir,"Figure1.",d,".png", sep="")
ggsave(outfile, plot = combplot, scale = 1, dpi = 600, limitsize = TRUE,
    width=8,height=2.5)
