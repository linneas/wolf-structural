#################
### written by LS
### Plotting Length distributions (incl for rejected, removed and uncurated sites)

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
#Read files and convert to tibbles
len_file="data/lengths/all_lengths.txt"

len_tib<-len_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(Type=V1, Set=V2, Length=V3) %>%
              mutate(Length=as.double(Length)) %>%
              mutate(Set=case_when(Set=="Removed" ~ "Rejected",
                                    Set=="Curated" ~ "Kept",
                                    TRUE ~ Set))

comb_tib <- bind_rows(len_tib)
comb_tib$Type<-factor(comb_tib$Type, levels=c("Deletions","Duplications","Inversions"))
comb_tib$Set<-factor(comb_tib$Set, levels=c("Rejected","Uncurated","Kept"))


################################################################################
# Plot only SVs below 7500bp (for visibility!)
acc_tib<-comb_tib %>% filter(Length<7500)

# Make main plot
mainplot<- ggplot() +
  geom_histogram(data=acc_tib, aes(x=Length, color=Set, fill=Set), alpha=0.5, bins=100) +
  scale_color_manual(name="",values=c("#8FD744", "#2C728E", "#440154")) +
  scale_fill_manual(name="",values=c("#8FD744", "#2C728E", "#440154")) +
  facet_wrap2(vars(Type), scales="free", ncol = 3, trim_blank = FALSE) +
  theme(panel.grid.major = element_line(colour = '#f5f4e6'),
        panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_blank(),
      legend.key.size = unit(0.4, 'cm'),
      legend.position="bottom",
      axis.text=element_text(size=8))+
      labs(x="SV length (bp)", y="Count")

# Inset plot 1, Deletions
p2<-ggplot(subset(acc_tib, Type=="Deletions" & Length<300), aes(x=Length, color=Set, fill=Set)) +
  geom_histogram(alpha=0.5, bins=30, show.legend = FALSE) +
  scale_color_manual(values=c("#8FD744", "#2C728E", "#440154")) +
  scale_fill_manual(values=c("#8FD744", "#2C728E", "#440154")) +
  theme(panel.grid.major = element_line(colour = '#FFFDF8'),
        panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(margin=margin(t=0), size=6),
        plot.margin=unit(c(0.5,0.5,0.5,0.5), "mm"))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c(0,100,200,300), limits=c(0,300))


# Inset plot 2, Duplications
p3<-ggplot(subset(acc_tib, Type=="Duplications" & Length<300), aes(x=Length, color=Set, fill=Set)) +
  geom_histogram(alpha=0.5, bins=30, show.legend = FALSE) +
  scale_color_manual(values=c("#8FD744", "#2C728E", "#440154")) +
  scale_fill_manual(values=c("#8FD744", "#2C728E", "#440154")) +
  theme(panel.grid.major = element_line(colour = '#FFFDF8'),
        panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(margin=margin(t=0), size=6),
        plot.margin=unit(c(0.5,0.5,0.5,0.5), "mm")) +
  scale_x_continuous(breaks=c(0,100,200,300), labels=c(0,100,200,300), limits=c(0,300))


# Inset plot 3, Inversions
p4<-ggplot(subset(acc_tib, Type=="Inversions" & Length<300), aes(x=Length, color=Set, fill=Set)) +
  geom_histogram(alpha=0.5, bins=30, show.legend = FALSE) +
  scale_color_manual(values=c("#8FD744", "#2C728E", "#440154")) +
  scale_fill_manual(values=c("#8FD744", "#2C728E", "#440154")) +
  theme(panel.grid.major = element_line(colour = '#FFFDF8'),
        panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(margin=margin(t=0), size=6),
        plot.margin=unit(c(0.5,0.5,0.5,0.5), "mm"))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c(0,100,200,300), limits=c(0,300))


# Combine inset plots into a data tibble
data.tb <- tibble(x = c(8000, 8000, 8000),
                  y = c(22500, 550, 1250),
                  Type = factor(c("Deletions", "Duplications", "Inversions")),
                  plot = list(p2, p3, p4))

# Combine main plot with inset plots
combplot<-mainplot +
  geom_plot(data=data.tb, aes(x, y, label = plot), vp.width = 0.5, vp.height = 0.6)
# Create output file
outfile=paste(plotdir,"FigureS3.",d,".jpg", sep="")
ggsave(outfile, plot = combplot, scale = 1, dpi = 300, limitsize = TRUE,
    width=8,height=2.5)
