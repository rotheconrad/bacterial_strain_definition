library(ggplot2)
library(ggExtra)
library(dplyr)
library(ggpubr)
library(cowplot)


setwd("/Users/djfeistel/Desktop/Projects/ANI_99.5_paper/MLST/UPDATED_DATA/Ecoli_check")
ANI_ecoli <- data.frame(read.table(
  file = 'Ecoli_fastANI_complete.ani_parsed.anisubset_MLSTs.ani', 
  sep = '\t', header = TRUE))

head(ANI_ecoli) 
nrow(ANI_ecoli)

mlst <- "^ST-11$"
ANI_ecoli_filt <- ANI_ecoli %>% filter(grepl(mlst, query_ST))
ANI_ecoli_filt <- ANI_ecoli_filt[ANI_ecoli_filt$ANI >= 99,]
tail(ANI_ecoli)
head(ANI_ecoli)

p1 <- ggplot(data = ANI_ecoli_filt, aes(x=ANI, 
                                        y=proportion, 
                                        color = ST,
                                        shape = ST)) +
  theme_bw() +
  geom_point(size = 4,
             alpha=0.5,
             show.legend = TRUE) +
  # color = "black",
  # fill = "grey") +
  xlab('ANI') + ylab('Shared / Total Fragments') + 
  #ggtitle(title_graph) +
  theme(plot.title = element_text(color="black", 
                                  size=25, 
                                  face="bold",
                                  hjust = 0.5,
                                  vjust=-1),
        axis.text.x=element_text(size=13, hjust=1, angle = 45),
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.position = "left",
        axis.title=element_text(size=15)
  ) + 
  scale_x_continuous(breaks = seq(99, 100, 0.2),
                     minor_breaks = seq(90, 110, 0.1),
                     limits = c(99, 100)) + 
  scale_y_continuous(breaks = seq(0.3, 1.1, 0.05),
                     minor_breaks = seq(0.2, 1, 0.01),
                     limits = c(0.7, 1))
#scale_color_manual(values = c("cyan3", "cyan3", "coral2", "grey40", "grey40","grey40")) +#coral2
#scale_shape_manual(values=c(5, 0, 1, 2))
#colors ST10 = lightseagreen, ST11 = mediumseagreen

ggMarginal(p1, type="density", fill = "gray69", groupColour=TRUE, groupFill=TRUE)

# 
# plot_grid(ST10, ST11, ST131, ST167,
#           labels=c('A','B','C','D'),
#           hjust = -13.75,
#           vjust = 2.5)

# 
# ggarrange(ST10, ST11, ST131, ST167,
#          labels = c("A", "B", "C", "D"),
#          ncol = 2,
#          nrow = 2,
#          hjust = -16,
#          vjust = 5.5,
#          align ="hv"
#          )


