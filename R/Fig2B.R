library(tidyverse)
library(viridis)

load("Fig2b_data.Rdata")

n_groups <- length(unique(data[["pheno"]]$group))

#Plot an MDS
MDS <- cmdscale(dist(t(data[["data_mds"]])), k = 2)
MDS <- data.frame(MDS)
MDS <- cbind(MDS, data[["pheno"]])
annotationBar <- plasma(n_groups)

p1 <- qplot(data = MDS , y = X1, x = X2 , colour = group, size = I(2))+
  geom_point(data = MDS, aes(y = X1, x = X2), size = I(2))+
  theme_bw()+
  scale_colour_manual(values = annotationBar) +
  theme(panel.border = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom",
        axis.text.x = element_text(colour = "black", size = 11), 
        axis.text.y = element_text(colour = "black", size = 11)) +
  ylab("Dim 1")+
  xlab("Dim 2")

print(p1)
