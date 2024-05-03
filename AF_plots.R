library(ggplot2)
library(data.table)
library(ggthemes)
library(gridExtra)
library(ggforce)
library(gridExtra)

# generate a AF plot of the haplotypes and add some segments to the plot to indicate which is a contious segment 
af = read.table("AF_HCGene_Hap_reps_togehter.txt", header = T)
head(af)

# merge the position with the gene name 
loc = fread("Location.txt")
head(loc, 20)
LOC = unique(loc[,c(1,4,6)])

la = merge(LOC, af, by.x="Gene_ID", by.y = "x1")
la = la[order(la$Chr, la$start),]


# test on 500.000 to 2 mio bp on chr 3h

la3 = la[la$Chr=="chr3H" & la$start < 2000000000 & la$start > 50000000,]

k1 = la3[,c("Gene_ID", "Chr", "start","F23P1K1")]
k2 = la3[,c("Gene_ID", "Chr", "start","F23P1K2")]
k1$set = "F23P1K1"
k2$set = "F23P1K2"
colnames(k2)[4] = colnames(k1)[4]
k = rbind(k1,k2)

#ggplot(k, aes(x=start, y=F23P1K1, color=set)) + geom_smooth()
#ggplot(k, aes(x=start, y=F23P1K1, color=set)) + geom_line()  + geom_smooth(span = 0.1, method = "loess") + theme_minimal()
p1 = ggplot(k, aes(x=start/1000000, y=F23P1K1, color=set, shape = set)) + theme_minimal() + geom_point(size=0.5) + geom_smooth(span = 0.08, method = "loess", se=F, size=2) +
  xlab("Genome Position chromosome 3H [Mb]") + ylab("Wild type allele frequency")
p2 = ggplot(k, aes(x=start/1000000, y=F23P1K1, color=set, shape = set)) + theme_minimal() + geom_point(size=0.5) + geom_smooth(span = 0.03, method = "loess", se=F, size=2) +
  xlab("Genome Position chromosome 3H [Mb]") + ylab("Wild type allele frequency")


ggplot(k, aes(x=start/1000000, y=F23P1K1, color=set, shape = set)) +  geom_point(alpha=0.1) +
  geom_rect(aes(xmin=40, xmax=420, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.01) + 
  theme_hc() + scale_color_gdocs() +  geom_smooth(aes(linetype = set),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
  xlab("") + ylab("Wild form allele frequency") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) +
  guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
  ylim(0,0.65) + theme(legend.position = "none", text = element_text(size = 25))

grid.arrange(p1,p2, ncol=1)

#https://stackoverflow.com/questions/71428312/is-there-a-way-to-plot-line-with-start-and-end-coordinate

# generate a dataframe where major genes are stored in , including their position and allele frequency.
MajorGenes = data.frame(Chr = c("chr1H", "chr1H", "chr2H", "chr4H", "chr5H", "chr5H", "chr7H", "chr2H"), start = c(500,550,60,630,550,600,45,450), 
                        Gene = c("PpdH2", "HvElf2", "HvPpdH1", "VrnH2", "HvCbf", "VrnH1", "VrnH3", "HvCen"), HAF=0.6, env="else",
                        Col = c("red", "yellow", "green", "blue", "orange", "gray", "black", "violet"))
MajorGenes$Gene = as.factor(MajorGenes$Gene)
### generate the above plot for all generations and chromosomes 
##########################################################################################################################
##########################################################################################################################
# create a list with all the population samples in it - so that we can build up on this. 
chrlist = unique(la$Chr)
plotlist = list()


for (CHR in chrlist){
  
                    la3 = la[la$Chr==CHR,]
                    la4 = as.data.frame(la3)
                    p1 = list()
                    
                    for (i in colnames(la)[4:ncol(la)]){
                      p1[[i]] = la4[,c("Gene_ID", "Chr", "start",i)]
                      p1[[i]]$set = i
                      colnames(p1[[i]])[4] = "HAF"
                    }
                    
                    # merge them together in one single dataframe
                    sink = data.frame(Gene_ID = character(), Chr = character(), start = integer(), HAF = numeric())
                    for (i in names(p1)){
                      sink = rbind(sink, p1[[i]])
                      
                    }
                    
                    # extract if it is organic, conventional or starting generation and which generation it is 
                    library(stringr)
                    sink$gen = str_split_fixed(sink$set, "P1", 2)[,1]
                    sink$env = str_split_fixed(sink$set, "P1", 2)[,2]
                    sink[grep("K", sink$env), "env"] = "Conventional"
                    sink[grep("Ö", sink$env), "env"] = "Organic"
                    sink[sink$gen=="F3", "env"] = "Origin"
                    
                    sink$gen = factor(sink$gen, levels = c("F3", "F12", "F16", "F22", "F23"))
                    sink$env = factor(sink$env, levels = c("Origin", "Organic", "Conventional"))
                    sink$Chr = as.factor(sink$Chr)
                    
                    #sink = sink[sample(1:nrow(sink), 5000),]
                    MG = MajorGenes[MajorGenes$Chr==CHR,]
                    
                    if (CHR == "chr1H"){
                      plotlist[[CHR]] = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) + #  geom_point(alpha=0.1) +
                                                 geom_rect(aes(xmin=40, xmax=420, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.01) + 
                                                 theme_classic() + scale_color_gdocs() +  geom_smooth(aes(linetype = env),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
                                                 #geom_point(data = MG, aes(x=start, y=HAF, color=Col), size=5) + 
                                                 xlab("") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + 
                                                 labs(title = "A", subtitle = "Wild form allele frequency [ratio]") +
                                                 guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), 
                                                 color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
                                                 ylim(0,0.65) + theme(legend.position = "none", text = element_text(size = 25), plot.subtitle = element_text(hjust = -0.03), 
                                                 plot.title = element_text(size = 30, face = "bold")) 

                    #  p = ggplot(sink[1:1000,], aes(x=start/1000000, y=HAF, color=gen, shape=env)) + #  geom_point(alpha=0.1) +
                    #                              geom_rect(aes(xmin=40, xmax=420, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.01) + 
                    #                              theme_classic() + scale_color_gdocs() +  geom_smooth(aes(linetype = env),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
                    #                              #geom_point(data = MG, aes(x=start, y=HAF, color=Col), size=5) + 
                    #                              xlab("") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + 
                    #                              labs(title = "A", subtitle = "Wild form allele frequency [ratio]") +
                    #                              guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), 
                    #                              color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
                    #                              ylim(0,0.65) + theme(legend.position = "none", text = element_text(size = 25), plot.subtitle = element_text(hjust = -0.03), 
                    #                              plot.title = element_text(size = 30, face = "bold")) 
                    # p                      
                      
                    } else if (CHR == "chr3H"){
                                                  plotlist[[CHR]] = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) +  # geom_point(alpha=0.1) +
                                                    #geom_point(data = MG, aes(x=start, y=HAF, color=Col, size = 15)) +
                                                    geom_rect(aes(xmin=35, xmax=450, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.01) + 
                                                    theme_classic() + scale_color_gdocs() +  geom_smooth(aes(linetype = env),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
                                                    xlab("") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + 
                                                    labs(title = "C", subtitle = "Wild form allele frequency [ratio]") +
                                                    guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), 
                                                    color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
                                                    ylim(0,0.65) + theme(legend.position = "none", text = element_text(size = 25), plot.subtitle = element_text(hjust = -0.03), 
                                                    plot.title = element_text(size = 30, face = "bold")) 


                                                  
                    } else if (CHR == "chr5H"){
                                                plotlist[[CHR]] = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) + # geom_point(alpha=0.1) +
                                                  #geom_point(data = MG, aes(x=start, y=HAF, color=Col, size = 15)) +
                                                  geom_rect(aes(xmin=25, xmax=500, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.01) + 
                                                  theme_classic() + scale_color_gdocs() +  geom_smooth(aes(linetype = env),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
                                                    xlab("") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + 
                                                    labs(title = "E", subtitle = "Wild form allele frequency [ratio]") +
                                                    guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), 
                                                    color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
                                                    ylim(0,0.65) + theme(legend.position = "none", text = element_text(size = 25), plot.subtitle = element_text(hjust = -0.03), 
                                                    plot.title = element_text(size = 30, face = "bold")) 

                      
                      
                      
                    } else if (CHR == "chr2H"){
                                                      plotlist[[CHR]] = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) + #  geom_point(alpha=0.1) +
                                                       # geom_point(data = MG, aes(x=start, y=HAF, color=Col, size = 15)) +
                                                        geom_rect(aes(xmin=160, xmax=500, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.01) + 
                                                        theme_classic() + scale_color_gdocs() +  geom_smooth(aes(linetype = env),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
                                                        xlab("") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + 
                                                        labs(title = "B") +
                                                        guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), 
                                                        color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
                                                        ylim(0,0.65) + theme(legend.position = "none", text = element_text(size = 25), plot.subtitle = element_text(hjust = -0.03), 
                                                        plot.title = element_text(size = 30, face = "bold"),  axis.text.y = element_blank()) 
                                                      
                    } else if (CHR == "chr4H"){
                                              plotlist[[CHR]] = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) + #  geom_point(alpha=0.1) +
                                               # geom_point(data = MG, aes(x=start, y=HAF, color=Col, size = 15)) +
                                                geom_rect(aes(xmin=110, xmax=420, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.01) + 
                                                theme_classic() + scale_color_gdocs() +  geom_smooth(aes(linetype = env),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
                                                xlab("") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + 
                                                labs(title = "D") +
                                                guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), 
                                                color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
                                                ylim(0,0.65) + theme(legend.position = "none", text = element_text(size = 25), plot.subtitle = element_text(hjust = -0.03), 
                                                plot.title = element_text(size = 30, face = "bold"),  axis.text.y = element_blank()) 
                                              
                    } else if (CHR == "chr6H" ){
                                              plotlist[[CHR]] = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) +  # geom_point(alpha=0.1) +
                                               # geom_point(data = MG, aes(x=start, y=HAF, color=Col, size = 15)) +
                                                geom_rect(aes(xmin=30, xmax=430, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.01) + 
                                                theme_classic() + scale_color_gdocs() +  geom_smooth(aes(linetype = env),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
                                                xlab("") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + 
                                                        labs(title = "F") +
                                                        guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), 
                                                        color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
                                                        ylim(0,0.65) + theme(legend.position = "none", text = element_text(size = 25), plot.subtitle = element_text(hjust = -0.03), 
                                                        plot.title = element_text(size = 30, face = "bold"),  axis.text.y = element_blank()) 
                                                          
                                  
                                  
                    } else if (CHR == "chrUn"){
                                    plotlist[[CHR]] = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) + # geom_point(alpha=0.1) +
                                      theme_classic() + scale_color_gdocs() +  geom_smooth(aes(linetype = env),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
                                      xlab("Genome position on chromosome [Mb]") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + 
                                      labs(title = "H") +
                                      guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
                                      ylim(0,0.65) + theme(legend.position = c(0.8, 0.7), legend.box = "horizontal", legend.background = element_rect(fill="lightyellow", 
                                      size=0.5, linetype="solid", colour ="darkblue"), text = element_text(size = 30), plot.title = element_text(size = 35, face = "bold"),
                                       axis.text.y = element_blank(), legend.text = element_text(size=25), legend.key.size = unit(3, 'cm'))
                                    
                    } else if (CHR == "chr7H"){
                                      plotlist[[CHR]] = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) + # geom_point(alpha=0.1) +
                                       # geom_point(data = MG, aes(x=start, y=HAF, color=Gene, size = 15)) +
                                        geom_rect(aes(xmin=65, xmax=480, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.01) + 
                                        theme_classic() + scale_color_gdocs() +   geom_smooth(aes(linetype = env),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
                                        xlab("Genome position on chromosome [Mb]") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + 
                                        labs(title = "G", subtitle = "Wild form allele frequency [ratio]") +
                                        guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
                                        ylim(0,0.65) + theme(legend.position = "none", text = element_text(size = 25),plot.subtitle = element_text(hjust = -0.03), 
                                                        plot.title = element_text(size = 30, face = "bold"))
                      
                    }
                    
                    ## generate two additional plots for the chr 3h and 5h 
                    if (CHR == "chr3H"){
                                                  plotlist[[paste0(CHR,"_2")]] = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) + # geom_point(alpha=0.1) +
                                               #     geom_point(data = MG, aes(x=start, y=HAF, color=Gene, size = 15)) +
                                                    geom_rect(aes(xmin=35, xmax=450, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.01) + 
                                                    theme_classic() + scale_color_gdocs() +   geom_smooth(aes(linetype = env),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
                                                    xlab("Genome position on chromosome [Mb]") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + 
                                                    guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
                                                    ylim(0,0.65) + theme(legend.position = "none", text = element_text(size = 14)) + 
                                                    labs(title = "", subtitle = "Wild form allele frequency [ratio]") 
                    }

                    if (CHR == "chr5H"){
                                                  plotlist[[paste0(CHR,"_2")]] = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) + # geom_point(alpha=0.1) +
                                               #     geom_point(data = MG, aes(x=start, y=HAF, color=Gene, size = 15)) +
                                                    geom_rect(aes(xmin=25, xmax=500, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.01) + 
                                                    theme_classic() + scale_color_gdocs() +   geom_smooth(aes(linetype = env),span = 0.03, method = "loess", se=F, size=1, alpha=0.7) +
                                                    xlab("Genome position on chromosome [Mb]") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + 
                                                    guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Farming system")) + 
                                                    ylim(0,0.65) + theme(legend.position = "none", text = element_text(size = 14)) + 
                                                    labs(title = "", subtitle = "Wild form allele frequency [ratio]") 
                    }
                    
                  
}

# 
# grid.arrange(plotlist[["chr1H"]], plotlist[["chr2H"]], plotlist[["chr3H"]], plotlist[["chr4H"]], 
#              plotlist[["chr5H"]], plotlist[["chr6H"]], plotlist[["chr7H"]], plotlist[["chrUn"]], ncol=2)


# print the legend only 
plegend = ggplot(MajorGenes, aes(x=start, y=HAF, color=Gene)) + geom_point() + facet_wrap(.~Chr) + theme(legend.position="bottom")

legend <- cowplot::get_legend(plegend)
grid::grid.newpage()
grid::grid.draw(legend) # plots the legend of the famous genes 

#########################################################################################################################
#########################################################################################################################
## make a zoom in plot for the chromosome 3H 30 mio to 160 mio 
#########################################################################################################################
#########################################################################################################################

la3 = la[la$Chr=="chr3H" & la$start < 160000000 & la$start > 30000000,]
la4 = as.data.frame(la3)
p1 = list()

for (i in colnames(la)[4:ncol(la)]){
  p1[[i]] = la4[,c("Gene_ID", "Chr", "start",i)]
  p1[[i]]$set = i
  colnames(p1[[i]])[4] = "HAF"
}

# merge them together in one single dataframe
sink = data.frame(Gene_ID = character(), Chr = character(), start = integer(), HAF = numeric())
for (i in names(p1)){
  sink = rbind(sink, p1[[i]])
  
}

# extract if it is organic, conventional or starting generation and which generation it is 
library(stringr)
sink$gen = str_split_fixed(sink$set, "P1", 2)[,1]
sink$env = str_split_fixed(sink$set, "P1", 2)[,2]
sink[grep("K", sink$env), "env"] = "Conventional"
sink[grep("Ö", sink$env), "env"] = "Organic"
sink[sink$gen=="F3", "env"] = "Origin"

sink$gen = factor(sink$gen, levels = c("F3", "F12", "F16", "F22", "F23"))
sink$env = factor(sink$env, levels = c("Origin", "Organic", "Conventional"))
sink$Chr = as.factor(sink$Chr)

p2 = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) +  theme_wsj(color = "gray") + scale_color_gdocs() +  geom_smooth(aes(linetype = env),span = 0.08, method = "loess", se=F, size=1, alpha=0.7) +
  xlab("") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) + #geom_point(alpha=0.1) +
  guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Agro-ecosystem")) + 
  ylim(0, 0.2) + theme(legend.position = c(0.45, - 0.1), legend.box = "horizontal", 
  legend.background = element_rect(fill="lightyellow", size=0.5, linetype="solid", colour ="darkblue"),
  text = element_text(size = 15), legend.text = element_text(size=20), legend.key.size = unit(2, 'cm')) + 
  annotate("segment", x = 45, y = 0.18, xend = 45, yend = 0.15, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x = 45, y = 0.19, label = "brt1/2", size = 8, color = "navy")


plotlist[["chr3H_3"]] = plotlist[["chr3H_2"]] + annotate("segment", x = 634, y = 0.58, xend = 634, yend = 0.55, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x = 634, y = 0.6, label = "Denso", size = 10, color = "navy") +
  geom_rect(aes(xmin = 30, xmax = 160, ymin = 0, ymax = 0.1), color = "black", alpha = 0, size = 2) +
  annotation_custom(ggplotGrob(p2), xmin = 0, xmax = 550, ymin = 0.25, ymax = 0.61) +
  geom_path(aes(x,y,group=grp, color="black", shape="none"), data=data.frame(x = c(30,0,160,550), y=c(0.1,0.25,0.1,0.25),grp=c(1,1,2,2)), linetype='dotdash', size=1.5, color="black") +
  theme(text = element_text(size = 30), plot.subtitle = element_text(hjust = -0.03))  



#########################################################################################################################
#########################################################################################################################
## make a zoom in plot for the chromosome 5H 350 mio to 550 mio 
#########################################################################################################################
#########################################################################################################################

la3 = la[la$Chr=="chr5H" & la$start < 520000000 & la$start > 350000000,]
la4 = as.data.frame(la3)
p1 = list()

for (i in colnames(la)[4:ncol(la)]){
  p1[[i]] = la4[,c("Gene_ID", "Chr", "start",i)]
  p1[[i]]$set = i
  colnames(p1[[i]])[4] = "HAF"
}

# merge them together in one single dataframe
sink = data.frame(Gene_ID = character(), Chr = character(), start = integer(), HAF = numeric())
for (i in names(p1)){
  sink = rbind(sink, p1[[i]])
  
}

# extract if it is organic, conventional or starting generation and which generation it is 
library(stringr)
sink$gen = str_split_fixed(sink$set, "P1", 2)[,1]
sink$env = str_split_fixed(sink$set, "P1", 2)[,2]
sink[grep("K", sink$env), "env"] = "Conventional"
sink[grep("Ö", sink$env), "env"] = "Organic"
sink[sink$gen=="F3", "env"] = "Origin"

sink$gen = factor(sink$gen, levels = c("F3", "F12", "F16", "F22", "F23"))
sink$env = factor(sink$env, levels = c("Origin", "Organic", "Conventional"))
sink$Chr = as.factor(sink$Chr)

p2 = ggplot(sink, aes(x=start/1000000, y=HAF, color=gen, shape=env)) +  theme_wsj(color = "gray") + scale_color_gdocs() +  geom_smooth(aes(linetype = env),span = 0.1, method = "loess", se=F, size=1, alpha=0.7) +
  xlab("") + ylab("") +  scale_linetype_manual(values=c("dotdash", "dashed", "solid")) +
  guides(shape = guide_legend(override.aes = list(alpha=1), title="Farming system"), color=guide_legend(title = "Generation"), linetype = guide_legend(title = "Agro-ecosystem")) + 
  ylim(0, 0.63) + theme(legend.position = c(0.45, - 0.1), legend.box = "horizontal", 
  legend.background = element_rect(fill="lightyellow", size=0.5, linetype="solid", colour ="darkblue"),
  text = element_text(size = 15), legend.text = element_text(size=20), legend.key.size = unit(2, 'cm'))  + 
  annotate("segment", x = 490, y = 0.48, xend = 490, yend = 0.35, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x = 490, y = 0.52, label = "SD1", size = 8, color = "navy") + 
  annotate("segment", x = 486, y = 0.59, xend = 486, yend = 0.505, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x = 486, y = 0.62, label = "QG[L/W]5H", size = 8, color = "navy")  

plotlist[["chr5H_3"]] = plotlist[["chr5H_2"]] + ylim(0, 0.85) + geom_rect(aes(xmin = 350, xmax = 520, ymin = 0, ymax = 0.45), color = "black", alpha = 0, size = 1) +
  geom_path(aes(x,y,group=grp, color="black", shape="none"), data=data.frame(x = c(350,180,520,650), y=c(0.45,0.5,0.45,0.5),grp=c(1,1,2,2)), linetype='dotdash', size=1, color="black") +
  annotation_custom(ggplotGrob(p2), xmin = 180, xmax = 650, ymin = 0.5, ymax = 0.9) +  theme(text = element_text(size = 30), plot.subtitle = element_text(hjust = -0.03))  


## save all plot to a file 

for (i in names(plotlist)[11:12]){
  png(paste0("F4_5_GenomePlot_", i, "3_adjusted.png"), width = 800, height = 350, units = "mm", res = 200)
  print(plotlist[[i]])
  dev.off()
}

png("F3_GenomePlot_allChr_pericentromere_noDots_adjusted_2.png", width = 1000, height = 1100, units = "mm", res = 250)
grid.arrange(plotlist[["chr1H"]], plotlist[["chr2H"]], plotlist[["chr3H"]], plotlist[["chr4H"]], 
             plotlist[["chr5H"]], plotlist[["chr6H"]], plotlist[["chr7H"]], plotlist[["chrUn"]], ncol=2)
dev.off()

png("F3_GenomePlot_allChr_pericentromere_Legend.png", width = 800, height = 50, units = "mm", res = 250)
grid::grid.draw(legend)
dev.off()

for (i in unique(la$Chr)){
  png(paste0("GenomePlot_", i, "3.png"), width = 500, height = 200, units = "mm", res = 200)
  print(plotlist[[i]])
  dev.off()
}




################################################################################
## counting crossing over events on the chromsomes by using a loess function
# test = loess(F23P1Ö1 ~ start, data = la[la$Chr=="chr3H",], span = 0.08)
# Test = la[la$Chr=="chr3H",]
# Test$loessF23o = predict(test, Test)
# 
# test = loess(F23P1Ö1 ~ start, data = la[la$Chr=="chr3H",], span = 0.01)
# Test$loess2F23o = predict(test, Test)
# 
# 
# test = loess(F23P1Ö1 ~ start, data = la[la$Chr=="chr3H",], span = 0.03)
# Test$loess4F23o = predict(test, Test)
# 
# a = ggplot(Test, aes(x=start/1000000, y=F23P1Ö1)) + theme_minimal() + geom_point(size=0.5)  +
#   xlab("Genome Position chromosome 3H [Mb]") + ylab("Wild type allele frequency")
# 
# b =  ggplot(Test, aes(x=start/1000000, y=loessF23o)) + theme_minimal() + geom_point(size=0.5)  +
#   xlab("Genome Position chromosome 3H [Mb]") + ylab("Wild type allele frequency")
# 
# c =  ggplot(Test, aes(x=start/1000000, y=loess2F23o)) + theme_minimal() + geom_point(size=0.5)  +
#   xlab("Genome Position chromosome 3H [Mb]") + ylab("Wild type allele frequency")
# 
# d = ggplot(Test, aes(x=start/1000000, y=loess4F23o)) + theme_minimal() + geom_point(size=0.5)  +
#   xlab("Genome Position chromosome 3H [Mb]") + ylab("Wild type allele frequency")
# 
# gridExtra::grid.arrange(a,b,c,d, ncol=1)


chrlist = unique(la$Chr)
droplist = list()
plotlist2 = list()

# first - convert all actual AF to loess corrected AF values
for (CHR in chrlist){
  la2 = la[la$Chr==CHR,]
  la2.1 = as.data.frame(la2)
  lA2.1 = apply(la2.1[,4:ncol(la2.1)],2, function(x) predict(loess(x ~ start, data = la2.1, span = 0.03), la2.1))
  ll2 = cbind(la2.1[,1:3], lA2.1)
  droplist[[CHR]] = ll2
  
}


# merge them back in a single file 
LA = data.frame(Gene_ID=character(),  Chr = character(), start=integer(), F3P1K1=numeric(), F3P1K2=numeric(), F12P1K1=numeric(),  F12P1K3=numeric(),
                F12P1Ö1=numeric(),  F12P1Ö2=numeric(),  F16P1K1=numeric(), F16P1K2=numeric(), F16P1Ö1=numeric(), F16P1Ö2=numeric(), F22P1KEP=numeric(), 
                F22P1K1=numeric(),  F22P1Ö1=numeric(), F22P1Ö2=numeric(), F23P1K1=numeric(), F23P1K2=numeric(), F23P1Ö1=numeric(),  F23P1Ö2=numeric())

for(i in chrlist){
  LA = rbind(LA, droplist[[i]])
}

## stack the LA df by measurements  
p1 = list()

for (i in colnames(la)[4:ncol(la)]){
  p1[[i]] = LA[,c("Gene_ID", "Chr", "start",i)]
  p1[[i]]$set = i
  colnames(p1[[i]])[4] = "HAF"
}

# merge them together in one single dataframe
sink = data.frame(Gene_ID = character(), Chr = character(), start = integer(), HAF = numeric())
for (i in names(p1)){
  sink = rbind(sink, p1[[i]])
  
}

# extract if it is organic, conventional or starting generation and which generation it is 
library(stringr)
sink$gen = str_split_fixed(sink$set, "P1", 2)[,1]
sink$env = str_split_fixed(sink$set, "P1", 2)[,2]
sink[grep("K", sink$env), "env"] = "Conventional"
sink[grep("Ö", sink$env), "env"] = "Organic"
sink[sink$gen=="F3", "env"] = "Origin"

sink$gen = factor(sink$gen, levels = c("F3", "F12", "F16", "F22", "F23"))
sink$env = factor(sink$env, levels = c("Origin", "Organic", "Conventional"))
sink$Chr = as.factor(sink$Chr)

## specifiy the centromeric regions 
CO = data.frame(Chr = c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H"), start = c(40,160,35,110,25,30,65), fin = c(420,500,450,420,500,420,480))

## candidate genes
MajorGenes = data.frame(Chr = c("chr1H", "chr1H", "chr2H", "chr4H", "chr5H", "chr5H", "chr7H", "chr2H"), start = c(500,550,60,630,550,600,45,450), 
                        Gene = c("PpdH2", "HvElf2", "HvPpdH1", "VrnH2", "HvCbf", "VrnH1", "VrnH3", "HvCen"), HAF=0.4, env="else",
                        Col = c("red", "yellow", "green", "blue", "orange", "gray", "black", "violet"))
MajorGenes$Gene = as.factor(MajorGenes$Gene)


#### Organic plot first                                           
# pick the 23rd generation to count the recombination events 
p = ggplot(sink[sink$gen=="F23" & sink$Chr!="chrUn" & sink$env=="Organic",]) + geom_rect(data = CO, aes(xmin=start, xmax=fin, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.8) +
  theme_minimal() + geom_point(aes(x=start/1000000, y=HAF), size=0.5, color="#E69F00")   + ggtitle("") + 
  xlab("Genome Position [Mb]") + ylab("") + facet_wrap(.~Chr, ncol = 1, strip.position = "bottom") + 
  theme(text = element_text(size = 25), axis.text.y = element_blank()) # + scale_x_continuous(breaks = seq(0, 800, 20))


### an table where to find the annotations 
## select all those position where a crossing over occured 
chr1 = c(12, 23, 35, 45, 330, 340, 355, 370, 450, 460, 470, 490, 520, 526, 530, 540, 545, 550)
chr2 = c(5, 15, 30, 48, 70, 85, 105, 125, 140, 600, 620, 630, 640, 665, 680, 695, 705, 720)
chr3 = c(8, 25, 30, 45, 60, 85, 100, 120, 135, 155, 200, 420, 430, 465, 500, 525, 545, 560, 580, 590, 610, 630, 650, 660, 670, 690)
chr4 = c(15, 20, 30, 60, 100, 340, 365, 520, 540, 565, 590, 600, 610, 615, 622, 630)
chr5 = c(10, 28, 75, 310, 330, 355, 400, 415, 432, 470, 491, 515, 540, 560, 575, 610, 615, 625, 650)
chr6 = c(10, 20, 25, 130, 155, 175, 305, 345, 375, 390, 420, 440, 470, 490, 510)
chr7 = c(7, 15, 23, 35, 65, 73, 80, 110, 410, 450, 465, 510, 545, 570, 590, 625, 630, 640, 650)


# create a new base file based on these values 
base = data.frame(Chr = c(rep("chr1H", length(chr1)), rep("chr2H", length(chr2)), rep("chr3H", length(chr3)), rep("chr4H", length(chr4)), rep("chr5H", length(chr5)), rep("chr6H", length(chr6)), 
                          rep("chr7H", length(chr7))),
                  start = c(chr1, chr2, chr3, chr4, chr5, chr6, chr7),
                  y = -0.05,
                  label = c(1:length(c(chr1, chr2, chr3, chr4, chr5, chr6, chr7)))
)

# write each second y value to 0.6
base[seq(2, nrow(base), by = 2),"y"] = 0.55


pOld = p + geom_text(data = base, size = 3.5, color = "blue", mapping = aes(x = start, y = y, label = label)) +  
        geom_point(data = MajorGenes, aes(x=start, y=HAF, color=Gene), size=5) + theme(legend.position = "right")


################################################################################  
### count the crossing overs in the conventional environment 
# pick the 23rd generation to count the recombination events 
pc = ggplot(sink[sink$gen=="F23" & sink$Chr!="chrUn" & sink$env=="Conventional",]) + 
    geom_rect(data = CO, aes(xmin=start, xmax=fin, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.8) +
    theme_minimal() + geom_point(aes(x=start/1000000, y=HAF), size=0.5, color="#999999")  +
    xlab("Genome Position [Mb]") + ylab("") + facet_wrap(.~Chr, ncol = 1, strip.position = "bottom") +  
    theme(text = element_text(size = 25), plot.subtitle = element_text(hjust = -0.1),  plot.title = element_text(size = 30, face = "bold")) + 
    labs(title = "B", subtitle = "Wild form haplotype allele frequency")

#+ scale_x_continuous(breaks = seq(0, 800, 20))


### an table where to find the annotations 
## select all those position where a crossing over occured 
chr1 = c(10, 30, 75, 335, 350, 366, 390, 420, 440, 455, 520)
chr2 = c(5, 15, 25, 640, 695)
chr3 = c(10, 20, 570, 580, 595, 628, 640, 650, 670, 685)
chr4 = c(10, 20, 25, 605, 610)
chr5 = c(10, 300, 330, 360, 495, 560, 570, 585, 600, 625, 640)
chr6 = c(5, 28, 130, 525)
chr7 = c(65, 585)


# create a new base file based on these values 
base2 = data.frame(Chr = c(rep("chr1H", length(chr1)), rep("chr2H", length(chr2)), rep("chr3H", length(chr3)), rep("chr4H", length(chr4)), rep("chr5H", length(chr5)), rep("chr6H", length(chr6)), 
                           rep("chr7H", length(chr7))),
                   start = c(chr1, chr2, chr3, chr4, chr5, chr6, chr7),
                   y = -0.05,
                   label = c(1:length(c(chr1, chr2, chr3, chr4, chr5, chr6, chr7)))
)

pCld = pc + geom_text(data  = base2, size = 3.5, color = "red", mapping = aes(x = start, y = y, label = label)) + 
        geom_point(data = MajorGenes, aes(x=start, y=HAF, color=Gene), size=5) + theme(legend.position = "none")

grid.arrange(pCld, pOld, ncol=2)   


# add a barplot with the number of crossing overs 
baseinfo = as.data.frame(table(base$Chr))
baseinfo$Env = "Organic"

Baseinfo = as.data.frame(table(base2$Chr))
Baseinfo$Env = "Conventional"

baseinfo = rbind(baseinfo, Baseinfo)
colnames(baseinfo)[3] = "Environment"

p3 = ggplot(baseinfo, aes(x=Var1, y=Freq, fill = Environment)) + geom_bar(stat = "identity", position=position_dodge(), color="black") + 
  scale_fill_manual(values=c('#999999','#E69F00')) + xlab("Chromosome") + ylab("") + theme_minimal() + guides(fill = guide_legend(title="Agro-ecosystem")) +
  geom_text(size = 5, aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +  
  theme(text = element_text(size = 25), legend.position=c(0.9, 0.9), plot.subtitle = element_text(hjust = -0.03), plot.title = element_text(size = 30, face = "bold"))  + 
  labs(title = "A", subtitle = "Crossing over [count]")
p3
##  merge the plots together                          

lay <- rbind(c(3,3),
             c(3,3),
             c(1,2),
             c(1,2),
             c(1,2),
             c(1,2),
             c(1,2))          

png("F6_COcount_adjusted.png", width = 600, height = 500, res = 200, units = "mm")
grid.arrange(pCld,pOld, p3, layout_matrix = lay)                              
dev.off()


  # make plots for all genertions 
 
    # F3
     po3 = ggplot(sink[sink$gen=="F3" & sink$Chr!="chrUn" & sink$env=="Origin",]) + 
              geom_hline(yintercept=0.125, linetype = 2) +
              geom_rect(data = CO, aes(xmin=start, xmax=fin, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.8) +
              theme_minimal() + geom_point(aes(x=start/1000000, y=HAF), size=0.5, color="#a14747") + ggtitle("F3 - origin send to the field") + 
              xlab("Genome Position [Mb]") + ylab("Wild form haplotype allele frequency") + facet_wrap(.~Chr, ncol = 1, strip.position = "bottom") +  
              theme(text = element_text(size = 25), plot.title = element_text(size=35)) + scale_y_continuous(breaks=seq(0,0.6,0.2)) +
              geom_point(data = MajorGenes, aes(x=start, y=HAF, color=Gene), size=5) + 
              geom_hline(yintercept=0.125)
    # F12
    po12 = ggplot(sink[sink$gen=="F12" & sink$Chr!="chrUn" & sink$env=="Organic",]) + geom_rect(data = CO, aes(xmin=start, xmax=fin, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.8) +
               geom_point(data = sink[sink$gen=="F3" & sink$Chr!="chrUn" & sink$env == "Origin",], aes(x=start/1000000, y=HAF), size=1.5, color="#a14747", alpha=0.05)  +
               geom_point(data = sink[sink$gen=="F12" & sink$Chr!="chrUn" & sink$env=="Conventional",], aes(x=start/1000000, y=HAF), size=1.5, color="#999999")  +
              theme_minimal() + geom_point(aes(x=start/1000000, y=HAF), size=1.5, color="#E69F00") + ggtitle("F12") + 
              xlab("Genome Position [Mb]") + ylab("Wild form haplotype allele frequency") + facet_wrap(.~Chr, ncol = 1, strip.position = "bottom") +  
              theme(text = element_text(size = 25), plot.title = element_text(size=35)) + scale_y_continuous(breaks=seq(0,0.6,0.2)) + 
              geom_point(data = MajorGenes, aes(x=start, y=HAF, color=Gene), size=5)
    # F16
    po16 = ggplot(sink[sink$gen=="F16" & sink$Chr!="chrUn" & sink$env=="Organic",]) + geom_rect(data = CO, aes(xmin=start, xmax=fin, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.8) +
               geom_point(data = sink[sink$gen=="F3" & sink$Chr!="chrUn" & sink$env == "Origin",], aes(x=start/1000000, y=HAF), size=1.5, color="#a14747", alpha=0.05)  +
                geom_point(data = sink[sink$gen=="F16" & sink$Chr!="chrUn" & sink$env=="Conventional",], aes(x=start/1000000, y=HAF), size=1.5, color="#999999")  +
              theme_minimal() + geom_point(aes(x=start/1000000, y=HAF), size=1.5, color="#E69F00") + ggtitle("F16") + 
              xlab("Genome Position [Mb]") + ylab("Wild form haplotype allele frequency") + facet_wrap(.~Chr, ncol = 1, strip.position = "bottom") +  
              theme(text = element_text(size = 25),   plot.title = element_text(size=35)) + scale_y_continuous(breaks=seq(0,0.6,0.2)) +
              geom_point(data = MajorGenes, aes(x=start, y=HAF, color=Gene), size=5)
    # F22
    po22 = ggplot(sink[sink$gen=="F22" & sink$Chr!="chrUn" & sink$env=="Organic",]) + geom_rect(data = CO, aes(xmin=start, xmax=fin, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.8) +
               geom_point(data = sink[sink$gen=="F3" & sink$Chr!="chrUn" & sink$env == "Origin",], aes(x=start/1000000, y=HAF), size=1.5, color="#a14747", alpha=0.05)  +
                geom_point(data = sink[sink$gen=="F22" & sink$Chr!="chrUn" & sink$env=="Conventional",], aes(x=start/1000000, y=HAF), size=1.5, color="#999999")  +
              theme_minimal() + geom_point(aes(x=start/1000000, y=HAF), size=1.5, color="#E69F00") + ggtitle("F22") + 
              xlab("Genome Position [Mb]") + ylab("Wild form haplotype allele frequency") + facet_wrap(.~Chr, ncol = 1, strip.position = "bottom") +  
              theme(text = element_text(size = 25), plot.title = element_text(size=35)) + scale_y_continuous(breaks=seq(0,0.6,0.2)) +
              geom_point(data = MajorGenes, aes(x=start, y=HAF, color=Gene), size=5)
    # F23
    po23 = ggplot(sink[sink$gen=="F23" & sink$Chr!="chrUn" & sink$env=="Organic",]) + geom_rect(data = CO, aes(xmin=start, xmax=fin, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.8) +
               geom_point(data = sink[sink$gen=="F3" & sink$Chr!="chrUn" & sink$env == "Origin",], aes(x=start/1000000, y=HAF), size=1.5, color="#a14747", alpha=0.05)  +
                geom_point(data = sink[sink$gen=="F23" & sink$Chr!="chrUn" & sink$env=="Conventional",], aes(x=start/1000000, y=HAF), size=1.5, color="#999999")  +
              theme_minimal() + geom_point(aes(x=start/1000000, y=HAF), size=1.5, color="#E69F00") + ggtitle("F23") + 
              xlab("Genome Position [Mb]") + ylab("Wild form haplotype allele frequency") + facet_wrap(.~Chr, ncol = 1, strip.position = "bottom") +  
              theme(text = element_text(size = 25), plot.title = element_text(size=35)) + scale_y_continuous(breaks=seq(0,0.6,0.2)) +
              geom_point(data = MajorGenes, aes(x=start, y=HAF, color=Gene), size=5)

    
  # bind the same year togehter

  png("F23_genomeplot.png", width = 700, height = 600, res = 200, units = "mm")
  po23
  dev.off()

  png("F22_genomeplot.png", width = 700, height = 600, res = 200, units = "mm")
  po22
  dev.off()

  png("F16_genomeplot.png", width = 700, height = 600, res = 200, units = "mm")
  po16
  dev.off()

  png("F12_genomeplot.png", width = 700, height = 600, res = 200, units = "mm")
  po12
  dev.off()

  png("F3_genomeplot.png", width = 700, height = 600, res = 200, units = "mm")
  po3
  dev.off()

##########################################################################################################################################################################################################################################
### create the plot from above and include the  QTLs for the 6 categories from plot 7 
# F3 in gestrichelt, F23 org organge, F23 con grau
# punkte dran für die 6 Kategorien (Farblich unterschiedlich)
# Form der punkte entsprechend gleich oder unterschiedlich selektiert. 


## specifiy the centromeric regions 
  CO = data.frame(Chr = c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H"), start = c(40,160,35,110,25,30,65), fin = c(420,500,450,420,500,420,480))

## candidate genes
  MajorGenes = data.frame(Chr = c("chr1H", "chr1H", "chr2H", "chr4H", "chr5H", "chr5H", "chr7H", "chr2H"), start = c(500,550,60,630,550,600,45,450), 
                        Gene = c("PpdH2", "HvElf2", "HvPpdH1", "VrnH2", "HvCbf", "VrnH1", "VrnH3", "HvCen"), HAF=0.4, env="else",
                        Col = c("red", "yellow", "green", "blue", "orange", "gray", "black", "violet"))
  MajorGenes$Gene = as.factor(MajorGenes$Gene)


################################################################################
      ay = read.csv("Region_wideAF.csv")
      reg = read.csv("QTL_Regio3.csv")
      
      row.names(ay) = ay[,1]
      ay = ay[,-1]
      ay = as.data.frame(t(ay))
      
      qwe = cbind(reg[,c(1:10)], ay)
      #qwe = qwe[,-c(8, 10, 12, 14)]
      
      qwe[grep("Root", qwe$Category),"Category"] = "Root"
      qwe[grep("Nutrient", qwe$Category),"Category"] = "Nutrient"
      qwe[grep("physiology", qwe$Category),"Category"] = "physiology"
      qwe[grep("drought", qwe$Category),"Category"] = "Drought tolerance"
      qwe[grep("Yield", qwe$Category),"Category"] = "Yield components"
      qwe[grep("Resistance", qwe$Category),"Category"] = "Biotic resistance"
      
      qwe[qwe$Category == "physiology","Category"] = "Yield physiology"
      
      ### add the positions 
      qwe$Chr = "."
      qwe$start = "."
      
      # add the chromosome 
      for (i in c("1H", "2H", "3H", "4H", "5H", "6H", "7H")){qwe[grep(i ,qwe$Position), "Chr"] = i}
      qwe$Chr = paste0("chr", qwe$Chr)
      
      # extract the physical position 
      qwe$start = as.matrix(c("60", "120cm", "130cm", "20cm", "30cm", "33cm", "48cm", "500", "93cm", "94cm", "97cm", "108cm", "106cm", "120cm", "718", "29", "50cm", "55", "8cm", "677", "80cm", "115cm", "105cm", "106cm", 
                  "109cm", "119", "122cm", "140cm", "138cm", "143cm", "580", "36cm", "39cm", "45cm", "48cm", "45cm", "45", "45", "51cm", "52cm", "55cm", "57cm", "63cm", "634", "655", "655", "667", "7cm", 
                  "8cm", "100cm", "111cm", "618", "613", "51cm", "65cm", "1", "57cm","130cm", "145cm", "25cm", "495", "48cm", "47cm", "580","510", "55cm", "647", "13", "16", "60cm", "50cm", "531", "538",
                  "9cm", "21cm", "134cm", "140cm", "642", "3cm", "77cm", "78cm", "110cm"))
      
      # convert the genomic positions to physical positions
      ref = as.data.frame(fread("Marker_reference_all.txt", header=T))
      
      # extract the genetic positions 
      qwe$genpos = 0
      qwe[grep("cm", qwe$start),"genpos"] = as.integer(str_split_fixed(qwe[grep("cm", qwe$start),"start"], "c", 2)[,1])
      
      # look up a range in the ref table around the cm locus and pick the mean physical position
      qwe$Pos = NA
      for(i in c(1:nrow(qwe))){
                                if(qwe$genpos[i] == 0){
                                                       qwe$Pos[i] = qwe$start[i] 
                                }else{
                                      # subset the chromsome
                                      ref2 = ref[ref$Chr==qwe$Chr[i],]
                                      qwe$Pos[i] = round(mean(ref2[ref2$Pos_genetic > qwe$genpos[i]-3 &ref2$Pos_genetic < qwe$genpos[i]+3,"Pos_phys"]) / 1000000)
                                }
        
      }
      
      qwe$Pos = as.integer(qwe$Pos)
      
      # select only relevant columns 
      candidateLoci = qwe[,c(1,4,6,8,9,10,20,23)]
      
      # add a column which indicates if it is equal or different in the af between the systems 
      candidateLoci$evolve = "."
      for (i in c(1:nrow(candidateLoci))){
                                          if(candidateLoci[i,"wtAF.org.F23"] > candidateLoci[i,"wtAF.conv.F23"] - 0.05 & candidateLoci[i,"wtAF.org.F23"] < candidateLoci[i,"wtAF.conv.F23"] + 0.05){
                                                                                                                                                                                                      candidateLoci[i,"evolve"] = "equal"
                                          }else{
                                                candidateLoci[i,"evolve"] = "segregating"
                                          }
      }
      
      # add a dummy allele frequency 
      candidateLoci$HAF = 0.6
      # convert hte categories to factors 
      candidateLoci$Category = factor(candidateLoci$Category, levels = unique(candidateLoci$Category))
      # adjust the seed domrancy and yield loci on chr 5H
      candidateLoci[65,"Pos"] = 495
      candidateLoci[61, "Pos"] = 485
      # some others are also not perfectly in line with the plot - we have to change this, too.
      candidateLoci = candidateLoci[order(candidateLoci$Chr, candidateLoci$Pos),]
      row.names(candidateLoci) = c(1:nrow(candidateLoci))

      candidateLoci[26:27,"Pos"] = 60 # 3h 
      candidateLoci[26:27,"evolve"] = "equal"
      candidateLoci[49,"Pos"] = 12 # 4h equal 
      candidateLoci[61,"evolve"] = "equal" 
      candidateLoci = candidateLoci[-1,]
           
    ##################
    ## making the plot 
      sink2 = sink 
      sink2$env = as.character(sink2$env)
      sink2[sink2$env == "Organic","env"] = "2 Organic"
      sink2[sink2$env == "Conventional","env"] = "1 Conventional"
      sink2[sink2$env == "Origin","env"] = "0 Origin"
      sink2$env = as.factor(sink2$env)
      
      pb = ggplot(sink2[sink2$gen=="F23" & sink2$Chr!="chrUn" ,]) + geom_rect(data = CO, aes(xmin=start, xmax=fin, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.8) +
        theme_minimal() + geom_point(aes(x=start/1000000, y=HAF, color=env), size=0.4)  + 
        scale_color_manual(values=c( '#999999', '#E69F00', 'chartreuse4', 'dodgerblue', 'darkturquoise', 'darkorange4', 'gold3', 'deeppink'), name=NULL) + 
        guides(colour = guide_legend(override.aes = list(size=8)), shape = guide_legend(override.aes = list(size=8)), fill = guide_legend(override.aes = list(size=8))) +
        xlab("Genome Position [Mb]") + ylab("") + facet_wrap(.~Chr, ncol = 1, strip.position = "bottom") +  
        theme(text = element_text(size = 25), legend.position = "bottom", plot.subtitle = element_text(hjust = -0.03), plot.title = element_text(size = 30, face = "bold"))  + 
        labs(title = "A", subtitle = "Wild form allele frequency [ratio]") 
      
      
      pB = pb + geom_point(data = candidateLoci, aes(x=Pos, y=HAF, color=Category, shape=evolve), size=5) + theme(legend.position = "top") + scale_shape_manual(name=NULL, values = c(15,17)) 
      pB      
      
      png("Candidate_QTLs_in_F23_adjusted.png", width = 600, height = 500, res = 200, units = "mm")
      pB
      dev.off()
      
###############################################################################################################################
      # Category plot as a sub plot - data created in "Phonotype_categories_Plot.jl !!! 
      xy = read.csv("Candidate_Groups_Pheno.csv")
      
      xy[grep("E1", xy$variable),c(9,10)] = ""
      xy[xy$varmin < 0, 8] = 0
      xy[xy$varmax > 35, 7] = 35
      xy[xy$p.signif == "ns",10] = ""

      # get the standard deviation
      xy$low = xy$value - sqrt(xy$variance)
      xy$high =  xy$value + sqrt(xy$variance)
      xy[xy$low < 0,"low"] = 0
      
      xy$gen = factor(xy$gen, levels = c("F3", "F12", "F16", "F22", "F23"))
            
      ## plot

      # p1 = ggplot(xy, aes(x=gen, y=value, shape=Environment, color=Environment)) + geom_point(size=4, stroke = 1.5) + geom_pointrange(aes(ymax = high, ymin =low))  +
      #   facet_wrap(~ay) + xlab("Generation") + ylab("Wild Allele Frequency [%]") + scale_shape_manual(values= c(8,23)) + theme_classic() +
      #   theme(text = element_text(size=20), axis.text.x = element_text(size=12), strip.background = element_rect(
      #     color="white", fill="gray94", size=1.5, linetype="solid"), strip.text.x = element_text(size = 12, color = "black", face = "bold.italic"), legend.position="top", plot.title = element_text(size=25)) +
      #   scale_color_manual(values=c('#999999', '#E69F00'))
      # 
      p2 = ggplot(xy, aes(x=gen, y=value/100, shape=Environment, color=Environment)) + geom_point(size=4, stroke = 1.5) + geom_pointrange(aes(ymax = high/100, ymin =low/100))  +
        geom_text(aes(label = p.signif), y=0.25, color = "black", size = 6) +
        facet_wrap(~ay) + xlab("Generation") + ylab("") + scale_shape_manual(values= c(8,23)) + theme_classic() +
        theme(text = element_text(size=28), axis.text.x = element_text(size=20), strip.background = element_rect(
          color="white", fill="gray94", size=1.5, linetype="solid"), 
          strip.text.x = element_text(size = 20, color = "black", face = "bold.italic"), 
          legend.position="right", axis.ticks.length=unit(.4, "cm"), 
          plot.subtitle = element_text(hjust = -0.08), plot.title = element_text(size = 30, face = "bold")) +
          scale_color_manual(values=c('#999999', '#E69F00'))  + 
        labs(title = "B", subtitle = "Wild form allele frequency [ratio]", color="Agro-ecosystem", shape = "Agro-ecosystem") 
    p2      
      # p3 = ggplot(xy, aes(x=gen, y=value, shape=Environment, color=Environment)) + geom_point(size=4, stroke = 1.5)  +
      #   geom_text(aes(label = p.signif), y=20, color = "black",  size = 6) +
      #   facet_wrap(~ay) + xlab("Generation") + ylab("Wild Allele Frequency [%]") + scale_shape_manual(values= c(8,23)) + theme_classic() +
      #   theme(text = element_text(size=20), axis.text.x = element_text(size=12), strip.background = element_rect(
      #     color="white", fill="gray94", size=1.5, linetype="solid"), strip.text.x = element_text(size = 12, color = "black", face = "bold.italic"), legend.position="top", plot.title = element_text(size=25)) +
      #   scale_color_manual(values=c('#999999', '#E69F00')) + ggtitle("B")
      

      lay <- rbind(c(1,1),
                   c(1,1),
                   c(2,2))
      
      # the other plot with another position of the figure
     pb = ggplot(sink2[sink2$gen=="F23" & sink2$Chr!="chrUn" ,]) + geom_rect(data = CO, aes(xmin=start, xmax=fin, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.8) +
        theme_minimal() + geom_point(aes(x=start/1000000, y=HAF, color=env), size=0.4)  + 
        scale_color_manual(values=c( '#999999', '#E69F00', 'chartreuse4', 'dodgerblue', 'darkturquoise', 'darkorange4', 'gold3', 'deeppink'), name=NULL) + 
        guides(colour = guide_legend(override.aes = list(size=8)), shape = guide_legend(override.aes = list(size=8)), fill = guide_legend(override.aes = list(size=8))) +
        xlab("Genome Position [Mb]") + ylab("") + facet_wrap(.~Chr, ncol = 1, strip.position = "bottom") +  
        theme(text = element_text(size = 25), legend.position = "bottom", plot.subtitle = element_text(hjust = -0.03), plot.title = element_text(size = 30, face = "bold"))  + 
        labs(title = "A", subtitle = "Wild form allele frequency [ratio]")      
      pB = pb + geom_point(data = candidateLoci, aes(x=Pos, y=HAF, color=Category, shape=evolve), size=5) + theme(legend.position = "top") + scale_shape_manual(name=NULL, values = c(15,17)) 
      pB      
      
      png("Candidate_QTLs_in_F23_&_Category_plot_adjusted.png", width = 500, height = 600, res = 200, units = "mm")
      grid.arrange(pB,p2, layout_matrix = lay)
      dev.off()
      
  #################################################################################################################################################################################################################
  ### Grafical summary submission
      
      
      # the other plot with another position of the figure
      pb = ggplot(sink2[(sink2$gen=="F23" | sink2$gen=="F3") & sink2$Chr!="chrUn" ,]) + geom_rect(data = CO, aes(xmin=start, xmax=fin, ymin=0, ymax=Inf), fill = "gray85", color="white", alpha=0.8) +
        theme_minimal() + geom_point(aes(x=start/1000000, y=HAF, color=env), size=0.4)  + 
        scale_color_manual(values=c("black", '#999999', '#E69F00', 'chartreuse4', 'dodgerblue', 'darkturquoise', 'darkorange4', 'gold3', 'deeppink'), name=NULL, 
                           labels = c("F3 - origin", "F23 - conventional", "F23 - organic", levels(candidateLoci$Category))) + 
        guides(colour = guide_legend(override.aes = list(size=8), ncol=3), shape = guide_legend(override.aes = list(size=8)), fill = guide_legend(override.aes = list(size=8))) +
        xlab("Genome Position [Mb]") + ylab("Wild form allele frequency") + facet_wrap(.~Chr, ncol = 1, strip.position = "bottom") +  theme(text = element_text(size = 25), legend.position = "bottom") + ggtitle("A")
      pB = pb + geom_point(data = candidateLoci, aes(x=Pos, y=HAF, color=Category, shape=evolve), size=3) + theme(legend.position = "top") + scale_shape_manual(name=NULL, values = c(15,17)) 
      
      
      png("Candidate_QTLs_in_F23_&_Category_plot_grafical_summary.png", width = 500, height = 600, res = 200, units = "mm")
      grid.arrange(pB,p2, layout_matrix = lay)
      dev.off()
      