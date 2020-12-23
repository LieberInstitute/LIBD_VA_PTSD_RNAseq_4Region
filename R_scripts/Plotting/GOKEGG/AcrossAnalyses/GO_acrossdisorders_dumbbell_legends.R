library("ggalt")
library("tidyr")
library(dplyr)
library(stringr)


##medialamyg
#Analysis/color legend

load("GO_BMedialAmyg_p.adjust0.1_plotinput.rda")

font <- "Helvetica"
text_color <- "#222222"
analysis_colors <- c("MDD" = "#CB8CAD", "PTSD" = "#000000","onlyPTSD" = "#3881B0","PTSDvsMDD" = "#E3712B")
x_colors<- c(rep("#C24A4A",11),"#222222",rep("#5A9C64",4))


all_medialamyg_sig2$point_sizes<- ifelse(all_medialamyg_sig2$GeneRatioDec < 0.01, NA,ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.01 & all_medialamyg_sig2$GeneRatioDec < 0.05, 4,
	 ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.05 & all_medialamyg_sig2$GeneRatioDec < 0.1, 6, ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.1 & all_medialamyg_sig2$GeneRatioDec < 0.15, 8,
	ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.15 & all_medialamyg_sig2$GeneRatioDec < 0.2, 10, NA)))))
point_sizes_vector <- all_medialamyg_sig2$point_sizes


wrap_MeA <- stringr::str_wrap(all_medialamyg_sig2$Description,55)
wrap_MeA[7] <- "transmembrane receptor protein serine/\nthreonine kinase signaling pathway"

medialamyg_dumbell <- ggplot(data=all_medialamyg_sig2, aes(x=0, xend=final.P.Adj,y = Description)) +
geom_dumbbell(colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) + 
geom_point(data=all_medialamyg_sig2,aes(y=all_medialamyg_sig2$Description, x=all_medialamyg_sig2$final.P.Adj,color=Analysis),size=4.5) +
scale_colour_manual(values=analysis_colors) +
annotate("text", y = all_medialamyg_sig2$Description, x=ifelse(all_medialamyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_medialamyg_sig2$final.P.Adj < 0, 0, 1), 
	label = wrap_MeA,vjust=0.5, lineheight = 0.85,color=text_color, family=font) +
ggtitle("GO MedialAmyg", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
scale_x_continuous(limits=c(-5.5,2),labels = c(seq(5.5, 0, by = -0.5),seq(0.5, 2, by = 0.5)), breaks=seq(-5.5, 2, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.77),
	plot.subtitle = element_text(family=font,size=20,margin=ggplot2::margin(t=0.25,r=0,b=9,l=0),hjust = 0.77),
	legend.background = element_blank(),
	legend.box.background = element_blank(), 
	legend.key=element_blank(),
	legend.title = element_text(family=font,size=18,color=text_color,face="bold"),
	legend.title.align = 0.5,
	legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold", margin = margin(t = 10, b = 5), hjust=0.79),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin = margin(r=-20)),
	axis.text.x = element_text(family=font, size=16,color=x_colors, face="bold"),
	axis.ticks.x=element_line(colour = x_colors, size=1.15),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
ggsave(plot=medialamyg_dumbell,file="png/GO_MedialAmyg_AcrossDisorders_dumbbell_colorlegend.png",width=10,height=8.5,dpi=300)

##BasoAmyg
#Size legend (gene ratio)

load("GO_BasoAmyg_p.adjust0.1_plotinput.rda")
x_colors<- c(rep("#C24A4A",5),"#222222",rep("#5A9C64",5))
font <- "Helvetica"
text_color <- "#222222"
analysis_colors <- c("onlyPTSD" = "#3881B0","PTSDvsMDD" = "#E3712B")
									
pointsizes <- c("0.05" = 4, "0.1" = 6,"0.15" = 8,"0.2" = 10)

all_basoamyg_sig2$point_sizes<- ifelse(all_basoamyg_sig2$GeneRatioDec < 0.01, NA,ifelse(all_basoamyg_sig2$GeneRatioDec >= 0.01 & all_basoamyg_sig2$GeneRatioDec < 0.05, 4,
	 ifelse(all_basoamyg_sig2$GeneRatioDec >= 0.05 & all_basoamyg_sig2$GeneRatioDec < 0.1, 6, ifelse(all_basoamyg_sig2$GeneRatioDec >= 0.1 & all_basoamyg_sig2$GeneRatioDec < 0.15, 8,
	ifelse(all_basoamyg_sig2$GeneRatioDec >= 0.15 & all_basoamyg_sig2$GeneRatioDec < 0.2, 10, NA)))))

all_basoamyg_sig2$Size <- ifelse(all_basoamyg_sig2$point_sizes == 4, 0.05, ifelse(all_basoamyg_sig2$point_sizes == 6, 0.1, 
	ifelse(all_basoamyg_sig2$point_sizes == 8, 0.15, ifelse(all_basoamyg_sig2$point_sizes == 10, 0.2, NA))))

all_basoamyg_sig2$Size <-factor(all_basoamyg_sig2$Size, levels = c("0.05","0.1","0.15","0.2"))

										
basoamyg_dumbell <- ggplot(data=all_basoamyg_sig2, aes(x=0, xend=final.P.Adj,y = Description)) +
geom_dumbbell(colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) + 
geom_dumbbell(data=subset(all_basoamyg_sig2,rownames(all_basoamyg_sig2) == "2567"), colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) +
geom_point(data=all_basoamyg_sig2,aes(y=all_basoamyg_sig2$Description, x=all_basoamyg_sig2$final.P.Adj,color=Analysis,size=Size)) +
scale_colour_manual(values=analysis_colors) +
scale_size_manual(values=pointsizes) +
annotate("text", y = all_basoamyg_sig2$Description, x=ifelse(all_basoamyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_basoamyg_sig2$final.P.Adj < 0, 0, 1), 
	label = stringr::str_wrap(all_basoamyg_sig2$Description,70),color=text_color, family=font, vjust=0.5,lineheight = 0.85) +
ggtitle("GO BasoAmyg", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
labs(size="Gene Ratio") +
scale_x_continuous(limits=c(-2.5,2.5),labels = c(seq(2.5, 0, by = -0.5),seq(0.5, 2.5, by = 0.5)), breaks=seq(-2.5, 2.5, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=20,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.background = element_blank(),
	legend.box.background = element_blank(), 
	legend.key=element_blank(),
	legend.title = element_text("Gene Ratio",family=font,size=18,color=text_color,face="bold"),
	legend.title.align = 0.5,
	legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold", margin = margin(t = 10)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=16,color=x_colors, face="bold"),
	axis.ticks.x=element_line(colour = x_colors, size=1.15),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
ggsave(plot=basoamyg_dumbell,file="png/GO_BasoAmyg_AcrossDisorders_dumbbell_sizelegend.png",width=10,height=8.5,dpi=300)
