# MinorityReport-ManhattanPlot.r
# 20160323 JAHorst 

# NOTE: you may need to install the packages that contain the libraries below. Google them for the latest.
# 1. uptate the INPUT & OUTPUT filenames
# 2. update the chromosome names
# 3. Run this file from the command line with the following command:
#        R CMD BATCH MinorityReport-ManhattanPlot.r

# load in file
input_data <- read.table("__INPUT_FILENAME__", header=TRUE, sep="\t")

library(ggplot2)
library("plyr")
library(scales)
library("cowplot")

# arrange chromosomes in proper order
chromosome_ouput_order <- levels(as.factor(input_data$chromosome))
## NOTE: upate chromosomes as you prefer for your organism
#chromosome_ouput_order <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","plastid","mito")
ordered_input_data<-arrange(transform(input_data,chromosome=factor(chromosome,levels=chromosome_ouput_order)),chromosome)

# set vertical axis limits symmetrically, with at least -3.9 to 3.9 range
extreme_y <- max( -1*min(ordered_input_data$log2_CNV,-3.9), max(ordered_input_data$log2_CNV,3.9))
lower_y = min(pretty(c(-1*extreme_y, extreme_y)))
upper_y = max(pretty(c(-1*extreme_y, extreme_y)))
most_reads <- max(ordered_input_data$mutant_tile_counts, ordered_input_data$wt_tile_counts)

# plot data
CNVs <- qplot() + geom_point( aes(x=start, y=log2_CNV, colour=probability), data=ordered_input_data, size=0.5) + guides(colour = guide_legend(override.aes = list(size=2))) + theme(legend.key=element_rect(fill='white'), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.background = element_rect(fill = 'grey95'), axis.line = element_blank() ) + facet_grid(~chromosome, scales="free_x", space="free_x")  + labs(x="", y="log2(CNV)")  + scale_x_continuous(breaks=seq(0,1e7,1e6), expand=c(0, 0))  + scale_colour_gradient(low="red", high="grey67", trans = "log", breaks=c(1e-3,1e-6,1e-9,1e-12,1e-15) ) + scale_y_continuous(breaks=pretty(c(-1*extreme_y, extreme_y)), limits=c(lower_y,upper_y), expand=c(0, 0))
mut_reads <- qplot() + geom_point( aes(x=start, y=ordered_input_data$mutant_tile_counts, color='Mutant'), data=ordered_input_data, size=0.3) + facet_grid(~chromosome, scales="free_x", space="free_x") + labs(x="", y="reads         ") + guides(colour = guide_legend(override.aes = list(size=2))) + theme(panel.background = element_rect(fill = 'grey95'), legend.key=element_rect(fill='white'), axis.text.x=element_blank(), axis.line = element_blank(), axis.ticks.x=element_blank(), strip.text.x = element_blank()) + scale_x_continuous(breaks=seq(0,1e7,1e6), expand=c(0, 0)) + scale_y_log10( limits=c(0.1,most_reads), labels = trans_format("log10", math_format(10^.x)) ) + labs(color="") + scale_colour_manual(values = c("blue"))
wt_reads  <- qplot() + geom_point( aes(x=start, y=ordered_input_data$wt_tile_counts,     color='Parent'), data=ordered_input_data, size=0.3) + facet_grid(~chromosome, scales="free_x", space="free_x") + labs(x="", y="")           + guides(colour = guide_legend(override.aes = list(size=2))) + theme(panel.background = element_rect(fill = 'grey95'), legend.key=element_rect(fill='white'), axis.text.x=element_blank(), axis.line = element_blank(), axis.ticks.x=element_blank(), strip.text.x = element_blank()) + scale_x_continuous(breaks=seq(0,1e7,1e6), expand=c(0, 0)) + scale_y_log10( limits=c(0.1,most_reads), labels = trans_format("log10", math_format(10^.x)) ) + labs(color="") + scale_colour_manual(values = c("grey20"))

# print it input_data together
pdf("__OUTPUT_FILENAME__.pdf", width = 12, height = 3, onefile=FALSE)
ggdraw() + draw_plot(CNVs, 0, 1/3, 1, 2/3) + draw_plot(mut_reads, 0.002725, .16, .98375, .31) + draw_plot(wt_reads, 0.002725, 0, .98375, .31) + panel_border(remove = TRUE)
