args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(reshape)



# Make scatter plot showing correlation of standard -log10 pvalues and slope -log10 pvalues
make_scatter_plot <- function(pvalue_df, output_file) {
    corr_coef <- cor(pvalue_df$standard, pvalue_df$slope)
    #PLOT!
    scatter <- ggplot(pvalue_df, aes(x = slope, y = standard)) + geom_point(alpha=.06,size=.00001) 
    scatter <- scatter + geom_smooth(method='lm')
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(x = "Expression Slope Dynamic QTL -log10(pvalue)", y = "Standard Dynamic QTL -log10(pvalue)", title=paste0("rho=",corr_coef))
    #scatter <- scatter + geom_abline() 

    ggsave(scatter, file=output_file, width=20, height=10.5, units="cm")
}







dynamic_qtl_results_comparison_file <- args[1]
output_root <- args[2]

data <- read.table(dynamic_qtl_results_comparison_file, header=TRUE, sep = "\t")

standard_pvalue <- -log10(data$standard_dynamic_qtl_pvalue)

slope_pvalue <- -log10(data$expression_slope_dynamic_qtl_pvalue)

pvalue_df <- data.frame(standard=standard_pvalue,slope=slope_pvalue)


######################
# Make scatter plot showing correlation of standard -log10 pvalues and slope -log10 pvalues
######################
output_file <- paste0(output_root, "_standard_dynamic_vs_expression_slope_dynamic_scatter.png")
make_scatter_plot(pvalue_df, output_file)

