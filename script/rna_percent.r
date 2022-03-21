# args <- commandArgs(T)

library(ggplot2)
library(readr)
library("optparse")
library(gridExtra)

option_list = list(
    make_option(c("-f","--file"), type = "character", default = NULL,
    help = "tsv file name", metavar = "character"),
    make_option(c("-o","--out"), type = "character", default = "output.pdf",
    help = "output pdf script in category name [default = %default]", metavar = "character"),
    make_option(c("-t","--title"), type = "character", default = NULL,
    help = "output script RNA for maintitle", metavar = "character"),
    make_option(c("-y","--ylab"), type = "character", default = NULL,
    help = "character using in ylab representing numerator", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
    print_help(opt_parser)
    stop ("At least one argument must be supplied (input file) .n", call.=FALSE)
}

rna <- read_tsv(opt$file, col_types = cols(group = col_character()), col_names = T)
ytitle <- paste(opt$ylab,"Bac-reads percent", sep = "/")

plot1 <- ggplot (data = rna, aes(x = group, y = ratio, group = group, fill = group)) +
geom_boxplot() + 
geom_jitter(color = 'black', alpha = 0.1, show.legend = FALSE) +
facet_wrap(~catgry) +
ggtitle(opt$title) +
labs(x = "Bacterial group", y = ytitle) +
# scale_y_continuous(limits = c(0, 5)) +
scale_fill_discrete(name="Category", 
                    breaks = c("1","2","3","4"),
                    labels = c("endo/epiphyte", "environment", "gut", "marine"))

alltitle <- paste0(opt$title,"_all")

plot2 <- ggplot (data = rna, aes(x = catgry, y = ratio, fill = catgry)) +
geom_boxplot() +
theme(legend.position = 'none') +
ggtitle(alltitle) +
labs(x = "reads group", y = ytitle) +
scale_fill_manual(values = colors()[6:8])

pdf(opt$out, width = 12, height = 3)
grid.arrange(plot1, plot2, ncol = 2, layout_matrix = rbind(c(1,1,1,2)))
dev.off()
