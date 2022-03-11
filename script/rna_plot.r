library(ggplot2)
library(readr)
library("optparse")

option_list = list(
    make_option(c("-f","--file"), type = "character", default = NULL,
    help = "tsv file name", metavar = "character"),
    make_option(c("-o","--out"), type = "character", default = "output.pdf",
    help = "output pdf script in category name [default = %default]", metavar = "character"),
    make_option(c("-t","--title"), type = "character", default = NULL,
    help = "figure main title", metavar = "character"),
    make_option(c("-y","--ylab"), type = "character", default = NULL,
    help = "ylab title", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
    print_help(opt_parser)
    stop ("At least one argument must be supplied (input file) .n", call.=FALSE)
}

num <- read_tsv(opt$file, col_types = cols(group = col_character()), col_names = T)

s <- ggplot (data = num, aes(x = group, y = num, group = group, fill = group)) +
geom_boxplot() +
facet_wrap(~catgry) +
labs(x = "bacterial group", y = opt$ylab) +
ggtitle(opt$title) +
scale_fill_discrete(name="Category", 
                    breaks = c("1","2","3","4"),
                    labels = c("endo/epiphyte", "environment", "gut", "marine"))

ggsave(s, file = opt$out, width = 8, height = 4)