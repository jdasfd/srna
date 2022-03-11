library(ggplot2)
library(readr)
library("optparse")

option_list = list(
    make_option(c("-d","--depth"), type = "character", default = NULL,
    help = "depth sum tsv", metavar = "character"),
    make_option(c("-o","--out"), type = "character", default = "output.pdf",
    help = "output pdf script in category name [default = %default]", metavar = "character"),
    make_option(c("-t","--title"), type = "character", default = NULL,
    help = "output script RNA for maintitle", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$depth)){
    print_help(opt_parser)
    stop ("At least one argument must be supplied (input depth file) .n", call.=FALSE)
}

line <- read_tsv(opt$depth, col_types = cols(group = col_character()), col_names = T)

pline <- ggplot(data = line, aes(x = pos, y = depth, group = group, color = group)) +
geom_line() +
labs(x = "relative position", y = "depth") +
ggtitle(opt$title) +
scale_x_continuous(breaks = seq (0, 10, 1)) +
theme(legend.position = "none")

ggsave(pline, file = opt$out, width = 6, height = 4)