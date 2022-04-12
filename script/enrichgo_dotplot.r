library(readr)
suppressMessages(library(clusterProfiler))
suppressMessages(library(biomaRt))
suppressMessages(library(org.At.tair.db))
library("optparse")

option_list = list(
    make_option(c("-f","--file"), type = "character", default = NULL,
    help = "gene list format", metavar = "character"),
    make_option(c("-o","--out"), type = "character", default = "output.pdf",
    help = "output dotplot pdf script [default = %default]", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
    print_help(opt_parser)
    stop ("At least one argument must be supplied (input file) .n", call.=FALSE)
}

gene <- read_tsv(opt$file, show_col_types = FALSE)
tair <- gene$TAIR

pmart <- useMart("plants_mart", "athaliana_eg_gene", host = "https://plants.ensembl.org")

genelist <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "go_id", "tair_locus"),
                filters = 'tair_locus',
                values = tair,
                mart = pmart)

# attributes: output ID formats, filters: our ID formats, values: our data, mart: chosen mart
# this step is for id transforming, tair_locus id to go_id for enrichGO

go <- enrichGO(gene = genelist$go_id,
                OrgDb = org.At.tair.db,
                ont = "ALL",
                keyType = 'GO',
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)

# extract values from data.frame.

pdf(file = opt$out)
dotplot(go, showCategory = 10)
dev.off()
