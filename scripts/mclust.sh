#!/usr/bin/env bash

usage () { echo "bash mlust.sh <count.tsv>" 1>&2; exit; }
[ $# -lt 1 ] && usage
TSV_FILE=$1

if [ ! -e ${TSV_FILE} ]; then
    echo 1>&2 "[${TSV_FILE}] is not a file" 1>&2; exit;
fi

name=${TSV_FILE%%.*}
catgry=${name##*_}

Rscript -e '
library(readr)
suppressMessages(library(mclust))
args <- commandArgs(T)
outB <- paste0(args[2], ".BIC.tmp.tsv")
outI <- paste0(args[2], ".ICL.tmp.tsv")

ratio <- read_tsv(args[1], show_col_types = FALSE)
dens <- Mclust(ratio$count)
BIC <- mclustBIC(ratio$count)
ICL <- mclustICL(ratio$count)
write.table(BIC, outB, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(ICL, outI, row.names = FALSE, col.names = TRUE, sep = "\t")
' ${TSV_FILE} ${catgry};

cat ${catgry}.BIC.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tBICE\n"; print "$i\t$a[1]\tBICV\n";
$i++;}' | tsv-filter --str-ne 2:NA > ${catgry}.BIC2.tmp.tsv;

cat ${catgry}.ICL.tmp.tsv | perl -n -e 'chomp;
@a = split/\t/, $_;if($_ =~ /"/){$i = 1;next;}else{
print "$i\t$a[0]\tICLE\n"; print "$i\t$a[1]\tICLV\n";
$i++;}' | tsv-filter --str-ne 2:NA > ${catgry}.ICL2.tmp.tsv;

cat ${catgry}.ICL2.tmp.tsv ${catgry}.BIC2.tmp.tsv | sed '1inum\tgrade\tmethod' > ${catgry}.mclust.tsv;

if [ -s ${catgry}.mclust.tsv ]; then
    echo "${catgry}.mclust.tsv is ok."
else
    echo 1>%2 "${catgry}.mclust.tsv is empty." 1>&2; exit;
fi

rm *.tmp.*;

Rscript -e '
library(readr)
library(ggplot2)
args <- commandArgs(T)
file <- paste0(args[1], ".mclust.tsv")
out <- paste0("../../figure/", args[1], ".mclust.pdf")

num <- read_tsv(file, show_col_types = FALSE)
p <- ggplot()+
geom_line(data = num, aes(x = num, y = grade, group = method, color = method)) +
scale_x_continuous(breaks = seq(0,9,1))
ggsave(p, file = out, width = 5, height = 4)
' ${catgry};

if [ -e "../../figure/${catgry}.mclust.pdf" ]; then
    echo "figure has been generated."
else
    echo "problem please check."
fi
