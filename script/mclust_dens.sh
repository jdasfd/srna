usage () { echo "bash mlust_dens.sh <count.tsv> <V_group>" 1>&2; exit; }
[ $# -lt 1 ] && usage
TSV_FILE=$1
clust=$2

if [ ! -e ${TSV_FILE} ]; then
    echo 1>&2 "[${TSV_FILE}] is not a file" 1>&2; exit;
fi

if [ "${clust}" -gt "9" ]; then
    echo 1>&2 "Make sure num < 9" 1>&2; exit;
fi

name=${TSV_FILE%%.*}
catgry=${name##*_}

Rscript -e '
library(readr)
library(ggplot2)
suppressMessages(library(mclust))
args <- commandArgs(T)
out <- paste0("../../figure/", args[3], ".density.pdf")

ratio <- read_tsv(args[1], show_col_types = FALSE)
dens <- densityMclust(ratio$count, G = args[2], model = "V")
br <- seq(min(ratio$count), max(ratio$count), length = max(ratio$count) - min(ratio$count))
x <- seq(min(ratio$count)-diff(range(ratio$count))/10, 
max(ratio$count)+diff(range(ratio$count))/10, length = 500)
cdens <- predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 3)

pdf(out, width = 6, height = 4)
plot(dens, what = "density", data = ratio$count, breaks = br)
matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 3)
dev.off()
' ${TSV_FILE} ${clust} ${catgry}

rm Rplots.pdf
