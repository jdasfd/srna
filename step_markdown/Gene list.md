* [Extracting gene list from bam files](#extracting-gene-list-from-bam-files)

    - [TAIR locus gene list from the 1mis bam files](#tair-locus-gene-list-from-the-1mis-bam-files)

        - [Convert gff to tsv](#convert-gff-to-tsv)

        - [Convert bam to tsv files](#convert-bam-to-tsv-files)

    - [clusterProfiler for GO enrichment analysis](#clusterprofiler-for-go-enrichment-analysis)

        - [Prepare clusterProfiler R packages](#prepare-clusterprofiler-r-packages)

        - [Extract 1mis gene list](#extract-1mis-gene-list)

        - [Gene list join with tier](#gene-list-join-with-tier)


# Extracting gene list from bam files

Because of the necessity to get the tair_locus id for enrichment analysis, the gene targeted by 1mis reads should be extracted according their positions.

## TAIR locus gene list from the 1mis bam files

### Convert gff to tsv

There were a few columns needed, including chrom, position, annotation.

```bash
mkdir -p /mnt/e/project/srna/output/gene/
cd /mnt/e/project/srna/output/gene/

cat /mnt/e/project/srna/annotation/plant/Atha/Atha.gff | \
grep -v '#' | tsv-filter --str-eq 3:gene | \
tsv-select -f 1,4,5,7,9 > Atha_genelist.tsv
# tsv-select get the specific columns
```

### Convert bam to tsv files

We only need those specific column for extraction.

```bash
mkdir -p /mnt/e/project/srna/output/gene/bam_tsv
cd /mnt/e/project/srna/output/bam/plant

for file in `ls SRR*_plant1mis.sort.bam | perl -p -e 's/\.sort\.bam$//'`
do
samtools view -@ 10 ${file}.sort.bam | \
perl -n -e 'chomp;if($_=~/XM:i:1/){print"$_\n";}else{next;}' | \
tsv-select -f 3,1,2,4,10 > ../../gene/bam_tsv/${file}.tsv
done

# 'XM:i:1' in sam/bam files represents 1mis alignment
```

## clusterProfiler for GO enrichment analysis

### Prepare clusterProfiler R packages

Before we started, we should install some packages for R

The raw R packages needed could be seen in <https://github.com/wang-q/dotfiles/blob/master/r/packages.R>. org.At.tair.db were added into path.

```bash
Rscript /mnt/e/project/srna/script/packages.R
```

### Extract 1mis gene list

There were only positions in bam files (which have already converted to tsv with our needed message). Using gff contained position <-> gene will help us getting gene list targeted by reads with 1 mismatch.

```bash
cd /mnt/e/project/srna/output/gene/bam_tsv

parallel -j 6 " \
perl /mnt/e/project/srna/script/extract_gene.pl \
-g ../Atha_genelist.tsv -i {}_plant1mis.tsv \
-o {}.gene_seq.tsv \
" ::: $(ls *_plant1mis.tsv | perl -p -e 's/_.+tsv$//')
```

After running this command, we would get tsv files contained gene list targeted by 1mis reads (col 1) and sequence (col 2).

```bash
mkdir -p /mnt/e/project/srna/output/gene/genelist
cd /mnt/e/project/srna/output/gene/bam_tsv

parallel -j 6 " \
cat {}.gene_seq.tsv | tsv-select -f 1 | sort | uniq | \
sed '1iTAIR' > ../genelist/{}.gene.tsv \
" ::: $(ls SRR*.gene_seq.tsv | perl -p -e 's/\.gene.+tsv$//')
```

There were tair_locus id we needed to transform them to GO_id and entrez_id for enrichment. The detailed transforming script could be seen in document script/enrichgo_dotplot.r

```bash
mkdir -p cd /mnt/e/project/srna/output/gene/plant_GO_figure
cd /mnt/e/project/srna/output/gene/genelist

for file in `ls SRR*.gene.tsv | perl -p -e 's/\.g.+tsv$//'`
do
Rscript /mnt/e/project/srna/script/enrichgo_dotplot.r \
-f ${file}.gene.tsv -o ../plant_GO_figure/${file}_GO.pdf;
done
```

There was a problem of dealing with Rscript using parallel, that is the output only saved the last command line using parallel. (The reason why )

### Gene list join with tier

Before this in [Statistical.md](https://github.com/jdasfd/srna/blob/main/step_markdown/Statistical.md), we seperated sequences into 3 tiers, which means the sequences frequently occured over 120 files. After this step, the GO enrichment would tell us gene lists targeted frequently by the 1mis reads (could not aligned to *A. tha* genome).

Meanwhile, the sequences from 3 different tiers allowed us using `tsv-join` to filter those uncommon reads.

```bash
parallel -j 6 " \
cat {}.gene_seq.tsv | tsv-select -f 2,1 | \
tsv-join --filter-file ../../tier/among/tier1.tsv --key-fields 1 | \
sort | uniq > ../tier/{}.gene_seq.tier1.tsv \
" ::: $(ls *.gene_seq.tsv | perl -p -e 's/\.ge.+tsv$//')

parallel -j 6 " \
cat {}.gene_seq.tsv | tsv-select -f 2,1 | \
tsv-join --filter-file ../../tier/among/tier2.tsv --key-fields 1 | \
sort | uniq > ../tier/{}.gene_seq.tier2.tsv \
" ::: $(ls *.gene_seq.tsv | perl -p -e 's/\.ge.+tsv$//')

parallel -j 6 " \
cat {}.gene_seq.tsv | tsv-select -f 2,1 | \
tsv-join --filter-file ../../tier/among/tier3.tsv --key-fields 1 | \
sort | uniq > ../tier/{}.gene_seq.tier3.tsv \
" ::: $(ls *.gene_seq.tsv | perl -p -e 's/\.ge.+tsv$//')
```

```bash
cd ../tier
cat SRR*.gene_seq.tier1.tsv | cut -f 2 > ../among/tier1.gene.tsv
cat SRR*.gene_seq.tier2.tsv | cut -f 2 > ../among/tier2.gene.tsv
cat SRR*.gene_seq.tier3.tsv | cut -f 2 > ../among/tier3.gene.tsv
```

```bash
cd ../among
cat tier1.gene.tsv | sort | uniq | sed '1iTAIR' > ../tier1.gene.tsv
cat tier2.gene.tsv | sort | uniq | sed '1iTAIR' > ../tier2.gene.tsv
cat tier3.gene.tsv | sort | uniq | sed '1iTAIR' > ../tier3.gene.tsv
```

```bash
Rscript /mnt/e/project/srna/script/enrichgo_dotplot.r -f tier1.gene.tsv -o plant_GO_figure/tier1_GO.pdf
Rscript /mnt/e/project/srna/script/enrichgo_dotplot.r -f tier2.gene.tsv -o plant_GO_figure/tier2_GO.pdf
Rscript /mnt/e/project/srna/script/enrichgo_dotplot.r -f tier3.gene.tsv -o plant_GO_figure/tier3_GO.pdf
```

## Gene list information

* Gene number of each file

```bash
cd /mnt/e/project/srna/output/gene/genelist

echo "file,gene_num" > ../gene_count.csv;
for file in `ls SRR*.gene.tsv | perl -p -e 's/\.ge.+tsv$//'`
do
num=`cat ${file}.gene.tsv | sed '1d' | sort | uniq | wc -l`;
echo "$file","$num" >> ../gene_count.csv;
done
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
gene <- read.csv(args[1])
p <- ggplot(gene, aes(x = file, y = gene_num)) +
geom_bar(stat = "identity") +
labs(x = "Seq files", y = "Gene number") +
theme(axis.text.x = element_text(size = 4, angle = 90))
ggsave(p, file = "/mnt/e/project/srna/output/figure/gene_count.pdf", width = 12, height = 4)
' gene_count.csv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
gene <- read.csv(args[1])
p <- ggplot(gene, aes(x = file, y = gene_num)) +
geom_bar(stat = "identity") +
labs(x = "Seq files", y = "Gene number") +
theme(axis.text.x = element_text(size = 4, angle = 90))
ggsave(p, file = "/mnt/e/project/srna/output/figure/gene_count_50.pdf", width = 12, height = 4)
' gene_count_50.csv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
library(dplyr)
library(forcats)
args <- commandArgs(T)
gene <- read_tsv(args[1], show_col_types = FALSE)
p <- gene %>%
mutate(gene = fct_reorder(gene, desc(gene_num))) %>%
ggplot(aes(x = gene, y = gene_num)) +
geom_bar(stat = "identity", position = "dodge") +
labs(x = "Gene", y = "File count") +
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(p, file = "/mnt/e/project/srna/output/figure/all_gene_count.pdf", width = 9, height = 4)
' all_gene_count.tsv
```

```bash
cat all_gene_count.tsv | tsv-filter -H --ge gene_num:100 > all_gene_count_100.tsv

cat all_gene_count_100.tsv | tsv-select -H -f gene | \
sed '1d' | sed '1iTAIR' > ../all_gene_count_100.tsv
```

```bash
Rscript /mnt/e/project/srna/script/enrichgo_dotplot.r -f all_gene_count_100.tsv -o plant_GO_figure/100_GO.pdf
```