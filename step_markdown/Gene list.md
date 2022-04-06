# Extracting gene list from bam files

Because of the necessity to get the tair_locus id for enrichment analysis, the gene targeted by 1mis reads should be extracted according their positions.

## Convert gff to tsv

There were a few columns needed, including chrom, position, annotation.

```bash
mkdir -p /mnt/e/project/srna/output/gene/
cd /mnt/e/project/srna/output/gene/

cat /mnt/e/project/srna/annotation/plant/Atha/Atha.gff | \
grep -v '#' | tsv-filter --str-eq 3:gene | \
tsv-select -f 1,4,5,7,9 > Atha_genelist.tsv
# tsv-select get the specific columns
```

## Convert bam to tsv files

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

## Extract 1mis gene list

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

cat SRR1042171.gene_seq.tsv | tsv-select -f 1 | sort | uniq | \
sed '1iTAIR' > ../genelist/SRR1042171.gene.tsv
```

## clusterProfiler for GO enrichment analysis

Before we started, we should install some packages for R

The raw R packages needed could be seen in <https://github.com/wang-q/dotfiles/blob/master/r/packages.R>. org.At.tair.db were added into path.

```bash
Rscript /mnt/e/project/srna/script/packages.R
```

There were tair_locus id we needed to transform them to GO_id and entrez_id for enrichment. The detailed transforming script could be seen in document script/enrichgo_dotplot.r

```bash
Rscript /mnt/e/project/srna/script/enrichgo_dotplot.r -f SRR1042171.gene.tsv
```