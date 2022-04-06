# Extracting gene list from bam files

Because of the necessity to get the 

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

## clusterProfiler for GO enrichment analysis

There were tair_locus id we needed to transform them to GO_id and entrez_id for enrichment. The detailed transforming script could be seen in document script/enrichgo_dotplot.r

```bash
