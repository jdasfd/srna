# Using chi-square test

## Chi-square test for testing tRNA reads between plant and bacteria

```bash
parallel --colsep '\t' -j 1 -k '
    echo "==> {1}"
    Rscript -e "
        x <- matrix(c({2},{4},{3},{5}), nrow=2)
        x
        chisq.test(x)
        "
'
```
### Mosdepth

Mosdepth for counting the coverage.

```bash
cd /mnt/e/project/srna/output/bam/bacteria

parallel -j 3 " \
mosdepth -t 4 {}_aliall {}_aliall.sort.bam \
" ::: $(ls *_aliall.sort.bam | perl -p -e 's/_aliall\.sort\.bam//')
# -t: threads, it has been said that 4 could reach the max speed

mv *_aliall.per-base.bed* *.txt ../../mosdepth
```

```bash
parallel -j 3 " \
mosdepth -t 4 {}_1mis {}_1mis.sort.bam \
" ::: $(ls *_1mis.sort.bam | perl -p -e 's/_1mis\.sort\.bam//')

mv *_1mis.per-base.bed* *.txt ../../mosdepth
```

```bash
parallel -j 3 " \
mosdepth -t 4 {}_unali {}_unali.sort.bam \
" ::: $(ls *_unali.sort.bam | perl -p -e 's/_unali\.sort\.bam//')

mv *_unali.per-base.bed* *.txt ../../mosdepth
```

Unzip all per-base.bed.gz.

```bash
cd /mnt/e/project/srna/output/mosdepth

gzip -d *.per-base.bed.gz
```

Convert bed format file to runlist files for better manipulating using a perl script. Each covered base will be counted.

```bash
parallel -j 10 " \
cat {}.per-base.bed | perl ../../script/bed2yml.pl > {}.yml \
" ::: $(ls *.per-base.bed | perl -p -e 's/\.per-base\.bed//')
```

### spanr  for tRNA coverage length

Use spanr from wang-q [intspan](https://github.com/wang-q/intspan).

Get all used bacteria genome size in .chr.sizes format.

```bash
mkdir -p /mnt/e/project/srna/output/opt
mkdir -p /mnt/e/project/srna/output/result
cd /mnt/e/project/srna/output/opt

faops size ../../genome/bacteria/bacteria.fna > bacteria.chr.sizes
```

Calculate the coverage.

The .csv file contains 4 columns, chr, chrLength, size and coverage. We need the column 2 ‘chrLength’ (representing genome length) and the column 3 ‘size’ (representing genome covered length). Column 4 coverage: column 3 (coverage length) / column 2 (chrLength), so it does not mean the reads cover.

```bash
cd /mnt/e/project/srna/output/mosdepth

parallel -j 10 " \
spanr stat ../opt/bacteria.chr.sizes {}.yml -o ../opt/{}.csv \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')
```

Extracte tRNA from .gff file.

```bash
cd /mnt/e/project/srna/annotation/bacteria

spanr gff bacteria.gff --tag tRNA > tRNA.yml
# --tag: selected gene name

spanr gff bacteria.gff --tag rRNA > rRNA.yml

spanr gff bacteria.gff --tag CDS > mRNA_all.yml
spanr compare --op diff mRNA_all.yml tRNA.yml -o mRNA1.yml
spanr compare --op diff mRNA1.yml rRNA.yml -o mRNA.yml
rm mRNA_all.yml mRNA1.yml
# mRNA region without tRNA and rRNA
# spanr compare could manipulate aggregation, including intersect (default), union, diff or xor
```

Compare two runlist files for their intersect part.

```bash
cd /mnt/e/project/srna/output/mosdepth

parallel -j 10 " \
spanr compare {}.yml ../../annotation/bacteria/tRNA.yml \
-o ../opt/{}.tRNA.intersect.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare --op diff {}.yml ../../annotation/bacteria/tRNA.yml \
-o ../opt/{}.tRNA.diff.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare {}.yml ../../annotation/bacteria/rRNA.yml \
-o ../opt/{}.rRNA.intersect.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare --op diff {}.yml ../../annotation/bacteria/rRNA.yml \
-o ../opt/{}.rRNA.diff.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare {}.yml ../../annotation/bacteria/mRNA.yml \
-o ../opt/{}.mRNA.intersect.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare --op diff {}.yml ../../annotation/bacteria/mRNA.yml \
-o ../opt/{}.mRNA.diff.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')
```

Get the same .csv results including tRNA length and tRNA covered length. rRNA and mRNA regions were used as control.

```bash
cd /mnt/e/project/srna/output/opt

spanr stat bacteria.chr.sizes ../../annotation/bacteria/tRNA.yml -o tRNA.csv
spanr stat bacteria.chr.sizes ../../annotation/bacteria/rRNA.yml -o rRNA.csv
spanr stat bacteria.chr.sizes ../../annotation/bacteria/mRNA.yml -o mRNA.csv
```

```bash
cd /mnt/e/project/srna/output/opt

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.tRNA.intersect.yml \
-o {}.tRNA.intersect.csv \
" ::: $(ls *.tRNA.intersect.yml | perl -p -e 's/\.tRNA.+yml$//')

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.tRNA.diff.yml \
-o {}.tRNA.diff.csv \
" ::: $(ls *.tRNA.diff.yml | perl -p -e 's/\.tRNA.+yml$//')

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.rRNA.intersect.yml \
-o {}.rRNA.intersect.csv \
" ::: $(ls *.rRNA.intersect.yml | perl -p -e 's/\.rRNA.+yml$//')

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.rRNA.diff.yml \
-o {}.rRNA.diff.csv \
" ::: $(ls *.rRNA.diff.yml | perl -p -e 's/\.rRNA.+yml$//')

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.mRNA.intersect.yml \
-o {}.mRNA.intersect.csv \
" ::: $(ls *.mRNA.intersect.yml | perl -p -e 's/\.mRNA.+yml$//')

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.mRNA.diff.yml \
-o {}.mRNA.diff.csv \
" ::: $(ls *.mRNA.diff.yml | perl -p -e 's/\.mRNA.+yml$//')
```

Convert .csv to .tsv format.

```bash
parallel -j 10 " \
cat {}.csv | csv2tsv -H > {}.tsv \
" ::: $(ls *.csv | perl -p -e 's/\.csv$//')

rm *.csv *.yml
```

Use tsv-utils join tsv together.

```bash
cd /mnt/e/project/srna/output/opt

parallel -j 3 " \
cat {}.tRNA.diff.tsv | cut -f 1,2,3 | \
tsv-join -H --filter-file tRNA.tsv --key-fields 1 --append-fields 3 | \
tsv-join -H --filter-file {}.tRNA.intersect.tsv --key-fields 1 --append-fields 3 | \
tsv-join --filter-file ../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 6 --sum 2 --sum 3 --sum 4 --sum 5 | \
tsv-join --filter-file ../../name.tsv --key-fields 1 --append-fields 2 \
> ../result/{}.tRNA.tsv \
" ::: $(ls *.tRNA.*.tsv | perl -p -e 's/\..*tsv$//' | uniq)

parallel -j 3 " \
cat {}.rRNA.diff.tsv | cut -f 1,2,3 | \
tsv-join -H --filter-file rRNA.tsv --key-fields 1 --append-fields 3 | \
tsv-join -H --filter-file {}.rRNA.intersect.tsv --key-fields 1 --append-fields 3 | \
tsv-join --filter-file ../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 6 --sum 2 --sum 3 --sum 4 --sum 5 | \
tsv-join --filter-file ../../name.tsv --key-fields 1 --append-fields 2 \
> ../result/{}.rRNA.tsv \
" ::: $(ls *.rRNA.*.tsv | perl -p -e 's/\..*tsv$//' | uniq)

parallel -j 3 " \
cat {}.mRNA.diff.tsv | cut -f 1,2,3 | \
tsv-join -H --filter-file mRNA.tsv --key-fields 1 --append-fields 3 | \
tsv-join -H --filter-file {}.mRNA.intersect.tsv --key-fields 1 --append-fields 3 | \
tsv-join --filter-file ../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 6 --sum 2 --sum 3 --sum 4 --sum 5 | \
tsv-join --filter-file ../../name.tsv --key-fields 1 --append-fields 2 \
> ../result/{}.mRNA.tsv \
" ::: $(ls *.mRNA.*.tsv | perl -p -e 's/\..*tsv$//' | uniq)
```

###  Chi-square test

Use chi-square tests performing chi-squared contingency  interspecies.

```bash
mkdir -p /mnt/e/project/srna/output/chi/inter
cd /mnt/e/project/srna/output/result

for tsv in `ls *.tsv`
do
cat $tsv | perl ../../script/interspecie.pl | \
parallel --colsep '\t' -j 1 -k '
    echo "==> {1}"
    Rscript -e "
        x <- matrix(c({2},{4},{3},{5}), nrow=2)
        x
        chisq.test(x)
        "
' > ../chi/inter/$tsv.inter.txt
done
```

Use a simple script for better looking.

```bash
cd /mnt/e/project/srna/output/chi/inter

for file in `ls *.inter.txt | perl -p -e 's/\.tsv.+txt$//'`
do
echo "${file}" >> ../result.tsv
cat ${file}.tsv.inter.txt | perl ../../../script/square.pl >> ../result.tsv
done
```


