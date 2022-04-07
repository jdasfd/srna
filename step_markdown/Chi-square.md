# Using chi-square test

Chi-square test could help us test whether the RNA region covered length was significantly differ to the whole genome covered length.


## Methods

For each sample's region: We tested the following (use the tRNA as an example, could be replaced by rRNA and mRNA)

Test 1:

<table>
    <tr>
        <td align="center">tRNA_covered length</td>
        <td align="center">chrom_covered legnth</td>
    </tr>
    <tr>
        <td align="center">tRNA length</td>
        <td align="center">chrom length</td>
    </tr>
</table>

Because the genome was way more longer than the tRNA-covered region, so we changed it to the ratio for chi-square test. R packages named chisq.test allowed three numbers input.

<table>
    <tr>
        <td align="center">tRNA_covered length</td>
        <td align="center">chrom_covered length</td> 
    </tr>
    <tr>
        <td colspan="2" align="center">tRNA/chrom length ratio</td>
    </tr>
</table>

Test 2:

For different bacterial groups within a sample, we also using chi-square test.

<table>
    <tr>
        <td align="center">tRNA_covered length group i</td>
        <td align="center">tRNA_covered length group j</td>
    </tr>
    <tr>
        <td align="center">tRNA length group i</td>
        <td align="center">tRNA length group j</td>
    </tr>
</table>

$i != j; i, j = 1,2,3,4$


## Mosdepth

Mosdepth is a fast BAM/CRAM depth calculation tool for WGS, exome, or targeted sequencing.

* Mosdepth for counting the coverage.

```bash
mkdir -p /mnt/e/project/srna/output/chi/mosdepth
cd /mnt/e/project/srna/output/chi/mosdepth
```

I do not know if it is my problem in reading help info. There is not any option for deciding output path. So be aware when using the mosdepth.

```bash
parallel -j 3 " \
mosdepth -t 4 {}_aliall \
../../bam/bacteria/{}_aliall.sort.bam \
" ::: $(ls ../../bam/bacteria/*_aliall.sort.bam | \
perl -p -e 's/^.*\/(.+)_.+bam$/$1/')

# -t: threads, it has been said that 4 could reach the max speed
# {}_aliall means prefix of mosdepth, details in help info
```

```bash
parallel -j 3 " \
mosdepth -t 4 {}_1mis \
../../bam/bacteria/{}_1mis.sort.bam \
" ::: $(ls ../../bam/bacteria/*_1mis.sort.bam | \
perl -p -e 's/^.*\/(.+)_.+bam$/$1/')
```

```bash
parallel -j 3 " \
mosdepth -t 4 {}_unali \
../../bam/bacteria/{}_unali.sort.bam \
" ::: $(ls ../../bam/bacteria/*_unali.sort.bam | \
perl -p -e 's/^.*\/(.+)_.+bam$/$1/')
```

* Convert bed to yml

Convert bed format file to runlist files for better manipulating using a perl script. Each covered base will be counted.

```bash
mkdir -p /mnt/e/project/srna/output/chi/yml
cd /mnt/e/project/srna/output/chi/mosdepth

parallel -j 10 " \
zcat {}.per-base.bed | \
perl ../../../script/bed2yml.pl > ../yml/{}.yml \
" ::: $(ls *.per-base.bed.gz | perl -p -e 's/\.per.+gz$//')
```


## spanr for RNA covered length

Use spanr from wang-q [intspan](https://github.com/wang-q/intspan).

* Get all used bacteria genome size in .chr.sizes format.

```bash
cd /mnt/e/project/srna/output/chi

faops size ../../genome/bacteria/bacteria.fna > bacteria.chr.sizes
```

* Get .csv with coverage

The .csv file contains 4 columns, chr, chrLength, size and coverage. We need the column 2 ‘chrLength’ (representing genome length) and the column 3 ‘size’ (representing genome covered length). Column 4 coverage: column 3 (coverage length) / column 2 (chrLength), so it does not mean the reads cover.

```bash
mkdir -p /mnt/e/project/srna/output/chi/opt
cd /mnt/e/project/srna/output/chi/yml

parallel -j 10 " \
spanr stat ../bacteria.chr.sizes {}.yml -o ../opt/{}.csv \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')
```

* Extracte tRNA from .gff file.

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

* Compare two runlist files to get intersect part.

```bash
cd /mnt/e/project/srna/output/chi/yml

parallel -j 10 " \
spanr compare {}.yml ../../../annotation/bacteria/tRNA.yml \
-o ../opt/{}.tRNA.intersect.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare --op diff {}.yml ../../../annotation/bacteria/tRNA.yml \
-o ../opt/{}.tRNA.diff.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare {}.yml ../../../annotation/bacteria/rRNA.yml \
-o ../opt/{}.rRNA.intersect.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare --op diff {}.yml ../../../annotation/bacteria/rRNA.yml \
-o ../opt/{}.rRNA.diff.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare {}.yml ../../../annotation/bacteria/mRNA.yml \
-o ../opt/{}.mRNA.intersect.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare --op diff {}.yml ../../../annotation/bacteria/mRNA.yml \
-o ../opt/{}.mRNA.diff.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')
```

* Get length results

Get the same .csv results including tRNA length and tRNA covered length. rRNA and mRNA regions were used as control.

```bash
mkdir -p /mnt/e/project/srna/output/chi/rna
cd /mnt/e/project/srna/output/chi

spanr stat bacteria.chr.sizes ../../annotation/bacteria/tRNA.yml -o tRNA.csv
spanr stat bacteria.chr.sizes ../../annotation/bacteria/rRNA.yml -o rRNA.csv
spanr stat bacteria.chr.sizes ../../annotation/bacteria/mRNA.yml -o mRNA.csv
```

```bash
cd /mnt/e/project/srna/output/chi/opt

parallel -j 10 " \
spanr stat ../bacteria.chr.sizes {}.tRNA.intersect.yml \
-o ../rna/{}.tRNA.intersect.csv \
" ::: $(ls *.tRNA.intersect.yml | perl -p -e 's/\.tRNA.+yml$//')

parallel -j 10 " \
spanr stat ../bacteria.chr.sizes {}.tRNA.diff.yml \
-o ../rna/{}.tRNA.diff.csv \
" ::: $(ls *.tRNA.diff.yml | perl -p -e 's/\.tRNA.+yml$//')

parallel -j 10 " \
spanr stat ../bacteria.chr.sizes {}.rRNA.intersect.yml \
-o ../rna/{}.rRNA.intersect.csv \
" ::: $(ls *.rRNA.intersect.yml | perl -p -e 's/\.rRNA.+yml$//')

parallel -j 10 " \
spanr stat ../bacteria.chr.sizes {}.rRNA.diff.yml \
-o ../rna/{}.rRNA.diff.csv \
" ::: $(ls *.rRNA.diff.yml | perl -p -e 's/\.rRNA.+yml$//')

parallel -j 10 " \
spanr stat ../bacteria.chr.sizes {}.mRNA.intersect.yml \
-o ../rna/{}.mRNA.intersect.csv \
" ::: $(ls *.mRNA.intersect.yml | perl -p -e 's/\.mRNA.+yml$//')

parallel -j 10 " \
spanr stat ../bacteria.chr.sizes {}.mRNA.diff.yml \
-o ../rna/{}.mRNA.diff.csv \
" ::: $(ls *.mRNA.diff.yml | perl -p -e 's/\.mRNA.+yml$//')
```

* Convert .csv to .tsv format.

```bash
cd /mnt/e/project/srna/output/chi/rna

parallel -j 10 " \
cat {}.csv | mlr --icsv --otsv cat > {}.tsv \
" ::: $(ls *.csv | perl -p -e 's/\.csv$//')

rm *.csv

cd ..

for file in `ls *.csv | perl -p -e 's/\.csv//'`
do
cat ${file}.csv | mlr --icsv --otsv cat > ${file}.tsv;
done

rm *.csv
```

* Use tsv-utils to join tsv together.

```bash
mkdir -p mnt/e/project/srna/output/chi/result
cd /mnt/e/project/srna/output/chi/rna

parallel -j 3 " \
cat {}.tRNA.diff.tsv | tsv-select -f 1,2,3 | \
tsv-join -H --filter-file ../tRNA.tsv --key-fields 1 --append-fields 3 | \
tsv-join -H --filter-file {}.tRNA.intersect.tsv \
--key-fields 1 --append-fields 3 | \
tsv-join --filter-file ../../../rawname.tsv \
--key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 6 --sum 2 --sum 3 --sum 4 --sum 5 | \
tsv-join --filter-file ../../../name.tsv \
--key-fields 1 --append-fields 2 | \
tsv-filter --ne 3:0 --ne 5:0 \
> ../result/{}.tRNA.tsv \
" ::: $(ls *.tRNA.diff.tsv | perl -p -e 's/\.t.+tsv$//')

parallel -j 3 " \
cat {}.rRNA.diff.tsv | tsv-select -f 1,2,3 | \
tsv-join -H --filter-file ../rRNA.tsv --key-fields 1 --append-fields 3 | \
tsv-join -H --filter-file {}.rRNA.intersect.tsv \
--key-fields 1 --append-fields 3 | \
tsv-join --filter-file ../../../rawname.tsv \
--key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 6 --sum 2 --sum 3 --sum 4 --sum 5 | \
tsv-join --filter-file ../../../name.tsv \
--key-fields 1 --append-fields 2 \
> ../result/{}.rRNA.tsv \
" ::: $(ls *.rRNA.diff.tsv | perl -p -e 's/\.r.+tsv$//')

parallel -j 3 " \
cat {}.mRNA.diff.tsv | tsv-select -f 1,2,3 | \
tsv-join -H --filter-file ../mRNA.tsv --key-fields 1 --append-fields 3 | \
tsv-join -H --filter-file {}.mRNA.intersect.tsv \
--key-fields 1 --append-fields 3 | \
tsv-join --filter-file ../../../rawname.tsv \
--key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 6 --sum 2 --sum 3 --sum 4 --sum 5 | \
tsv-join --filter-file ../../../name.tsv \
--key-fields 1 --append-fields 2 \
> ../result/{}.mRNA.tsv \
" ::: $(ls *.mRNA.diff.tsv | perl -p -e 's/\.m.+tsv$//')
```

##  Chi-square test

* Chi-square tests of Test 1

```bash
mkdir -p /mnt/e/project/srna/output/chi/test1
cd /mnt/e/project/srna/output/chi/result

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