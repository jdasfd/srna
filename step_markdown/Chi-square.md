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

For different bacterial groups of the same RNA regions within a sample, we also using chi-square test.

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


For diffrent RNA regions of the same group within a sample, we also using chi-square test

<table>
    <tr>
        <td align="center">group i tRNA_covered length</td>
        <td align="center">group i rRNA_covered length</td>
    </tr>
    <tr>
        <td align="center">group i tRNA length</td>
        <td align="center">group i rRNA length</td>
</table>

$i = 1,2,3,4$

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

spanr gff bac_trna.gff > tRNA.yml

spanr gff bac_rrna.gff > rRNA.yml

spanr gff bac_mrna.gff > mRNA.yml
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
--key-fields 1 --append-fields 2 \
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

### Chi-square tests of Test 1

* Use goodness-of-fit tests performing genome covered length and tRNA covered length consistent with expectation.


```bash
cd /mnt/e/project/srna/output/chi/result

for tsv in `ls *.tsv`
do
name=${tsv%.*};
cat ${tsv} | tsv-summarize --sum 2,3,4,5 | \
awk -v name=$name '{print name"\t"$1"\t"$2"\t"$3"\t"$4}' \
>> ../test1.all.tsv
done
```

```bash
cat test1.all.tsv | head -n 3 | tsv-select -f 1,2,4 | \
perl -n -e '$_=~s/SRR.*_1mis\.(.+)/$1/;print"$_";' | mlr --itsv --omd cat
```

| mRNA | 620294455 | 524734478 |
| --- | --- | --- |
| rRNA | 602358469 | 3641469 |
| tRNA | 608199399 | 3614083 |


```bash
for tsv in `ls *.tsv`
do
name=${tsv%.*};
cat ${tsv} | tsv-summarize --group-by 6 --sum 2,3,4,5 | \
awk -v name=$name '{print name"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}' \
>> ../test1.group.tsv
done
```


```bash
cd ..

cat test1.all.tsv | tsv-filter --ne 3:0 --ne 5:0 | \
parallel --colsep '\t' -j 1 -k '
        echo "==> {1}"
        Rscript -e "
            x <- c({5}, {3})
            a <- {4}/{2}
            b <- 1-({4}/{2})
            y <- matrix(c(x, a, b), nrow=2)
            old.warn <- options()$warn
            options(warn = -1)
            y
            chisq.test(x, p = c(a, b))
        "
' > test1.chi.all.txt

cat test1.group.tsv | tsv-filter --ne 4:0 --ne 6:0 | \
parallel --colsep '\t' -j 1 -k '
        echo "==> {1} of group{2}"
        Rscript -e "
            x <- c({6}, {4})
            a <- {5}/{3}
            b <- 1-({5}/{3})
            y <- matrix(c(x, a, b), nrow=2)
            old.warn <- options()$warn
            options(warn = -1)
            y
            chisq.test(x, p = c(a, b))
        "
' > test1.chi.group.txt
```

### Chi-square tests of Test 2

* Use chi-square tests for the same RNA region

```bash
cd /mnt/e/project/srna/output/chi/result

for file in `ls *.tsv | perl -p -e 's/\.tsv//'`
do
name=${file%.*};
rna=${file#*.};
cat ${file}.tsv | tsv-summarize --group-by 6 --sum 2,3,4,5 | \
tsv-select -f 1,4,5 | sort -nk 1 | \
perl -n -e '{chomp;print "$_\t"}' | \
perl -n -e '{$_=~s/\t$//;print"$_\n";}' | \
awk -v name=$name -v rna=$rna '{print name"\t"rna"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' \
>> ../test2.rna.tsv;
done
```

```bash
cd ..

cat test2.rna.tsv | \
parallel --colsep '\t' -j 1 -k '
    echo "==> {1}.{2} {3}vs{6}"
    Rscript -e "
        old.warn <- options()$warn
        options(warn = -1)
        x <- matrix(c({5},{4},{8},{7}), nrow=2)
        x
        chisq.test(x)
        "
    echo "==> {1}.{2} {3}vs{9}"
    Rscript -e "
        old.warn <- options()$warn
        options(warn = -1)
        x <- matrix(c({5},{4},{11},{10}), nrow=2)
        x
        chisq.test(x)
        "
    echo "==> {1}.{2} {3}vs{12}"
    Rscript -e "
        old.warn <- options()$warn
        options(warn = -1)
        x <- matrix(c({5},{4},{14},{13}), nrow=2)
        x
        chisq.test(x)
        "
' > test2.chi.rna.txt
```

* Use chi-square tests for the same bacterial group

```bash
for file in `ls *.tsv | perl -p -e 's/\.tsv//'`
do
name=${file%.*};
rna=${file#*.};
cat ${file}.tsv | tsv-summarize --group-by 6 --sum 2,3,4,5 | \
tsv-select -f 1,4,5 | sort -nk 1 | \
awk -v name=$name -v rna=$rna \
'{print name"\t"rna"\t"$1"\t"$2"\t"$3}'>> ../test2.rna_all.tsv;
done

cd ..

cat test2.rna_all.tsv | \
perl ../../scripts/extract_chi_group.pl | \
sed '1,4d' > test2.group.tsv

rm test2.rna_all.tsv
```

```bash
cat test2.group.tsv | \
parallel --colsep '\t' -j 1 -k '
    echo "==> {1}.group{2} {9}vs{3}"
    Rscript -e "
        old.warn <- options()$warn
        options(warn = -1)
        x <- matrix(c({11},{10},{5},{4}), nrow=2)
        x
        chisq.test(x)
        "
    echo "==> {1}.group{2} {9}vs{6}"
    Rscript -e "
        old.warn <- options()$warn
        options(warn = -1)
        x <- matrix(c({11},{10},{8},{7}), nrow=2)
        x
        chisq.test(x)
        "
    echo "==> {1}.group{2} {3}vs{6}"
    Rscript -e "
        old.warn <- options()$warn
        options(warn = -1)
        x <- matrix(c({5},{4},{8},{7}), nrow=2)
        x
        chisq.test(x)
        "    
' > test2.chi.group.txt
```

## Chi-square value extraction and plot

* Use a simple script extracting chi-square value from the results

```bash
cat test1.chi.all.txt | perl ../../scripts/chi_square_extraction.pl | \
sed '1ifile\trna\tchi' > test1.chi.all.tsv

cat test1.chi.group.txt | perl ../../scripts/chi_square_extraction.pl | 
perl -n -e 'chomp;
@a = split/\t/,$_;
$a[1]=~/^(.+)\sof/;$rna=$1;
$a[1]=~/\sof(.+)$/;$group=$1;
print "$a[0]\t$group\t$rna\t$a[2]\n";
' | sed '1ifile\tgroup\trna\tchi' > test1.chi.group.tsv

cat test2.chi.rna.txt | perl ../../scripts/chi_square_extraction.pl | \
perl -n -e 'chomp;
@a = split/\t/,$_;
$a[1]=~/^(.+)\s/;$rna=$1;
$a[1]=~/\s(.+)$/;$group=$1;
print "$a[0]\t$group\t$rna\t$a[2]\n";
' | sed '1ifile\tgroup\trna\tchi' > test2.chi.rna.tsv

cat test2.chi.group.txt | perl ../../scripts/chi_square_extraction.pl | \
perl -n -e 'chomp;
@a = split/\t/,$_;
$a[1]=~/^(.+)\s/;$group=$1;
$a[1]=~/\s(.+)$/;$rna=$1;
print "$a[0]\t$rna\t$group\t$a[2]\n";
' | sed '1ifile\trna\tgroup\tchi' > test2.chi.group.tsv
```

```bash
cat test1.chi.all.tsv | \
tsv-summarize -H --group-by rna --mean chi | \
mlr --itsv --omd cat
```

| rna | chi_mean |
| --- | --- |
| mRNA | 735650.497818 |
| rRNA | 16293788.7036 |
| tRNA | 487963.950754 |

```bash
cat test1.chi.group.tsv | \
tsv-summarize -H --group-by rna,group --mean chi | mlr --itsv --omd cat
```

| rna | group | chi_mean |
| --- | --- | --- |
| mRNA |  group3 | 216156.605921 |
| mRNA |  group1 | 270605.336612 |
| mRNA |  group2 | 289235.392047 |
| mRNA |  group4 | 1016.08996373 |
| rRNA |  group3 | 4018554.83301 |
| rRNA |  group1 | 5952724.08608 |
| rRNA |  group2 | 6759234.38281 |
| rRNA |  group4 | 58027.6115616 |
| tRNA |  group3 | 531174.634131 |
| tRNA |  group1 | 947827.700424 |
| tRNA |  group2 | 78130.1847635 |
| tRNA |  group4 | 2187.89569302 |


* Plot

```bash
Rscript -e '
library(ggplot2)
library(readr)

args <- commandArgs(T)

chi <- read_tsv(args[1], show_col_types = FALSE)

plot1 <- ggplot (data = chi, aes(x = rna, y = chi)) +
geom_boxplot() + 
geom_jitter(color = "black", alpha = 0.1, show.legend = FALSE) +
labs(x = "RNA region", y = "Chi square (RNA/genome)") +
scale_y_continuous(limits = c(0,5000000))

pdf("test1.chi.all.pdf")
plot(plot1)
dev.off()
' test1.chi.all.tsv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)

args <- commandArgs(T)

chi <- read_tsv(args[1], show_col_types = FALSE)

plot1 <- ggplot (data = chi, aes(x = rna, y = chi, group = rna, fill = rna)) +
geom_boxplot() + 
geom_jitter(color = "black", alpha = 0.1, show.legend = FALSE) +
facet_wrap(~group) +
labs(x = "RNA region", y = "Chi square (RNA/genome)") +
scale_y_continuous(limits = c(0,2000000))

pdf("test1.chi.group.pdf")
plot(plot1)
dev.off()
' test1.chi.group.tsv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)

args <- commandArgs(T)

chi <- read_tsv(args[1], show_col_types = FALSE)

plot1 <- ggplot (data = chi, aes(x = group, y = chi, group = group, fill = group)) +
geom_boxplot() + 
geom_jitter(color = "black", alpha = 0.1, show.legend = FALSE) +
facet_wrap(~rna) +
labs(x = "group", y = "Chi square") +
scale_y_continuous(limits = c(0,5000000))

pdf("test2.chi.group.pdf", width = 10, height = 5)
plot(plot1)
dev.off()
' test2.chi.group.tsv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)

args <- commandArgs(T)

chi <- read_tsv(args[1], show_col_types = FALSE)

plot1 <- ggplot (data = chi, aes(x = rna, y = chi, group = rna, fill = rna)) +
geom_boxplot() + 
geom_jitter(color = "black", alpha = 0.1, show.legend = FALSE) +
facet_wrap(~group) +
labs(x = "RNA region", y = "Chi square") +
scale_y_continuous(limits = c(0,25000))

pdf("test2.chi.rna.pdf", width = 12, height = 5)
plot(plot1)
dev.off()
' test2.chi.rna.tsv
```

```bash
cat test2.rna.tsv | head -n 3 | \
tsv-select -f 2,3,4,6,7,9,10,12,13 | mlr --itsv --omd cat
```

| mRNA | 1 | 163640953 | 2 | 198309617 | 3 | 150084354 | 4 | 12699554 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rRNA | 1 | 964059 | 2 | 1277029 | 3 | 1350924 | 4 | 49457 |
| tRNA | 1 | 202353 | 2 | 3088334 | 3 | 303062 | 4 | 20334 |

```bash
cat test2.rna.tsv | tsv-summarize --group-by 2 --mean 5,8,11,14 | sed
'1iRNA\t1\t2\t3\t4' | mlr --itsv --omd cat
```

| RNA | 1 | 2 | 3 | 4 |
| --- | --- | --- | --- | --- |
| mRNA | 80354.2861111 | 109970.054167 | 98965.7236111 | 10804.4680556 |
| rRNA | 62917.7111111 | 75304.1861111 | 62095.9972222 | 1080.06944444 |
| tRNA | 9638.07083333 | 13973.6916667 | 9185.92083333 | 103.230555556 |