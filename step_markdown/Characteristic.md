# sRNA reads characteristics

The sRNA sequences has some specific characteristics. We wanted to find if there were patterns of the sRNA distribution.

## Ratio of sRNA reads / bac-reads from different RNA regions of 4 groups

- Using bed of rna to extract mapping reads from different RNA regions

```bash
mkdir -p /mnt/e/project/srna/output/bam/rna
cd /mnt/e/project/srna/output/bam/bacteria

parallel -j 6 " \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/trna.bed \
{}_aliall.sort.bam > ../rna/{}_aliall.trna.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/trna.bed \
{}_mis.sort.bam > ../rna/{}_mis.trna.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/trna.bed \
{}_unali.sort.bam > ../rna/{}_unali.trna.sort.bam; \
" ::: $(cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d')

parallel -j 6 " \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/rrna.bed \
{}_aliall.sort.bam > ../rna/{}_aliall.rrna.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/rrna.bed \
{}_mis.sort.bam > ../rna/{}_mis.rrna.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/rrna.bed \
{}_unali.sort.bam > ../rna/{}_unali.rrna.sort.bam; \
" ::: $(cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d')

parallel -j 6 " \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/mrna.bed \
{}_aliall.sort.bam > ../rna/{}_aliall.mrna.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/mrna.bed \
{}_mis.sort.bam > ../rna/{}_mis.mrna.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/mrna.bed \
{}_unali.sort.bam > ../rna/{}_unali.mrna.sort.bam; \
" ::: $(cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d')
```

- Use CIGAR to filter all bac-reads from different region

Three different RNA regions `.bam` files from the previous step.

```bash
mkdir -p /mnt/e/project/srna/output/bam/rna_tsv
cd /mnt/e/project/srna/output/bam/rna

for file in `ls SRR*.bam | perl -p -e 's/\.sort\.bam$//'`
do
samtools view -@ 10 ${file}.sort.bam | \
tsv-filter --not-iregex 6:X | \
tsv-select -f 1,2,3,4,6,10,11 > ../rna_tsv/${file}.tsv;
done
```

```bash
cd /mnt/e/project/srna/output/count/
mkdir all rna
cd /mnt/e/project/srna/output/bam/bac_tsv

for file in `ls *.tsv | perl -p -e 's/\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/all/${file}.all.tsv;
done

cd /mnt/e/project/srna/output/bam/rna_tsv

for file in `ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.trna.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/rna/${file}.trna.tsv;
done

for file in `ls *.rrna.tsv | perl -p -e 's/\.rrna\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.rrna.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/rna/${file}.rrna.tsv;
done

for file in `ls *.mrna.tsv | perl -p -e 's/\.mrna\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.mrna.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/rna/${file}.mrna.tsv;
done
```

- RNA reads ratio in all bac-reads (group)

```bash
cd /mnt/e/project/srna/output/count/rna

for file in `ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//'`
do
for rna in "trna" "rrna" "mrna"
do
cat ${file}.${rna}.tsv | tsv-join --filter-file ../all/${file}.all.tsv \
--key-fields 2 --append-fields 3 >> ../bac_reads_group.${rna}.tsv;
done
done
```

```bash
cd ..

for rna in "trna" "rrna" "mrna"
do
cat bac_reads_group.${rna}.tsv | \
tsv-filter --ne 5:0 | \
perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[4];
printf"%s\t%s\t%.3f\t%s\n",$a[0],$a[1],$b,$a[3];
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > bac_ratio_group.${rna}.tsv
done
```

- Plot

```bash
Rscript /mnt/e/project/srna/scripts/rna_percent.r \
-f bac_ratio_group.trna.tsv -t tRNA_region -y "Bac-reads in tRNA" -o ../figure/trna_reads.pdf

Rscript /mnt/e/project/srna/scripts/rna_percent.r \
-f bac_ratio_group.rrna.tsv -t rRNA_region -y "Bac-reads in rRNA" -o ../figure/rrna_reads.pdf

Rscript /mnt/e/project/srna/scripts/rna_percent.r \
-f bac_ratio_group.mrna.tsv -t mRNA_region -y "Bac-reads in mRNA" -o ../figure/mrna_reads.pdf
```

- Summary

| group | catgry | ratio_mean     | ratio_median |
| ----- | ------ | -------------- | ------------ |
| 1     | aliall | 12.9401138614  | 8.463        |
| 2     | aliall | 10.7372512315  | 6.531        |
| 3     | aliall | 7.22172772277  | 4.5755       |
| 4     | aliall | 0.161110169492 | 0.072        |
| 1     | mis    | 11.1145114943  | 8.39         |
| 2     | mis    | 9.25039408867  | 6.641        |
| 3     | mis    | 15.660540107   | 13.418       |
| 4     | mis    | 1.85691011236  | 1.307        |
| 1     | unali  | 19.3068579235  | 19.266       |
| 2     | unali  | 15.6841386139  | 15.868       |
| 3     | unali  | 14.3077142857  | 16.3145      |
| 4     | unali  | 7.67607462687  | 5.882        |

## Sequence among all files

So the main reason I do this step is to select those frequently occurred among all sequence files, that is, the most likely sRNA appeared among *A. tha* sRNA-seq files.

- Extract all reads and count their frequencies occured among files

```bash
mkdir -p /mnt/e/project/srna/output/tier/all
cd /mnt/e/project/srna/output/bam/bac_tsv

parallel -j 4 " \
cat {}_aliall.tsv | tsv-summarize --group-by 3,10 --count \
> ../../tier/all/{}_aliall.tsv; \
cat {}_mis.tsv | tsv-summarize --group-by 3,10 --count \
> ../../tier/all/{}_mis.tsv; \
cat {}_unali.tsv | tsv-summarize --group-by 3,10 --count \
> ../../tier/all/{}_unali.tsv; \
" ::: $(cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d')

cd ../../tier/all

if [[ -e ../among/reads_*.all.tsv ]]
then
rm ../among/reads_*.all.tsv;
else
echo "OK for next step"
fi

for file in `ls *_aliall.tsv | perl -p -e 's/_.+tsv$//'`
do
cat ${file}_aliall.tsv | tsv-select -f 2 | \
tsv-summarize --group-by 1 --count >> ../among/reads_aliall.all.tsv;
cat ${file}_mis.tsv | tsv-select -f 2 | \
tsv-summarize --group-by 1 --count >> ../among/reads_mis.all.tsv;
cat ${file}_unali.tsv | tsv-select -f 2 | \
tsv-summarize --group-by 1 --count >> ../among/reads_unali.all.tsv;
done
```

- Count and sort

```bash
cd /mnt/e/project/srna/output/tier/among

for catgry in "aliall" "mis"
do
cat reads_${catgry}.all.tsv | \
tsv-summarize --group-by 1 --count | \
tsv-filter --ge 2:5 | \
sort -nk 2 -r | sed '1iseq\tcount' \
> all_file_${catgry}.all.tsv;
done

cat reads_unali.all.tsv | \
tsv-summarize --group-by 1 --count | \
tsv-filter --ge 2:11 | \
sort -nk 2 -r | sed '1iseq\tcount' \
> all_file_unali.all.tsv
```

- Plot

```bash
bash ../../../scripts/mclust.sh all_file_mis.all.tsv
bash ../../../scripts/mclust.sh all_file_aliall.all.tsv
bash ../../../scripts/mclust.sh all_file_unali.all.tsv

# choose rignt n according to the previous step figure

bash ../../../scripts/mclust_dens.sh all_file_mis.all.tsv 3
bash ../../../scripts/mclust_dens.sh all_file_aliall.all.tsv 3
bash ../../../scripts/mclust_dens.sh all_file_unali.all.tsv 3
```

- Summary

According to the results generated by `mclust`, the frequencies of reads occured among all 203 files in different catagories were counted and generated as plots.

| catgry | tier1 | tier2 | tier3  |
| ------ | ----- | ----- | ------ |
| aliall | 5-8   | 9-28  | 29-193 |
| mis    | 5-6   | 7-16  | 17-154 |
| unali  | 11-13 | 14-24 | 25-149 |

```bash
cat all_file_aliall.all.tsv | tsv-filter -H --ge count:5 --le count:8 | sed '1d' > aliall.tier3.all.tsv
cat all_file_aliall.all.tsv | tsv-filter -H --ge count:9 --le count:28 | sed '1d' > aliall.tier2.all.tsv
cat all_file_aliall.all.tsv | tsv-filter -H --ge count:29 --le count:193 | sed '1d' > aliall.tier1.all.tsv
cat all_file_mis.all.tsv | tsv-filter -H --ge count:5 --le count:6 | sed '1d' > mis.3.all.tsv
cat all_file_mis.all.tsv | tsv-filter -H --ge count:7 --le count:16 | sed '1d' > mis.tier2.all.tsv
cat all_file_mis.all.tsv | tsv-filter -H --ge count:17 --le count:154 | sed '1d' > mis.tier1.all.tsv
cat all_file_unali.all.tsv | tsv-filter -H --ge count:11 --le count:13 | sed '1d' > unali.tier3.all.tsv
cat all_file_unali.all.tsv | tsv-filter -H --ge count:14 --le count:24 | sed '1d' > unali.tier2.all.tsv
cat all_file_unali.all.tsv | tsv-filter -H --ge count:25 --le count:149 | sed '1d' > unali.tier1.all.tsv
```

### Extract reads in tiers and ratio count

- Extract reads

```bash
cd /mnt/e/project/srna/output/tier/
mkdir tier1 tier2 tier3
cd /mnt/e/project/srna/output/bam/bac_tsv

for file in `cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d'`
do
for catgry in "aliall" "mis" "unali"
do
cat ${file}_${catgry}.tsv | tsv-select -f 10,3 | \
tsv-join --filter-file ../../tier/among/${catgry}.tier1.all.tsv --key-fields 1 \
> ../../tier/tier1/${file}_${catgry}.tsv;
cat ${file}_${catgry}.tsv | tsv-select -f 10,3 | \
tsv-join --filter-file ../../tier/among/${catgry}.tier2.all.tsv --key-fields 1 \
> ../../tier/tier2/${file}_${catgry}.tsv;
cat ${file}_${catgry}.tsv | tsv-select -f 10,3 | \
tsv-join --filter-file ../../tier/among/${catgry}.tier3.all.tsv --key-fields 1 \
> ../../tier/tier3/${file}_${catgry}.tsv;
done
done
```

- Ratio count

```bash
mkdir -p /mnt/e/project/srna/output/count/tier
cd /mnt/e/project/srna/output/tier

for dir in "tier1" "tier2" "tier3"
do
cd ${dir};
for file in `ls *.tsv | perl -p -e 's/\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.tsv | \
tsv-summarize --group-by 2 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/tier/${file}.${dir}.tsv;
done
cd ..;
done
```

```bash
cd /mnt/e/project/srna/output/count/tier

for file in `ls *.tier1.tsv | perl -p -e 's/\.tier1\.tsv$//'`
do
for tier in "tier1" "tier2" "tier3"
do
cat ${file}.${tier}.tsv | tsv-join --filter-file ../all/${file}.all.tsv \
--key-fields 2 --append-fields 3 >> ../bac_reads_group.${tier}.tsv;
done
done
```

```bash
cd ..

for tier in "tier1" "tier2" "tier3"
do
cat bac_reads_group.${tier}.tsv | \
tsv-filter --ne 5:0 | \
perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[4];
printf"%s\t%s\t%.3f\t%s\n",$a[0],$a[1],$b,$a[3];
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > bac_ratio_group.${tier}.tsv
done
```

- Plot

```bash
Rscript /mnt/e/project/srna/scripts/rna_percent.r \
-f bac_ratio_group.tier1.tsv -t tier1 -y "Bac-reads in tRNA" -o ../figure/tier1_reads.pdf

Rscript /mnt/e/project/srna/scripts/rna_percent.r \
-f bac_ratio_group.tier2.tsv -t tier2 -y "Bac-reads in rRNA" -o ../figure/tier2_reads.pdf

Rscript /mnt/e/project/srna/scripts/rna_percent.r \
-f bac_ratio_group.tier3.tsv -t tier3 -y "Bac-reads in mRNA" -o ../figure/tier3_reads.pdf
```

- Summary

```bash
for tier in "tier1" "tier2" "tier3"
do
cat bac_ratio_group.${tier}.tsv | \
tsv-summarize -H --group-by group,catgry --mean ratio --median ratio | \
mlr --itsv --omd cat
done
```

tier1:
| group | catgry | ratio_mean    | ratio_median |
| ----- | ------ | ------------- | ------------ |
| 1     | aliall | 84.4243497537 | 91.416       |
| 2     | aliall | 82.778044335  | 91.692       |
| 3     | aliall | 81.3253366337 | 89.0585      |
| 4     | aliall | 87.3819310345 | 95.408       |
| 1     | mis    | 54.7773399015 | 57.961       |
| 2     | mis    | 52.8013891626 | 59.662       |
| 3     | mis    | 46.4350246305 | 47.479       |
| 4     | mis    | 36.43272      | 32.52        |
| 1     | unali  | 29.832729064  | 29.483       |
| 2     | unali  | 36.5572660099 | 31.955       |
| 3     | unali  | 17.9435320197 | 17.045       |
| 4     | unali  | 6.56031666667 | 3.852        |

tier2:
| group | catgry | ratio_mean    | ratio_median |
| ----- | ------ | ------------- | ------------ |
| 1     | aliall | 6.95022167488 | 5.239        |
| 2     | aliall | 7.15395073892 | 4.413        |
| 3     | aliall | 7.67451485149 | 5.8685       |
| 4     | aliall | 3.95294059406 | 1.931        |
| 1     | mis    | 14.8673399015 | 13.889       |
| 2     | mis    | 15.150364532  | 15.302       |
| 3     | mis    | 15.2956157635 | 14.085       |
| 4     | mis    | 16.8176060606 | 14.8055      |
| 1     | unali  | 19.5638275862 | 15.517       |
| 2     | unali  | 13.7571231527 | 13.56        |
| 3     | unali  | 15.9602413793 | 14.13        |
| 4     | unali  | 17.0745890411 | 12.903       |

tier3:
| group | catgry | ratio_mean    | ratio_median |
| ----- | ------ | ------------- | ------------ |
| 1     | aliall | 1.98148275862 | 0.87         |
| 2     | aliall | 2.39063054187 | 0.844        |
| 3     | aliall | 2.41216336634 | 1.2985       |
| 4     | aliall | 1.57330150754 | 0.408        |
| 1     | mis    | 54.7773399015 | 57.961       |
| 2     | mis    | 52.8013891626 | 59.662       |
| 3     | mis    | 46.4350246305 | 47.479       |
| 4     | mis    | 36.43272      | 32.52        |
| 1     | unali  | 5.39365174129 | 5.473        |
| 2     | unali  | 4.92791089109 | 4.7175       |
| 3     | unali  | 6.06823383085 | 6            |
| 4     | unali  | 10.1163176471 | 5.392        |

### Length distribution

Count different reads length in each tier

- Length distribution

```bash
cd /mnt/e/project/srna/output/tier/among

parallel -j 6 " \
cat {}.all.tsv | tsv-select -f 1 | \
perl ../../../script/seq_len.pl | sort -nk 1 | \
sed '1ilen\tnum' > ../../length/{}.len.tsv
" ::: $(ls *.tier*.all.tsv | perl -p -e 's/\.all\.tsv$//')
```

- Plot

```bash
cd /mnt/e/project/srna/output/length

for file in `ls *.len.tsv | perl -p -e 's/\.len\.tsv$//'`
do
Rscript -e '
library(readr)
library(ggplot2)
args <- commandArgs(T)
input <- paste0(args[1], ".len.tsv")
out <- paste0("../figure/length_", args[1], ".pdf")

length <- read_tsv(input, show_col_types = FALSE)
nmax <- max(length$len)
nmin <- min(length$len)
p <- ggplot(length, aes(x = len, y = num)) +
geom_line() +
scale_x_continuous(limits = c(nmin, nmax), breaks = seq(nmin, nmax, by = 1))

ggsave(p, file = out, width = 9, height = 4)
' ${file}
done
```

## tRF region

tRF3/5 regions and other_tRNA regions were extracted from the tRNA.bed file according to tRF characteristics.

- Using bed of trna to extract mapping reads from different RNA regions

```bash
faops count ../../genome/bacteria/bacteria.fna | tsv-select -f 1,2 | sed '1d' > bac.genome

bedtools flank -i bac_trna.gff -b 50 -g bac.genome > bac_trf1.gff
cat bac_trf1.gff | convert2bed --input=gff --output=bed > trf1.bed

cat trna.bed | perl -e 'while(<>){
    chomp;
    @a = split/\t/,$_;
    $end1 = $a[1] + 25;
    $end2 = $a[2] - 25;
    print"$a[0]\t$a[1]\t$end1\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\n";
    print"$a[0]\t$end2\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\n";
}' > trf3_5.bed

bedtools subtract -a trna.bed -b trf3_5.bed > other_trf.bed
```

```bash
mkdir -p /mnt/e/project/srna/output/bam/trf
cd /mnt/e/project/srna/output/bam/bacteria

parallel -j 6 " \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/trf3_5.bed \
{}_aliall.sort.bam > ../trf/{}_aliall.trf3_5.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/trf3_5.bed \
{}_mis.sort.bam > ../trf/{}_mis.trf3_5.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/trf3_5.bed \
{}_unali.sort.bam > ../trf/{}_unali.trf3_5.sort.bam; \
" ::: $(cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d')

parallel -j 6 " \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/trf1.bed \
{}_aliall.sort.bam > ../trf/{}_aliall.trf1.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/trf1.bed \
{}_mis.sort.bam > ../trf/{}_mis.trf1.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/trf1.bed \
{}_unali.sort.bam > ../trf/{}_unali.trf1.sort.bam; \
" ::: $(cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d')

parallel -j 6 " \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/other_trf.bed \
{}_aliall.sort.bam > ../trf/{}_aliall.other.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/other_trf.bed \
{}_mis.sort.bam > ../trf/{}_mis.other.sort.bam; \
samtools view -@ 2 -F 4 -bh -L ../../../annotation/bacteria/other_trf.bed \
{}_unali.sort.bam > ../trf/{}_unali.other.sort.bam; \
" ::: $(cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d')
```

- Use CIGAR to filter all bac-reads from different tRF region

Three different RNA regions `.bam` files from the previous step.

```bash
mkdir -p /mnt/e/project/srna/output/bam/trf_tsv
cd /mnt/e/project/srna/output/bam/trf

for file in `ls SRR*.bam | perl -p -e 's/\.sort\.bam$//'`
do
samtools view -@ 10 ${file}.sort.bam | \
tsv-filter --not-iregex 6:X > ../trf_tsv/${file}.tsv;
done
```

```bash
cd /mnt/e/project/srna/output/count/
mkdir all rna
cd /mnt/e/project/srna/output/bam/bac_tsv

for file in `ls *.tsv | perl -p -e 's/\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/all/${file}.all.tsv;
done

cd /mnt/e/project/srna/output/bam/rna_tsv

for file in `ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.trna.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/rna/${file}.trna.tsv;
done

for file in `ls *.rrna.tsv | perl -p -e 's/\.rrna\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.rrna.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/rna/${file}.rrna.tsv;
done

for file in `ls *.mrna.tsv | perl -p -e 's/\.mrna\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.mrna.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/rna/${file}.mrna.tsv;
done
```

## Extract tRF sequence from different bacteria

### Get the *A. thaliana* transcripts

```bash
mkdir -p /mnt/e/project/srna/genome/plant_CDS/Atha
cd /mnt/e/project/srna/genome/plant_CDS/Atha
gffread ../../../annotation/plant/Atha/Atha.gff -g ../../plant/Atha/Atha.fna \
-w Atha_transcript.fa -x Atha_CDS.fa
```

### Extract sequence from the tRF3/5 region and high frequency sequence

```bash
mkdir -p /mnt/e/project/srna/output/seq/trf
cd /mnt/e/project/srna/output/seq/trf
mkdir 1 2 3
cd /mnt/e/project/srna/output/bam/trf

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:1 \
> ../../seq/trf/1/{}.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:2 \
> ../../seq/trf/2/{}.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:3 \
> ../../seq/trf/3/{}.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
mkdir -p /mnt/e/project/srna/output/seq/other
cd /mnt/e/project/srna/output/seq/other
mkdir 1 2 3
cd /mnt/e/project/srna/output/bam/trf

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:1 \
> ../../seq/other/1/{}.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:2 \
> ../../seq/other/2/{}.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:3 \
> ../../seq/other/3/{}.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/seq/trf/1

cat *_1mis.trf3_5.tsv | tsv-select -f 2,4,5 | \
perl -n -e 'chomp;@a=split/\t/,$_;$l=length($a[2]);print"$a[0]\t$a[1]\t$l\t$a[2]\n";' | \
tsv-summarize --group-by 1,2,3,4 --count | \
perl -n -e 'chomp;@a = split/\t/,$_;print">$a[0]_$a[1]_$a[2]\n";print"$a[3]\n"' \
> ../1mis.1.fasta
```

---

## Reads cover depth and position in different RNA region

### Reads depth and position distribution in different RNA

Divide tRNA regions to 10 separate domains. Every domains reads were counted as relative depth of RNA regions.

```bash
cd /mnt/e/project/srna/output/depth/trna

cat *_1mis.trna.txt >> ../1mis.trna.txt
cat *_unali.trna.txt >> ../unali.trna.txt
cat *_aliall.trna.txt >> ../aliall.trna.txt

cd ..

for file in `ls *.txt | perl -p -e 's/\.trna\.txt//'`
do
cat ${file}.trna.txt | tsv-summarize --group-by 1,2 --sum 3 \
> ${file}.trna.tsv
done

parallel -j 3 " \
perl ../../script/depth.pl -b ../../annotation/bacteria/bac_trna.bed \
-t {}.trna.tsv -o {}.trna.depth.tsv \
" ::: $(ls *.trna.tsv | perl -p -e 's/\.trna\.tsv//')

rm *.txt
```

```bash
parallel -j 3 " \
cat {}.depth.tsv | \
tsv-join --filter-file ../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 1,2,4 --sum 3 | tsv-select -f 3,2,4 | \
tsv-join --filter-file ../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 4,2 --sum 3 | sort -nk 1,1 -nk 2,2 | \
sed '1i\group\tpos\tdepth' > {}.sum.tsv \
" ::: $(ls *.depth.tsv | perl -p -e 's/\.dep.+\.tsv$//')
```

```bash
parallel -j 3 " \
Rscript ../../script/line_cover.r -d {}.sum.tsv \
-o {}.pdf -t {} \
" ::: $(ls *.sum.tsv | perl -p -e 's/\.sum\.tsv$//')
```

### Bacteria with tRNA coverage

```bash
cd /mnt/e/project/srna/output/count/trna

for file in `ls *_1mis.tsv | perl -p -e 's/\.tsv//'`
do
cat ${file}.tsv | \
tsv-filter --ne 3:0 | \
tsv-select -f 1,3 \
>> ../1mis_tRNA.name.tsv
done

cd ..

cat 1mis_tRNA.name.tsv | tsv-summarize --group-by 1 --sum 2 | \
sort -r -nk2 | tsv-join --filter-file ../../name.tsv --key-fields 1 --append-fields 2 \
> bac_count.trna.tsv

cat bac_count.trna.tsv | head -n 50 | sed '1i\name\tcount\tgroup'> top50.trna.tsv

cp top50.trna.tsv /mnt/c/Users/59717/Documents/
```

```R
library(ggplot2)
library(readr)
library(dplyr)
library(forcats)
library(ggforce)

t50t <- read_tsv("top50.trna.tsv")

t50t$group <- as.character(t50t$group)

tplot <- t50t %>%
mutate(name = fct_reorder(name, desc(count))) %>%
ggplot(aes(x = name, y = count, fill = group)) +
geom_bar(stat="identity") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4)) +
geom_col() +
facet_zoom(ylim = c(0, 30000))
# as.character could change the numeric variables to character variables
# mutate in forcats could sort for the bar plot, desc could sort reversely
# facet_zoom could seperate the bar plot into 2 different y axis resolution
# facet_zoom function was from ggforce
```

### tRNA reads statistical information

Reads extracted from tRNA regions were considered as our target sRNAs. First, we used a script to check file

```bash
mkdir -p /mnt/e/project/srna/output/kegg
cd /mnt/e/project/srna/output/target/trna

for file in `ls *.fasta | perl -p -e 's/\.fasta//'`
do
cat ${file}.fasta | perl ../../../scripts/seq.pl | sort -nk 2 -r > ${file}_seq.tsv;
done

touch all_trna.tsv
cat *_seq.tsv | cut -f 1 >> all_trna.tsv

cat all_trna.tsv | perl -n -e 'while(<>){chomp; my $seq = $_; $hash{$seq}++}
END{
for my $seq (sort keys %hash){print "$seq\t","$hash{$seq}\n";}
}
' | tsv-filter --gt 2:20 | sort -nk 2 -r > ../../kegg/target_trna.tsv


cat *_rrna_seq.tsv | cut -f 1 >> all_rrna.tsv


cat all_rrna.tsv | perl -n -e 'while(<>){chomp; my $seq = $_; $hash{$seq}++}
END{
for my $seq (sort keys %hash){print "$seq\t","$hash{$seq}\n";}
}
' | tsv-filter --gt 2:20 | sort -nk 2 -r > ../kegg/target_rrna.tsv

rm all_trna.tsv all_rrna.tsv
```

```bash
cat target_trna.tsv | cut -f 1 | perl -n -e 'while(<>){chomp; my $seq = $_; $hash{$seq} = length($seq)}
END{
for my $seq (sort keys %hash){print "$seq\t","$hash{$seq}\n";}
}
' | cut -f 2 | perl -n -e 'while(<>){chomp; my $num = $_; $count{$num}++;}
END{
for my $num (sort keys %count){print "$num\t","$count{$num}\n";}
}
' | sed '1i length\tcount'> length.tsv
```

```R
library(ggplot2)
library(readr)

length <- read_tsv("length.tsv")
l <- ggplot(data = length, aes(x = length, y = count))+
geom_point()+
geom_line()+
scale_x_continuous(breaks=seq(23,35,1))

plot(l)
```

```bash
for file in `ls *.fasta | perl -p -e 's/\.fasta$//'`
do
cat ${file}.fasta | perl -n -e '
while(<>){chomp;
if($_ =~ /^>(.+?)$/)
{print"$1\n";}
else{next;}
}' > ../../sequence/${file}.tsv
done
```

```bash
cd /mnt/e/project/srna/output/bam/rna

parallel -j 10 " \
samtools view {}.trna.sort.bam | \
tsv-select -f 1,3 > ../../sequence/{}.bam.tsv \
" ::: $(ls *_1mis.trna.sort.bam | perl -p -e 's/\.trna.+?bam$//')

cd /mnt/e/project/srna/output/sequence

parallel -j 10 " \
tsv-join {}_1mis.tsv --filter-file {}_1mis.bam.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 2 --count | \
tsv-join --filter-file ../../rawname.tsv --key-fields 1 --append-fields 2 \
> {}.tsv \
" ::: $(ls *.tsv | perl -p -e 's/\_.+?tsv$//' | uniq)

rm *_1mis*.tsv
```

---

```R
bar <- ggplot(data = name, mapping = aes(x = 'category', y = count, fill = cate)) +
geom_bar(stat = 'identity', position = 'stack') +
coord_polar(theta = 'y') +
labs(x = '', y = '', title = '') +

```

```bash
cat plot.count | perl -n -e 'chomp;
my $p,$chr;
@a=split(/\t/,$_);
if($a[0]==$chr){
$p++;
print"$a[0]\t$p\t$a[2]\n";
}
else{
$chr=$a[0];
$p=1;
print"$a[0]\t$p\t$a[2]\n";
}
' > chr.count
```

```R
library(readr)
library(ggplot2)

trnad <- read_tsv("1mis.trna.depth")
names(trnad) = c("bac", "pos", "depth")

depth <- ggplot(data = trnad, aes(x = pos, y = depth, group = bac))+
geom_point()+
geom_line()+
theme(legend.position = 'none')

plot(depth)
```

```bash
Rscript -e '
library(readr)
library(ggplot2)
library(dplyr)
library(forcats)

count <- read_tsv("all_seq.count.tsv")
count$num <- as.character(count$num)

tplot <- count %>%
mutate(seq = fct_reorder(seq, desc(count))) %>%
ggplot(aes(x = seq, y = num)) +
geom_bar(stat="identity") +
theme(axis.text.x = NULL)
ggsave(tplot, file = "freq.pdf", width = 7, height = 4)
'
```

## Reads TPM among different groups (waiting for update)

The goal is to filter out those reads among one file.

```bash
cd /mnt/e/project/srna/annotation/bacteria

cat trna.bed | \
perl -e 'my $sum = 0;
while(<>){chomp;@a = split/\t/,$_;
$length = $a[2]-$a[1];
$sum = $sum+$length;
}
print "$sum";
' | awk '{print $1"\ttrna"}'> ../../output/cover/length.tsv

cat rrna.bed | \
perl -e 'my $sum = 0;
while(<>){chomp;@a = split/\t/,$_;
$length = $a[2]-$a[1];
$sum = $sum+$length;
}
print "$sum";
' | awk '{print $1"\trrna"}'>> ../../output/cover/length.tsv

cat mrna.bed | \
perl -e 'my $sum = 0;
while(<>){chomp;@a = split/\t/,$_;
$length = $a[2]-$a[1];
$sum = $sum+$length;
}
print "$sum";
' | awk '{print $1"\tmrna"}'>> ../../output/cover/length.tsv
```

```bash
cd /mnt/e/project/srna/output/cover

perl ../../scripts/tpm.pl | sed '1ifile\tgroup\tcatgry\tttpm\trtpm\tmtpm' > tpm.tsv

tsv-select -H -f file,group,catgry,ttpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > ttpm.tsv
tsv-select -H -f file,group,catgry,rtpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > rtpm.tsv
tsv-select -H -f file,group,catgry,mtpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > mtpm.tsv

Rscript /mnt/e/project/srna/scripts/rna_plot.r -f ttpm.tsv -o ../figure/ttpm.pdf -t trna -y TPM
Rscript /mnt/e/project/srna/scripts/rna_plot.r -f rtpm.tsv -o ../figure/rtpm.pdf -t rrna -y TPM
Rscript /mnt/e/project/srna/scripts/rna_plot.r -f mtpm.tsv -o ../figure/mtpm.pdf -t mrna -y TPM
```
