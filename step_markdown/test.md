# sequence length distribution

cat tier1.tsv | perl -e 'my %len_count;while(<>){chomp;@a=split/\t/,$_;$len=length($a[0]);$len_count{$len}++;}END{for my $len(sort {$len_count{$a} <=> $len_count{$b}} keys %len_count){print"$len\t$len_count{$len}\n";}}' | sed '1ilength\tcount' > seq_length_tier1.tsv

Rscript -e '
library(readr)
library(ggplot2)
args <- commandArgs(T)
count <- read_tsv(args[1])

tplot <- ggplot(count, aes(x = length, y = count)) +
geom_line()
ggsave(tplot, file = "seq_length_tier1.pdf", width = 7, height = 4)
' seq_length_tier1.tsv

# TEST: Extract bam genelist with 1mis hit:

```bash
cat SRR1042171_plant1mis.tsv | \
perl -n -e 'chomp;if($_=~/XM:i:1/){print"$_\n";}else{next;}' | \
tsv-select -f 3,1,2,4,10 > ../SRR1042171_plant1mis.tsv
# 'XM:i:1' in sam/bam files represents 1mis alignment
```



```bash
parallel -j 6 " \
perl /mnt/e/project/srna/script/extract_gene.pl \
-g Atha_genelist.tsv -i SRR1042171_plant1mis.tsv \
-o {}.genelist.tsv \
" ::: $(ls )
```

## Only in tRF region

Remove the intersection of trf3_5 and other_trf.

```bash
mkdir -p /mnt/e/project/srna/output/bam/trf_tsv
cd /mnt/e/project/srna/output/bam/trf

parallel -j 6 " \
samtools view -@ 2 {}.sort.bam > ../trf_tsv/{}.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
parallel -j 6 " \
samtools view -@ 2 {}.sort.bam > ../trf_tsv/{}.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
parallel -j 6 " \
cat {}.trf3_5.tsv | tsv-join --filter-file {}.other_trf.tsv --exclude \
> {}.only_trf.tsv \
" ::: $(ls *.trf3_5.tsv | perl -p -e 's/\.tr.+tsv$//')
```

```bash
mkdir -p /mnt/e/project/srna/output/trf/tsv
cd /mnt/e/project/srna/output/bam/trf_tsv

parallel -j 6 " \
cat {}.only_trf.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-select -f 3,2 | perl ../../../script/filter.pl | \
tsv-join --filter-file ../../count/all/{}.all.tsv --key-fields 1 --append-fields 2 \
> ../../trf/tsv/{}.tsv
" ::: $(ls *.only_trf.tsv | perl -p -e 's/\.on.+tsv$//')
```

```bash
cd /mnt/e/project/srna/output/trf/tsv
bash ../../../script/group_rna_count.sh > ../name_count.only_trf.tsv
```

```bash
cd ..

cat name_count.only_trf.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.only_trf.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.only_trf.tsv -t "only_trf" -y "Bac-reads in tRNA" -o ../figure/only_trf_percent.pdf
```

## 