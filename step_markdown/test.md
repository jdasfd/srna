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

