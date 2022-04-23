- [Bacterial related reads statistical information](#bacterial-related-reads-statistical-information)

    - [Basic percentages of reads aligned to plant or bacteria](#basic-percentages-of-reads-aligned-to-plant-or-bacteria)

    - [Filter by ratio of aligning to plant](#filter-by-ratio-of-aligning-to-plant)

# Bacterial related reads statistical information

In this markdown, I recorded all reads aligned to bacteria and their characristics.

## Ratio of reads aligned to bacteria / all non-plant reads among categories

* All 240 sRNA-seq files

```bash
cd /mnt/e/project/srna/output/bam/bacteria

parallel -j 4 " \
samtools idxstats -@ 3 {}.sort.bam | grep '*' | \
tsv-select -f 4 > ../../count/all/{}.unknown.tsv \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools idxstats -@ 3 {}.sort.bam | grep -v '*' | \
tsv-select -f 1,3 | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 > ../../count/all/{}.bac.tsv \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/count/all

parallel -j 6 " \
cat {}.bac.tsv | tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 > {}.all.tsv \
" ::: $(ls *.bac.tsv | perl -p -e 's/\.bac\.tsv$//')

parallel -j 12 " \
perl ../../../script/count.pl -a {}.all.tsv \
-r {}.unknown.tsv -o {}.tsv \
" ::: $(ls *.all.tsv | perl -p -e 's/\.all\.tsv$//')

rm *.all.tsv *.unknown.tsv
```

```bash
cd /mnt/e/project/srna/output/count/all

rm ../bacreads_group.tsv
# because next step will use >> so clear first

for file in `ls *.bac.tsv | perl -p -e 's/\.bac\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.tsv | awk -v name=$name -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
>> ../bacreads_group.tsv
done

cd ..
sed -i '1i\name\tgroup\tratio\tcatgry' bacreads_group.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
 -n 5 -f bacreads_group.tsv -t RNA_in_group -y "Bac-reads" -o ../figure/bacreads_group.pdf
# add the scale_y_continuous limits
```

#### After filtering

```bash
cat bacreads_group.tsv | \
tsv-join -H --filter-file plant_50.tsv --key-fields name \
> bacreads_group_50.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-n 5 -f bacreads_group_50.tsv -t RNA_in_group -y "Bac-reads" -o ../figure/bacreads_group_50.pdf
```

```bash
cat bacreads_group_50.tsv | \
tsv-summarize -H --group-by group,catgry --mean ratio | \
mlr --itsv --omd cat
```

| group | catgry | ratio_mean |
| --- | --- | --- |
| 3 | 1mis | 0.634710151378 |
| 1 | 1mis | 1.02873775248 |
| 4 | 1mis | 0.00620741048405 |
| 2 | 1mis | 1.50359081252 |
| 1 | aliall | 1.48569685888 |
| 4 | aliall | 0.208781526894 |
| 3 | aliall | 1.16217552165 |
| 2 | aliall | 2.14279758871 |
| 1 | unali | 3.40977514859 |
| 4 | unali | 0.0586593955001 |
| 2 | unali | 4.55701684177 |
| 3 | unali | 7.21581767507 |

```bash
cat bacreads_group_50.tsv | \
tsv-summarize -H --group-by group,catgry --median ratio | \
mlr --itsv --omd cat
```

| group | catgry | ratio_median |
| --- | --- | --- |
| 3 | 1mis | 0.212392565465 |
| 1 | 1mis | 0.425748152575 |
| 4 | 1mis | 0.0033455682172 |
| 2 | 1mis | 0.208726262452 |
| 1 | aliall | 0.937687844361 |
| 4 | aliall | 0.112857058589 |
| 3 | aliall | 0.750917485459 |
| 2 | aliall | 1.45728697456 |
| 1 | unali | 2.02597994901 |
| 4 | unali | 0.0434075553298 |
| 2 | unali | 3.1837129376 |
| 3 | unali | 1.20582662298 |

