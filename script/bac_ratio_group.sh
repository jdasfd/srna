cd /mnt/e/project/srna/output/bam/bac_tsv
for file in `ls *.tsv | perl -p -e 's/\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
all=`faops filter -l 0 ../../fastq/${name}_plant${catgry}.fq.gz stdout | grep '^>' | wc -l`;
cat ${file}.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/all/${file}.all.tsv;
echo -e "$name\tall\t$all\t$catgry" >> ../../count/all/${file}.all.tsv;
done