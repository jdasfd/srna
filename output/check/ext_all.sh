cd /mnt/e/project/srna/output/test/fastq
for file in `cat ../group1_trna_50.tsv | cut -f 1 | sed '1d'`
do
faops filter -l 0 ${file}_aliall.fq.gz stdout | grep '^>' | \
perl -p -e 's/^>//' > ../name/${file}_aliall.name.tsv
samtools view -@ 10 ../plant/${file}_plantall1mis.sort.bam | \
tsv-join --filter-file ../name/${file}_aliall.name.tsv --key-fields 1 | \
tsv-select -f 1,3,4,6,10 > ../all/${file}_bam.tsv
done