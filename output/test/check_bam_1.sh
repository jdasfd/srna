echo -e "name\tbam"
cd /mnt/e/project/srna/output/test/plant
for file in `ls *_plantall.sort.bam | perl -p -e 's/_.+bam$//'`
do
bam=`samtools view -c -@ 10 ${file}_plantall.sort.bam`
echo -e "${file}\t${bam}"
done