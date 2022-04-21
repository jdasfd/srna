echo "name,count,group";
cd /mnt/e/project/srna/output/bam/plant
for file in `ls *_plantall1mis.sort.bam | perl -p -e 's/_.+bam$//'`
do
plant=`samtools view -F 4 -c -@ 10 ${file}_plantall1mis.sort.bam`;
un=`samtools view -f 4 -c -@ 10 ${file}_plantall1mis.sort.bam`;
echo "${file},${plant},plant";
echo "${file},${un},unknown";
done