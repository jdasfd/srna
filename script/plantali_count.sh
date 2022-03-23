echo "name,count,group";
cd /mnt/e/project/srna/output/bam/plantrna
for file in `ls SRR*.bam | perl -p -e 's/\.bam$//'`
do
name=${file%%*_};
catgry=${file#*.};
count=`samtools view ${file}.bam | wc -l`;
echo "${name},${count},${catgry}";
done