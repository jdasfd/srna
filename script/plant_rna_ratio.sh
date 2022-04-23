echo "name,count,group";
cd /mnt/e/project/srna/output/bam/plant
for file in `ls SRR*_plantall.sort.bam | perl -p -e 's/_plant.+bam$//'`
do
array=(trna rrna mrna mirna snrna snorna lncrna ncrna)
for rna in ${array[@]}
do
count=`samtools view -@ 10 -L ../../../annotation/plant/Atha/Atha_${rna}.bed -F 4 -c ${file}_plantall.sort.bam`;
echo "${file},${count},${rna}";
done
done