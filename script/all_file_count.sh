echo "name,count,group";
cd /mnt/e/project/trf/output/bam/plant
for file in `ls SRR*.sort.bam | perl -p -e 's/_.+bam$//'`
do
cd ..;
alip=`samtools view -c -F 4 -@ 10 ./plant/${file}_plantall.sort.bam`;
alib1=`samtools view -c -F 4 -@ 10 ./bacteria/${file}_unali.sort.bam`;
un1=`samtools view -c -f 4 -@ 10 ./bacteria/${file}_unali.sort.bam`;
alib2=`samtools view -c -F 4 -@ 10 ./bacteria/${file}_1mis.sort.bam`;
un2=`samtools view -c -f 4 -@ 10 ./bacteria/${file}_1mis.sort.bam`;
alib=${alib1}+${alib2};
un=${un1}+${un2};
echo "${file},${alip},plant";
echo "${file},${alib},bacteria";
echo "${file},${un},unknown";
cd ./plant;
done