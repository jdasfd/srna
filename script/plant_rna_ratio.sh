echo "name,count,catgry";
cd /mnt/e/project/srna/output/bam/plant

for file in `ls SRR*_plant.sort.bam | perl -p -e 's/_plant.+bam$//'`
do
array=(trna rrna mrna mirna snrna snorna lncrna ncrna)
for rna in ${array[@]}
do
count=`samtools view -@ 10 -L ../../../annotation/plant/Atha/Atha_${rna}.bed \
${file}_plant.sort.bam | perl -n -e 'chomp;if($_=~/^@/){print "$_\n";}else{
@array = split/\t/, $_;if($array[5] =~ /^((?!X).)*$/ && $array[1] != 4){print "$_\n";}else{next;}}
' | wc -l`;
echo "${file},${count},${rna}";
done
done