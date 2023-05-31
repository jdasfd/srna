echo "name,count,catgry";
cd /mnt/e/project/srna/output/bam/bacteria
for file in `ls SRR*_mis.sort.bam | perl -p -e 's/_.+bam$//'`
do
plant=`samtools view -c -@ 10 ${file}_aliall.sort.bam`;
mis=`samtools view -c -@ 10 ${file}_mis.sort.bam`;
bacmis=`cat ../bac_tsv/${file}_mis.tsv | wc -l`;
un=`samtools view -c -@ 10 ${file}_unali.sort.bam`;
bacun=`cat ../bac_tsv/${file}_unali.tsv | wc -l`;
let no1=$mis-$bacmis;
let no2=$un-$bacun;
let no=$no1+$no2;
let bac=$bacmis+$bacun;
echo "${file},$plant,plant";
echo "${file},$bac,bacteria";
echo "${file},$no,unknown";
done