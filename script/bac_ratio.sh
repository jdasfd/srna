echo "name,ratio,catgry"
cd /mnt/e/project/srna/output/fastq
for file in `ls *_plantall.fq.gz | perl -p -e 's/_.+gz$//'`
do
plant=`faops filter -l 0 ${file}_plantall.fq.gz stdout | grep '^>' | wc -l`;
mis=`faops filter -l 0 ${file}_plantmis.fq.gz stdout | grep '^>' | wc -l`;
un=`faops filter -l 0 ${file}_plantunali.fq.gz stdout | grep '^>' | wc -l`;
bacp=`cat ../bam/bac_tsv/${file}_aliall.tsv | wc -l`;
bacm=`cat ../bam/bac_tsv/${file}_mis.tsv | wc -l`;
bacu=`cat ../bam/bac_tsv/${file}_unali.tsv | wc -l`;
rp=`echo "scale=4;$bacp*100/$plant" | bc | awk '{printf "%.4f", $0}'`;
rm=`echo "scale=4;$bacm*100/$mis" | bc | awk '{printf "%.4f", $0}'`;
ru=`echo "scale=4;$bacu*100/$un" | bc | awk '{printf "%.4f", $0}'`;
echo "${file},${rp},aliall";
echo "${file},${rm},mis";
echo "${file},${ru},unali";
done

# using bc for float calculationz