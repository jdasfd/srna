cd /mnt/e/project/srna/output/fastq
for file in `ls *_plantall.fq.gz | perl -p -e 's/_.+gz$//'`
do
plant=`faops filter -l 0 ${file}_plantall.fq.gz stdout | grep '^>' | wc -l`;
1mis=`faops filter -l 0 ${file}_plant1mis.fq.gz stdout | grep '^>' | wc -l`;
un=`faops filter -l 0 ${file}_plantunali.fq.gz stdout | grep '^>' | wc -l`;
let no=${1mis}+${un};
echo "${file},${plant},plant";
echo "${file},${no},unknown";
done