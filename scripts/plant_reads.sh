echo "name,count,catgry"
cd /mnt/e/project/srna/output/fastq
for file in `ls *_plantall.fq.gz | perl -p -e 's/_.+gz$//'`
do
plant=`faops filter -l 0 ${file}_plantall.fq.gz stdout | grep '^>' | wc -l`;
mis=`faops filter -l 0 ${file}_plantmis.fq.gz stdout | grep '^>' | wc -l`;
un=`faops filter -l 0 ${file}_plantunali.fq.gz stdout | grep '^>' | wc -l`;
let no=${mis}+${un};
echo "${file},${plant},plant";
echo "${file},${no},unknown";
# echo "${file},${mis},target-plant";
# echo "${file},${un},unknown";
done

# variable 1mis was not available