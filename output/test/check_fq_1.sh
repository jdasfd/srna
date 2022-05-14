echo -e "name\taliall\tmis\tunali\tsum"
cd /mnt/e/project/srna/output/test/fastq
for file in `ls SRR*_1mis.fq.gz | perl -p -e 's/_.+gz$//'`
do
a=`faops filter -l 0 ${file}_aliall.fq.gz stdout | wc -l`
let aliall=$a/2
m=`faops filter -l 0 ${file}_1mis.fq.gz stdout | wc -l`
let mis=$m/2
u=`faops filter -l 0 ${file}_unali.fq.gz stdout | wc -l`
let unali=$u/2
let sum=$aliall+$mis+$unali
echo -e "${file}\t${aliall}\t${mis}\t${unali}\t${sum}"
done