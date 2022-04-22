cd /mnt/e/project/srna/output/fastq

for file in `ls SRR*_plantall.fq.gz | perl -p -e 's/_pl.+gz$//'`
do
cat <(seqtk seq -A ${file}_plantall1mis.fq.gz | grep '>' | sed 's/>//') | \
tsv-join --filter-file <(seqtk seq -A ${file}_plantall.fq.gz | grep '>' | sed 's/>//') \
--key-fields 1 -e | perl -n -e 'chomp;$_ =~ /^(.+)\s/;print"$1\n";' > ${file}.list;
seqtk subseq ${file}_plantall1mis.fq.gz ${file}.list | gzip > ${file}_plant1mis.fq.gz;
list=`cat ${file}.list | wc -l`;
snum=`seqtk seq ${file}_plant1mis.fq.gz | grep '^@' | wc -l`;
if [[ $list = $snum ]]; then
    echo "==> extract ${file} reads OK";
    rm ${file}.list;
else
    echo "==> ${file} list has problem, please check";
fi
done