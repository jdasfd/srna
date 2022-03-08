echo "name,num,group";
for file in `ls SRR*.sort.bam | perl -p -e 's/_.+bam$//' | uniq`
do
ab=`samtools view --count -F 4 -@ 10 ${file}_aliall.sort.bam`;
ap=`samtools view --count -@ 10 ${file}_aliall.sort.bam`;
bb=`samtools view --count -F 4 -@ 10 ${file}_1mis.sort.bam`;
bp=`samtools view --count -@ 10 ${file}_1mis.sort.bam`;
cb=`samtools view --count -F 4 -@ 10 ${file}_unali.sort.bam`;
cp=`samtools view --count -@ 10 ${file}_unali.sort.bam`;
aliall=`echo "scale=4;$ab*100/$ap" | bc | awk '{printf "%.4f", $0}'`;
mis=`echo "scale=4;$bb*100/$bp" | bc | awk '{printf "%.4f", $0}'`;
unali=`echo "scale=4;$cb*100/$cp" | bc | awk '{printf "%.4f", $0}'`;
echo "${file},${aliall},aliall";
echo "${file},${mis},mis1";
echo "${file},${unali},unali";
done

# bc is the shell caculator for floating-point calculation
# -c, --count: Print only the count of matching records
# awk can modify numbers because bc won't show you the 0 before the decimal point
# *100: percentage