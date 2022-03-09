echo -e "name\tgroup\tall\trna\tcatgry";
# echo -e: characters preceded by a slash will be escaped characters
for file in `ls *.tsv | perl -p -e 's/\.tsv//'`
do
name=${file%%_*};
catgry=${file#*_};
# ${%%},${#} were changing variable in shell
cat ${file}.tsv | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 4 --sum 2,3 | sort | \
tsv-filter --ne 3:0 | \
tsv-filter --not-empty 1 | \
awk -v name=$name -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"$3"\t"catgry}'
# awk -v: pass external variables to awk command, otherwise there will be mistakes
done