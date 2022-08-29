# from https://unix.stackexchange.com/a/490660/190054

indir=$1
for i in ${indir}*.tsv
do
  for j in ${indir}*.tsv
  do
    if [ "$i" \< "$j" ]
    then
     echo -e "$i\t$j"
    fi
  done
done
