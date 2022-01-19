scrpath=$(dirname "$0")
dstpath=$(dirname "$1")
name=$(basename $1 .gtf)

$scrpath/exonRefine -p $dstpath"/"$name.r $1 
$scrpath/groupGenes $dstpath"/"$name.r.gtf $dstpath"/"$name.rg.gtf 

rm $dstpath"/"$name.r.gtf
