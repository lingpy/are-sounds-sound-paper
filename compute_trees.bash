for f in `ls cognate_classes_phylip`
do
    db=`basename $f .phy`
    raxml-ng --msa cognate_classes_phylip/$f --model BIN --prefix cognate_classes_mltree/$db.phy --redo
done


for f in `ls correspondences_phylip`
do
    db=`basename $f .phy`
    raxml-ng --msa correspondences_phylip/$f --model BIN --prefix correspondences_mltree/$db.phy --redo
done


