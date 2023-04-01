for f in `ls cognate_classes_phylip`
do
    db=`basename $f .phy`
    raxml-ng --msa cognate_classes_phylip/$f --model BIN --prefix cognate_classes_mltree/$db.phy --redo
    cp cognate_classes_mltree/$db.phy.raxml.bestTree cognate_classes_mltree/$db"_ml.tre"
done


for f in `ls correspondences_phylip`
do
    db=`basename $f .phy`
    raxml-ng --msa correspondences_phylip/$f --model BIN --prefix correspondences_mltree/$db.phy --redo
    cp correspondences_mltree/$db.phy.raxml.bestTree correspondences_mltree/$db"_ml.tre"
done


