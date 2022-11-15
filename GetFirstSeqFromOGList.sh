#script by Tomas Carrasco
Folder="/srv/public/users/eksramos/turtles/analyses/new_curation/after_TOGA/11_orthology/Results_26jul22/Results_Jul26/Orthogroup_Sequences"
while read OG; do
sed 's\>\>'"${OG}"';\' $Folder/$OG.fa > Mod_$OG.fa
head -2 Mod_$OG.fa
rm Mod_*.fa
done < /srv/public/users/eksramos/turtles/analyses/new_curation/after_TOGA/13_cafe2/cafe5/results_all_09nov22/Base_branch_probabilities_OG_rCheMyd1_significant.tab
