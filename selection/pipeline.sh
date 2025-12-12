#!/bin/bash
#SBATCH --mem=150gb
#SBATCH --ntasks=1
#SBATCH --job-name=codeml_pipeline
#SBATCH --mail-type=END
#SBATCH --mail-user=slq4@kent.ac.uk
#SBATCH --output=/home/slq4/codeML.rerun/codeml_pipeline._%j.out
#SBATCH -p biosoc2

source /home/slq4/anaconda3/bin/activate trf-env # environment with translatorx, raxml, and codeML installed

echo "start of pipeline shell script with ${1}.fasta.txt"

mkdir -p ${1}.alignments
cd ${1}.alignments
cp ../species.tree species.tree


/home/slq4/codeML.rerun/selection/translatorX.pl -i ../${1}.fasta.txt -p F -o mafft_translatorx #align fasta sequences 

/home/slq4/codeML.rerun/selection/pal2nal.pl mafft_translatorx.aa_ali.fasta mafft_translatorx.nt_ali.fasta -output fasta > ${1}_pal2nal_mafft_out.fasta 

/home/slq4/codeML.rerun/selection/one_line_fasta.pl ${1}_pal2nal_mafft_out.fasta

mv mafft_translatorx.nt_ali.fasta ${1}_nuc_mafft_aln.fasta
mv mafft_translatorx.aa_ali.fasta ${1}_prot_mafft_aln.fasta

raxmlHPC -f a -m GTRGAMMA -p 12345 -# 100 -x 12345 -# 500 -s ${1}_nuc_mafft_aln.fasta -n ${1}_newick

printf "11  1\n" > ${1}_gene.tree
sed 's/\:[0-9]*\.[0-9]*//g' RAxML_bestTree.${1}_newick >> ${1}_gene.tree

num=$( grep '>' ${1}_nuc_mafft_aln.fasta | wc -l )
len=$( sed -n '2,2p' ${1}_nuc_mafft_aln.fasta | sed 's/\r//' | sed 's/\n//' | wc --m )
perl /home/slq4/codeML.rerun/selection/FASTAtoPHYL.pl ${1}_nuc_mafft_aln.fasta $num $len
mv ${1}_nuc_mafft_aln.phy ${1}_mafft.phy

# single species
cp species.tree ${1}_branch_human.tree 
cp species.tree ${1}_branch_baboon.tree 
cp species.tree ${1}_branch_mouse.tree 
cp species.tree ${1}_branch_rat.tree 
cp species.tree ${1}_branch_sheep.tree 
cp species.tree ${1}_branch_cow.tree 
cp species.tree ${1}_branch_pig.tree 
cp species.tree ${1}_branch_chicken.tree 
cp species.tree ${1}_branch_opossum.tree 
cp species.tree ${1}_branch_chimp.tree
cp species.tree ${1}_branch_gibbon.tree

#sister pairs 
cp species.tree ${1}_branch_humanbaboon.tree 
cp species.tree ${1}_branch_mouserat.tree 
cp species.tree ${1}_branch_cowsheep.tree
cp species.tree ${1}_branch_cowsheeppig.tree 

# primate and mammal
cp species.tree ${1}_branch_primate.tree  
cp species.tree ${1}_branch_mammal.tree  
cp species.tree ${1}_branch_placental.tree




# Add tags 
sed -i 's/human/human#1/' ${1}_branch_human.tree 
sed -i 's/baboon/baboon#1/' ${1}_branch_baboon.tree 
sed -i 's/mouse/mouse#1/' ${1}_branch_mouse.tree 
sed -i 's/rat/rat#1/' ${1}_branch_rat.tree 
sed -i 's/sheep/sheep#1/' ${1}_branch_sheep.tree 
sed -i 's/cow/cow#1/' ${1}_branch_cow.tree 
sed -i 's/pig/pig#1/' ${1}_branch_pig.tree 
sed -i 's/chicken/chicken#1/' ${1}_branch_chicken.tree 
sed -i 's/opossum/opossum#1/' ${1}_branch_opossum.tree 
sed -i 's/chimp/chimp#1/' ${1}_branch_chimp.tree
sed -i 's/gibbon/gibbon#1/' ${1}_branch_gibbon.tree

sed -i 's/human/human#1/' ${1}_branch_humanbaboon.tree 
sed -i 's/baboon/baboon#1/' ${1}_branch_humanbaboon.tree 
sed -i 's/mouse/mouse#1/' ${1}_branch_mouserat.tree  
sed -i 's/rat/rat#1/' ${1}_branch_mouserat.tree 
sed -i 's/cow/cow#1/' ${1}_branch_cowsheep.tree
sed -i 's/sheep/sheep#1/' ${1}_branch_cowsheep.tree
sed -i 's/cow/cow#1/' ${1}_branch_cowsheeppig.tree
sed -i 's/sheep/sheep#1/' ${1}_branch_cowsheeppig.tree
sed -i 's/pig/pig#1/' ${1}_branch_cowsheeppig.tree

sed -i 's/human/human#1/' ${1}_branch_primate.tree 
sed -i 's/chimp/chimp#1/' ${1}_branch_primate.tree 
sed -i 's/baboon/baboon#1/' ${1}_branch_primate.tree
sed -i 's/gibbon/gibbon#1/' ${1}_branch_primate.tree

sed -i 's/human/human#1/' ${1}_branch_mammal.tree 
sed -i 's/chimp/chimp#1/' ${1}_branch_mammal.tree 
sed -i 's/gibbon/gibbon#1/' ${1}_branch_mammal.tree 
sed -i 's/baboon/baboon#1/' ${1}_branch_mammal.tree
sed -i 's/mouse/mouse#1/' ${1}_branch_mammal.tree 
sed -i 's/pig/pig#1/' ${1}_branch_mammal.tree
sed -i 's/cow/cow#1/' ${1}_branch_mammal.tree 
sed -i 's/opossum/opossum#1/' ${1}_branch_mammal.tree
sed -i 's/sheep/sheep#1/' ${1}_branch_mammal.tree
sed -i 's/rat/rat#1/' ${1}_branch_mammal.tree


sed -i 's/human/human#1/' ${1}_branch_placental.tree
sed -i 's/chimp/chimp#1/' ${1}_branch_placental.tree
sed -i 's/gibbon/gibbon#1/' ${1}_branch_placental.tree
sed -i 's/baboon/baboon#1/' ${1}_branch_placental.tree
sed -i 's/mouse/mouse#1/' ${1}_branch_placental.tree
sed -i 's/pig/pig#1/' ${1}_branch_placental.tree
sed -i 's/cow/cow#1/' ${1}_branch_placental.tree
sed -i 's/sheep/sheep#1/' ${1}_branch_placental.tree
sed -i 's/rat/rat#1/' ${1}_branch_placental.tree



for species in human baboon mouse rat sheep cow pig chicken opossum chimp gibbon; do
    mkdir -p branch_model/${species}/${1}_CODEML/
    cp ../codeML_branch_ctrl.txt branch_model/${species}/${1}_CODEML/
    cp ${1}_mafft.phy branch_model/${species}/${1}_CODEML/mafft.phy
    cp ${1}_branch_${species}.tree branch_model/${species}/${1}_CODEML/branch${species}.tree
    sed -i "s/<parameter1>/mafft.phy/" branch_model/${species}/${1}_CODEML/codeML_branch_ctrl.txt
    sed -i "s/<parameter2>/branch${species}.tree/" branch_model/${species}/${1}_CODEML/codeML_branch_ctrl.txt
    sed -i "s/<parameter3>/codeML.out/" branch_model/${species}/${1}_CODEML/codeML_branch_ctrl.txt
done

mkdir -p null_model/${1}_CODEML/
mkdir -p branch_model/human_baboon/${1}_CODEML/
mkdir -p branch_model/mouse_rat/${1}_CODEML/
mkdir -p branch_model/cow_sheep/${1}_CODEML/
mkdir -p branch_model/cow_sheep_pig/${1}_CODEML/
mkdir -p branch_model/primate/${1}_CODEML/
mkdir -p branch_model/mammal/${1}_CODEML/
mkdir -p branch_model/placental/${1}_CODEML/

cp ../codeML_M0_ctrl.txt null_model/${1}_CODEML/  
cp ../codeML_branch_ctrl.txt branch_model/human_baboon/${1}_CODEML/  
cp ../codeML_branch_ctrl.txt branch_model/mouse_rat/${1}_CODEML/
cp ../codeML_branch_ctrl.txt branch_model/cow_sheep/${1}_CODEML/
cp ../codeML_branch_ctrl.txt branch_model/cow_sheep_pig/${1}_CODEML/ 
cp ../codeML_branch_ctrl.txt branch_model/primate/${1}_CODEML/ 
cp ../codeML_branch_ctrl.txt branch_model/mammal/${1}_CODEML/ 
cp ../codeML_branch_ctrl.txt branch_model/placental/${1}_CODEML/

cp ${1}_mafft.phy null_model/${1}_CODEML/mafft.phy
cp ${1}_mafft.phy branch_model/human_baboon/${1}_CODEML/mafft.phy
cp ${1}_mafft.phy branch_model/mouse_rat/${1}_CODEML/mafft.phy
cp ${1}_mafft.phy branch_model/cow_sheep/${1}_CODEML/mafft.phy
cp ${1}_mafft.phy branch_model/cow_sheep_pig/${1}_CODEML/mafft.phy
cp ${1}_mafft.phy branch_model/primate/${1}_CODEML/mafft.phy
cp ${1}_mafft.phy branch_model/mammal/${1}_CODEML/mafft.phy
cp ${1}_mafft.phy branch_model/placental/${1}_CODEML/mafft.phy

cp species.tree null_model/${1}_CODEML/species.tree  
cp ${1}_branch_humanbaboon.tree  branch_model/human_baboon/${1}_CODEML/branchhumanbaboon.tree 
cp ${1}_branch_mouserat.tree branch_model/mouse_rat/${1}_CODEML/branchmouserat.tree 
cp ${1}_branch_cowsheep.tree  branch_model/cow_sheep/${1}_CODEML/branchcowsheep.tree 
cp ${1}_branch_cowsheeppig.tree branch_model/cow_sheep_pig/${1}_CODEML/branchcowsheeppig.tree
cp ${1}_branch_primate.tree branch_model/primate/${1}_CODEML/branchprimate.tree #test for all primate positive selection
cp ${1}_branch_mammal.tree branch_model/mammal/${1}_CODEML/branchmammal.tree #test for all primate positive selection
cp ${1}_branch_placental.tree branch_model/placental/${1}_CODEML/branchplacental.tree

sed -i 's/<parameter1>/mafft.phy/' null_model/${1}_CODEML/codeML_M0_ctrl.txt 
sed -i 's/<parameter2>/species.tree/' null_model/${1}_CODEML/codeML_M0_ctrl.txt 
sed -i 's/<parameter3>/codeML.out/' null_model/${1}_CODEML/codeML_M0_ctrl.txt 

sed -i 's/<parameter1>/mafft.phy/' branch_model/human_baboon/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter2>/branchhumanbaboon.tree/' branch_model/human_baboon/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter3>/codeML.out/' branch_model/human_baboon/${1}_CODEML/codeML_branch_ctrl.txt

sed -i 's/<parameter1>/mafft.phy/' branch_model/mouse_rat/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter2>/branchmouserat.tree/' branch_model/mouse_rat/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter3>/codeML.out/' branch_model/mouse_rat/${1}_CODEML/codeML_branch_ctrl.txt

sed -i 's/<parameter1>/mafft.phy/' branch_model/cow_sheep/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter2>/branchcowsheep.tree/' branch_model/cow_sheep/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter3>/codeML.out/' branch_model/cow_sheep/${1}_CODEML/codeML_branch_ctrl.txt

sed -i 's/<parameter1>/mafft.phy/' branch_model/cow_sheep_pig/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter2>/branchcowsheeppig.tree/' branch_model/cow_sheep_pig/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter3>/codeML.out/' branch_model/cow_sheep_pig/${1}_CODEML/codeML_branch_ctrl.txt

sed -i 's/<parameter1>/mafft.phy/' branch_model/primate/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter2>/branchprimate.tree/' branch_model/primate/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter3>/codeML.out/' branch_model/primate/${1}_CODEML/codeML_branch_ctrl.txt

sed -i 's/<parameter1>/mafft.phy/' branch_model/mammal/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter2>/branchmammal.tree/' branch_model/mammal/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter3>/codeML.out/' branch_model/mammal/${1}_CODEML/codeML_branch_ctrl.txt

sed -i 's/<parameter1>/mafft.phy/' branch_model/placental/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter2>/branchplacental.tree/' branch_model/placental/${1}_CODEML/codeML_branch_ctrl.txt
sed -i 's/<parameter3>/codeML.out/' branch_model/placental/${1}_CODEML/codeML_branch_ctrl.txt



cd null_model/${1}_CODEML/
codeml codeML_M0_ctrl.txt
cd ../../

branches=("human" "baboon" "mouse" "rat" "sheep" "cow" "pig" "chicken" "chimp" "gibbon" "opossum" "human_baboon" "mouse_rat" "cow_sheep" "cow_sheep_pig"  "primate" "mammal" "placental")
for branch in "${branches[@]}"; do
    cd branch_model/$branch/${1}_CODEML/
    codeml codeML_branch_ctrl.txt
    cd ../../../
done


grep 'lnL' null_model/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' > ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/human/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/baboon/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/mouse/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/rat/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/sheep/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/cow/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/pig/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/chicken/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/chimp/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/gibbon/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/opossum/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/human_baboon/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/mouse_rat/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/cow_sheep/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/cow_sheep_pig/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/primate/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/mammal/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt
grep 'lnL' branch_model/placental/${1}_CODEML/codeML.out | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> ${1}_lnL_branch_mods.txt

#head -2 ${1}_nuc_mafft_aln.fasta | tail -n 1 > mafft.alignment
#head -2 ${1}_one_line.fa | tail -n 1 > pal2nal.alignment

w_h_back=$( grep 'w (dN/dS) for branches'  branch_model/human/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_h_for=$( grep 'w (dN/dS) for branches'  branch_model/human/${1}_CODEML/codeML.out | sed 's/..* //')

w_b_back=$( grep 'w (dN/dS) for branches'  branch_model/baboon/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_b_for=$( grep 'w (dN/dS) for branches'  branch_model/baboon/${1}_CODEML/codeML.out | sed 's/..* //')

w_mou_back=$( grep 'w (dN/dS) for branches'  branch_model/mouse/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_mou_for=$( grep 'w (dN/dS) for branches'  branch_model/mouse/${1}_CODEML/codeML.out | sed 's/..* //')

w_r_back=$( grep 'w (dN/dS) for branches'  branch_model/rat/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_r_for=$( grep 'w (dN/dS) for branches'  branch_model/rat/${1}_CODEML/codeML.out | sed 's/..* //')

w_s_back=$( grep 'w (dN/dS) for branches'  branch_model/sheep/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_s_for=$( grep 'w (dN/dS) for branches'  branch_model/sheep/${1}_CODEML/codeML.out | sed 's/..* //')

w_c_back=$( grep 'w (dN/dS) for branches'  branch_model/cow/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_c_for=$( grep 'w (dN/dS) for branches'  branch_model/cow/${1}_CODEML/codeML.out | sed 's/..* //')

w_pig_back=$( grep 'w (dN/dS) for branches'  branch_model/pig/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_pig_for=$( grep 'w (dN/dS) for branches'  branch_model/pig/${1}_CODEML/codeML.out | sed 's/..* //')

w_ch_back=$( grep 'w (dN/dS) for branches'  branch_model/chicken/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_ch_for=$( grep 'w (dN/dS) for branches'  branch_model/chicken/${1}_CODEML/codeML.out | sed 's/..* //')

w_chimp_back=$( grep 'w (dN/dS) for branches'  branch_model/chimp/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_chimp_for=$( grep 'w (dN/dS) for branches'  branch_model/chimp/${1}_CODEML/codeML.out | sed 's/..* //')

w_g_back=$( grep 'w (dN/dS) for branches'  branch_model/gibbon/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_g_for=$( grep 'w (dN/dS) for branches'  branch_model/gibbon/${1}_CODEML/codeML.out | sed 's/..* //')

w_o_back=$( grep 'w (dN/dS) for branches'  branch_model/opossum/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_o_for=$( grep 'w (dN/dS) for branches'  branch_model/opossum/${1}_CODEML/codeML.out | sed 's/..* //')

w_hb_back=$( grep 'w (dN/dS) for branches'  branch_model/human_baboon/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_hb_for=$( grep 'w (dN/dS) for branches'  branch_model/human_baboon/${1}_CODEML/codeML.out | sed 's/..* //')

w_mr_back=$( grep 'w (dN/dS) for branches'  branch_model/mouse_rat/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_mr_for=$( grep 'w (dN/dS) for branches'  branch_model/mouse_rat/${1}_CODEML/codeML.out | sed 's/..* //')

w_cs_back=$( grep 'w (dN/dS) for branches'  branch_model/cow_sheep/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_cs_for=$( grep 'w (dN/dS) for branches'  branch_model/cow_sheep/${1}_CODEML/codeML.out | sed 's/..* //')

w_csp_back=$( grep 'w (dN/dS) for branches'  branch_model/cow_sheep_pig/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_csp_for=$( grep 'w (dN/dS) for branches'  branch_model/cow_sheep_pig/${1}_CODEML/codeML.out | sed 's/..* //')

w_p_back=$( grep 'w (dN/dS) for branches'  branch_model/primate/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_p_for=$( grep 'w (dN/dS) for branches'  branch_model/primate/${1}_CODEML/codeML.out | sed 's/..* //')

w_m_back=$( grep 'w (dN/dS) for branches'  branch_model/mammal/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_m_for=$( grep 'w (dN/dS) for branches'  branch_model/mammal/${1}_CODEML/codeML.out | sed 's/..* //')

w_pl_back=$( grep 'w (dN/dS) for branches'  branch_model/placental/${1}_CODEML/codeML.out | sed 's/..*: *//' | sed 's/ ..*//' )
w_pl_for=$( grep 'w (dN/dS) for branches'  branch_model/placental/${1}_CODEML/codeML.out | sed 's/..* //')


printf "w_back_human\tw_fore_human\tw_back_baboon\tw_fore_baboon\tw_back_mouse\tw_fore_mouse\tw_back_rat\tw_fore_rat\tw_back_sheep\tw_fore_sheep\tw_back_cow\tw_fore_cow\tw_back_pig\tw_fore_pig\tw_back_chicken\tw_fore_chicken\tw_back_chimp\tw_fore_chimp\tw_back_gibbon\tw_fore_gibbon\tw_back_opossum\tw_fore_opossum\tw_back_humanbaboon\tw_fore_humanbaboon\tw_back_mouserat\tw_fore_mouserat\tw_back_cowsheep\tw_fore_cowsheep\tw_back_cowsheeppig\tw_fore_cowsheeppig\tw_back_primate\tw_fore_primate\tw_back_mammal\tw_fore_mammal\tw_back_placental\tw_fore_placental\n" > ${1}_w_est_branches.tsv 

printf $w_h_back"\t"$w_h_for"\t"$w_b_back"\t"$w_b_for"\t"$w_mou_back"\t"$w_mou_for"\t"$w_r_back"\t"$w_r_for"\t"$w_s_back"\t"$w_s_for"\t"$w_c_back"\t"$w_c_for"\t"$w_pig_back"\t"$w_pig_for"\t"$w_ch_back"\t"$w_ch_for"\t"$w_chimp_back"\t"$w_chimp_for"\t"$w_g_back"\t"$w_g_for"\t"$w_o_back"\t"$w_o_for"\t"$w_hb_back"\t"$w_hb_for"\t"$w_mr_back"\t"$w_mr_for"\t"$w_cs_back"\t"$w_cs_for"\t"$w_csp_back"\t"$w_csp_for"\t"$w_p_back"\t"$w_p_for"\t"$w_m_back"\t"$w_m_for"\t"$w_pl_back"\t"$w_pl_for"\n" >> ${1}_w_est_branches.tsv 

Rscript --save ../stats.calc.2.R $1 ${1}_lnL_branch_mods.txt ${1}_w_est_branches.tsv 

cd ../

tail -n 1 ${1}.alignments/${1}.stats.table >> stats.table

mkdir -p results/${1}.output
mkdir -p trees/


cp ${1}.alignments/null_model/${1}_CODEML/codeML.out results/${1}.output/${1}.null.codeml.out
cp ${1}.alignments/null_model/${1}_CODEML/codeML_M0_ctrl.txt results/${1}.output

branches=("human" "baboon" "mouse" "rat" "sheep" "cow" "pig" "chicken" "chimp" "gibbon" "opossum" "human_baboon" "mouse_rat" "cow_sheep" "cow_sheep_pig"  "primate" "mammal" "placental")
for branch in "${branches[@]}"; do
	cp ${1}.alignments/branch_model/$branch/${1}_CODEML/codeML.out results/${1}.output/${1}.${branch}.codeml.out
	cp ${1}.alignments/branch_model/$branch/${1}_CODEML/codeML_branch_ctrl.txt results/${1}.output/${1}.${branch}.ctrl.out

done


cp ${1}.alignments/${1}_mafft.phy results/${1}.output/
cp ${1}.alignments/${1}_branch_human.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_baboon.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_mouse.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_rat.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_sheep.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_cow.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_pig.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_chicken.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_chimp.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_gibbon.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_opossum.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_humanbaboon.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_mouserat.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_cowsheep.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_cowsheeppig.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_primate.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_mammal.tree results/${1}.output/
cp ${1}.alignments/${1}_branch_placental.tree results/${1}.output/
cp ${1}.alignments/${1}.stats.table results/${1}.output/

cp ${1}.alignments/${1}_gene.tree trees/

rm ${1}.alignments/ 

