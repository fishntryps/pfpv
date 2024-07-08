
### this script was run on Broad's UGER (Univa GridEngine for Research) cluster using the Linux distribution RedHat 7.7 x86_64 

### UGER was run interactively with access to 2 cores (-pe smp 2 -binding linear:2) and 4g memory (h_vmem=4g)

### commands run in seconds to minutes

##############################################################################################################
### the script tabulates frequency of P. falciparum >0.50 IBD sharing within and between epidemiological zones ####
##############################################################################################################

cd /seq/plasmodium/pschwabl/Pf_IBD_g50_v4/

cut -f1-10 ../Pf_IBD_g90/Pf_test_for_R_only_p025_mono_no_ge10clones.txt > myIBD.txt

##### requisites

cp ../Pf_IBD_g50_v3/CLEAN* .
cp ../Pf_IBD_g50_v3/degr* .
sed -re 's/^(Venez.*)Cristinas_Border(.*)/\1Venezuela\2/g' ../Pf_IBD_g50_v3/COORDS_grouped_v3.tab > COORDS_grouped_v4.tab


#cut -f4 myIBD.txt | while read i; do awk '$1=="'$i'"' COORDS_grouped_v4.tab; done > s1_groups_v4.txt
#cut -f5 myIBD.txt | while read i; do awk '$1=="'$i'"' COORDS_grouped_v4.tab; done > s2_groups_v4.txt


paste s1_groups_v4.txt s2_groups_v4.txt | cut -f6,7,15,16 | paste myIBD.txt - > myIBD_grouped_v4.txt



cut -f12,14 myIBD_grouped_v4.txt  | perl -pe 's/\t/\n/g' | sort | uniq > Pf_groups_v4.txt




threshold=0.50
threshasperc=50

cat Pf_groups_v4.txt | while read i; do awk '$7=="'$i'"' COORDS_grouped_v4.tab | awk '{sum1+=$2; sum2+=$3} END {print  "'$i'\t"sum1/NR"\t"sum2/NR}' ; done > groups_v4_centroids.txt

cat Pf_groups_v4.txt | while read i; do cat Pf_groups_v4.txt | while read j; do awk '($12=="'$i'" && $14=="'$j'") || ($14=="'$i'" && $12=="'$j'")' myIBD_grouped_v4.txt > "$i"_"$j"_grocomp_v4.txt; awk '$3>"'$threshold'"'  "$i"_"$j"_grocomp_v4.txt | wc -l | awk '{print "'$i'\t'$j'\t"$0}' | paste - <( cat "$i"_"$j"_grocomp_v4.txt | wc -l); done; done | awk ' $4!=0 {print $0"\t"$3/$4}' | while read line; do echo "$line" | cut -f1,2 | perl -pe 's/\t/\n/g' | sort -V | perl -pe 's/\n/\t/g' | paste - <(echo "$line" | cut -f3-); done |  sort | uniq | sort -nk5,5 > perc_g"$threshasperc"_amoung_groups_v4.txt


cut -f1,2 perc_g"$threshasperc"_amoung_groups_v4.txt | while read i; do g1=$(echo "$i" | cut -f1); g2=$(echo "$i" | cut -f2); paste <(awk '$1=="'$g1'"' groups_v4_centroids.txt) <(awk '$1=="'$g2'"' groups_v4_centroids.txt) | awk '{print sqrt((($2-$5)^2) + (($3-$6)^2))}'; done | paste perc_g"$threshasperc"_amoung_groups_v4.txt - > perc_g"$threshasperc"_amoung_groups_v4_with_distances.txt


cat perc_g"$threshasperc"_amoung_groups_v4_with_distances.txt | while read i; do L1=$(echo "$i" | cut -f1); L2=$(echo "$i" | cut -f2); awk '$3>"'$threshold'"' myIBD_grouped_v4.txt | awk '($12=="'$L1'" && $14=="'$L2'") || ($12=="'$L2'" && $14=="'$L1'")' | awk '{print $1"xxxx"$12"\t"$2"xxxx"$14}' | perl -pe 's/\t/\n/g' | sort | uniq > temp2_"$L1"_"$L2".txt; paste <(egrep -c "xxxx$L1$" temp2_"$L1"_"$L2".txt) <(egrep -c "xxxx$L2$" temp2_"$L1"_"$L2".txt) | paste <(echo "$i") -; done > perc_g"$threshasperc"_amoung_groups_v4_with_distances_and_indcount.txt


use R-4.0

cat perc_g"$threshasperc"_amoung_groups_v4_with_distances_and_indcount.txt | while read i; do L1=$(echo "$i" | cut -f1); L2=$(echo "$i" | cut -f2); awk '$3>"'$threshold'"' myIBD_grouped_v4.txt | awk '($12=="'$L1'" && $14=="'$L2'") || ($12=="'$L2'" && $14=="'$L1'")' | awk '{print $1"\t"$2}' | perl -pe 's/\t/\n/g' | sort | uniq | perl -pe 's/\n/"||\$1=="/g' | sed 's/||\$1\=\=\"$/\)/g' | sed 's/^/\(\$1==\"/g' | awk '{print $0"\n"$0}' | sed '2s/\$1/\$2/g' | perl -pe 's/\n/ && /g' | sed 's/ && $//g'| sed -re 's/(.*)/awk \x27\1\x27 myIBD_grouped_v4.txt/g' | bash | awk '$3>0.9 {print $1"\t"$2"\t"$12"\t"$14}' > both_"$L1"_"$L2"_for_igraph.txt

cat <(echo -e 'hold\thold\thold\thold')  both_"$L1"_"$L2"_for_igraph.txt | awk '($3=="'$L1'" && $4=="'$L1'") || $1=="hold"' | cut -f1,2 > "$L1"_for_igraph.txt
cat <(echo -e 'hold\thold\thold\thold') both_"$L1"_"$L2"_for_igraph.txt | awk '($3=="'$L2'" && $4=="'$L2'") || $1=="hold"' | cut -f1,2 > "$L2"_for_igraph.txt

echo "library(igraph)
l1 <- read.table(\""$L1"_for_igraph.txt\", header=F)
l1g <- graph_from_edgelist(as.matrix(l1),directed=F)
l2 <- read.table(\""$L2"_for_igraph.txt\", header=F)
l2g <- graph_from_edgelist(as.matrix(l2),directed=F)
write.table(cbind(components(l1g)\$no-1, components(l2g)\$no-1), file=\"both_"$L1"_"$L2"_numcc.txt\", row.names=F, col.names=F, sep=\"\\t\")" > command_both_"$L1"_"$L2".txt
Rscript command_both_"$L1"_"$L2".txt
paste <(echo "$i") both_"$L1"_"$L2"_numcc.txt; done > Pf_perc_g"$threshasperc"_amoung_groups_v4_with_distances_and_indcount_and_cccount.txt



cat Pf_groups_v4.txt | while read i; do awk '{print $1"\t"$12"xxx"$2"\t"$14}' myIBD_grouped_v4.txt | perl -pe 's/xxx/\n/g' | sort | uniq | egrep -c "[[:space:]]$i$" | awk '{print "'$i'\t"$0}'; done > Pf_sampsizes_grouped_v4_myIBD.txt

cut -f1 Pf_sampsizes_grouped_v4_myIBD.txt | while read i; do awk '$3 > "'$threshold'" && $12=="'$i'" && $14=="'$i'" {print $1"xxx"$2}' myIBD_grouped_v4.txt | perl -pe 's/xxx/\n/g' | sort | uniq | wc -l; done | paste Pf_sampsizes_grouped_v4_myIBD.txt - | sed '1iloc\tnumsamps\tnumlocallyclonal' > Pf_sampsizes_grouped_v4_myIBD_with_locallyclonalnums.txt


