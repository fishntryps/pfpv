
##################################################################################################################
### this script prepares analysis input files from a joint call VCF containing all sites with 1+ minor alleles ###
##################################################################################################################

### this script was run on Broad's UGER (Univa GridEngine for Research) cluster using the Linux distribution RedHat 7.7 x86_64 

### UGER was run with 12-hour access to 6 cores (-pe smp 6 -binding linear:6), 16g memory (h_vmem=16g)

### most commands run in seconds to minutes except hmmIBD which requires multiple hours

### the script uses bcftools v1.20 and vcftools v0.1.15 via dotkit and hmmIBD v2.0.4 and THEREALMcCOIL v2 via github

### this script relates to P. falciparum, but the same procedures was used to process a P. vivax joint call

### first we find INDELS and set mask plus minus 5nt around them

vcftools --vcf /gsap/garage-protistvector/pschwabl/Pf_lists_for_June_jointcall/Pf_Guyana_SWGA2020_Jan22_Venez.snp.indel.recalibrated.filtered.ALL.Unsorted.mac1Pass.recode.vcf --remove-indels --recode --recode-INFO-all --out Guyanaplus_Pf

vcftools --vcf /gsap/garage-protistvector/pschwabl/Pf_lists_for_June_jointcall/Pf_Guyana_SWGA2020_Jan22_Venez.snp.indel.recalibrated.filtered.ALL.Unsorted.mac1Pass.recode.vcf --keep-only-indels --recode --recode-INFO-all --out Guyanaplus_Pf_INDELS
egrep -v '#' Guyanaplus_Pf_INDELS.recode.vcf | awk '{print $1"\t"$2-5"\t"$2+5}' | sed '1ichrom\tchromStart\tchromEnd'> indel_mask.txt

vcftools --vcf Guyanaplus_Pf.recode.vcf --exclude-bed indel_mask.txt --recode --recode-INFO-all --out Guyanaplus_Pf_5ntINDEL

### then exclude noncore regions

vcftools --vcf Guyanaplus_Pf_5ntINDEL.recode.vcf --bed /seq/plasmodium/data/bed/Miles_2016_Pf3D7_core_only.bed --recode --recode-INFO-all --out Guyanaplus_Pf_5ntINDEL_CORE

### identify and keep only good samples with max 50 missing sites
vcftools --vcf Guyanaplus_Pf_5ntINDEL_CORE.recode.vcf --maf 0.02 --out Guyanaplus_Pf_5ntINDEL_CORE_maf02 --missing-indv
sort -nk5,5 Guyanaplus_Pf_5ntINDEL_CORE_maf02.imiss | awk '$5<=0.5' | cut -f1 > good_samples_ie_max50miss_applying_maf02_to_Guyanaplus_Pf_5ntINDEL_CORE.txt
vcftools --vcf Guyanaplus_Pf_5ntINDEL_CORE.recode.vcf --keep good_samples_ie_max50miss_applying_maf02_to_Guyanaplus_Pf_5ntINDEL_CORE.txt --recode --recode-INFO-all --out Guyanaplus_Pf_5ntINDEL_CORE_GOOD

### identify and remove sites with >7.5 percent heterozygosity

egrep -v '#' Guyanaplus_Pf_5ntINDEL_CORE_GOOD.recode.vcf | while read i; do paste <(echo "$i" | egrep -o '0/[1-3]\:|1/[2-3]\:|2/3\:' | wc -l) <(echo "$i" | egrep -o '[0-3]/[0-3]\:' | wc -l) | paste <(echo "$i" | cut -f1,2) -; done > hetfrac_Guyanaplus_Pf_5ntINDEL_CORE_GOOD.txt

awk '$4!=0 && ($3/$4)>0.075 {print $1"\t"$2}' hetfrac_Guyanaplus_Pf_5ntINDEL_CORE_GOOD.txt > het_075_to_remove_from_Guyanaplus_Pf_5ntINDEL_CORE_GOOD.sites

vcftools --vcf Guyanaplus_Pf_5ntINDEL_CORE_GOOD.recode.vcf --mac 1 --exclude-positions het_075_to_remove_from_Guyanaplus_Pf_5ntINDEL_CORE_GOOD.sites --recode --recode-INFO-all --out Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets

### for hmmIBD using Guyana and small Venezuelan comparator set 

egrep CHROM Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets.recode.vcf | perl -pe 's/\t/\n/g' | egrep '^[G,C][0-9][A-Z]|CEM|Venez|SPT26229|PW0065\-C' > only_guyZ_Pf.txt

# MAF 0.02
vcftools --vcf Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets.recode.vcf  --keep only_guyZ_Pf.txt --maf 0.02 --recode --recode-INFO-all --out Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02

# get bi allelic
bcftools view --max-alleles 2 Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02.recode.vcf > Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl.vcf

# remove heterozygous
bcftools filter -S . -e 'GT=="het"' Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl.vcf > Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss.vcf

# remove monomorphic
bcftools view -c 1 Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss.vcf > Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed.vcf

# remove missing alleles by 30%
bcftools view -i 'F_MISSING < 0.3'  Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed.vcf > Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30.vcf

# make sure individuals are missing no more than 50% sites
vcftools --vcf Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30.vcf --missing-indv --out Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30
vcftools --vcf Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30.vcf --remove <(sed 1d Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30.imiss | awk '$5>0.5 {print $1}') --recode --recode-INFO-all --out Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30_indmax50

# remove missing alleles by 30%
bcftools view -i 'F_MISSING < 0.3'  Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30_indmax50.recode.vcf > Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30_indmax50_less30.vcf

### make input for hmmIBD (https://github.com/glipsnort/hmmIBD)

egrep -v '#' Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30_indmax50_less30.vcf | cut -f1,2,10- | sed -re 's/\:[^/]*([[:space:]].)/:\1/g' | sed -re 's/\:[^/]*$/:/g' | sed 's-1/1:-1-g' | sed 's-0/0:-0-g' | sed 's-2/2:-2-g' | sed 's_\./\.:_-1_g' | sed 's_\.:_-1_g' | sed -re 's/Pf3D7_([0-9]*)_v3/\1/g' | sed 's/^0//g' | sed 's/^M76611/15/g' | sed 's/^PFC10_API_IRAB/16/g' | sort -Vk1,2 | cat <(egrep CHROM Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30_indmax50_less30.vcf | cut -f1,2,10- | sed 's/#CHROM/chrom/g' | sed 's/POS/pos/g') - > Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30_indmax50_less30_hmmIBD_input.txt

### run hmmIBD in a screen
/gsap/garage-protistvector/pschwabl/hmmIBD-master/hmmIBD -i Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30_indmax50_less30_hmmIBD_input.txt -o Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf02_biAl_hets2miss_noFixed_less30_indmax50_less30_hmmIBD

### for polyclonality analysis

vcftools --vcf Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets.recode.vcf  --keep only_guyZ_Pf.txt --maf 0.1 --recode --recode-INFO-all --out Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10

bcftools view --max-alleles 2 Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10.recode.vcf > Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl.vcf

vcftools --vcf Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl.vcf --thin 50000 --recode --recode-INFO-all --out Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl_thin50

bcftools view -i 'F_MISSING < 0.10'  Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl_thin50.recode.vcf > Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl_thin50_less10_forCOIL.vcf

vcftools --vcf Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl_thin50_less10_forCOIL.vcf --missing-indv --out Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl_thin50_less10_forCOIL

vcftools --vcf Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl_thin50_less10_forCOIL.vcf --remove <(sed 1d Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl_thin50_less10_forCOIL.imiss | awk '$5>0.5 {print $1}') --recode --recode-INFO-all --out Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl_thin50_less10_forCOIL_indmax50

bcftools view -i 'F_MISSING < 0.10'  Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl_thin50_less10_forCOIL_indmax50.recode.vcf > Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl_thin50_less10_forCOIL_indmax50_less10.vcf

### make input for THEREALMcCOIL v2 and run with R-4.2.2 following the categorical method using instructions from github below

### https://github.com/EPPIcenter/THEREALMcCOIL/tree/master/categorical_method

egrep -v '##' Guyanaplus_Pf_5ntINDEL_CORE_GOOD_cut075Hets_onlyGuyZ_maf10_biAl_thin50_less10_forCOIL_indmax50_less10.vcf | cut -f10- | sed -re 's/\:[^/]*([[:space:]].)/:\1/g' | sed -re 's/\:[^/]*$/:/g' | sed 's-1/1:-1-g' | sed 's-0/0:-0-g' | sed 's-2/2:-1-g' | sed 's-3/3:-1-g' | sed 's_\./\.:_-1_g' | sed 's_\.:_-1_g' | sed 's_0/[0-3]:_0.5_g' | awk '{print "site"NR-1"\t"$0}' | sed '1s/site0/ind_name/' | datamash transpose > input_test_thin50.txt
