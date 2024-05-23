#!/bin/bash

# ==================================
# === install conda environments ===
# ==================================

for env in default phesant; do
	if [ ! -d "envs/${env}" ]; then
		mamba env create --file envs/${env}.yml -p envs/${env}
	fi
done

# =============================
# === prepare genotype data ===
# =============================

# set working directory and load conda environment
cd /slow/projects/ukb_faah
conda activate envs/default

# create genotypes folder and download utility 'ukbgene'
mkdir -p data/genotypes
wget -P data/genotypes biobank.ndph.ox.ac.uk/showcase/util/ukbgene; chmod 770 data/genotypes/ukbgene

# download meta and relatedness files
keyfile=$(ls data/basket/20240307_4017567/data/*.key)
mkdir -p data/genotypes/meta
wget -P data/genotypes/meta biobank.ndph.ox.ac.uk/ukb/ukb/docs/ukb_genetic_data_description.txt
wget -P data/genotypes/meta biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_bim.tar
wget -P data/genotypes/meta biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_qc.txt
./data/genotypes/ukbgene rel -a./${keyfile}; \mv *rel* data/genotypes/meta/ # get relatedness

# get genotype .bed .bim .fam
mkdir -p data/genotypes/chr1
./data/genotypes/ukbgene cal -c1 -a./${keyfile}; \mv *.bed data/genotypes/chr1/ # get genotype calls (.bed file)
./data/genotypes/ukbgene cal -c1 -a./${keyfile} -m; \mv *.fam data/genotypes/chr1/ukb_cal_chr1_v2.fam # get individual ids (.fam file)
tar -xvf data/genotypes/meta/ukb_snp_bim.tar ukb_snp_chr1_v2.bim; \mv *.bim data/genotypes/chr1/ukb_cal_chr1_v2.bim

# extract faah snp
plink \
--bfile data/genotypes/chr1/ukb_cal_chr1_v2 \
--extract <(echo rs324420) \
--make-bed \
--recodeA \
--freq \
--missing \
--out data/genotypes/rs324420

# ============================
# === run sample filtering ===
# ============================

# set working directory and load conda environment
cd /slow/projects/ukb_faah
conda activate envs/default

# run sample filtering (unrelated, no reported vs. genetic sex mismatch, no sex aneuploidy, no outliers in heterozygosity and missingness, etc.)
basketFile="data/basket/20240307_4017567/data/ukb678162.RData"
kinshipFile="data/genetics/2024/meta/ukb42032_rel_s487950.dat"
panFile="data/basket/20210327_2008261/data/Files for retman/all_pops_non_eur_pruned_within_pop_pc_covs.tsv"
bridgeFile="data/basket/20210327_2008261/data/ukb42032bridge31063.txt"
targetDir="results/sample/"
Rscript code/sampleFiltering.R "${basketFile}" "${kinshipFile}" "${panFile}" "${bridgeFile}" "${targetDir}"

# prepare phenotype file in ukb basket (white-British only)
cd data/basket/20240307_4017567/data
idFile='results/sample/iid.discovery.txt'
basketFile="data/basket/20240307_4017567/data/ukb678162.csv"
outFile="results/sample/ukb678162.phesant.csv"
awk -F',' 'NR==FNR { id[$1]; next } FNR==1 { print } $1 in id { print }' OFS='\t' <(awk 'NR > 1 { print "\""$1"\"" }' "${idFile}") "${basketFile}" > "${outFile}"
awk -F',' 'NR==1 { gsub(/-/,"."); gsub(/[.]/,"_"); gsub(/,"/,",\"x"); print; next} NR > 1 { print }' "${outFile}" > "${outFile}.tmp"
\mv "${outFile}.tmp" "${outFile}"
chmod 770 "${outFile}"

head -1 ukb678162.phesant.csv | awk -F',' '{for(i=1; i<=NF; i++) {if($i~/x6138_0_0/) print i}}' 

# prepare pain phenotype file
idFile='results/sample/iid.discovery.txt'
rapFile='data/rap/pain_participant.csv'
outFile="results/sample/pain.phesant.csv"
awk -F',' 'NR==FNR { id[$1]; next } FNR==1 { print } "\""$1"\"" in id { print }'  <(awk 'NR > 1 { print "\""$1"\"" }' "${idFile}") "${rapFile}" > "${outFile}"
awk -F',' 'NR==1 { gsub(/-/,"."); gsub(/[.]/,"_"); gsub(/,/,"\",\"x") } 
		   NR >1 { gsub(",","\",\"") }
		   { print "\""$0"\"" }' "${outFile}" > "${outFile}.tmp"
awk -F',' '{out=$1; for(i=2; i<=NF; i++) {if(i!=41) out=out","$i}; print out}' "${outFile}.tmp" > "${outFile}"
rm -f "${outFile}.tmp"

# extract sample genotypes (all individuals)
idFile='results/sample/r2024.vars.gz'
genotypeFile='data/genotypes/rs324420.raw'
outFile="results/sample/r2024.genotypes.txt"
awk 'NR==1 { next }
	 NR==FNR { id[$1]; next} 
	 FNR==1 { print $1, $2, $NF }
	 $1 in id { print $1, $2, $NF }' OFS='\t' <(zcat "${idFile}") "${genotypeFile}" > "${outFile}"
chmod 770 "${outFile}"

# ============================================
# === run cross-trait association analysis ===
# ============================================

# set working directory and activate conda environment
cd /slow/projects/ukb_faah # cd /home/groups/markett/ukb_faah
conda activate envs/phesant

# settings
trait="rs324420_A"
traitFile="results/sample/r2024.genotypes.txt"
zcat results/sample/r2024.vars.gz > results/sample/r2024.vars.txt; covsFile="results/sample/r2024.vars.txt"
covs=$(echo "sex,age,age2,ageXsex,age2Xsex,array,$(echo ac{1..21} | sed 's/ /,/g' ),$(echo PC{1..20} | sed 's/ /,/g' )")
covs=$(echo "sex,age,age2,ageXsex,age2Xsex,array,$(echo ac{1..21} | sed 's/ /,/g' ),$(echo PC{1..10} | sed 's/ /,/g' )")
covs="sex,age,age2,ageXsex,age2Xsex"
targetDir="results/phesant"
phenoFile="results/sample/ukb678162.phesant.csv"
phesantDir="/home/groups/markett/software/PHESANT" # phesantDir="/fast/software/PHESANT"
nparts=25
standardise="FALSE"
genetic="TRUE"
sensitivity="TRUE"
taskset -c 0-24 ./code/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}" "${phesantDir}" "${nparts}" "${standardise}" "${genetic}" "${sensitivity}"

# run pain analysis
trait="rs324420_A"
traitFile="results/sample/r2024.genotypes.txt"
zcat results/sample/r2024.vars.gz > results/sample/r2024.vars.txt; covsFile="results/sample/r2024.vars.txt"
covs=$(echo "sex,age,age2,ageXsex,age2Xsex,array,$(echo ac{1..21} | sed 's/ /,/g' ),$(echo PC{1..20} | sed 's/ /,/g' )")
covs=$(echo "sex,age,age2,ageXsex,age2Xsex,array,$(echo ac{1..21} | sed 's/ /,/g' ),$(echo PC{1..10} | sed 's/ /,/g' )")
covs="sex,age,age2,ageXsex,age2Xsex"
targetDir="results/phesant.pain"
phenoFile="results/sample/pain.phesant.csv"
phesantDir="/home/groups/markett/software/PHESANT" # phesantDir="/fast/software/PHESANT" # 
nparts=1
standardise="FALSE"
genetic="TRUE"
sensitivity="TRUE"
taskset -c 0-28 ./code/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}" "${phesantDir}" "${nparts}" "${standardise}" "${genetic}" "${sensitivity}"

# combine 
mkdir -p results/combined
awk 'NR==1 { print } FNR==1 { next } { print }' results/phesant/phesant.output/results-combined.txt results/phesant.pain/phesant.output/results-combined.txt > results/combined/phesant.combined.txt

# create phewas plot and get result summary
trait="rs324420_A"
phesantResults="results/combined/phesant.combined.txt"
ncovs=34
targetDir="results/combined/"
imaging="TRUE"
multipleTesting="bonferroni"
ylim=8
ybreaks=2
width=12.7
height=3.8
repel_nudge_y=10
desat=FALSE
code/phesant.output.R "${trait}" "${phesantResults}" "${ncovs}" "${targetDir}" "${imaging}" "${multipleTesting}" "${ylim}" "${ybreaks}" "${width}" "${height}" "${repel_nudge_y}" "${desat}" 

# create phewas supplementum table
conda activate envs/default
traits="rs324420_A"
phesantSummary="results/combined/phesant.summary.txt"
phesantPlot=""
outputFile="results/combined/phewas"
code/phesant.combine.R "${traits}" "${phesantSummary}" "${phesantPlot}" "${outputFile}"

# create qq plots of phewas categories
inputFile="results/combined/phesant.summary.txt"
outFile="results/combined/phewas.qq.png"
categoryCol="custom_category"
categories="Blood and urine assays,Cancer register,Cognitive function,Diet by 24-hour recall,Digestive health,Family history and early life factors,Hospital Inpatient - Administration,Imaging,Lifestyle and environment,Maternity and sex-specific factors,Medical history and conditions,Medications,Mental health,Operations,Physical activity,Physical measures,Sociodemographics,Work environment"
categoriesRename="Blood and urine,Cancer register,Cognitive function,Diet,Digestive health,Family history,Hospital Inpatient,Imaging,Lifestyle,Maternity,Medical history,Medications,Mental health,Operations,Physical activity,Physical measures,Sociodemographics,Work environment"
pCol="pvalue"
ylim=6
ysteps=2
xlim=4.5
xsteps=2
patchCols=6
width=7
height=4.5
code/phesant.qq.R "${inputFile}" "${outFile}" "${categoryCol}" "${categories}" "${categoriesRename}" "${pCol}" "${ylim}" "${ysteps}" "${xlim}" "${xsteps}" "${patchCols}" "${width}" "${height}"

# draw qq plot for category 'pain and anxiety'
outFile="results/combined/phewas.anxiety.png"
categories="Anxiety and Pain"
categoriesRename="Anxiety and Pain"
ylim=6
ysteps=2
xlim=3
xsteps=2
width=2
height=3
code/phesant.qq.R "${inputFile}" "${outFile}" "${categoryCol}" "${categories}" "${categoriesRename}" "${pCol}" "${ylim}" "${ysteps}" "${xlim}" "${xsteps}" "${patchCols}" "${width}" "${height}"

# make surface plots
matlab -r "\
	workingDir = '$(pwd)';\
	ENIGMAtoolboxPath = '/home/groups/markett/software/ENIGMA/matlab';\
	surfCorr = 'results/combined/phesant.summary.txt';\
	mappingFileSurfArea = 'code/surfplot.mapping.surfarea.txt';\
	mappingFileThickAvg = 'code/surfplot.mapping.thickavg.txt';\
	mappingFileGrayVol = 'code/surfplot.mapping.grayvol.txt';\
	mappingFileSubcortical = 'code/surfplot.mapping.subcortical.txt';\
	rCol = 'rho';\
	pCol = 'pvalue';\
	barTitle = 'Correlation (r)';\
	outFile = 'results/combined/surfplot.png';\
	colorBar = 'horizontal';\
	cblim = '0.03';\
	run code/surfplot.m;\
	exit"
