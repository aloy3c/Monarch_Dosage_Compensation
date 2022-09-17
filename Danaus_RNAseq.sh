###QC
mkdir FastQC_raw
fastqc -t 16 -o FastQC_raw Reads_raw/*gz
multiqc FastQC_raw -o multiqc/FastQC_raw

###DE analysis with DEseq
module load R/3.3.2
R
install.packages("lattice")
source("http://bioconductor.org/biocLite.R")
biocLite('edgeR')
biocLite('limma')
biocLite('DESeq2')
biocLite('ctc')
biocLite('Biobase')
install.packages('gplots')
install.packages('ape')
q()

/nfs/apps/7/arch/generic/Trinity/2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix genes.counts.matrix --method voom --samples_file sample_rep_info.txt --min_reps_min_cpm 1,1 --output voom 2>&1 | tee -a DE.log

##Copy Trinity directory and modify "run_DE_analysis.pl" so it takes "0,0" for "--min_reps_min_cpm"
cp -r /nfs/apps/7/arch/generic/Trinity/2.4.0 .
cp /nfs/apps/7/arch/generic/Trinity/2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl .
cp run_DE_analysis.pl 2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl_
2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl_ \
--matrix genes.counts.matrix --method voom --samples_file sample_rep_info.txt --min_reps_min_cpm 0,0 --output voom_ 2>&1 | tee -a DE.log

#format result file for BETA input -p argument
#extract logFC/PValue
awk -F '[\t_]' 'NR>1 {print $1,$5,$7}' OFS='\t' genes.counts.matrix.HdF_vs_HdM.voom.DE_results > genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf
#tidy up non-DE genes & remove gene name duplicates (or can be removed by above script field separator '_')
sed -i -r 's/-0.278556176774529/0/;s/0.528696052790396/1/' genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf

#build annotation file for '-r' following the format of "refseqID	chrom	strand	start	end	gname2"
grep 'gene' ../danaus_plexippus_v3_core_32_85_1.gff | awk -F '[\t=;]' '{print $10,$1,$7,$4-1,$5,$10}' OFS='\t' > RefSeqGene.dp3.BETA.bed

#modify scripts under "/nfs/apps/7/arch/generic/build/BETA/BETA_1.0.7/BETA/ and/or?? /nfs/apps/7/arch/generic/build/BETA/BETA_1.0.7/build/lib/BETA
#fileformat_check.py on line23,24,91;
#PScore.py on line8 "chroms" value; 
#expr_combine.py on line259 "chrom" value;

#Update build
module load BETA/1.0.7
cd /nfs/apps/7/arch/generic/build/BETA/BETA_1.0.7
python setup.py clean
python setup.py build
python setup.py install --prefix=/nfs/apps/7/arch/generic/BETA/1.0.7

	BETA plus -p trim.HD-F-327-W200-G600-FDR0.05-island_c3.bed \
	-e genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf -k BSF --gs ../Dp_genome_v3.fasta -r RefSeqGene.dp3.BETA.bed -n F327 --df 0.05 -c 0.01 2>&1 | tee -a BETA.log

#make 3-column bed file for BETA -e argument
for i in 327 39 416 420; do
	mkdir $i
	cp ../SICER/$i/trim.HD-F-$i-vs-trim.HD-M-$i-W200-G600-E-union.island $i
	BETA plus -p $i/trim.HD-F-$i-vs-trim.HD-M-$i-W200-G600-E-union.island -e genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf \
	-k BSF --gs ../Dp_genome_v3.fasta -r RefSeqGene.dp3.BETA.bed -n $i-union_fdr0.01 -o $i --df 0.01 -c 0.01 2>&1 | tee -a BETA.log
	
	awk -F '[\t]' '{print $1,$2,$3}' OFS='\t' ../SICER/$i/trim.HD-F-$i-W200-G600-increased-islands-summary-FDR0.05 \
	> $i/trim.HD-F-$i-W200-G600-increased-islands-summary-FDR0.05
	BETA plus -p $i/trim.HD-F-"$i"-W200-G600-increased-islands-summary-FDR0.05_c3.bed -e genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf \
	-k BSF --gs ../Dp_genome_v3.fasta -r RefSeqGene.dp3.BETA.bed -n $i-up_fdr0.01 -o $i --df 0.01 -c 0.01 2>&1 | tee -a BETA.log
	
	awk -F '[\t]' '{print $1,$2,$3}' OFS='\t' ../SICER/$i/trim.HD-F-$i-W200-G600-decreased-islands-summary-FDR0.05 \
	> $i/trim.HD-F-$i-W200-G600-decreased-islands-summary-FDR0.05
	BETA plus -p $i/trim.HD-F-$i-W200-G600-decreased-islands-summary-FDR0.05_c3.bed -e genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf \
	-k BSF --gs ../Dp_genome_v3.fasta -r RefSeqGene.dp3.BETA.bed -n $i-de_fdr0.01 -o $i --df 0.01 -c 0.01 2>&1 | tee -a BETA.log
done

#####According to tests on parameters, use FDR (--df) 0.01, for -c no difference between 0.01 and 0.001 (default) 
###None detected when feeding "increased" data set for 327, even under most strict parameter
	
BETA plus -p ../SICER/H3K9me3/trim.HD-F-39-vs-trim.HD-M-39-W200-G600-E-union.island \
-e genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf -k BSF --gs ../Dp_genome_v3.fasta -r RefSeqGene.dp3.BETA.bed -n 39 --df 0.05 -c 0.01 2>&1 | tee -a BETA.log
awk '{print $1,$2,$3}' OFS='\t' ../SICER/H3K9me3/trim.HD-F-39-W200-G600-increased-islands-summary-FDR0.05 > trim.HD-F-39-W200-G600-increased-islands-summary-FDR0.05_c3.bed
awk '{print $1,$2,$3}' OFS='\t' ../SICER/H3K9me3/trim.HD-F-39-W200-G600-decreased-islands-summary-FDR0.05 > trim.HD-F-39-W200-G600-decreased-islands-summary-FDR0.05_c3.bed
BETA plus -p trim.HD-F-39-W200-G600-increased-islands-summary-FDR0.05_c3.bed \
-e genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf -k BSF --gs ../Dp_genome_v3.fasta -r RefSeqGene.dp3.BETA.bed -n 39+ --df 0.05 -c 0.01 2>&1 | tee -a BETA.log
BETA plus -p trim.HD-F-39-W200-G600-decreased-islands-summary-FDR0.05_c3.bed \
-e genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf -k BSF --gs ../Dp_genome_v3.fasta -r RefSeqGene.dp3.BETA.bed -n 39- --df 0.05 -c 0.01 2>&1 | tee -a BETA.log

#####Try Beta minus
BETA plus -p $i/trim.HD-F-$i-vs-trim.HD-M-$i-W200-G600-E-union.island -e genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf \
	-k BSF --gs ../Dp_genome_v3.fasta -r RefSeqGene.dp3.BETA.bed -n $i-union_fdr0.01 -o $i --df 0.01 -c 0.01 2>&1 | tee -a BETA.log

BETA minus -p $i/trim.HD-F-416-W200-G600-increased-islands-summary-FDR0.05_c3.bed -r RefSeqGene.dp3.BETA.bed -n $i-increased -d 1000
BETA minus -p $i/trim.HD-F-416-W200-G600-decreased-islands-summary-FDR0.05_c3.bed -r RefSeqGene.dp3.BETA.bed -n $i-decreased -d 1000
###Reduce distance reduced number of target genes (2k vs 13K) but also reduces score.

# Add "#" to header to result file so that can be viewed in IGV 
awk -F '[\t]' '{print $1,$2,$3,$5}' OFS='\t' BETA_OUTPUT/327-_downtarget.txt > 327-_downtarget_c4.txt
sort -k1,1 -k2,2n 327-_downtarget_c4.txt > 327-_downtarget_c4_sorted.txt

###Make Scaff-Gene-Linkage file
#Correct DPSCF300403
awk '/DPSCF300403/{$3="A";$4="NA"}1' OFS="\t" dp3_scaff_gene_linkage.bed > dp3_scaff_gene_linkage_corrected.bed 

join -1 1 -2 2 -t $'\t' <(sort ../BETA/RefSeqGene.dp3.BETA.bed) dp3_scaff_gene_linkage_corrected.bed | cut -f1-5,7-9 > RefSeqGene_Scaff_Cord_Linkage.bed

grep 'Z' RefSeqGene_Scaff_Cord_Linkage.bed > RefSeqGene_Scaff_Cord_Linkage_Z.bed
grep 'Z_anc' RefSeqGene_Scaff_Cord_Linkage.bed > RefSeqGene_Scaff_Cord_Linkage_Z_anc.bed
#624 RefSeqGene_Scaff_Cord_Linkage_Z_anc.bed

grep 'Z_neo' RefSeqGene_Scaff_Cord_Linkage.bed > RefSeqGene_Scaff_Cord_Linkage_Z_neo.bed
#472 RefSeqGene_Scaff_Cord_Linkage_Z_neo.bed

grep -w 'A' RefSeqGene_Scaff_Cord_Linkage.bed > RefSeqGene_Scaff_Cord_Linkage_A_all.bed
#14034 RefSeqGene_Scaff_Cord_Linkage_A_all.bed

#Autosomal only as assigned to an autosome
awk '$7=="A"&&$8!="NA"' RefSeqGene_Scaff_Cord_Linkage.bed > RefSeqGene_Scaff_Cord_Linkage_A.bed
12974 RefSeqGene_Scaff_Cord_Linkage_A.bed

#Total number of Scaffolds that have a gene
cut -f2,6-8 RefSeqGene_Scaff_Cord_Linkage.bed | sort -u > Scaff_with_Gene_List.txt
872

#Number of linkage assigned AS scaffolds
grep -wv '1\|NA' Scaff_with_Gene_List.txt | wc -l
441

###Make Scaff-Chimera-Break-Bed file
#Mark chimera scaffolds
awk 'BEGIN{OFS="\t"}/DPSCF300001/||/DPSCF300028/{$1=$1"-1";print;getline;$1=$1"-2";print;getline;$1=$1"-3";print;next}
/DPSCF300044/{$1=$1"-1";print;getline;$1=$1"-2";print;next}1' dp3_chimera_break.bed > dp3_chimera_break_marked.bed

#Match using chimera-marked scaffold then delete the chimera mark
grep -Ff <(cut -f6 RefSeqGene_Scaff_Cord_Linkage_A.bed | sort -u) dp3_chimera_break_marked.bed | sed 's/-.//g' > dp3_chimera_break_A.bed 
grep -Ff <(cut -f6 RefSeqGene_Scaff_Cord_Linkage_Z.bed | sort -u) dp3_chimera_break_marked.bed | sed 's/-.//g' > dp3_chimera_break_Z.bed 

###Get gene lists by linkage
for i in Z A; do grep -Ff <(cut -f1 ../BEDFILE/RefSeqGene_Scaff_Cord_Linkage_$i.bed) \
genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf > genes.counts.matrix.HdF_vs_HdM.voom.DE_results_$i.bsf; done

###Make RefSeq bed for BETA
awk '{$6=$1;NF=6;print}' OFS='\t' ../BEDFILE/RefSeqGene_Scaff_Cord_Linkage_Z.bed > RefSeqGene.dp3.BETA.Z.bed
awk '{$6=$1;NF=6;print}' OFS='\t' ../BEDFILE/RefSeqGene_Scaff_Cord_Linkage_A.bed > RefSeqGene.dp3.BETA.A.bed

for i in 327 39 416 420; do
	mkdir $i
	cp ../SICER/$i/trim.HD-F-$i-vs-trim.HD-M-$i-W200-G600-E-union.island $i
	BETA plus -p $i/trim.HD-F-$i-vs-trim.HD-M-$i-W200-G600-E-union.island -e genes.counts.matrix.HdF_vs_HdM.voom.DE_results.bsf \
	-k BSF --gs ../Dp_genome_v3.fasta -r RefSeqGene.dp3.BETA.bed -n $i-union_fdr0.01 -o $i --df 0.01 -c 0.01 2>&1 | tee -a BETA.log
done
	
module load BEDTools/2.26.0 BETA/1.0.7
for l in Z A; do 
	awk 'NF=3' OFS='\t' ../SICER/$i/trim.HD-F-$i-W200-G600-increased-islands-summary-FDR0.05 \
> $i/trim.HD-F-$i-W200-G600-increased-islands-summary-FDR0.05.bed

l=A
	intersectBed -a ../BEDFILE/dp3_chimera_break_$l.bed -b $i/trim.HD-F-$i-W200-G600-increased-islands-summary-FDR0.05.bed \
> $i/trim.HD-F-$i-W200-G600-increased-islands-summary-FDR0.05.$l.bed
	BETA basic -p $i/trim.HD-F-"$i"-W200-G600-increased-islands-summary-FDR0.05.$l.bed -e genes.counts.matrix.HdF_vs_HdM.voom.DE_results_$l.bsf \
-k BSF -r RefSeqGene.dp3.BETA.$l.bed -n $i-up_$l -o $i --df 0.05 -c 1 2>&1 | tee -a BETA_linkage_break.log

	awk 'NF=3' OFS='\t' ../SICER/$i/trim.HD-F-$i-W200-G600-decreased-islands-summary-FDR0.05 \
> $i/trim.HD-F-$i-W200-G600-decreased-islands-summary-FDR0.05.bed

l=A
	intersectBed -a ../BEDFILE/dp3_chimera_break_$l.bed -b $i/trim.HD-F-$i-W200-G600-decreased-islands-summary-FDR0.05.bed \
> $i/trim.HD-F-$i-W200-G600-decreased-islands-summary-FDR0.05.$l.bed
	BETA basic -p $i/trim.HD-F-"$i"-W200-G600-decreased-islands-summary-FDR0.05.$l.bed -e genes.counts.matrix.HdF_vs_HdM.voom.DE_results_$l.bsf \
-k BSF -r RefSeqGene.dp3.BETA.$l.bed -n $i-de_$l -o $i --df 0.01 -c 0.01 2>&1 | tee -a BETA_linkage_break.log	
done

##### GATK
# during Variant Calling
##### ERROR MESSAGE: VCFHeaderLine: ID cannot contain an equals sign
#samtools view -h $DIRNAME/split.bam | sed 's/HRSCAF=/HRSCAF/g' | samtools view -b -o $DIRNAME/split_.bam -@ 16
#samtools index $DIRNAME/split_.bam -@ 16
sed 's/HRSCAF=/HRSCAF/g' ../Danaus_HiRise/Dapl_Zhan_v3_HiC.RN.fasta > Dapl_Zhan_v3_HiC.RN.fasta
REF=Dapl_Zhan_v3_HiC.RN.fasta
#gatk-launch CreateSequenceDictionary -R ref.fasta
java -jar /nfs/apps/7/arch/generic/picard/2.8.1-Java-1.8.0_102/picard.jar CreateSequenceDictionary R=$REF O=$(basename $REF .fasta).dict
samtools faidx $REF

for i in STAR_out/Hd?[1-3]/Aligned.sortedByCoord.out.bam; do
DIRNAME=$(dirname ${i/STAR_out/GATK})
mkdir -p $DIRNAME
samtools view -h $i | sed 's/HRSCAF=/HRSCAF/g' | samtools view -b -o ${i/STAR_out/GATK} -@ 16
# 2. Add read groups, sort, mark duplicates, and create index
# change "sample" name, important!!!
java -jar /nfs/apps/7/arch/generic/picard/2.8.1-Java-1.8.0_102/picard.jar AddOrReplaceReadGroups I=${i/STAR_out/GATK} O=$DIRNAME/rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=$(basename $DIRNAME)
java -jar /nfs/apps/7/arch/generic/picard/2.8.1-Java-1.8.0_102/picard.jar MarkDuplicates I=$DIRNAME/rg_added_sorted.bam O=$DIRNAME/dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics
# 3. Split'N'Trim and reassign mapping qualities
java -jar /nfs/apps/7/arch/generic/GATK/3.7-Java-1.8.0_102/GenomeAnalysisTK.jar -T SplitNCigarReads -R $REF -I $DIRNAME/dedupped.bam -o $DIRNAME/split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
##### ERROR MESSAGE: VCFHeaderLine: ID cannot contain an equals sign
# 6. Variant calling
# java -jar /nfs/apps/7/arch/generic/GATK/3.7-Java-1.8.0_102/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF -I $DIRNAME/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $DIRNAME/pass1_raw.vcf
java -jar /nfs/apps/7/arch/generic/GATK/3.7-Java-1.8.0_102/GenomeAnalysisTK.jar -T HaplotypeCaller -ERC GVCF -nct 16 -R $REF -I $DIRNAME/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $DIRNAME/pass1_raw.g.vcf
# need to change RG name since first run did not do it
SA=$(basename $DIRNAME)
sed -i "s/FORMAT\tsample/FORMAT\t$SA/" $DIRNAME/pass1_raw.g.vcf
done
# -nt 16 flag does not work; --javaOptions "-Xmx50g" work for GATK4 not 3.
java -Xmx90g -jar /nfs/apps/7/arch/generic/GATK/3.7-Java-1.8.0_102/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $REF \
    --variant_index_type LINEAR --variant_index_parameter 128000 \
    $(echo GATK/*/pass1_raw.g.vcf  | sed 's/GATK/-V GATK/g') -o GATK/HdFM123.vcf
# -G StandardAnnotation -G AS_StandardAnnotation

# Now run again to generate F|M vcf
for i in GATK/HdF?/pass1_raw.g.vcf; do
sed 's/HdF./HdF/' $i > ${i/pass/HdF.pass}
done
for i in GATK/HdM?/pass1_raw.g.vcf; do
sed 's/HdM./HdM/' $i > ${i/pass/HdM.pass}
done
java -Xmx90g -jar /nfs/apps/7/arch/generic/GATK/3.7-Java-1.8.0_102/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $REF \
    --variant_index_type LINEAR --variant_index_parameter 128000 \
    $(echo GATK/*/Hd?.pass1_raw.g.vcf  | sed 's/GATK/-V GATK/g') -o GATK/HdFM.vcf

for i in GATK/HdFM.vcf GATK/HdFM123.vcf; do
# 7. Variant filtering
java -jar /nfs/apps/7/arch/generic/GATK/3.7-Java-1.8.0_102/GenomeAnalysisTK.jar -T VariantFiltration -R $REF \
    -V $i -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${i/.vcf/_filtered.vcf}
# 8. Filter out indels and scaffolds
vcftools --vcf ${i/vcf/_filtered.vcf} --remove-indels --remove-filtered-all --recode --recode-INFO-all \
    $(echo chr{1..30} | sed 's/chr/--chr chr/g') --stdout | sed '/contig=<ID=ScaoKnI/d' | gzip -c > ${i/.vcf/}_SNPs_chr.vcf.gz
done

# 5. Base Recalibration: use 1 pass to do calibration
# java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $DIRNAME/split_.bam -L 20 -knownSites dbsnp.vcf -knownSites gold_indels.vcf -o recal_data.table

grep -vc "#" GATK/HdFM123.AS_filtered.vcf
174088
# Number filtered
grep -v "#" GATK/HdFM123.AS_filtered.vcf | grep -c "PASS"
126645
# After filtering, kept 122619 out of a possible 174088 Sites

# Filter out indels and scaffolds
vcftools --vcf GATK/HdFM123.AS_filtered.vcf --remove-indels --remove-filtered-all --recode --recode-INFO-all \
    $(echo chr{1..30} | sed 's/chr/--chr chr/g') --out GATK/HdFM123.AS_SNPs_chr
# After filtering, kept 122619 out of a possible 174088 Sites

# Calculate total transcript length per chr
$ Danaus_HiRise

awk '/exon/{print $1,$4,$5}' Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.exon.gtf
grep -w exon Zhan_v3_HiC.deeptools.gtf | awk 'BEGIN{getline;i=$1;a=$4;b=$5}/$1==i/{getline;a+=$4;b+=$5}{print i,b-a;i=$1;a=$4;b=$5;next}END{print i,b-a}'

### SRA FTP upload; male: ?'; female: 9-11-2019
$ /rfs/jwalters/EXPORT/Monarch_RNA_ChIP
lftp -u subftp,w4pYB9VQ -e 'cd uploads/lg356@cornell.edu_RfJE9S4l;mkdir Dp_RNA_male_tissues;mput -P 3 *gz' ftp://ftp-private.ncbi.nlm.nih.gov
lftp -u subftp,w4pYB9VQ -e 'cd uploads/lg356@cornell.edu_RfJE9S4l;mkdir Dp_Chicago_HiC;cd Dp_Chicago_HiC;mput -P 4 *gz' ftp://ftp-private.ncbi.nlm.nih.gov

### TMM
/nfs/apps/7/arch/generic/Trinity/2.4.0/util/abundance_estimates_to_matrix.pl --est_method eXpress eXpress_out/*/results.xprs \
    --name_sample_by_basedir --out_prefix Dapl_HD

method=edgeR
# self calculated "count" file not correct! Should use "effective count"
matrix_file=Dapl_xprs.hd_counts.matrix
samples_file=Dapl_hd_rep.txt
outdir=edgeR_out

module load R/3.4.2
/nfs/apps/7/arch/generic/Trinity/2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --method $method \
    --matrix $matrix_file --samples_file $samples_file --output $outdir --min_reps_min_cpm 2,1




