#!/bin/sh

srun -p intel -N 1 -n 1 -c 16 -t 72:00:00 --mem 100G --pty /bin/bash

##### STAR temporary directory must be absent at time of running
module load STAR/2.5.2b
ts=$(date +"%Y%m%d")

mkdir STAR_ref.Dapl
STAR --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_ref.Dapl \
    --genomeFastaFiles ../Danaus_HiRise/Dapl_Zhan_v3_HiC.RN.fasta 2>&1 | tee STAR.$ts.log

for i in Reads_trim/*R1.P*
do
j=$(basename $i)
outDir=STAR_out/${j%%_*}
mkdir -p $outDir
STAR --runMode alignReads --runThreadN 16 --genomeDir STAR_ref.Dapl \
    --readFilesIn $i ${i/R1/R2} --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 16 \
    --outFileNamePrefix $outDir/ --outTmpDir /tmp/l338g110 \
    --outFilterMultimapNmax 10 --alignMatesGapMax 700 2>&1 | tee -a $outDir/STAR.$ts.log
done

multijoin1 () {
f1=$1; f2=$2; shift 2
if [ $# -gt 0 ]; then
    join -t "|" $f1 $f2 | multijoin1 - $@
else
    join -t "|" $f1 $f2
fi
}

multijoin1 STAR_out/*/Log.final.out | awk 'BEGIN{FS="|\t ";OFS="|\t ";printf "%+47s |\t", "Sample name";print "F327" ,"F39" ,"F416" ,"F420" ,"M327" ,"M39" ,"M416" ,"M420" ,"M"}NR>3{$1=$1;print}' > STAR_out/Log.final.out.matrix

### Merge BAMs M samples: M327 M39 M416 M420; plus F420
# samtools merge greatly reduce the size of two files combined!
$ Danaus_ChIP
mkdir bamFile
#for j in 327 39 416 420 ""; do samtools merge BAM/M$j.star.merge.bam STAR_out/M$j/*bam ../ChIP/STAR_out/M$j/*bam; done
for i in STAR_out/M?* STAR_out/F420; do
    samtools merge ${i/STAR_out/bamFile}.star.merge.bam $i/*bam ../ChIP/$i/*bam
done

# for F327 F39 F416 plus M, use 2018; for F, use 2016
for i in STAR_out/F* STAR_out/M; do
    mv $i/*bam ${i/STAR_out/bamFile}.star.2018.bam
done

mv ../ChIP/STAR_out/F/*bam bamFile/F.star.2016.bam

for i in bamFile/*; do samtools index $i; done


##### calculate reads mapping ratios
for
for s in F M; do
    echo $s
    printf Z-neo"\t"
    samtools view bamFile/$s$m.star*bam chr1:0-5685560 -c -@ 16
    printf Z-anc"\t"
    samtools view bamFile/$s$m.star*bam chr1:5685561-15616146 -c -@ 16
for i in {2..30}; do
    printf chr$i"\t"
    samtools view bamFile/$s$m.star*bam chr$i -c -@ 16
done; done > $m.mapped.reads.count.per.chr.txt


for s in M; do
echo $s
printf Z-neo"\t"
samtools view ../ChIP/STAR_out/M/Aligned.sortedByCoord.out.bam ScaoKnI_86;HRSCAF=140:0-5685560 -c -@ 16
printf Z-anc"\t"
samtools view ../ChIP/STAR_out/M/Aligned.sortedByCoord.out.bam ScaoKnI_86;HRSCAF=140:5685561-15616146 -c -@ 16
for i in {2..30}; do
printf chr$i"\t"
samtools view ../ChIP/STAR_out/M/Aligned.sortedByCoord.out.bam chr$i -c -@ 16
done; done > M.2016.mapped.reads.count.per.chr.txt

for bam in *bam; do
echo $bam
samtools sort $bam -@16 -o ${bam/bam/sorted.bam}
samtools index $bam -@16
printf Z-neo"\t"
samtools view $bam chr1:0-5685560 -c -@ 16
printf Z-anc"\t"
samtools view $bam chr1:5685561-15616146 -c -@ 16
printf chr$i"\t"
samtools view $bam chr$i -c -@ 16
done; done > M.2016.mapped.reads.count.per.chr.txt

##### remove duplicates: silly, 'cause bamCompare remove duplicates automatically!
## Plus, results almost same!
mkdir bamFile_rmdup
for i in bamFile/*416*bam; do
samtools view $i -M -u -@ 16 | samtools sort -@ 16 -o ${i/bamFile/bamFile_rmdup}
done

samtools rmdup $i ${i/bamFile/bamFile_rmdup} 2>&1 | tee -a bamFile_rmdup/rmdup.log

##### Coverage analysis
genomeCoverageBed -ibam F327.star.2018.bam -bga > F327.star.2018.cov.bedgraph


##### deepTools
module load SAMtools/1.6
# for i in STAR_out/*/*bam; do samtools index $i; done

module load deeptools/2.5.6-Python-2.7.12
mkdir deepTools
cp ../Danaus_ATACseq/deepTools_Dapl/bl.chr1.chr30.bed deepTools

### plotCoverage'
dir=plotCoverage
mkdir $dir
j=416
for s in F M; do
plotCoverage -b ../bamFile//$s$j*bam -o $dir/plotCoverage_$s$j.Zn.pdf -p 16 \
    -l $(echo STAR_out/$s*/*bam | sed 's/STAR_out\/\([FM0-9]\+\)\/[ABCa-z.]\+/\1/g') \
    -T "Sample Coverage of "$s" Samples" --ignoreDuplicates --minMappingQuality 255 \
    --skipZeros --outRawCounts $dir/plotCoverage_$s$j.Zn.txt -r chr1:1:5685560 2>&1 | tee $dir/plotCoverage_$s$j.Zn.log
done


### multiBamSummary
$ deepTools
dirIn=../bamFile
dirOut=multiBamSummary
mkdir $dirOut

for i in {1..30}; do
multiBamSummary bins --bamfiles $dirIn/[F,M]416*bam $dirIn/[F,M].*bam \
    --labels H4K16ac-F H4K16ac-M input-F input-M \
    -out $dirOut/FM416_input.chr$i.npz --outRawCounts $dirOut/FM416_input.chr$i.tab -p 8 --region $i
plotCorrelation -in $dirOut/FM416_input.chr$i.npz \
    --corMethod pearson --skipZeros --removeOutliers \
    --plotTitle "Pearson Correlation of H4K16ac ChIP and Input samples: chr$i" \
    --whatToPlot scatterplot \
    -o $dirOut/scatterplot_PearsonCorr_FM416_input.chr$i.pdf \
    --outFileCorMatrix $dirOut/PearsonCorr_FM416_input.chr$i.tab
plotCorrelation -in $dirOut/FM416_input.chr$i.npz \
    --corMethod pearson --skipZeros --removeOutliers \
    --plotTitle "Pearson Correlation of H4K16ac ChIP and Input samples" \
    --whatToPlot heatmap --plotNumbers \
    -o $dirOut/heatmap_PearsonCorr_FM416_input.chr$i.pdf
done

--minMappingQuality 30 \

plotCorrelation \
-in readCounts.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o heatmap_SpearmanCorr_readCounts.png   \
--outFileCorMatrix SpearmanCorr_readCounts.tab

# --labels uses labels defined i n multiBamSummary

### bamCompare
# bamCompare to get scale factors using chrA
cd deepTools
ts=$(date +"%Y%m%d")
dirIn=../bamFile
dirOut=bamCompare
mkdir $dirOut

# little difference with Scale Factor calculated with using and not using bl50
# for female samples, follow two steps: first calculate scale factor, then using SF to make bigwig
### re-run w/o --skipNAs flag
for j in 327 39 416 420; do
echo "bamCompare -b1 $dirIn/F$j.*bam -b2 $dirIn/F.*bam \
    -o $dirOut/F$j.SES.A.blZW.bigwig -of bigwig  \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactorsMethod SES --ratio log2 \
    -bs 50 --smoothLength 150 -p 16 -e --centerReads \
    --maxFragmentLength 700 -bl bl.chr1.chr30.bed -v \
    &>> $dirOut/F.A.blZW.$ts.log " | tee -a $dirOut/F.A.blZW.$ts.log | bash
done

grep "scale" $dirOut/F.A.blZW.$ts.log > $dirOut/F.A.blZW.scaleFactor.txt

# w/o --skipNAs does not change scale factors
bamCompare -b1 $dirIn/F327.*bam -b2 $dirIn/F.*bam \
    -o $dirOut/F327.SES.bigwig -of bigwig \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactors 0.43377435:1 --ratio log2 \
    -bs 50 --smoothLength 150 -p 16 -e --centerReads \
    --maxFragmentLength 700
bamCompare -b1 $dirIn/F39.*bam -b2 $dirIn/F.*bam \
    -o $dirOut/F39.SES.bigwig -of bigwig \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactors 1:0.65761043 --ratio log2 \
    -bs 50 --smoothLength 150 -p 16 -e --centerReads \
    --maxFragmentLength 700
bamCompare -b1 $dirIn/F416.*bam -b2 $dirIn/F.*bam \
    -o $dirOut/F416.SES.bigwig -of bigwig \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactors 0.37349834:1 --ratio log2 \
    -bs 50 --smoothLength 150 -p 16 -e --centerReads \
    --maxFragmentLength 700
bamCompare -b1 $dirIn/F420.*bam -b2 $dirIn/F.*bam \
    -o $dirOut/F420.SES.bigwig -of bigwig \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactors 1:0.28393517 --ratio log2 \
    -bs 50 --smoothLength 150 -p 16 -e --centerReads \
    --maxFragmentLength 700

###! can be done with --ignoreForNormalization chr1 chr30 flag!

# for male samples, tested that including vs not including Z does not change scale factor much
### re-run w/o --skipNAs flag
for j in 327 39 416 420; do
echo "bamCompare -b1 $dirIn/M$j.*bam -b2 $dirIn/M.*bam \
    -o $dirOut/M$j.SES.bigwig -of bigwig \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactorsMethod SES --ratio log2 \
    -bs 50 --smoothLength 150 -p 16 -e --centerReads \
    --maxFragmentLength 700 -v &>> $dirOut/M.$ts.log " | tee -a $dirOut/M.$ts.log | bash
done

# -bs 10 for TSS.b500
bamCompare -b1 $dirIn/F416.*bam -b2 $dirIn/F.*bam \
    -o $dirOut/F416.SES.bs10.bigwig -of bigwig \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactors 0.37349834:1 --ratio log2 -bs 10 \
    -p 16 -e --centerReads --maxFragmentLength 700
bamCompare -b1 $dirIn/M$j.*bam -b2 $dirIn/M.*bam \
    -o $dirOut/M$j.SES.bs10.bigwig -of bigwig  \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactorsMethod SES --ratio log2 -bs 10 \
    -p 16 -e --centerReads --maxFragmentLength 700

# use ratio instead of log2
for s in F M; do
bamCompare -b1 $dirIn/"$s"416.*bam -b2 $dirIn/$s.*bam \
    -o $dirOut/"$s"416.ratio.bigwig -of bigwig \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactorsMethod SES --ratio ratio \
    -bs 50 --smoothLength 150 -p 16 -e --centerReads \
    --maxFragmentLength 700
done

# bedgraph for plotting, increase bin/smoothlength 10x
# -bs 50 useless, too large even to load in IGV
j=416
j=420
bamCompare -b1 $dirIn/F416.*bam -b2 $dirIn/F.*bam \
    -o $dirOut/F416.SES.bs500.bedgraph -of bedgraph \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactors 0.37349834:1 --ratio log2 \
    -bs 500 --smoothLength 1500 -p 16 -e --centerReads \
    --maxFragmentLength 700
bamCompare -b1 $dirIn/M$j.*bam -b2 $dirIn/M.*bam \
    -o $dirOut/M$j.SES.bs500.bedgraph -of bedgraph  \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactorsMethod SES --ratio log2 \
    -bs 500 --smoothLength 1500 -p 16 -e --centerReads \
    --maxFragmentLength 700


# bedgraph file contains scaffolds and chr names need to be chanaged.

find -name ?416.SES*bedgraph | xargs sed -i -e '/chr/!d;s/chr/dp/'

bedtools makewindows -g ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chr.seq_length.txt -w 10000 > Dapl_Zhan_v3_HiC.chr.10k.bed
# require bed4 format, fourth column can't be duplicate
awk '{$4=NR}1' Dapl_Zhan_v3_HiC.chr.10k.bed > Dapl_Zhan_v3_HiC.chr.10k.bed4

module load kentUtils
m=SES
#m=ratio
j=416
j=420
for s in F M; do
bigWigAverageOverBed -bedOut="$s"$j.$m.bs50.10k.bed bamCompare/"$s"$j.$m.bigwig Dapl_Zhan_v3_HiC.chr.10k.bed4 "$s"$j.$m.bs50.mean.tab
cut -f1-3,5 "$s"$j.$m.bs50.10k.bed > "$s"$j.$m.bs50.10k.bedgraph
rm "$s"$j.$m.bs50.10k.bed
sed -i 's/chr/dp/' "$s"$j.$m.bs50.10k.bedgraph
done

cp ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chr.bed Dapl_Zhan_v3_HiC.chr.bed
awk '{$4=NR}1' Dapl_Zhan_v3_HiC.chr.bed > Dapl_Zhan_v3_HiC.chr.bed4
m=ratio
for s in F M; do
bigWigAverageOverBed -bedOut="$s"$j.$m.bs50.chr.bed bamCompare/"$s"$j.$m.bs50.bigwig Dapl_Zhan_v3_HiC.chr.bed4 "$s"$j.$m.bs50.mean.tab
done

# F
chr1    0       5685560 1.6911
chr1    5685561 15616146 1.49697
# A mean
awk 'NR>2{a+=$5}END{print a/(NR-2)}' F416.$m.bs50.chr.bed
1.29536

echo "1.6911/1.29536" | bc -l
1.30550580533596837944
echo "1.49697/1.29536" | bc -l
1.15564013092885375494

# M
chr1    0       5685560 1.17851
chr1    5685561 15616146 1.12804
awk 'NR>2{a+=$5}END{print a/(NR-2)}' M416.$m.bs50.chr.bed
1.20352

echo "1.17851/1.20352" | bc -l
.97921929008242488699
echo "1.12804/1.20352" | bc -l
.93728396703004520074

### bamCoverage to get absolute value tracks
# scaling between inputs
ts=$(date +"%Y%m%d")
dirIn=../bamFile
dirOut=bamCompare
echo "bamCompare -b1 $dirIn/F.*bam -b2 $dirIn/M.*bam \
    -o $dirOut/FM.SES.bigwig -of bigwig \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactorsMethod SES --ratio log2 \
    -bs 50 --smoothLength 150 -p 16 -e --centerReads \
    --maxFragmentLength 700 --ignoreForNormalization chr1 chr30 \
    -v &>> $dirOut/FM.input.$ts.log " | tee $dirOut/FM.input.$ts.log | bash


### plotProfile FM all marks in one
cp ../../Danaus_ATACseq/*/*gtf .

mkdir computeMatrix plotProfile
### re-run w/o --skipNAs flag, add -bs 50
# metagene
for j in 327 39 416 420; do
computeMatrix scale-regions -S bamCompare/?$j.SES.bigwig \
    -R Zhan_v3_HiC.Z-neo.gtf Zhan_v3_HiC.Z-anc.gtf Zhan_v3_HiC.A.gtf  \
    -out computeMatrix/FM.$j.SES.metagene.gz --metagene --skipZeros \
    -b 3000 -m 5000 -a 3000 -p 16 -bs 50
plotProfile -m computeMatrix/FM.$j.SES.metagene.gz -out plotProfile/FM.$j.SES.metagene.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel F_$j M_$j --plotType se -T metagene \
    --colors olivedrab darkorange black --yAxisLabel "log2(ChIP:Input)"
# intergenic
computeMatrix scale-regions -S bamCompare/?$j.SES.bigwig -R Zhan_v3_HiC_intergenic.bed6 \
    -out computeMatrix/FM.$j.SES.intergenic.gz --skipZeros -p 16 -bs 50
plotProfile -m computeMatrix/FM.$j.SES.intergenic.gz -out plotProfile/FM.$j.SES.intergenic.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel F_$j M_$j --plotType se -T intergenic \
    --colors olivedrab darkorange black --yAxisLabel "log2(ChIP:Input)" --startLabel TES --endLabel TSS
# intron
computeMatrix scale-regions -S bamCompare/?$j.SES.bigwig -R Zhan_v3_HiC_intron.bed6 \
    -out computeMatrix/FM.$j.SES.intron.gz --skipZeros -p 16 -bs 50
plotProfile -m computeMatrix/FM.$j.SES.intron.gz -out plotProfile/FM.$j.SES.intron.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel F_$j M_$j --plotType se -T intron \
    --colors olivedrab darkorange black --yAxisLabel "log2(ChIP:Input)" --startLabel "" --endLabel ""
done


# by group
# metagene
for j in 327 39 416 420; do
plotProfile -m computeMatrix/FM.$j.SES.metagene.gz -out plotProfile/FM.$j.SES.metagene.bySex.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel F M --plotType se \
    --colors deeppink blue --yAxisLabel "log2(ChIP:Input)" --numPlotsPerRow 1 --perGroup -T ""
# intergenic
plotProfile -m computeMatrix/FM.$j.SES.intergenic.gz -out plotProfile/FM.$j.SES.intergenic.bySex.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel F M --plotType se --startLabel TES --endLabel TSS \
	--colors deeppink blue --yAxisLabel "log2(ChIP:Input)" --numPlotsPerRow 1 --perGroup -T ""
# intron
plotProfile -m computeMatrix/FM.$j.SES.intron.gz -out plotProfile/FM.$j.SES.intron.bySex.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel F M --plotType se --startLabel "" --endLabel "" \
	--colors deeppink blue --yAxisLabel "log2(ChIP:Input)" --numPlotsPerRow 1 --perGroup -T ""
done


# Z only unsorted
# computeMatrixOperations sort -m
j=416
computeMatrix scale-regions -S bamCompare/?$j.SES.bigwig \
    -R Zhan_v3_HiC.Z-neo.gtf Zhan_v3_HiC.Z-anc.gtf --sortRegions keep \
    -out computeMatrix/FM.$j.SES.metagene.Zs.gz --metagene --skipZeros \
    -b 3000 -m 5000 -a 3000 -p 16 -bs 50
# plotProfile -m computeMatrix/FM.$j.SES.metagene.Zs.gz -out plotProfile/FM.$j.SES.metagene.Zs.pdf \
    --regionsLabel Z-neo Z-anc --samplesLabel F_$j M_$j --plotType se -T metagene \
    --colors olivedrab darkorange --yAxisLabel "log2(ChIP:Input)"
basename=FM.$j.SES.metagene.Zs
basename=FM.$j.SES.metagene.Zs.TSS
computeMatrix reference-point -S bamCompare/?$j.SES.bigwig \
    -R Zhan_v3_HiC.Z-neo.gtf Zhan_v3_HiC.Z-anc.gtf \
    -out computeMatrix/$basename.gz --metagene --skipZeros \
    --referencePoint TSS -b 500 -a 0 --nanAfterEnd -p 16 -bs 10
plotHeatmap -m computeMatrix/$basename.gz -out plotHeatmap/$basename.pdf \
	--outFileSortedRegions plotHeatmap/$basename.bed \
    --outFileNameMatrix plotHeatmap/$basename.tab \
    --regionsLabel "neo-Z" "anc-Z" --samplesLabel F$j M$j -T "" -T "" \
    --yAxisLabel "log2(ChIP:Input)" --xAxisLabel "gene distance" --sortRegions no \
    --boxAroundHeatmaps no --heatmapHeight 14 --heatmapWidth 6 \
    --colorList black,gold black,gold --alpha 0.8

# Get Z gene list ordered in coordinates
awk -F '[\t"]' '$1=="chr1"&&$3=="gene"{print $10}' Zhan_v3_HiC.deeptools.gtf > Z.gene.order.list.txt


### reference-point
# metagene
# average lengths of intron and intergenic regions
awk '{l+=$3-$2}END{print l/NR}' Zhan_v3_HiC_intergenic.bed6
26644.5
awk '{l+=$3-$2}END{print l/NR}' Zhan_v3_HiC_intron.bed6
799.605

for j in 327 39 416 420; do
computeMatrix reference-point -S bamCompare/?$j.SES.bigwig \
    -R Zhan_v3_HiC.Z-neo.gtf Zhan_v3_HiC.Z-anc.gtf Zhan_v3_HiC.A.gtf  \
    -out computeMatrix/FM.$j.SES.metagene.TSS.gz --metagene --skipZeros \
    --referencePoint TSS -b 500 -a 500 --nanAfterEnd -p 16 -bs 10
plotProfile -m computeMatrix/FM.$j.SES.metagene.TSS.gz -out plotProfile/FM.$j.SES.metagene.TSS.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel Female Male --plotType se -T metagene \
    --colors olivedrab darkorange black --yAxisLabel "log2(ChIP:Input)" \
    --refPointLabel TSS
basename=FM.$j.SES.metagene.TSS
plotHeatmap -m computeMatrix/$basename.gz -out plotHeatmap/$basename.pdf \
	--boxAroundHeatmaps no --samplesLabel Female Male --xAxisLabel "" \
	--yAxisLabel "log2(ChIP:Input)" --regionsLabel Z-neo Z-anc A -T "" \
	--colorList mediumblue,lightyellow --alpha 0.8 --heatmapWidth 8 \
    --whatToShow "heatmap and colorbar" --refPointLabel TSS
# intergenic: !can't have --nanAfterEnd
computeMatrix reference-point -S bamCompare/?$j.SES.bigwig -R Zhan_v3_HiC_intergenic.bed6 \
    -out computeMatrix/FM.$j.SES.intergenic.center.1.5k.gz --skipZeros \
    --referencePoint center -b 1500 -a 1500 -p 16 -bs 50
plotProfile -m computeMatrix/FM.$j.SES.intergenic.center.1.5k.gz -out plotProfile/FM.$j.SES.intergenic.center.1.5k.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel Female Male --plotType se -T intergenic \
    --colors olivedrab darkorange black --yAxisLabel "log2(ChIP:Input)" \
    --refPointLabel center
# intron
computeMatrix reference-point center -S bamCompare/?$j.SES.bigwig -R Zhan_v3_HiC_intron.bed6 \
    -out computeMatrix/FM.$j.SES.intron.center.gz --skipZeros \
    --referencePoint center -b 300 -a 300 --nanAfterEnd -p 16 -bs 10
plotProfile -m computeMatrix/FM.$j.SES.intron.center.gz -out plotProfile/FM.$j.SES.intron.center.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel Female Male --plotType se -T intron \
    --colors olivedrab darkorange black --yAxisLabel "log2(ChIP:Input)" \
    --refPointLabel center
done

### Calculate 500bp upstream of TSS and get data table
j=416
basename=FM$j.TSS.b500
# little difference with keeping Zeros, only slightly increase autosomal genes
basename=FM$j.TSS.b500.keepZeros
basename=FM$j.bs10.TSS.b500
computeMatrix reference-point -S bamCompare/?$j.SES.bs10.bigwig \
    -R Zhan_v3_HiC.Z-neo.gtf Zhan_v3_HiC.Z-anc.gtf Zhan_v3_HiC.A.gtf  \
    -out computeMatrix/$basename.gz --metagene --skipZeros \
    --referencePoint TSS -b 500 -a 0 --nanAfterEnd -p 16 -bs 10
plotProfile -m computeMatrix/$basename.gz -out plotProfile/$basename.pdf \
    --outFileSortedRegions plotProfile/$basename.bed \
    --outFileNameData plotProfile/$basename.tab \
    --regionsLabel Z-neo Z-anc A --samplesLabel Female Male --plotType se -T metagene \
    --colors olivedrab darkorange black --yAxisLabel "log2(ChIP:Input)" \
    --refPointLabel ""
plotHeatmap -m computeMatrix/$basename.gz -out plotHeatmap/$basename.pdf \
	--outFileSortedRegions plotHeatmap/$basename.bed \
    --outFileNameMatrix plotHeatmap/$basename.tab \
    --boxAroundHeatmaps no --samplesLabel Female Male --xAxisLabel "" \
	--yAxisLabel "log2(ChIP:Input)" --regionsLabel Z-neo Z-anc A -T "" \
	--colorList mediumblue,lightyellow --alpha 0.8 --heatmapWidth 8 \
    --whatToShow "heatmap and colorbar" --refPointLabel TSS

# in plotHeatmap/*.tab file, 3rd line 103, 4th and forward all 100, corresponding to 50 bins for each sex
Z-neo:473       Z-anc:597       A:13773 Female ...
# get gene id from plotHeatmap/*.bed file 4th field> plotHeatmap/FM416.TSS.b500.IDs.txt
paste <(awk 'NR>1{sub(/-TA/,"");print $4}' plotHeatmap/FM416.TSS.b500.bed) <(awk 'NR>3' plotHeatmap/FM416.TSS.b500.tab) > plotHeatmap/FM416.TSS.b500.matrix
# total 14843 genes


# gtf by chr
mkdir chrA_gtf
for i in {2..30}; do
	grep -w chr$i Zhan_v3_HiC.A.gtf > chrA_gtf/Zhan_v3_HiC.chr$i.gtf
done

#for j in 327 39 416 420; do
basename=FM$j.SES.chr
basename=FM$j.SES.bs50.chr
computeMatrix scale-regions -S bamCompare/[F,M]$j.SES.bigwig \
    -R Zhan_v3_HiC.Z-neo.gtf Zhan_v3_HiC.Z-anc.gtf GTF_chrA/Zhan_v3_HiC.chr{2..30}.gtf \
    -out computeMatrix/$basename.gz --metagene --skipZeros \
    -b 3000 -m 5000 -a 3000 -p 16 -bs 50
# all in one
plotProfile -m computeMatrix/$basename.gz -out plotProfile/$basename.pdf \
    --samplesLabel Female Male -T "" --plotType se --yAxisLabel "log2(ChIP:Input)" \
    --regionsLabel Z-neo Z-anc chr{2..30} --outFileNameData plotProfile/$basename.tab
sed -i '1,2d' plotProfile/$basename.tab

# by chr
plotProfile -m computeMatrix/$basename.gz -out plotProfile/$basename.perGroup.pdf \
    --samplesLabel Female Male -T "" --plotType se --yAxisLabel "log2(ChIP:Input)" \
    --regionsLabel Z-neo Z-anc chr{2..30} \
    --perGroup --colors deeppink blue --numPlotsPerRow 4


### Exclude chr30
grep -v chr30 Zhan_v3_HiC.A.gtf > Zhan_v3_HiC.A.xchr30.gtf

j=416
basename=FM$j.SES.xchr30
computeMatrix scale-regions -S bamCompare/[F,M]$j.SES.bigwig \
    -R Zhan_v3_HiC.Z-neo.gtf Zhan_v3_HiC.Z-anc.gtf Zhan_v3_HiC.A.xchr30.gtf  \
    -out computeMatrix/$basename.gz --metagene --skipZeros \
    -b 3000 -m 5000 -a 3000 -p 16
plotProfile -m computeMatrix/$basename.gz -out plotProfile/$basename.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel F$j M$j -T metagene \
    --colors olivedrab darkorange black --plotType se \
    --yAxisLabel "log2(sample:input)"


### plotHeatmap sorted by expression
mkdir plotHeatmap
# expression sorted by F or M, plotting all marks in each: Not good!
for s in F M; do
basename=expr$s
computeMatrix scale-regions -S bamCompare/$s*.SES.bigwig \
	-R Zhan_v3_HiC.expr$s.Z-neo.gtf Zhan_v3_HiC.expr$s.Z-anc.gtf Zhan_v3_HiC.expr$s.A.gtf \
	-out computeMatrix/$basename.gz --metagene -b 3000 -m 5000 -a 3000 -p 16
plotHeatmap -m computeMatrix/$basename.gz -out plotHeatmap/$basename.pdf \
	--boxAroundHeatmaps no --sortRegions no --samplesLabel H3K27me3 H3k9me3 H4k16ac H4k20me1 \
	--yAxisLabel "log2(sample:input)" --regionsLabel Z-neo Z-anc A -T $s" sorted by fpkm" \
	--xAxisLabel "gene distance" --colorList mediumblue,lightyellow --alpha 0.8 --heatmapWidth 8 \
done

# expression sorted by F or M, plotting signal of both in one according to each.
#for j in 327 39 416 420; do
j=416
for s in F M; do
basename=expr.sorted.$s$j
computeMatrix scale-regions -S bamCompare/?$j.SES.bigwig \
	-R Zhan_v3_HiC.deeptools.expr$s.gtf \
	-out computeMatrix/$basename.gz --metagene -b 3000 -m 5000 -a 3000 -p 16 --skipZeros
plotHeatmap -m computeMatrix/$basename.gz -out plotHeatmap/$basename.pdf \
	--boxAroundHeatmaps no --sortRegions no --samplesLabel "F$j sorted by $s" "M$j sorted by $s" \
	--yAxisLabel "log2(sample:input)" --regionsLabel Z-neo Z-anc A -T "" \
	--xAxisLabel "gene distance" --colorList purple,lightyellow mediumblue,lightyellow \
	--alpha 0.8 --heatmapWidth 8
done; done

	--outFileNameMatrix computeMatrix/$basename.tab --outFileSortedRegions computeMatrix/$basename.bed

### plotHeatmap sorted by kcluster
s=F
s=M
for l in A Z-neo Z-anc; do
basename=416.$l
computeMatrix scale-regions -S bamCompare/?416.SES.bigwig \
	-R Zhan_v3_HiC.$l.gtf -out computeMatrix/$basename.gz --missingDataAsZero \
	--metagene -b 3000 -m 5000 -a 3000 -p 16
plotHeatmap -m computeMatrix/$basename.gz -out plotHeatmap/$basename.pdf \
	--outFileSortedRegions plotHeatmap/$basename.bed \
	--samplesLabel F416 M416 -T "" --kmeans 3 \
	--yAxisLabel "log2(sample:input)" --xAxisLabel "gene distance" \
	--boxAroundHeatmaps no --heatmapWidth 8 --colorList mediumblue,lightyellow --alpha 0.8
done

--sortUsingSamples 1

for i in plotHeatmap/?416.*.bed; do cut -f4,13 $i > ${i/bed/txt}; done


### filtering non-expressed genes
cd /scratch/l338g110/Danaus_ChIP/deepTools

mkdir GTF_nonzero
grep -Fv -f <(awk '$8==0{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.-F0.gtf
grep -Fv -f <(awk '$9==0{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.-M0.gtf
grep -Fv -f <(awk '$8==0||$9==0{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.-FM0.gtf
for s in F M FM; do
awk '/chr1\>/&&$5<5685560' Zhan_v3_HiC.deeptools.-"$s"0.gtf > Zhan_v3_HiC.-"$s"0.Z-neo.gtf
awk '/chr1\>/&&$5>5685560' Zhan_v3_HiC.deeptools.-"$s"0.gtf > Zhan_v3_HiC.-"$s"0.Z-anc.gtf
awk '!/chr1\>/' Zhan_v3_HiC.deeptools.-"$s"0.gtf > Zhan_v3_HiC.-"$s"0.A.gtf
done

mkdir GTF_nonzero
mv *0*gtf GTF_nonzero

j=416
dir=GTF_nonzero
for s in F M; do
gtf=-"$s"0
computeMatrix scale-regions -S bamCompare/$s$j.SES.bigwig \
    -R $dir/Zhan_v3_HiC.$gtf.Z-neo.gtf $dir/Zhan_v3_HiC.$gtf.Z-anc.gtf $dir/Zhan_v3_HiC.$gtf.A.gtf  \
    -out computeMatrix/$s$j.SES.$gtf.Zs.A.gz --metagene --skipZeros \
    -b 3000 -m 5000 -a 3000 -p 16
plotProfile -m computeMatrix/$s$j.SES.$gtf.Zs.A.gz -out plotProfile/$s$j.SES.$gtf.Zs.A.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel $sF$j -T metagene \
    --colors olivedrab darkorange black --plotType se \
    --yAxisLabel "log2(sample:input)"
done

# use FPKM=1 cutoff, output sorted
sort -Vr -k8 ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix | awk '$8<=1.13{sub(/-TA/,"");print $1}' | xargs -i grep {} Zhan_v3_HiC.deeptools.gtf > $dir/Zhan_v3_HiC.deeptools.-F1.gtf
sort -Vr -k8 ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix | awk '$8>1.13{sub(/-TA/,"");print $1}' | xargs -i grep {} Zhan_v3_HiC.deeptools.gtf > $dir/Zhan_v3_HiC.deeptools.+F1.gtf
sort -Vr -k9 ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix | awk '$9<=1.14{sub(/-TA/,"");print $1}' | xargs -i grep {} Zhan_v3_HiC.deeptools.gtf > $dir/Zhan_v3_HiC.deeptools.-M1.gtf
sort -Vr -k9 ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix | awk '$9>1.14{sub(/-TA/,"");print $1}' | xargs -i grep {} Zhan_v3_HiC.deeptools.gtf > $dir/Zhan_v3_HiC.deeptools.+M1.gtf

for i in - +; do
for s in F M; do
awk '/chr1\>/&&$5<5685560' $dir/Zhan_v3_HiC.deeptools.$i"$s"1.gtf > $dir/Zhan_v3_HiC.$i"$s"1.Z-neo.gtf
awk '/chr1\>/&&$5>5685560' $dir/Zhan_v3_HiC.deeptools.$i"$s"1.gtf > $dir/Zhan_v3_HiC.$i"$s"1.Z-anc.gtf
awk '!/chr1\>/' $dir/Zhan_v3_HiC.deeptools.$i"$s"1.gtf > $dir/Zhan_v3_HiC.$i"$s"1.A.gtf
done; done

# heatmap & profile on A sorted by FPKM
# setting -bs 50 in computeMatrix skips 8 regions that is smaller than bin size.
# use --missingDataAsZero
dir=GTF_nonzero
j=416
j=420
# use
for s in F M; do
computeMatrix scale-regions -S bamCompare/[F,M]$j.SES.bigwig \
    -R $dir/Zhan_v3_HiC.+"$s"1.A.gtf $dir/Zhan_v3_HiC.-"$s"1.A.gtf \
    -out computeMatrix/$s$j.SES.1.A.bs50.gz --missingDataAsZero \
    --metagene -b 3000 -m 5000 -a 3000 -p 16 -bs 50
plotHeatmap -m computeMatrix/$s$j.SES.1.A.bs50.gz -out plotHeatmap/$s$j.SES.1.A.bs50.pdf \
    --regionsLabel "FPKM>1" "FPKM<1" --samplesLabel F$j M$j -T "" -T "" \
    --yAxisLabel "log2(ChIP:Input)" --xAxisLabel "gene distance" --sortRegions no \
    --boxAroundHeatmaps no --heatmapHeight 14 --heatmapWidth 6 \
    --colorList black,gold black,gold --alpha 0.8
plotProfile -m computeMatrix/$s$j.SES.1.A.bs50.gz -out plotProfile/$s$j.SES.1.A.bs50.pdf \
    --outFileNameData plotProfile/$s$j.SES.1.A.bs50.tab \
    --regionsLabel "FPKM>1" "FPKM<1" --samplesLabel F$j M$j -T "sorted by $s FPKM" \
    --colors gold black --plotType se --yAxisLabel "log2(ChIP:Input)" \
    --plotHeight 6 --plotWidth 8 --yMin -0.3 --yMax 0.9
done


for s in F M; do
computeMatrix scale-regions -S bamCompare/[F,M]$j.SES.bigwig \
    -R $dir/Zhan_v3_HiC.+"$s"1.A.gtf $dir/Zhan_v3_HiC.-"$s"1.A.gtf \
    -out computeMatrix/$s$j.SES.1.A.bs50.gz --missingDataAsZero \
    --metagene -b 3000 -m 5000 -a 3000 -p 16 -bs 50
plotProfile -m computeMatrix/$s$j.SES.1.A.bs50.gz -out plotProfile/$s$j.SES.1.A.bs50.pdf \
    --outFileNameData plotProfile/$s$j.SES.1.A.bs50.tab \
    --regionsLabel "FPKM>1" "FPKM<1" --samplesLabel F$j M$j -T "sorted by $s FPKM" \
    --colors gold black --plotType se --yAxisLabel "log2(ChIP:Input)" \
    --plotHeight 6 --plotWidth 8 --yMin -0.3 --yMax 0.9
sed -i '1,2d' plotProfile/$s$j.SES.1.A.bs50.tab
done

#
sed -ie '1,2d' F416.SES.1.A.bs50.tab
sed -ie '1,2d' M416.SES.1.A.bs50.tab

### filtering sex-biased genes
mkdir GTF_sexbias
grep -F -f <(awk '$4>0&&$7<0.05{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/edgeR_out/Dapl_xprs.hd_counts.matrix.HdF_vs_HdM.edgeR.DE_results) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.Fbias.gtf
grep -F -f <(awk '$4<0&&$7<0.05{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/edgeR_out/Dapl_xprs.hd_counts.matrix.HdF_vs_HdM.edgeR.DE_results) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.Mbias.gtf
grep -Fv -f <(awk '$7<0.05{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/edgeR_out/Dapl_xprs.hd_counts.matrix.HdF_vs_HdM.edgeR.DE_results) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.nbias.gtf
for s in F M n; do
awk '/chr1\>/&&$5<5685560' Zhan_v3_HiC.deeptools."$s"bias.gtf > Zhan_v3_HiC."$s"bias.Z-neo.gtf
awk '/chr1\>/&&$5>5685560' Zhan_v3_HiC.deeptools."$s"bias.gtf > Zhan_v3_HiC."$s"bias.Z-anc.gtf
awk '!/chr1\>/' Zhan_v3_HiC.deeptools."$s"bias.gtf > Zhan_v3_HiC."$s"bias.A.gtf
done

mv *bias*gtf GTF_sexbias

j=416
dir=GTF_sexbias
for s in F M n; do
gtf="$s"bias
computeMatrix scale-regions -S bamCompare/[F,M]$j.SES.bigwig \
    -R $dir/Zhan_v3_HiC.$gtf.Z-neo.gtf $dir/Zhan_v3_HiC.$gtf.Z-anc.gtf $dir/Zhan_v3_HiC.$gtf.A.gtf  \
    -out computeMatrix/FM$j.SES.$gtf.Zs.A.gz --metagene --skipZeros \
    -b 3000 -m 5000 -a 3000 -p 16
plotProfile -m computeMatrix/FM$j.SES.$gtf.Zs.A.gz -out plotProfile/FM$j.SES.$gtf.Zs.A.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel F M -T metagene \
    --colors olivedrab darkorange black --plotType se \
    --yAxisLabel "log2(sample:input)"
done


awk '$4>0&&$7<0.05' edgeR_out/Dapl_xprs.hd_counts.matrix.HdF_vs_HdM.edgeR.DE_results | sort -nr -k4,4 | less

DPOGS206682    histone acetyltransferase Tip60 [Culex quinquefasciatus]    0.0    68.11%    NCBI RefSeq
DPOGS216125    PREDICTED: histone acetyltransferase KAT2A-like [Nasonia vitripennis]    0.0    64.04%    NCBI nr
DPOGS210871    putative MYST histone acetyltransferase [Bombyx mori]    7e-138    97.12%    NCBI nr
DPOGS204995    histone acetyltransferase type B catalytic subunit, putative n=2 Tax=Arthropoda RepID=B7PSE2_IXOSC    4e-114    49.25%    EBI UniRef50
DPOGS206458    PREDICTED: similar to histone acetyltransferase, putative [Nasonia vitripennis]    7e-63    87.88%    NCBI RefSeq
DPOGS216126    PREDICTED: histone acetyltransferase KAT2A-like [Bombus impatiens]    5e-46    79.80%    NCBI nr
DPOGS210870    putative MYST histone acetyltransferase [Bombyx mori]    2e-27    95.16%    NCBI nr

# no candiate HAT genes show a pronoucesd sex-biased pattern
grep -F -f <(cut -f1 hat.txt) edgeR_out/Dapl_xprs.hd_counts.matrix.HdF_vs_HdM.edgeR.DE_results

DPOGS216125-TA  HdF     HdM     -0.521376924931141      5.76099517627635        0.199170632690332       0.843142836955632
DPOGS206682-TA  HdF     HdM     0.251051875007094       5.20801003670929        0.201059282303535       0.843406159058517
DPOGS204995-TA  HdF     HdM     -0.185776238620895      2.54876893975311        0.638132977832314       0.99630320679398
DPOGS206458-TA  HdF     HdM     0.0944621061490936      3.29578495271598        0.74854201993504        1
DPOGS210871-TA  HdF     HdM     -0.00308379519398053    4.22162753260144        1       1


### filtering based on (global) quantiles
mkdir GTF_quantile
grep -F -f <(awk '$8<=0.01{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.F.q0.gtf
grep -F -f <(awk '$8>0.01&&$8<=3.028175{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.F.q1.gtf
grep -F -f <(awk '$8>3.028175&&$8<=14.832{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.F.q2.gtf
grep -F -f <(awk '$8>14.832&&$8<=41.38729{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.F.q3.gtf
grep -F -f <(awk '$8>41.38729{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.F.q4.gtf

grep -F -f <(awk '$9<=0.01{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.M.q0.gtf
grep -F -f <(awk '$9>0.01&&$9<=2.867975{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.M.q1.gtf
grep -F -f <(awk '$9>2.867975&&$9<=13.98499{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.M.q2.gtf
grep -F -f <(awk '$9>13.98499&&$9<=39.14146{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.M.q3.gtf
grep -F -f <(awk '$9>39.14146{sub(/-TA/,"");print $1}' ../../Danaus_RNAseq/Dapl_xprs.mean.fpkm.matrix) Zhan_v3_HiC.deeptools.gtf > Zhan_v3_HiC.deeptools.M.q4.gtf

for q in {0..4}; do
for s in F M; do
gtf=$s.q$q
awk '/chr1\>/&&$5<5685560' Zhan_v3_HiC.deeptools.$gtf.gtf > Zhan_v3_HiC.$gtf.Z-neo.gtf
awk '/chr1\>/&&$5>5685560' Zhan_v3_HiC.deeptools.$gtf.gtf > Zhan_v3_HiC.$gtf.Z-anc.gtf
awk '!/chr1\>/' Zhan_v3_HiC.deeptools.$gtf.gtf > Zhan_v3_HiC.$gtf.A.gtf
done; done

mv *q*gtf GTF_

# by quantile
j=416
dir=GTF_quantile
mkdir plotProfile/quantile
for s in F M; do
for l in Z-neo Z-anc A; do
computeMatrix scale-regions -S bamCompare/$s$j.SES.bigwig \
    -R $dir/Zhan_v3_HiC.$s.q0.$l.gtf $dir/Zhan_v3_HiC.$s.q1.$l.gtf \
    $dir/Zhan_v3_HiC.$s.q2.$l.gtf $dir/Zhan_v3_HiC.$s.q3.$l.gtf $dir/Zhan_v3_HiC.$s.q4.$l.gtf \
    -out computeMatrix/$s$j.SES.$l.quantile.gz --metagene --skipZeros \
    -b 3000 -m 5000 -a 3000 -p 16 -bs 50
plotProfile -m computeMatrix/$s$j.SES.$l.quantile.gz -out plotProfile/quantile/$s$j.SES.$l.quantile.pdf \
    --regionsLabel q0 q1 q2 q3 q4 --samplesLabel $s$j.$l -T "" --plotType se --yAxisLabel "log2(ChIP:Input)"
done; done

# by linkage
j=416
dir=GTF_quantile
for q in {0..4}; do
for s in F M; do
gtf=$s.q$q
computeMatrix scale-regions -S bamCompare/$s$j.SES.bigwig \
    -R $dir/Zhan_v3_HiC.$gtf.Z-neo.gtf $dir/Zhan_v3_HiC.$gtf.Z-anc.gtf $dir/Zhan_v3_HiC.$gtf.A.gtf \
    -out computeMatrix/$j.SES.$gtf.Zs.A.gz --metagene --skipZeros \
    -b 3000 -m 5000 -a 3000 -p 16 -bs 50
plotProfile -m computeMatrix/$j.SES.$gtf.Zs.A.gz -out plotProfile/quantile/$j.SES.$gtf.Zs.A.pdf \
    --regionsLabel Z-neo Z-anc A --samplesLabel $gtf -T "" --colors olivedrab darkorange black \
    --plotType se --yAxisLabel "log2(ChIP:Input)"
done; done


##### MACS2
$ Danaus_ChIP

# alignments includes scaffold fragments!
mkdir bamFile.A
for i in bamFile/*bam; do
	samtools view -h -@ 16 $i -q 255 | awk '/chr/&&$3!="chr1"' | samtools view -b -@ 16 > ${i/bamFile/bamFile.A}
done
mkdir bamFile.Z
for i in bamFile/*bam; do
	samtools view -b -@ 16 $i chr1 -q 255 > ${i/bamFile/bamFile.Z}
done

module load MACS2

j=416
mkdir MACS2.$j
for s in F M; do
	macs2 callpeak -c bamFile.A/$s.*bam -t bamFile.A/$s416*bam -n $s.A --outdir MACS2.$j -g 198015805 \
-f BAMPE --bw 147 -B --keep-dup 1 --call-summits 2>&1 | tee -a MACS2.$j/macs2.$s.A.$(date +"%Y%m%d").log
	macs2 callpeak -c bamFile.Z/$s.*bam -t bamFile.Z/$s416*bam -n $s.Z --outdir MACS2.$j -g 15360736 \
-f BAMPE --bw 147 -B --keep-dup 1 --call-summits 2>&1 | tee -a MACS2.$j/macs2.$s.Z.$(date +"%Y%m%d").log
done


wc -l MACS2.*/F.?_peaks.narrowPeak

wc -l MACS2.*/M.?_peaks.narrowPeak


for i in bamFile/*bam; do
	samtools view -h -@ 16 $i -q 255 | awk '/chr/&&$3!="chr1"' | samtools view -Sb -@ 16 > ${i/bamFile/bamFile.A}
done

samtools view -s 0.2 -@ 16 -b bamFile.Z/M.star.2018.bam > bamFile.Z/M.star.2018.0.2.bam
samtools view -s 0.5 -@ 16 -b bamFile.A/F416.star.2018.bam > bamFile.A/F416.star.2018.0.5.bam
samtools view -s 0.3 -@ 16 -b bamFile.A/M.star.2018.bam > bamFile.A/M.star.2018.0.3.bam

j=416
mkdir MACS2.$j.sub.narrow
	macs2 callpeak -c bamFile.A/F.star.2016.bam -t bamFile.A/F416.star.2018.0.5.bam -n F.A --outdir MACS2.$j.sub -g 198015805 \
-f BAMPE --bw 147 -B --keep-dup 1 --call-summits 2>&1 | tee -a MACS2.$j.sub/macs2.$s.A.$(date +"%Y%m%d").log
	macs2 callpeak -c bamFile.A/M.star.2018.0.3.bam -t bamFile.A/M416.star.merge.bam -n M.A --outdir MACS2.$j.sub -g 198015805 \
-f BAMPE --bw 147 -B --keep-dup 1 --call-summits 2>&1 | tee -a MACS2.$j.sub/macs2.$s.A.$(date +"%Y%m%d").log

	macs2 callpeak -c bamFile.Z/F.star.2016.bam -t bamFile.Z/F416.star.2018.bam -n F.Z --outdir MACS2.$j.sub -g 15360736 \
-f BAMPE --bw 147 -B --keep-dup 1 --call-summits 2>&1 | tee -a MACS2.$j.sub/macs2.$s.Z.$(date +"%Y%m%d").log
	macs2 callpeak -c bamFile.Z/M.star.2018.0.2.bam -t bamFile.Z/M416.star.merge.bam -n M.Z --outdir MACS2.$j.sub -g 15360736 \
-f BAMPE --bw 147 -B --keep-dup 1 --call-summits 2>&1 | tee -a MACS2.$j.sub/macs2.$s.Z.$(date +"%Y%m%d").log

# indeed downsampling (smaller library) results in more peaks called!

for s in F M; do
for l in A Z; do
	macs2 bdgcmp -c $dir/$s."$l"_control_lambda.bdg -t $dir/$s."$l"_treat_pileup.bdg -p 1 -m logFE --outdir $dir --o-prefix bdgcmp_$s.$l
	macs2 bdgbroadcall -i $dir/bdgcmp_$s."$l"_logFE.bdg -g 36 --outdir $dir --o-prefix bdgbroadcall_$s.$l
done; done
# bdgbroadcall does not produce any broadpeaks!

j=416
dir=MACS2.$j.sub
mkdir $dir
	macs2 callpeak --broad -c bamFile.A/F.star.2016.bam -t bamFile.A/F416.star.2018.0.5.bam -n F.A --outdir MACS2.$j.sub -g 198015805 \
-f BAMPE --bw 147 -B --keep-dup 1 2>&1 | tee -a MACS2.$j.sub/macs2.F.A.$(date +"%Y%m%d").log
	macs2 callpeak --broad -c bamFile.A/M.star.2018.0.3.bam -t bamFile.A/M416.star.merge.bam -n M.A --outdir MACS2.$j.sub -g 198015805 \
-f BAMPE --bw 147 -B --keep-dup 1 2>&1 | tee -a MACS2.$j.sub/macs2.M.A.$(date +"%Y%m%d").log

	macs2 callpeak --broad -c bamFile.Z/F.star.2016.bam -t bamFile.Z/F416.star.2018.bam -n F.Z --outdir MACS2.$j.sub -g 15360736 \
-f BAMPE --bw 147 -B --keep-dup 1 2>&1 | tee -a MACS2.$j.sub/macs2.F.Z.$(date +"%Y%m%d").log
	macs2 callpeak --broad -c bamFile.Z/M.star.2018.0.2.bam -t bamFile.Z/M416.star.merge.bam -n M.Z --outdir MACS2.$j.sub -g 15360736 \
-f BAMPE --bw 147 -B --keep-dup 1 2>&1 | tee -a MACS2.$j.sub/macs2.M.Z.$(date +"%Y%m%d").log

wc -l MACS2.416.sub/*broadPeak
  10375 MACS2.416.sub/F.A_peaks.broadPeak
    922 MACS2.416.sub/F.Z_peaks.broadPeak
  12824 MACS2.416.sub/M.A_peaks.broadPeak
   1130 MACS2.416.sub/M.Z_peaks.broadPeak
  25251 total

# make bigwig score tracks
module load kentUtils/359
# parallel command not available??
# find $dir -name '*.bdg' | parallel "sort -k1,1 -k2,2n {} > {.}.sort.bdg"

find $dir -name '*.bdg' | xargs -P 16 -i sort -k1,1 -k2,2n {} > {}.sorted
find $dir -name '*.sorted' | xargs -P 16 -i bedGraphToBigWig {} ../Danaus_HiRise/Dapl_Zhan_v3_HiC.chr.bed {}.bigwig

cat $dir/F.Z_control_lambda.bdg $dir/F.A_control_lambda.bdg | sort -k1,1 -k2,2n > $dir/F_control_lambda.bdg
cat $dir/M.Z_control_lambda.bdg $dir/M.A_control_lambda.bdg | sort -k1,1 -k2,2n > $dir/M_control_lambda.bdg

cat $dir/F.Z_treat_pileup.bdg $dir/F.A_treat_pileup.bdg | sort -k1,1 -k2,2n > $dir/F_treat_pileup.bdg
cat $dir/M.Z_treat_pileup.bdg $dir/M.A_treat_pileup.bdg | sort -k1,1 -k2,2n > $dir/M_treat_pileup.bdg

find $dir -name ?_*bdg | xargs -P 16 -i bedGraphToBigWig {} ../Danaus_HiRise/Dapl_Zhan_v3_HiC.chr.seq_length.txt {}.bigwig

cd $dir
cat F.Z_peaks.broadPeak F.A_peaks.broadPeak > F_peaks.broadPeak
cat M.Z_peaks.broadPeak M.A_peaks.broadPeak > M_peaks.broadPeak

cat F.Z_peaks.gappedPeak F.A_peaks.gappedPeak > F_peaks.gappedPeak
cat M.Z_peaks.gappedPeak M.A_peaks.gappedPeak > M_peaks.gappedPeak

for i in $dir/?_*broadPeak; do cut -f1-3 $i | sort -V | awk 'BEGIN{getline;a=$1}$1!=a{print "#\n"$0;a=$1;next}1' > $i.bed; done

# breakpoint
5685560

##### visualized peaks in deepTools
bamCompare -b1 $dirIn/F420.*bam -b2 $dirIn/F.*bam \
    -o $dirOut/F420.SES.bigwig -of bigwig  \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactors 1:0.28393517 --ratio log2 --skipNAs \
    -bs 50 --smoothLength 150 -p 16 -e --centerReads \
    --maxFragmentLength 700
echo "bamCompare -b1 $dirIn/M$j.*bam -b2 $dirIn/M.*bam \
    -o $dirOut/M$j.SES.bigwig -of bigwig  \
    --ignoreDuplicates --minMappingQuality 255 \
    --scaleFactorsMethod SES --ratio log2 --skipNAs \
    -bs 50 --smoothLength 150 -p 16 -e --centerReads \
    --maxFragmentLength 700 -v &>> $dirOut/M.$ts.log " | tee -a $dirOut/M.$ts.log | bash

-bs 10

j=416
for s in F M; do
basename=$s$j.broadPeak
computeMatrix reference-point -S bamCompare/$s$j.SES.bigwig \
    -R ../MACS2.416.sub/"$s"_peaks.broadPeak.bed -out computeMatrix/$basename.gz \
    --referencePoint center -a 1000 -b 1000 -p 16
plotHeatmap -m computeMatrix/$basename.gz -out plotHeatmap/$basename.pdf \
    --samplesLabel $s --refPointLabel 0 -T "" --legendLocation none \
    --regionsLabel Z-neo Z-anc 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 \
    --xAxisLabel "distance from peak center" --yAxisLabel "log2(sample:input)" \
    --boxAroundHeatmaps no --heatmapWidth 6 --heatmapHeight 42 \
    --colorList lightblue,lemonchiffon --alpha 0.8
done

    --zMin -1.8 --zMax 3.5 --yMin 0 --yMax 2.8	\

j=416
for f in exon intron intergenic; do
# total length for features
awk 'BEGIN{getline;a=$1;sum1=$2;sum2=$3}$1==a{sum1+=$2;sum2+=$3}
	$1!=a{print a"\t",sum2-sum1;sum1=sum2=0;a=$1;sum1+=$2;sum2+=$3}
	END{print a"\t",sum2-sum1}' deepTools/Zhan_v3_HiC_$f.bed6 > Zhan_v3_HiC_$f.length.txt
for s in F M; do
# peak length sum per each chr, for A
intersectBed -a deepTools/Zhan_v3_HiC_$f.bed6 -b MACS2.416.sub/$s.A_peaks.broadPeak -wo | cut -f1,16 | \
	awk 'BEGIN{getline;a=$1;sum=$2}$1==a{sum+=$2}$1!=a{print a"\t",sum;sum=0;a=$1;sum+=$2}END{print a"\t",sum}' \
	> broadPeakCovSum.$s$j.A.$f.txt
# do the same for Z
awk '{print}/#/{exit}' deepTools/Zhan_v3_HiC_$f.bed6 | intersectBed -a stdin -b MACS2.416.sub/$s.Z_peaks.broadPeak -wo | cut -f1,16 | \
	awk '{a=$1;sum+=$2}END{print a"\t",sum}' \
	> broadPeakCovSum.$s$j.Z-neo.$f.txt
awk '/#/{f=1;next}/chr2/{f=0}f' deepTools/Zhan_v3_HiC_$f.bed6 | intersectBed -a stdin -b MACS2.416.sub/$s.Z_peaks.broadPeak -wo | cut -f1,16 | \
	awk '{a=$1;sum+=$2}END{print a"\t",sum}' \
	> broadPeakCovSum.$s$j.Z-anc.$f.txt
# join three, order is important
cat broadPeakCovSum.$s$j.Z-neo.$f.txt broadPeakCovSum.$s$j.Z-anc.$f.txt broadPeakCovSum.$s$j.A.$f.txt > broadPeakCovSum.$s$j.$f.txt
done
# join with length file, don't use join: two chr1 confuses it
#join <(sed '/#/d' Zhan_v3_HiC_$f.length.txt) broadPeakCovSum.$s$j.$f.txt | sed '2d;3d' > $f.length_broadPeakCovSum.$s$j.txt
paste <(sed '/#/d' Zhan_v3_HiC_$f.length.txt) broadPeakCovSum.F$j.$f.txt broadPeakCovSum.M$j.$f.txt | \
	cut -f1,2,4,6 | sed 's/ //g' > $f.length_broadPeakCovSum.FM$j.txt
done

mkdir PeakIntersect
mv *txt PeakIntersect


### DIFF
j=416
dir1=MACS2.$j.sub
dir2=MACS2.$j.sub.diff
mkdir $dir2
# A
macs2 bdgdiff --c1 $dir1/F.A_control_lambda.bdg --c2 $dir1/M.A_control_lambda.bdg \
	--t1 $dir1/F.A_treat_pileup.bdg --t2 $dir1/M.A_treat_pileup.bdg \
	--d1 16584650 --d2 11282014 --outdir $dir2 --o-prefix Diff.$j.$l
# Z
macs2 bdgdiff --c1 $dir1/F."$l"_control_lambda.bdg --c2 $dir1/M."$l"_control_lambda.bdg \
	--t1 $dir1/F."$l"_treat_pileup.bdg --t2 $dir1/M."$l"_treat_pileup.bdg \
	--d1 1712922 --d2 944581 --outdir $dir2 --o-prefix Diff.$j.$l
done


##### Correlation between gene density and peak density
# 100Kb bins
bedtools makewindows -g ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chr.seq_length.txt -w 100000 > Dapl_Zhan_v3_HiC.chr.100k.bed
intersectBed -a Dapl_Zhan_v3_HiC.chr.100k.bed -b <(grep -w gene Zhan_v3_HiC.deeptools.gtf) -c > Zhan_v3_HiC.chr.gene.density.100k.txt
intersectBed -a Dapl_Zhan_v3_HiC.chr.100k.bed -b ../MACS2.416.sub/F_peaks.broadPeak.bed -c > ../MACS2.416.sub/F_peaks.broadPeak.100k.txt
intersectBed -a Dapl_Zhan_v3_HiC.chr.100k.bed -b ../MACS2.416.sub/M_peaks.broadPeak.bed -c > ../MACS2.416.sub/M_peaks.broadPeak.100k.txt

# 10Kb bins: does not work well for circos! Counts too low!
bedtools makewindows -g ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chr.seq_length.txt -w 10000 > Dapl_Zhan_v3_HiC.chr.10k.bed
intersectBed -a Dapl_Zhan_v3_HiC.chr.10k.bed -b <(grep -w gene Zhan_v3_HiC.deeptools.gtf) -c > Zhan_v3_HiC.chr.gene.density.10k.txt


##### BETA
mkdir BETA

wget http://cistrome.org/BETA/src/BETA_1.0.7.zip BETA
unzip BETA_1.0.2.zip
cd BETA_1.0.2
python setup.py install --prefix=/scratch/l338g110/Danaus_ChIP/BETA/BETA_1.0.7
#Modify PYTHONPATH if necessary by adding the following two lines in your .bashrc (or .bash_profile of mac) file in home directory.
export PATH=/scratch/l338g110/Danaus_ChIP/BETA/BETA_1.0.7/bin:$PATH
export PYTHONPATH=//scratch/l338g110/Danaus_ChIP/BETA/BETA_1.0.7/:$PYTHONPATH
BETA

# for -r genome annotation
awk -F '[\t\"]' '$3=="transcript"{print $10,$1,$7,$4,$5,$10}' OFS="\t" deepTools/Zhan_v3_HiC.deeptools.gtf > BETA/Zhan_v3_HiC.BETA.bed
# for -p peak file, take 7th field (fold-change)
for i in MACS2.416.sub/?_peaks.broadPeak; do cut -f1-4,7 $i > BETA/$(basename $i).bed; done
for i in MACS2.416.sub/?.?_peaks.broadPeak; do cut -f1-4,7 $i > BETA/$(basename $i).bed; done
# for -e DE file, take 6th field (pvalue) or 7th (FDR)
cut -f1,4,6 ../Danaus_RNAseq/edgeR_out/Dapl_xprs.hd_counts.matrix.HdF_vs_HdM.edgeR.DE_results | sed '1d' > BETA/Dapl_xprs.HdF_vs_HdM.edgeR.bsf


cd /nfs/apps/7/arch/generic/build/BETA/BETA_1.0.7
python setup.py clean
python setup.py build
python setup.py install --prefix=/nfs/apps/7/arch/generic/BETA/1.0.7

# note to remove first line of peak file! Otherwise cause error!
cd BETA

# BETA basic does not identify supression/activation fucntion
# BETA basic -p $dir2/Diff.416.A_c3.0_cond1.bed -e BETA/Dapl_xprs.HdF_vs_HdM.edgeR.bsf -k BSF -r Zhan_v3_HiC.BETA.bed \
	--da 500 -o F416.A.broad -d 10000 -c 0.05

for s in F M; do
cut -f1-5 ../MACS2.416.sub/"$s"_peaks.broadPeak > "$s"_peaks.broadPeak.5c.bed
BETA minus -p "$s"_peaks.broadPeak.5c.bed -r Zhan_v3_HiC.BETA.bed \
	-o "$s"416.broad -n MACS2.416.sub.$s.broadPeak
done


tblastn -db ../Danaus_HiRise/BLAST_DB/Dapl_cds -query mof_dm.prot.fa -out mof_dm.tblastn.out -evalue 1e-5
NP_5110511malesabsentonthefirstDrosophilamelanogaster    DPOGS210871    1e-88    63.95%    245    732    MOF protein [Bombyx mori]    Uncharacterized protein n=42 Tax=Metazoa RepID=F1QX50_DANRE


Database: Dapl_cds
15,006 sequences; 20,869,128 total letters

Query= NP_511051.1 males absent on the first [Drosophila melanogaster]

Length=827
Score     E
Sequences producing significant alignments:                          (Bits)  Value

DPOGS210871-TA                                                      323     4e-105
DPOGS206682-TA                                                      311     6e-97
DPOGS213223-TA                                                      285     4e-81
DPOGS202817-TA                                                      148     2e-37
DPOGS202818-TA                                                      123     1e-32

grep -f <(awk '/e-/&&/-TA/{print $1}' mof_dm.tblastn.out) ../Danaus_RNAseq/edgeR_out/matrix.counts.matrix.HdF_vs_HdM.edgeR.DE_results
DPOGS206682-TA  HdF     HdM     0.256166891610559       5.34495571733485        0.202152976435368       0.953714650849684
DPOGS202818-TA  HdF     HdM     -0.628037255189841      0.821301192420635       0.309114407706765       1
DPOGS202817-TA  HdF     HdM     0.273087450521546       1.96873422714799        0.555286856770402       1
DPOGS213223-TA  HdF     HdM     -0.0686697259630548     6.64403319055375        0.718144860726004       1
DPOGS210871-TA  HdF     HdM     0.0683022420744969      4.83155813118392        0.735965809713101       1



# Dmel clamp isoform A 561aa B 566 aa
tblastn -db ../Danaus_HiRise/BLAST_DB/Dapl_cds -query clamp_A_dm.prot.fa -out clamp_A_dm.tblastn.out -evalue 1e-5
# same as B and blastp, both highest hit: DPOGS213107; zinc finger
tblastn -db ../Danaus_HiRise/BLAST_DB/Dapl_cds -query clamp_B_dm.prot.fa -out clamp_B_dm.tblastn.out -evalue 1e-5
NP_6101371Chromatin-linkedadaptorforMSLproteinsisoformADrosophilamelanogaster    DPOGS213107    7e-131    53.86%    435    1602    zinc finger protein [Culex quinquefasciatus]    Zinc finger protein n=3 Tax=Endopterygota RepID=G6CSE7_DANPL

tblastn -db ../Danaus_HiRise/BLAST_DB/Dapl_cds -query clamp_A_dm.prot.fa -out clamp_A_dm.tblastn.out -evalue 1e-5


# Dmel MSL1 (isoform A&B identical)
tblastn -db ../Danaus_HiRise/BLAST_DB/Dapl_cds -query msl1_dm.prot.fa -out msl1_dm.tblastn.out -evalue 1e-5
NP_4768961male-specificlethal1isoformADrosophilamelanogaster    DPOGS201387    7e-10    31.34%    128    2301    MSL1 protein [Bombyx mori]    MSL1 protein n=1 Tax=Bombyx mori RepID=A5JPL7_BOMMO
AAI189981MSL1proteinHomosapiens    DPOGS201387    2e-13    41.13%    128    2301    MSL1 protein [Bombyx mori]    MSL1 protein n=1 Tax=Bombyx mori RepID=A5JPL7_BOMMO

NP_5234671male-specificlethal2isoformADrosophilamelanogaster    DPOGS215626    1e-09    38.57%    60    1287    MSL2 protein [Bombyx mori]    MSL2 protein n=3 Tax=Neoptera RepID=A5JPL8_BOMMO

# Dmel mle is 1293 aa
NP_4766411malelessisoformADrosophilamelanogaster    DPOGS206766    0.0    56.94%    133    2415    MLE protein [Bombyx mori]    Dosage compensation regulator n=2 Tax=Coelomata RepID=F4W6K4_ACREC

# Bmori masc: > DPOGS203394
tblastn -db ../Danaus_HiRise/BLAST_DB/Dapl_cds -query masc_bm.prot.fa -out masc_bm.tblastn.out -evalue 1e-5
# DPOGS203394 is not in the hit-list of CLAMP blastp

# C. elegans sir-2.1
tblastn -db ../Danaus_HiRise/BLAST_DB/Dapl_cds -query SIR-2.1.fa -out sir_ce.tblastn.out -evalue 1e-5
Sequences producing significant alignments:                          (Bits)  Value

DPOGS204184-TA                                                      208     1e-61
DPOGS208880-TA                                                      106     9e-26
DPOGS206876-TA                                                      91.7    1e-20
DPOGS205396-TA                                                      61.2    4e-10

grep -f <(awk '/e-/&&/-TA/{print $1}' sir_ce.tblastn.out) ../Danaus_RNAseq/edgeR_out/matrix.counts.matrix.HdF_vs_HdM.edgeR.DE_results

DPOGS206876-TA  HdF     HdM     -0.258423630418048      4.67958584466369        0.275911808543543       1
DPOGS204184-TA  HdF     HdM     0.088923910171518       4.74752908308835        0.673141296131861       1
DPOGS205396-TA  HdF     HdM     -0.0817102855595828     5.03068711468357        0.757441377212563       1
DPOGS208880-TA  HdF     HdM     0.0288408988400014      5.12280179765927        0.892433166576232       1


<plot>

file  = data/Zhan_v3_HiC.chr.repeat.density.100k.txt
color   = vdgrey
orientation = in
r1    = 0.85r
r0    = 0.75r
#min   = 0
#max   = 25
fill_color = vdgrey_a2

</plot>

### GATK
java‐jarpicard.jarMarkDuplicates\
I=bams/exp_design/NA12878_rnaseq_20.bam\
O=sandbox/NA12878_rnaseq_20_dedup.bam\
CREATE_INDEX=true\
M=output.metrics


##### Bowtie
module load FastQC multiqc
for i in Reads_raw/*; do mv $i ${i/6440_7149_299??_HT7L3BGXX_/}; done

for i in Reads_raw/*R1*
do
out=${i/raw/trim}
java -jar /nfs/apps/7/arch/generic/trimmomatic/0.36/trimmomatic-0.36.jar PE -threads 16 -phred33 \
    $i ${i/R1/R2} ${out/R1/R1.P} ${out/R1/R1.U} ${out/R1/R2.P} ${out/R1/R2.U} \
    ILLUMINACLIP:/nfs/apps/7/arch/generic/trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10:8:true \
LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:5 2>&1 | tee -a Reads_trim/trimmomatic.$ts.log
done
# Reads dropped all very low
for i in Reads_trim/HD*; do mv $i ${i/HD-/}; done
for i in Reads_trim/?-*; do mv $i ${i/-/}; done

fastqc -t 16 Reads_trim/*T* -o FastQC_trim
# based on report, at 50m already reach saturation: ~8% dup level with ~15mi; but ~15% dup with 50mi
multiqc FastQC_trim/*P* -o MultiQC/FastQC_trim_P_all
# Mapping with Bowtie, combine two sets of sequencing reads
# Note to set -X 700!
mkdir BAM
module load Bowtie2 SAMtools BEDTools R
ts=$(date +"%Y%m%d")
for i in F M F416 M416 F420 M420; do
j=$(echo Reads_trim/${i}_*R1.P* | sed 's/ /,/g')
echo "bowtie2 -x ../Danaus_HiRise/BowtieIndex/Dapl.chr -1 $j -2 ${j//R1/R2} --very-sensitive-local --no-unal -X 1000 -p 16 2>> bowtie2.$ts.log \
    | samtools sort -m 5G -@ 16 -o BAM/$i.bam" | tee -a bowtie2.$ts.log | bash
done
# bowtie log run on 29th with "-X 1000"
sed -i '/Warning/d' bowtie2.20190829.log
# >90 overall alignment rate; >70 unique mapping

samtools view -c -@ 15 BAM/F.bam
30228717

16321018 reads; of these:
12118856 (74.25%) aligned concordantly exactly 1 time
804423 (20.84%) aligned exactly 1 time
2413319

22873196


for i in BAM/*.bam; do
samtools view -@ 16 $i | grep -v "XS:" | cat <(samtools view -H $i) - | samtools view -b - -@ 16 > ${i/bam/unique.bam}
done

for i in BAM/*.bam; do samtools index $i -@ 16; done


dir=deepToopls_Bowtie2
mkdir $dir

#bamPEFragmentSize
for s in F M; do
bamPEFragmentSize -b BAM/$s*bam -hist $dir/bamPEFragmentSize_$s.pdf -p 16 \
    --samplesLabel $(echo BAM/$s*bam | sed 's/BAM\///g;s/.bam//g') \
    -T "Sample Fragment Size" 2>&1 | tee -a $dir/bamPEFragmentSize_$s.log
done

#### MACS
module load bioepic/0.1.25-Python-2.7.12 Jellyfish/2.2.6

epic-effective --read-length=36 --nb-cpu=16 ../Danaus_HiRise/Dapl_Zhan_v3_HiC.chr.fasta
Genome length:  245185305
Number unique 36-mers: 213161022 (File: effective_genome_size, Log level: INFO, Time: Thu, 29 Aug 2019 16:21:13 )
Effective genome size: 0.869387429234 (File: effective_genome_size, Log level: INFO, Time: Thu, 29 Aug 2019 16:21:13 )

module load MACS2
j=416
dir=MACS2.$j
dir=MACS2.$j.cutoff
mkdir $dir
for s in F M; do
	macs2 callpeak -c BAM/$s.unique.bam -t BAM/$s$j.unique.bam -n $s$j --outdir $dir -g 213161022 \
-f BAMPE --bw 147 -B --keep-dup 1 --broad --cutoff-analysis 2>&1 | tee -a $dir/$s$j.$(date +"%Y%m%d").log
done

# --cutoff-analysis option has no effect on output results


total fragments in treatment: 45199050
total fragments in control: 10985069

for s in F M; do
	macs2 callpeak -c bamFile.A/$s.*bam -t bamFile.A/$s416*bam -n $s.A --outdir MACS2.$j -g 198015805 \
-f BAMPE --bw 147 -B --keep-dup 1 -B --call-summits 2>&1 | tee -a MACS2.$j/macs2.$s.A.$(date +"%Y%m%d").log
	macs2 callpeak -c bamFile.Z/$s.*bam -t bamFile.Z/$s416*bam -n $s.Z --outdir MACS2.$j -g 15360736 \
-f BAMPE --bw 147 -B --keep-dup 1 --call-summits 2>&1 | tee -a MACS2.$j/macs2.$s.Z.$(date +"%Y%m%d").log
done

grep -v "#" MACS2.416/F416_peaks.xls | awk 'NR>2{ a[$1] += $4; n+= 1} END { for (i in a) print i, n, a[i] }' OFS="\t"
grep -v "#" MACS2.416/F416_peaks.xls | cut -f1 | uniq -c

# Get top ranking peaks; F4:peak_name; F7:FC; F8:-log10(pvalue); F9: -log10(qvalue)
for i in F M; do
awk '$1=="chr1"&&$3<5685560' ${i}416_peaks.broadPeak | sort -rk7,7 > ${i}416_peaks.broadPeak.Zn
awk '$1=="chr1"&&$3>=5685560' ${i}416_peaks.broadPeak | sort -rk7,7 > ${i}416_peaks.broadPeak.Za
awk '$1!="chr1"' ${i}416_peaks.broadPeak | sort -rk7,7 > ${i}416_peaks.broadPeak.A
done

for i in ?416_peaks.broadPeak.*; do
fastaFromBed  -name -fo $i.fa -fi ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chr.fasta -bed $i
done

module load ProSampler
mkdir ProSampler
# Top 10% peaks
ProSampler -i <(head -n 44 F416_peaks.broadPeak.Zn.fa) -b ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chrZ.fasta -o ProSampler/F416_peaks.broadPeak.Zn
ProSampler -i <(head -n 64 F416_peaks.broadPeak.Za.fa) -b ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chrZ.fasta -o ProSampler/F416_peaks.broadPeak.Za
ProSampler -i <(head -n 2034 F416_peaks.broadPeak.A.fa) -b ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chrA.fasta -o ProSampler/F416_peaks.broadPeak.A

ProSampler -i <(head -n 92 M416_peaks.broadPeak.Zn.fa) -b ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chrZ.fasta -o ProSampler/M416_peaks.broadPeak.Zn
ProSampler -i <(head -n 138 M416_peaks.broadPeak.Za.fa) -b ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chrZ.fasta -o ProSampler/M416_peaks.broadPeak.Za
ProSampler -i <(head -n 2464 M416_peaks.broadPeak.A.fa) -b ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chrA.fasta -o ProSampler/M416_peaks.broadPeak.A

i=F
for i in F M; do
n=$(wc -l ${i}416_peaks.broadPeak.Zn.fa)
ProSampler -i ${i}416_peaks.broadPeak.Zn.fa -b ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chrZ.fasta -o ${i}416_peaks.broadPeak.Zn
ProSampler -i ${i}416_peaks.broadPeak.Za.fa -b ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chrZ.fasta -o ${i}416_peaks.broadPeak.Za
ProSampler -i ${i}416_peaks.broadPeak.A.fa -b ../../Danaus_HiRise/Dapl_Zhan_v3_HiC.chrA.fasta -o ${i}416_peaks.broadPeak.A
done

# Make TSS+-500 bed files: - strand use $5, + strand use $4
grep -w 'transcript' ../../Danaus_HiRise/Zhan_v3_HiC.deeptools.gtf | awk '$7=="-"?a=$5"\t"$5+500:a=$4-500"\t"$4{print $1,a,$7}' OFS="\t" > tss.+500.bed
grep -w 'transcript' ../../Danaus_HiRise/Zhan_v3_HiC.deeptools.gtf | awk '$7=="-"?a=$5-500"\t"$5:a=$4"\t"$4+500{print $1,a,$7}' OFS="\t" > tss.-500.bed
grep -w 'transcript' ../../Danaus_HiRise/Zhan_v3_HiC.deeptools.gtf | awk '$7=="-"?a=$5-500"\t"$5+500:a=$4-500"\t"$4+500{print $1,a,$7}' OFS="\t" > tss.+-500.bed

#grep -w 'transcript' ../../Danaus_HiRise/Zhan_v3_HiC.deeptools.gtf | awk '$7=="-"?a=$5"\t"$5+1000:a=$4-1000"\t"$4{print $1,a,$7}' OFS="\t" > tss.+1000.bed
#grep -w 'transcript' ../../Danaus_HiRise/Zhan_v3_HiC.deeptools.gtf | awk '$7=="-"?a=$5-1000"\t"$5:a=$4"\t"$4+1000{print $1,a,$7}' OFS="\t" > tss.-1000.bed

intersectBed -a tss.+500.bed -b F416_peaks.broadPeak -wo | less


wc -l ?416_peaks.broadPeak
10700 F416_peaks.broadPeak
13455 M416_peaks.broadPeak

# Over a third of peaks are within TSS+-500 region
intersectBed -a M416_peaks.broadPeak -b tss.+500.bed tss.-500.bed -u | wc -l
4601
intersectBed -a M416_peaks.broadPeak -b tss.+-500.bed -u | wc -l
4601
intersectBed -a F416_peaks.broadPeak -b tss.+-500.bed -u | wc -l
3754

intersectBed -a <(head -n 22 F416_peaks.broadPeak) -b tss.+-500.bed -u | wc -l


### SRA deposition
$ /panfs/panfs.ittc.ku.edu/scratch/l338g110/Danaus_ChIP/Reads_raw
lftp -u subftp,w4pYB9VQ -e 'cd uploads/lg356@cornell.edu_RfJE9S4l;mkdir Dp_ChIP;cd Dp_ChIP;mput -P 6 *gz' ftp://ftp-private.ncbi.nlm.nih.gov
lftp -u subftp,w4pYB9VQ -e 'cd uploads/lg356@cornell.edu_RfJE9S4l/Dp_ChIP;mput -P 6 F4* M4* M_*' ftp://ftp-private.ncbi.nlm.nih.gov



