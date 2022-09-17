setwd("~/Documents/Danaus/RNAseq_DC")

# read fpkm matrix
ms_fpkm = read.table("manduca_fpkm_head.txt", header = T)
dp_fpkm = read.table("danaus.head.fpkm.tmm.matrix.txt", header = T)

rbh_DpMs = read.table("rbh_DpMs.txt")

# calculate mean
ms_fpkm$HdF = rowMeans(ms_fpkm[,2:5])
ms_fpkm$HdM = rowMeans(ms_fpkm[,6:9])
ms_fpkm$Ratio = ms_fpkm[,10]/ms_fpkm[,11]

#
sp1.fpkm = dp_fpkm[match(rbh_DpMs$V1, dp_fpkm$Gene),c(1,2,7,11,12)]
sp2.fpkm = ms_fpkm[match(rbh_DpMs$V2, ms_fpkm$gene_id), c(1,10:12)]

rbh.fpkm = cbind(sp1.fpkm, sp2.fpkm)
colnames(rbh.fpkm)[c(1,3:5,6:9)] = c("Gene_sp1", "HdF.sp1", "HdM.sp1", "Ratio.sp1", 
                                     "Gene_Sp2", "HdF.sp2", "HdM.sp2", "Ratio.sp2")

rbh.fpkm1 = subset(rbh.fpkm, HdF.sp1&HdF.sp2>0)
rbh.fpkm2 = subset(rbh.fpkm, HdM.sp1&HdM.sp2>0)

# factor IR.F
#x.IR.F = median(rbh.fpkm1[rbh.fpkm1$Linkage=="A",3])/median(rbh.fpkm1[rbh.fpkm1$Linkage=="A",7])
x.IR.F = with(rbh.fpkm1[rbh.fpkm1$Linkage=="A",], median(HdF.sp1)/median(HdF.sp2))
# factor IR.M
x.IR.M = median(rbh.fpkm2[rbh.fpkm2$Linkage=="A",4])/median(rbh.fpkm2[rbh.fpkm2$Linkage=="A",8])

rbh.fpkm1$IR.F = with(rbh.fpkm1, HdF.sp1/HdF.sp2/x.IR.F)
rbh.fpkm2$IR.M = with(rbh.fpkm2, HdM.sp1/HdM.sp2/x.IR.M)

gdo = c("grey", "darkorange", "olivedrab")
gdo1 = c("grey40", "darkorange4", "darkgreen")
border = gdo1

pdf("DpMs_ratio_1.pdf", width = 5, height = 8, pointsize = 14)
#pdf("DpMs_HmMs_ratio.pdf", width = 12, height = 8, pointsize = 14)
#par(mfrow = c(1,2))
boxplot(IR.F ~ Linkage, rbh.fpkm1, boxwex = 0.25, outline = F, notch = T, 
        col = c("grey60", "darkorange", "olivedrab"), xaxt = "n", las = 1,
        xlim = c(0, 3), ylim = c(0, 4.25), main = "D.plexippus:M.sexta",
        ylab = "interspecies FPKM ratio", cex.lab = 1.1, at = c(0.25, 1.25, 2.25), 
        range =1, whisklty = 5, staplewex = 0.25, boxlty = 0)
text(c(0.25, 1.25, 2.25), rep(4,3), paste("n=", tapply(rbh.fpkm1$IR.F, rbh.fpkm1$Linkage, length), sep=""), cex = 0.75)
boxplot(IR.M ~ Linkage, rbh.fpkm2, boxwex = 0.25, outline = F,
        notch = T, col = c("grey60", "darkorange", "olivedrab"),
        add = T, at = c(0.75, 1.75, 2.75), axes = F,
        range =1, whisklty = 5, staplewex = 0.25, boxlty = 0)
text(c(0.75, 1.75, 2.75), rep(4.25,3), paste("n=", tapply(rbh.fpkm2$IR.M, rbh.fpkm2$Linkage, length), sep=""), cex = 0.75)
abline(h = 1, lty = 5, col = "red4")
axis(side = 1, at = seq(0.25, 2.75, by=0.5), labels = rep(c("F", "M"),3))
axis(side = 1, at = c(0.5, 1.5, 2.5), tick = F, labels = c("A", "Z-anc", "Z-neo"), 
     pos = -0.5,  cex.axis = 1.25)
dev.off()

rbind(with(rbh.fpkm1, tapply(IR.F, Linkage, median)), 
      with(rbh.fpkm2, tapply(IR.M, Linkage, median)))
A     Z-anc     Z-neo
[1,] 1.018771 0.7624407 0.8950314
[2,] 1.027424 0.7944142 1.0082098

mwu = function (i, df) {
  a = wilcox.test(df[[i]][df$Linkage=="A"], df[[i]][df$Linkage=="Z-anc"])
  b = wilcox.test(df[[i]][df$Linkage=="A"], df[[i]][df$Linkage=="Z-neo"])
  c = wilcox.test(df[[i]][df$Linkage=="Z-anc"], df[[i]][df$Linkage=="Z-neo"])
  print(c("A~Za:", a$p.value, "A~Zn:", b$p.value, "Za~Zn:", c$p.value))
}

# female
mwu(10, rbh.fpkm1)
[1] "A~Za:"                "1.26239767947249e-07" "A~Zn:"                "0.131547062649483"   
[5] "Za~Zn:"               "0.0100089002349225"  

pairwise.wilcox.test(rbh.fpkm1$IR.F, rbh.fpkm1$Linkage, p.adjust.method = "BH")
A       Z-anc
Z-anc 3.8e-07 -    
  Z-neo 0.132   0.015
pairwise.wilcox.test(rbh.fpkm1$IR.F, rbh.fpkm1$Linkage, p.adjust.method = "bonferroni")
A       Z-anc
Z-anc 3.8e-07 -    
  Z-neo 0.39    0.03 

# male
mwu(10, rbh.fpkm2)
[1] "A~Za:"                "2.76647265359157e-05" "A~Zn:"                "0.472103560878504"   
[5] "Za~Zn:"               "0.0154104657295191" 

pairwise.wilcox.test(rbh.fpkm2$IR.M, rbh.fpkm2$Linkage, p.adjust.method = "BH")
A       Z-anc
Z-anc 8.3e-05 -    
  Z-neo 0.472   0.023
pairwise.wilcox.test(rbh.fpkm2$IR.M, rbh.fpkm2$Linkage, p.adjust.method = "bonferroni")
A       Z-anc
Z-anc 8.3e-05 -    
  Z-neo 1.000   0.046

# Individual chromosome
library("naturalsort")

linkage = read.table('../chrom_linkage.txt',header=F,sep="")
colnames(linkage) = c("linkage", "chrom", "chrom1", "Gene_sp1")

rbh.fpkmF = merge(linkage, rbh.fpkm1, by = "Gene_sp1")
rbh.fpkmF$chrom1 = factor(rbh.fpkmF$chrom1, levels(rbh.fpkmF$chrom1)[c(31,30,naturalorder(levels(rbh.fpkmF$chrom1)[1:29]))])

rbh.fpkmM = merge(linkage, rbh.fpkm2, by = "Gene_sp1")
rbh.fpkmM$chrom1 = factor(rbh.fpkmM$chrom1, levels(rbh.fpkmM$chrom1)[c(31,30,naturalorder(levels(rbh.fpkmM$chrom1)[1:29]))])

pdf("Boxplot_by_chr_DpMs.pdf", width = 14, height = 12, pointsize = 12)
par(mfrow=c(2,1))
boxplot(log(IR.F,2) ~ chrom1, rbh.fpkmF, boxwex = 0.5, outline = F, notch = F, 
        col = c("olivedrab","darkorange",rep("grey60", 29)), xaxt = "n", las = 1,
        xlim = c(1, 31), ylim = c(-3.75, 5), cex.lab = 1.25,
        main = "D.plexippus:M.sexta in Female (WZ)", ylab = "interspecies FPKM ratio", 
        cex.lab = 1.1, range = 1, whisklty = 5, staplewex = 0.25, boxlty = 0)
axis(1, at = 1:31, labels = c("Zn","Za",2:30))
abline(h = 0, lty = 5, col = "red4")
text(1:31, c(4.75,4.25), paste("n=", tapply(rbh.fpkmF$IR.F, rbh.fpkmF$chrom1, length), sep=""), cex = 0.8)
boxplot(log(IR.M,2) ~ chrom1, rbh.fpkmM, boxwex = 0.5, outline = F, notch = F, 
        col = c("olivedrab","darkorange",rep("grey60", 29)), xaxt = "n", las = 1,
        xlim = c(1, 31), ylim = c(-3.75, 5), cex.lab = 1.25,
        main = "D.plexippus:M.sexta in Male (ZZ)", ylab = "interspecies FPKM ratio", 
        cex.lab = 1.1, range =1, whisklty = 5, staplewex = 0.25, boxlty = 0)
axis(1, at = 1:31, labels = c("Zn","Za",2:30))
abline(h = 0, lty = 5, col = "red4")
text(1:31, c(4.75,4.25), paste("n=", tapply(rbh.fpkmM$IR.M, rbh.fpkmM$chrom1, length), sep=""), cex = 0.8)
dev.off()

## Hm~Ms
df = rbh.fpkm0.HmMs
# female
wilcox.test(df[[10]][df$Chrom!="cZ"], df[[10]][df$Chrom=="cZ"])
W = 1248600, p-value = 0.08627
# male
wilcox.test(df[[11]][df$Chrom!="cZ"], df[[11]][df$Chrom=="cZ"])
W = 1135400, p-value = 0.1896

# upper 0.75 quantile
pdf("DpMs_ratio.pdf", width = 5, height = 8, pointsize = 14)
df = subset(rbh.fpkm1, IR.F>quantile(IR.F, 0.75))
boxplot(log(IR.F, 2) ~ Linkage, df, boxwex = 0.25, outline = F, 
        notch = T, col = c("grey", "darkorange", "olivedrab"),
        xlim = c(0, 3), xaxt = "n", main = "D.plexippus:M.sexta",
        ylab = expression("log"[2]*" (FPKM ratio)"), cex.lab = 1.1,
        ylim = c(1, 6), at = c(0.25, 1.25, 2.25))
text(c(0.25, 1.25, 2.25), rep(5.75,3), paste("n=", tapply(df$IR.F, df$Linkage, length), sep=""), cex = 0.75)

df = subset(rbh.fpkm2, IR.M>quantile(IR.M, 0.75))
boxplot(log(IR.M, 2) ~ Linkage, df, boxwex = 0.25, outline = F, 
        notch = T, col = c("grey", "darkorange", "olivedrab"),
        add = T, at = c(0.75, 1.75, 2.75), axes = F)
text(c(0.75, 1.75, 2.75), rep(5.5,3), paste("n=", tapply(df$IR.M, df$Linkage, length), sep=""), cex = 0.75)
abline(h = 0, lty = 2)
axis(side = 1, at = seq(0.25, 2.75, by=0.5), labels = rep(c("F", "M"),3))
axis(side = 1, at = c(0.5, 1.5, 2.5), tick = F, labels = c("A", "Z-anc", "Z-neo"), 
     pos = -6.5,  cex.axis = 1.25)
dev.off()

# continuous threshold test
quantG = function (q, i, df) {
  df = df[df[[i]]>quantile(df[[i]], q/100),]
  rA = df[[10]][df$Linkage=="A"]
  ra = df[[10]][df$Linkage=="Z-anc"]
  rn = df[[10]][df$Linkage=="Z-neo"]
  return(c(length(rA), length(ra), length(rn), median(rA), median(ra), median(rn)))
}

x = data.frame(matrix(ncol = 6, nrow = 0))
colnames(x) = c("obs_A", "obs_Za", "obs_Zn", "RA", "Ra", "Rn" )
F.quantG = x
M.quantG = x
for (q in 1:80) {
  F.quantG[q,] = quantG(q, 3, rbh.fpkm1)
  M.quantG[q,] = quantG(q, 4, rbh.fpkm2)
}

pdf("DpMS_ratio_quantG.pdf", width = 8, height = 8, pointsize = 20)
plot(F.quantG$RA,, type = "l", xlab = "Quantile cut-off threshold", ylab = "Dp:Ms ratio", 
     ylim = c(0.7, 1.6), main = "Global Quantile", lwd = 3, col = "black")
lines(F.quantG$Ra, lwd = 3, col = "darkorange")
lines(F.quantG$Rn, lwd = 3, col = "olivedrab")

lines(M.quantG$RA, lwd = 3, col = "black", lty = 2)
lines(M.quantG$Ra, lwd = 3, col = "darkorange4", lty = 2)
lines(M.quantG$Rn, lwd = 3, col = "darkblue", lty = 2)
legend (5, 1.3, c("Z-neo:female", "Z-anc:female", "Z-neo:male", "Z-anc:male"), 
        col = c("darkorange", "olivedrab", "darkorange4", "darkblue"), 
        bty = "n", cex = 0.75, lty = c(1,1,2,2), lwd = 3)
dev.off()

# only require one sex >0 in both species: result very similar
rbh.fpkm.f0 = subset(rbh.fpkm, HdF.sp1&HdF.sp2>0)
rbh.fpkm.f0$IR.F = with(rbh.fpkm.f0, HdF.sp1/HdF.sp2)

rbh.fpkm.m0 = subset(rbh.fpkm.f0, HdM.sp1&HdM.sp2>0)
rbh.fpkm.m0$IR.M = with(rbh.fpkm.m0, HdM.sp1/HdM.sp2)

boxplot(log(IR.F, 2) ~ Linkage, rbh.fpkm.f0, boxwex = 0.25, col = c("grey", "darkorange", "olivedrab"))
boxplot(log(IR.M, 2) ~ Linkage, rbh.fpkm.m0, boxwex = 0.25, col = c("grey", "darkorange", "olivedrab"),
        add = T, at = c(1.3, 2.3, 3.3), axes = F)

### Danaus
pdf("Dp_expression_head.pdf", width = 8, height = 8, pointsize = 10)
par(mfrow = c(1,2), mar=c(5.1,6,4.1,3))
boxplot(log(HdF, 2) ~ Linkage, dp_fpkm, subset = HdF>0, notch = T, boxwex = 0.25, 
        col = c("grey", "darkorange", "olivedrab"), xlim = c(0.35, 4), xaxt = "n", 
        main = "D. plexippus head (F|M > 0)", ylab = expression("log"[2]*" (FPKM)"))
boxplot(log(HdM, 2) ~ Linkage, dp_fpkm, subset = HdM>0, notch = T, boxwex = 0.25, 
        col = c("grey", "darkorange", "olivedrab"), at = c(1.4, 2.4, 3.4), add = T, axes = F)
# add male A median as reference
abline(h = log(median(with(subset(dp_fpkm, Linkage=="A"&HdM>0), HdM)), 2), lty = 2)
axis(side = 1, at = c(1, 1.4, 2, 2.4, 3, 3.4), labels = rep(c("F", "M"),3))
axis(side = 1, at = c(1.2, 2.2, 3.2), tick = F, labels = c("A", "Z-anc", "Z-neo"), 
     pos = -11,  cex.axis = 1.25)
# require both sex >0
boxplot(log(HdF, 2) ~ Linkage, dp_fpkm, subset = HdF&HdM>0, notch = T, boxwex = 0.25, 
        col = c("grey", "darkorange", "olivedrab"), xlim = c(0.35, 4), xaxt = "n", 
        main = "D. plexippus head (F&M >0)", ylab = expression("log"[2]*" (FPKM)"))
boxplot(log(HdM, 2) ~ Linkage, dp_fpkm, subset = HdF&HdM>0, notch = T, boxwex = 0.25, 
        col = c("grey", "darkorange", "olivedrab"), at = c(1.4, 2.4, 3.4), add = T, axes = F)
abline(h = log(median(with(subset(dp_fpkm, Linkage=="A"&HdF&HdM>0), HdM)), 2), lty = 2)
axis(side = 1, at = c(1, 1.4, 2, 2.4, 3, 3.4), labels = rep(c("F", "M"),3))
axis(side = 1, at = c(1.2, 2.2, 3.2), tick = F, labels = c("A", "Z-anc", "Z-neo"), 
     pos = -11,  cex.axis = 1.25)
dev.off()

### Ms~Hm
hm_fpkm = read.table("RNAseq_DC/heliconius.head.ave.tmm.txt", header = T)
hm_fpkm$Ratio = with(hm_fpkm, HD.F/HD.M)
#hm_fpkm = na.omit(hm_fpkm)
levels(hm_fpkm$Chrom)[levels(hm_fpkm$Chrom)!="cZ"] <- "A"

rbh_HmMs = read.table("rbh_HmMs.txt")

sp1.fpkm = hm_fpkm[match(rbh_HmMs$V1, hm_fpkm$HM_ID),c(2,1,3:5)]
sp2.fpkm = ms_fpkm[match(rbh_HmMs$V2, ms_fpkm$gene_id), c(1,10:12)]

rbh.fpkm.HmMs = cbind(sp1.fpkm, sp2.fpkm)

colnames(rbh.fpkm.HmMs)[c(1:5,6:9)] = c("Gene_sp1", "Linkage", "HdF.sp1", "HdM.sp1", "Ratio.sp1", 
                                     "Gene_Sp2", "HdF.sp2", "HdM.sp2", "Ratio.sp2")

rbh.fpkm.HmMs.f0 = subset(rbh.fpkm.HmMs, HdF.sp1&HdF.sp2>0)
rbh.fpkm.HmMs.m0 = subset(rbh.fpkm.HmMs, HdM.sp1&HdM.sp2>0)

# factor IR.F
x.IR.F = with(subset(rbh.fpkm.HmMs.f0, Linkage=="A"), median(HdF.sp1)/median(HdF.sp2))
# factor IR.M
x.IR.M = with(subset(rbh.fpkm.HmMs.m0, Linkage=="A"), median(HdM.sp1)/median(HdM.sp2))

rbh.fpkm.HmMs.f0$IR.F = with(rbh.fpkm.HmMs.f0, HdF.sp1/HdF.sp2/x.IR.F)
rbh.fpkm.HmMs.m0$IR.M = with(rbh.fpkm.HmMs.m0, HdM.sp1/HdM.sp2/x.IR.M)

levels(t$Chrom)[levels(t$Chrom)!="cZ"] <- "A"

t1=rbh.fpkm.HmMs.f0
t2=rbh.fpkm.HmMs.m0

pdf("HmMs_ratio_1.pdf", width = 4, height = 8, pointsize = 14)
boxplot(IR.F ~ Linkage, t1,  boxwex = 0.25, outline = F, notch = T, xaxt = "n", las = 1,
        col = c("grey60", "darkorange"), xlim = c(0.75, 2.65), ylim = c(0, 5.5), 
        main = "H.melpomene:M.sexta", cex.lab = 1.1,
        range = 1, whisklty = 5, staplewex = 0.25, boxlty = 0)
boxplot(IR.M ~ Linkage, t2, boxwex = 0.25, outline = F, notch = T, 
        col = c("grey60", "darkorange"), at = c(1.4, 2.4), add = T, axes = F,
        range = 1, whisklty = 5, staplewex = 0.25, boxlty = 0)
text(c(1, 2, 1.4, 2.4), c(5, 5, 5.5, 5.5), cex = 0.75,
     paste("n=", c(tapply(t1$IR.F, t1$Linkage, length), tapply(t2$IR.M, t2$Linkage, length)), sep=""))
abline(h = 1, lty = 5, col = "red4")
axis(side = 1, at = c(1, 1.4, 2, 2.4), labels = rep(c("F", "M"),2))
axis(side = 1, at = c(1.2, 2.2), tick = F, labels = c("A", "Z"), 
     pos = -0.6,  cex.axis = 1.25)
dev.off()

rbind(tapply(t1$IR.F, t1$Linkage, median), tapply(t2$IR.M, t2$Linkage, median))
A        cZ
[1,] 1.005109 0.9358769
[2,] 0.979488 1.0461904

wilcox.test(t1$IR.F[t1$Linkage=="A"], t1$IR.F[t1$Linkage=="cZ"])
wilcox.test(t2$IR.M[t1$Linkage=="A"], t2$IR.M[t1$Linkage=="cZ"])

pdf("HmMs_ratio_log2.pdf", width = 10, height = 8, pointsize = 15)
par(mfrow = c(1,3), mar = c(5.1,5,4.1,2.1))

boxplot(log(IR.F, 2) ~ Chrom, t,  boxwex = 0.25, outline = F, notch = T, 
        col = c("grey", "darkorange"), ylim = c(-5, 5),
        xlim = c(0.35, 3), xaxt = "n", main = "H.melpomene:M.sexta", 
        ylab = expression("log"[2]*" (FPKM ratio)"), cex.lab = 1.1)
boxplot(log(IR.M, 2) ~ Chrom, t, boxwex = 0.25, outline = F, notch = T, 
        col = c("grey", "darkorange"), at = c(1.4, 2.4), add = T, axes = F)
abline(h = 0, lty = 2)
axis(side = 1, at = c(1, 1.4, 2, 2.4), labels = rep(c("F", "M"),2))
axis(side = 1, at = c(1.2, 2.2), tick = F, labels = c("A", "Z"), 
     pos = -6,  cex.axis = 1.25)

boxplot(log(HdF.sp2, 2) ~ Chrom, t,  boxwex = 0.25, outline = F, notch = T, 
        col = c("grey", "darkorange", "olivedrab"), ylim = c(-4.5, 10),
        xlim = c(0.35, 3), xaxt = "n", main = "M. sexta",
        ylab = expression("log"[2]*" (FPKM)"), cex.lab = 1.1)
boxplot(log(HdM.sp2, 2) ~ Chrom, t, boxwex = 0.25, outline = F, 
        notch = T, col = c("grey", "darkorange", "olivedrab"),
        at = c(1.4, 2.4), add = T, axes = F)
abline(h = median(log(t$HdM.sp2[t$HdM.sp2>0],2)), lty = 2)
axis(side = 1, at = c(1, 1.4, 2, 2.4), labels = rep(c("F", "M"),2))
axis(side = 1, at = c(1.2, 2.2), tick = F, labels = c("A", "Z"), 
     pos = -6,  cex.axis = 1.25)

boxplot(log(HdF.sp1, 2) ~ Chrom, t,  boxwex = 0.25, outline = F, notch = T, 
        col = c("grey", "darkorange"), ylim = c(-4.3, 9),
        xlim = c(0.35, 3), xaxt = "n", main = "H. melpomene",
        ylab = expression("log"[2]*" (FPKM)"), cex.lab = 1.1)
boxplot(log(HdM.sp1, 2) ~ Chrom, t, boxwex = 0.25, outline = F, 
        notch = T, col = c("grey", "darkorange"),
        at = c(1.4, 2.4), add = T, axes = F)
abline(h = median(log(t$HdM.sp1[t$HdM.sp1>0],2)), lty = 2)
axis(side = 1, at = c(1, 1.4, 2, 2.4), labels = rep(c("F", "M"),2))
axis(side = 1, at = c(1.2, 2.2), tick = F, labels = c("A", "Z"), 
     pos = -5.7,  cex.axis = 1.25)

dev.off()
