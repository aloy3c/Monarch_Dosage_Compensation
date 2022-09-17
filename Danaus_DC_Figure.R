setwd("~/Documents/Danaus_RNAseq/Figure")

# Fig 1A
pdf("Fig_1A_Dapl_ScaffoldMeanRatios-Windows.pdf", width = 6, height = 5, pointsize = 18, useDingbats = FALSE)
par(mar=c(4.1,5.1,1.1,1.1))
with(subset(ratiosFrame, lengths > 50000), {
  plot(log10(lengths), ratios.windows, pch = 19, cex = 0.75, ylim = c(-1, 1.5), las = 1,
       ylab = expression("log"[2]*"(Mean M:F Read Counts)"), 
       xlab = expression("log"[10]*"(Scaffold Length (bp))"), col = "grey50")
  abline(h = 0, col = "black", lty = 5)
  abline(h = 1, col = "royalblue", lty = 5)
  points(log10(lengths)[1], ratios.windows[1], col = "royalblue", pch = 19, cex = 1.25)
  points(log10(lengths)[1], ratios.windows[1], col = zcol, pch = 1, cex = 2)
  text(7.15, 1.3, "Z", col = zcol, font = 2, cex = 1.25)
})
dev.off()

library("randomcoloR")
# Fig 1B
#chr.col = randomColor(30)
chr.col = distinctColorPalette(30)
chr.col1 = chr.col

# labels
dp.labs = c("1(Z)",2:30)
sl.labs = c(1:30, "31(Z)")

# Dp cordinates
ideo.dp = read.table("../# Danaus_HiRise/Dapl_Zhan_v3_HiC.chr.seq_length_dp.txt", header = F, sep = "", as.is = T)
a = 0
for (i in 1:30) {
  ideo.dp$Start[i] = a 
  ideo.dp$End[i] = ideo.dp$Start[i] + ideo.dp[i,2]
  a = ideo.dp$End[i] + 2000000
}
ideo.dp$col = chr.col

# Sl cordinates
ideo.sl = read.table("../# Danaus_HiRise/karyotype.sl.txt", header = F, sep = "", as.is = T)
sl.order = rev(c(3,13,19,10,4,26,12,24,18,8,5,11,28,1,17,23,6,25,22,7,20,14,21,27,9,15,16,29,2,31,30))
ideo.sl.r = ideo.sl[sl.order,c(3,6)]
a = 0
for (i in 1:31) {
  ideo.sl.r$Start[i] = a 
  ideo.sl.r$End[i] = ideo.sl.r$Start[i] + ideo.sl.r[i,2]
  a = ideo.sl.r$End[i] + 2000000
}

# Synteny
#ideo.sl$V6 = ideo.sl$V6*0.6
Sl.Dp = read.table("../# Danaus_HiRise/synteny.SlDp.0.8.txt", header = F, sep = "", as.is = T)
#plot(0, xlim = c(0, 325000000), ylim = c(0, 50), type = "n", yaxt = "n", xaxt = "n", xlab = "", ylab = "")

Sl.Dp_1 = subset(Sl.Dp, V4!="dp1")
Sl.Dp_Z = subset(Sl.Dp, V4=="dp1")

# re-position Dp: 
# Dp_space = (450000000-245185305-29*2000000)/2
Dp_space = 80000000
  
ideo.dp$Start_r = ideo.dp$Start + Dp_space
ideo.dp$End_r = ideo.dp$End + Dp_space

pdf("Fig_1B_Synteny_DpSl_r.pdf", width = 12, height = 5, pointsize = 16)
par(mar=c(0, 0, 0, 1.1))
plot(0, xlim = c(0, 450000000), ylim = c(10, 50), type = "n", axes = F, xlab = "", ylab = "")
#yaxt = "n", xaxt = "n"
segments(0+Dp_space, 40, 5685560+Dp_space, 40, lwd = 5, lend = 1, col = "olivedrab")
segments(5685561+Dp_space, 40, ideo.dp[1,7], 40, lwd = 5, lend = 1, col = "darkorange")
# ideograms:Dp
for (i in 2:30) {
  segments(ideo.dp[i,6], 40, ideo.dp[i,7], 40, lwd = 5, lend = 1, col = chr.col[i])
}
# labels
for (i in 1:30) {
  text((ideo.dp[i,6]+ideo.dp[i,7])/2, 42, dp.labs[i], col = "black", cex = 0.5)
}

# ideograms:Sl
for (i in 1:31) {
  segments(ideo.sl.r[i,3], 20, ideo.sl.r[i,4], 20, lwd = 5, lend = 1, col = "grey50")
}
for (i in 1:31) {
  text((ideo.sl.r[i,3]+ideo.sl.r[i,4])/2, 18, rownames(ideo.sl.r)[i], col = "black", cex = 0.5)
}

# synteny: autosomes
for (i in 1:nrow(Sl.Dp_1)) {
  segments(ideo.dp$Start_r[ideo.dp$V1==Sl.Dp_1[i,4]]+Sl.Dp_1[i,5], 39.5, 
         ideo.sl.r$Start[ideo.sl.r$V3==Sl.Dp_1[i,1]]+Sl.Dp_1[i,2], 20.5, lwd = 0.2, lend = 1, 
         col = ideo.dp$col[ideo.dp$V1==Sl.Dp_1[i,4]])
}
# synteny: Z
for (i in 1:nrow(Sl.Dp_Z)) {
  if (Sl.Dp_Z[i,5]<5685560) {
    b = "olivedrab"
  } else {
    b = "darkorange"
  }
  segments(ideo.dp$Start_r[ideo.dp$V1==Sl.Dp_Z[i,4]]+Sl.Dp_Z[i,5], 39.5, 
           ideo.sl.r$Start[ideo.sl.r$V3==Sl.Dp_Z[i,1]]+Sl.Dp_Z[i,2], 20.5, lwd = 0.2, lend = 1, col = b)
}

# species title
text(ideo.dp[30,4]/2 + Dp_space, 45, "Danaus plexippus", font = 3, cex = 0.8)
text(ideo.sl.r[31,4]/2, 15, "Spodoptera litura", font = 3, cex = 0.8)

dev.off()

### Scaled
# scaling factor
245185305/438940887
[1] 0.5585839
s_fac = 0.6
ideo.sl.r$Start_s[1] =  ideo.sl.r$Start[1]
ideo.sl.r$End_s[1] =  ideo.sl.r$End[1]*s_fac
for (i in 2:31) {
  ideo.sl.r$Start_s[i] =  ideo.sl.r$End_s[i-1] + 2000000
  ideo.sl.r$End_s[i] =  ideo.sl.r$Start_s[i] + ideo.sl.r[i,2]*s_fac
}

pdf("Fig_1B_Synteny_DpSl_s.pdf", width = 12, height = 5, pointsize = 16)
par(mar=c(0, 0, 0, 1.1))
plot(0, xlim = c(0, 300000000), ylim = c(10, 50), type = "n", axes = F, xlab = "", ylab = "")
segments(0, 40, 5685560, 40, lwd = 5, lend = 1, col = "olivedrab")
segments(5685561, 40, ideo.dp[1,4], 40, lwd = 5, lend = 1, col = "darkorange")
# ideograms:Dp
for (i in 2:30) {
  segments(ideo.dp[i,3], 40, ideo.dp[i,4], 40, lwd = 5, lend = 1, col = chr.col[i])
}
# labels
dp.labs = c("1(Z)",2:30)
for (i in 1:30) {
  text((ideo.dp[i,3]+ideo.dp[i,4])/2, 42, dp.labs[i], col = "black", cex = 0.5)
}

# ideograms:Sl
for (i in 1:31) {
  segments(ideo.sl.r[i,5], 20, ideo.sl.r[i,6], 20, lwd = 5, lend = 1, col = "grey50")
}
# labels
#sl.labs = c(1:30, "31(Z)")
for (i in 1:31) {
  text((ideo.sl.r[i,5]+ideo.sl.r[i,6])/2, 18, rownames(ideo.sl.r)[i], col = "black", cex = 0.5)
}

# synteny: autosomes
Sl.Dp_1 = subset(Sl.Dp, V4!="dp1")
for (i in 1:nrow(Sl.Dp_1)) {
  segments(ideo.dp$Start[ideo.dp$V1==Sl.Dp_1[i,4]]+Sl.Dp_1[i,5], 39.5, 
           ideo.sl.r$Start_s[ideo.sl.r$V3==Sl.Dp_1[i,1]]+Sl.Dp_1[i,2]*s_fac, 20.5, lwd = 0.2, lend = 1, 
           col = ideo.dp$col[ideo.dp$V1==Sl.Dp_1[i,4]])
}
# synteny: Z
Sl.Dp_Z = subset(Sl.Dp, V4=="dp1")
for (i in 1:nrow(Sl.Dp_Z)) {
  if (Sl.Dp_Z[i,5]<5685560) {
    b = "olivedrab"
  } else {
    b = "darkorange"
  }
  segments(ideo.dp$Start[ideo.dp$V1==Sl.Dp_Z[i,4]]+Sl.Dp_Z[i,5], 39.5, 
           ideo.sl.r$Start_s[ideo.sl.r$V3==Sl.Dp_Z[i,1]]+Sl.Dp_Z[i,2]*s_fac, 20.5, lwd = 0.2, lend = 1, col = b)
}

# species
text(ideo.dp[30,4]/2, 45, "Danaus plexippus", font = 3, cex = 0.8)
text(ideo.sl.r[31,6]/2, 15, "Spodoptera litura", font = 3, cex = 0.8)

dev.off()

# Fig 2A
library("vioplot")
pdf("Fig2A_Vioplot_by_linkage_0.01.pdf", width = 10, height = 7.5, pointsize = 16, useDingbats = FALSE)
#df1 = fpkm.mean[fpkm.mean$HdF>0.01,]
#df2 = fpkm.mean[fpkm.mean$HdM>0.01,]
#par(mfrow=c(1,2), mar=c(2.1,3.1,3.1,0), oma = c(2,2.5,2,2))
par(mfrow=c(1,2), mar=c(3.1,3.1,3.1,0), oma = c(2,2.5,1.5,2))
#par(mar = c(5,5,4,2) + 0.1, mfrow=c(1,2), mgp = c(2.5, 0.75, 0))
# female
plot(1, 1, xlim = c(0, 4), ylim = range(log(df1$HdF,2), log(df2$HdM,2)), 
     type = 'n', xlab = '', ylab = expression("log"[2]*"(FPKM)"), xaxt = 'n', main = "", las = 1)
mtext(expression("log"[2]*"(FPKM)"), side = 2, line = 3)
vioplot(log(df1$HdF[df1$linkage=="A"], 2), at = 0.75, add =T, border = "grey40", 
        rectCol = "grey40", col = "grey")
vioplot(log(df1$HdF[df1$linkage=="Z-anc"], 2), at = 2, add = T, border = "darkorange3", 
        rectCol = "darkorange3", col = "orange")
vioplot(log(df1$HdF[df1$linkage=="Z-neo"], 2), at = 3.25, add = T, border = "olivedrab", 
        rectCol = "olivedrab", col = "olivedrab3")
axis(3, c(0.75, 2, 3.25), paste("n=", tapply(df1$HdF, df1$linkage, length), sep=""), cex = 0.75, tick = F, line = -0.75)
axis(1, c(0.75, 2, 3.25), c("Aut", "anc-Z", "neo-Z"), font = 2)
abline(h=log(median(df1$HdF[df1$linkage=="A"]),2), lty = 5)
text(c(2, 3.25), c(11, 13), c("***", "p=0.94"), font = 2)
mtext("Female (WZ)", 3, line = 2, cex = 1.2, font = 2)
axis(1, c(0, 2, 3.25), c("median Z:A ratio", "0.57", "1.09"), tick = F, line = 1.5)
# male
plot(1, 1, xlim = c(0, 4), ylim = range(log(df1$HdF,2), log(df2$HdM,2)), 
     type = 'n', xlab = '', ylab = expression("log"[2]*"(FPKM)"), xaxt = 'n', main = "", yaxt = 'n')
axis(2, at = seq(-5, 15, by = 5), labels = F)
vioplot(log(df2$HdM[df2$linkage=="A"], 2), at = 0.75, add =T, border = "grey40", 
        rectCol = "grey40", col = "grey")
vioplot(log(df2$HdM[df2$linkage=="Z-anc"], 2), at = 2, add = T, border = "darkorange3", 
        rectCol = "darkorange3", col = "orange")
vioplot(log(df2$HdM[df2$linkage=="Z-neo"], 2), at = 3.25, add = T, border = "olivedrab", 
        rectCol = "olivedrab", col = "olivedrab3")
axis(3, c(0.75, 2, 3.25), paste("n=", tapply(df2$HdM, df2$linkage, length), sep=""), cex = 0.75, tick = F, line = -0.75)
axis(1, c(0.75, 2, 3.25), c("Aut", "anc-Z", "neo-Z"), font = 2)
abline(h=log(median(df2$HdM[df2$linkage=="A"]),2), lty = 5)
text(c(2, 3.25), c(11.5, 13), c("***", "p=0.43"), font = 2)
mtext("Male (ZZ)", 3, line = 2, cex = 1.2, font = 2)
axis(1, c(0, 2, 3.25), c("median Z:A ratio", "0.56", "1.10"), tick = F, line = 1.5)
dev.off()

# Fig 2B
pdf("Fig2B_FM_FPKM_scatterplot_TMM.pdf", width = 10, height = 5.5, pointsize = 16)
par(mar = c(5,5,4,2) + 0.1, mfrow=c(1,2), mgp = c(2.5, 0.75, 0))
#df3 = subset(fpkm.mean, HdF>0.01&HdM>0.01)
with(subset(df3, linkage=="A"), {
  smoothScatter(log(HdM,2), log(HdF,2), nbin = 256, 0.25, 
                colramp = colorRampPalette(c("white", "darkblue")), nrpoints = length(HdF)/100, 
                pch = 20, col = "darkgrey", cex = 0.5, xlab = expression("log"[2]*"(FPKM) \u2642"), 
                ylab = expression("log"[2]*"(FPKM) \u2640"), main = "Autosomes", 
                las = 1, xaxs = "i", yaxs = "i", cex.lab = 1, cex.main = 1)
  abline(lm(log(HdF,2)~log(HdM,2)))
  lm(log(HdF,2)~log(HdM,2))
})
#title("Autosomes", adj = 0, cex = 1.5)
with(subset(df3, linkage=="Z-anc"), {
  plot(log(HdM,2), log(HdF,2), pch = 20, col = adjustcolor("darkorange", alpha.f = 0.2),
       xlab = expression("log"[2]*"(FPKM) \u2642"), ylab = expression("log"[2]*"(FPKM) \u2640"),
       main = "Z", las = 1, cex.lab = 1, cex.main = 1)
  abline(lm(log(HdF,2)~log(HdM,2)), col = "darkorange", lwd = 2)
  lm(log(HdF,2)~log(HdM,2))
})
with(subset(df3, linkage=="Z-neo"), {
  points(log(HdM,2), log(HdF,2), pch = 20, col = adjustcolor("olivedrab", alpha.f = 0.2))
  abline(lm(log(HdF,2)~log(HdM,2)), col = "olivedrab", lwd = 2)
  lm(log(HdF,2)~log(HdM,2))
})
legend (-6, 10.5, c("anc-Z", "neo-Z"), col = c("darkorange", "olivedrab"), pch = 20, bty = "n")
#legend (4.5, -1.5, c("anc-Z", "neo-Z"), col = c("darkorange", "olivedrab"), pch = 20, bty = "n")
dev.off()

Coefficients:
(Intercept)  log(HdM, 2)  
-0.002805     0.987001  

Coefficients:
  (Intercept)  log(HdM, 2)  
-0.1785       1.0085  

Coefficients:
  (Intercept)  log(HdM, 2)  
-0.06745      0.98615  

# Fig 3A&B
#range(rbh.fpkm1$IR.F, rbh.fpkm2$IR.M)
pdf("Fig3_DpMs_HmMs_ratio.pdf", width = 8, height = 8, pointsize = 22)
layout(matrix(c(1,1,1,1,2,2,2,1,1,1,1,2,2,2), 2, 7, byrow = TRUE))
#layout.show(2)
par(mar=c(4.1,4.1,3.1,0), oma = c(1,1,1,2))
boxplot(IR.F ~ Linkage, rbh.fpkm1, boxwex = 0.25, outline = F, notch = T, 
        col = c("grey50", "darkorange", "olivedrab"), xaxt = "n", las = 1,
        xlim = c(0, 3), ylim = c(0, 4.5), main = "D. plexippus : M. sexta",
        ylab = "Interspecies FPKM ratio", cex.lab = 1.1, at = c(0.25, 1.25, 2.25), 
        range =1, whisklty = 5, staplewex = 0.25, boxlty = 0)
text(c(0.25, 1.25, 2.25), rep(4.25,3), paste("n=", tapply(rbh.fpkm1$IR.F, rbh.fpkm1$Linkage, length), sep=""), cex = 0.75)
#axis(3, c(0.75, 2, 3.25), paste("n=", tapply(df1$HdF, df1$linkage, length), sep=""), cex = 0.75, tick = F, line = -0.75)
boxplot(IR.M ~ Linkage, rbh.fpkm2, boxwex = 0.25, outline = F,
        notch = T, col = c("grey50", "darkorange", "olivedrab"),
        add = T, at = c(0.75, 1.75, 2.75), axes = F,
        range =1, whisklty = 5, staplewex = 0.25, boxlty = 0)
text(c(0.75, 1.75, 2.75), rep(4.5,3), paste("n=", tapply(rbh.fpkm2$IR.M, rbh.fpkm2$Linkage, length), sep=""), cex = 0.75)
text(seq(1.25, 2.75, by = 0.5), 3.75, c("***", "***", "ns", "ns"))
axis(side = 1, at = seq(0.25, 2.75, by=0.5), labels = rep(c("F", "M"),3))
axis(side = 1, at = c(0.5, 1.5, 2.5), tick = F, labels = c("Aut", "anc-Z", "neo-Z"), 
     pos = -0.6,  cex.axis = 1.25, font = 2)
abline(h = 1, lty = 5, col = "red4")

boxplot(IR.F ~ Linkage, t1,  boxwex = 0.25, outline = F, notch = T, xaxt = "n", las = 1,
        col = c("grey50", "darkorange"), xlim = c(0.75, 2.65), ylim = c(0, 5.75), 
        main = "H. melpomene : M.sexta", cex.lab = 1.1,
        range = 1, whisklty = 5, staplewex = 0.25, boxlty = 0)
boxplot(IR.M ~ Linkage, t2, boxwex = 0.25, outline = F, notch = T, 
        col = c("grey50", "darkorange"), at = c(1.4, 2.4), add = T, axes = F,
        range = 1, whisklty = 5, staplewex = 0.25, boxlty = 0)
text(c(1, 2, 1.4, 2.4), c(5.5, 5.5, 5.75, 5.75), cex = 0.75,
     paste("n=", c(tapply(t1$IR.F, t1$Linkage, length), tapply(t2$IR.M, t2$Linkage, length)), sep=""))
text(c(2, 2.4), 4.75, c("ns", "ns"))
axis(side = 1, at = c(1, 1.4, 2, 2.4), labels = rep(c("F", "M"),2))
axis(side = 1, at = c(1.2, 2.2), tick = F, labels = c("Aut", "Z"), 
     pos = -0.75,  cex.axis = 1.25, font = 2)
abline(h = 1, lty = 5, col = "red4")
dev.off()

pairwise.wilcox.test(rbh.fpkm1$IR.F, rbh.fpkm1$Linkage, p.adjust.method = "BH")
A       Z-anc
Z-anc 3.8e-07 -    
  Z-neo 0.132   0.015
pairwise.wilcox.test(rbh.fpkm1$IR.F, rbh.fpkm1$Linkage, p.adjust.method = "bonferroni")
A       Z-anc
Z-anc 3.8e-07 -    
  Z-neo 0.39    0.03 

pairwise.wilcox.test(rbh.fpkm2$IR.M, rbh.fpkm2$Linkage, p.adjust.method = "BH")
A       Z-anc
Z-anc 8.3e-05 -    
  Z-neo 0.472   0.023
pairwise.wilcox.test(rbh.fpkm2$IR.M, rbh.fpkm2$Linkage, p.adjust.method = "bonferroni")
A       Z-anc
Z-anc 8.3e-05 -    
  Z-neo 1.000   0.046


gdo = c("grey50", "darkorange", "olivedrab")
# Fig 4C
# previous use of "tF" is wrong!)

# t: all na.omit
pairwise.wilcox.test(t$Female, t$linkage, p.adjust.method = "BH")
data:  t$Female and t$linkage 

A       Z-anc  
Z-anc 0.44    -      
  Z-neo 2.9e-13 2.1e-07

pairwise.wilcox.test(t$Male, t$linkage, p.adjust.method = "BH")
data:  t$Male and t$linkage 

A       Z-anc  
Z-anc 4.2e-07 -      
  Z-neo 0.072   1.4e-06

# tss.fpkm.F and tss.fpkm.M are actually FPKM>0.01
pairwise.wilcox.test(tss.fpkm.F$Female, tss.fpkm.F$linkage, p.adjust.method = "none")
A       Z-anc  
Z-anc 0.48    -      
  Z-neo 1.4e-10 8.7e-06
P value adjustment method: none      

pairwise.wilcox.test(tss.fpkm.F$Female, tss.fpkm.F$linkage, p.adjust.method = "BH")
A       Z-anc  
Z-anc 0.48    -      
  Z-neo 4.2e-10 1.3e-05
P value adjustment method: BH 

pairwise.wilcox.test(tss.fpkm.F$Female, tss.fpkm.F$linkage, p.adjust.method = "bonferroni")
A       Z-anc  
Z-anc 1       -      
  Z-neo 4.2e-10 2.6e-05
P value adjustment method: bonferroni 

pairwise.wilcox.test(tss.fpkm.M$Male, tss.fpkm.M$linkage, p.adjust.method = "none")
A       Z-anc  
Z-anc 2.4e-08 -      
  Z-neo 0.57    2.2e-05
P value adjustment method: none 

pairwise.wilcox.test(tss.fpkm.M$Male, tss.fpkm.M$linkage, p.adjust.method = "BH")
A       Z-anc  
Z-anc 7.2e-08 -      
  Z-neo 0.57    3.3e-05
P value adjustment method: BH 

pairwise.wilcox.test(tss.fpkm.M$Male, tss.fpkm.M$linkage, p.adjust.method = "bonferroni")
A       Z-anc  
Z-anc 7.2e-08 -      
  Z-neo 1       6.5e-05
P value adjustment method: bonferroni 

pdf("Fig_4C_TSS_Boxplot_by_linkage_0.01.pdf", width = 8, height = 7, pointsize = 16)
par(mfrow=c(1, 2), mar=c(3.1, 3.1, 4.1, 0), oma = c(0, 2, 0, 2))

boxplot(Female~linkage, tss.fpkm.F, col = gdo, boxwex = 0.5, 
        main = "Female (WZ)", ylim = c(-1.3, 2.6), notch = T, xaxt = "n", las = 1, 
        outline = F, range = 1, whisklty = 5, staplewex = 0.5, boxlty = 0)
mtext(expression("log"[2]*"(ChIP:Input)"), side = 2, line = 3, cex = 1.25)
axis(1, at = 1:3, c("Aut","anc-Z","neo-Z"))
axis(3, at = 1:3, paste("n=", tapply(tss.fpkm.F$Female,tss.fpkm.F$linkage,length), sep=""), 
     cex.axis = 0.75, tick = F, line = -0.75)
text(c(2, 3), c(2.45, 2.6), c("p=0.48", "***"), font = c(1, 2))
abline(h = median(tss.fpkm.F$Female[tss.fpkm.F$linkage=="A"]), lty = 5, col = "red4")

boxplot(Male~linkage, tss.fpkm.M, col = gdo, boxwex = 0.5, yaxt = "n", 
        main = "Male (ZZ)", ylim = c(-1.3, 2.6), notch = T, xaxt = "n", las = 1,
        outline = F, range = 1, whisklty = 5, staplewex = 0.5, boxlty = 0)
axis(side = 2, at = seq(-1, 3.0, by = 1), labels = F, las =1)
axis(1, at = 1:3, c("Aut","anc-Z","neo-Z"))
axis(3, 1:3, paste("n=", tapply(tss.fpkm.M$Male,tss.fpkm.M$linkage,length), sep=""), 
     cex.axis = 0.75, tick = F, line = -0.75)
text(c(2, 3), c(2, 2.3), c("***", "p=0.57"), font = c(2, 1))
abline(h = median(tss.fpkm.M$Male[tss.fpkm.M$linkage=="A"]), lty = 5, col = "red4")

dev.off()

#####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")
browseVignettes("DECIPHER")
#####


##### Supplementary
# Fig.S3.A
pdf("Fig.S3.A_Fraction_non_expressed_by_chr.pdf", width = 12, height = 6, pointsize = 16)
par(mfrow = c(1, 2), mar=c(5.1,1.1,4.1,1.1), oma=c(0,6,0,1))

plot(0, type = "l", col = "olivedrab", lwd = 3, las =1,
     xlim = c(0, 500), ylim = c(0.05, 0.5), xaxt = "n", 
     ylab = "fraction of genes (FPKM < cut-off)",
     xlab = "FPKM cut-off", main = "Female (WZ)")
mtext("fraction of genes (FPKM < cut-off)", side = 2, line = 4, cex = 1.25)
axis(side = 1, at = seq(0, 500, by = 100), seq(0, 5, by = 1))
for (i in 3:31) { points(noexprFchr[,i], type = "l", col = "grey60") }
points(noexprF[,1], type = "l", col = "grey20", lwd = 3)
points(noexprFchr[,1], type = "l", col = "olivedrab", lwd = 3)
points(noexprFchr[,2], type = "l", col = "darkorange", lwd = 3)

plot(0, type = "l", col = "olivedrab", lwd = 3, las =1,
     xlim = c(0, 500), ylim = c(0.05, 0.5), xaxt = "n", yaxt = "n", 
     ylab = "", xlab = "FPKM cut-off", main = "Male (ZZ)")
axis(side = 1, at = seq(0, 500, by = 100), seq(0, 5, by = 1))
for (i in 3:31) { points(noexprFchr[,i], type = "l", col = "grey60") }
points(noexprM[,1], type = "l", col = "grey20", lwd = 3)
points(noexprMchr[,1], type = "l", col = "olivedrab", lwd = 3)
points(noexprMchr[,2], type = "l", col = "darkorange", lwd = 3)

legend (300, 0.2, c("aut", "aut mean", "anc-Z", "neo-Z"), 
        col = c("grey60", "grey20", "olivedrab", "darkorange"), 
        bty = "n", cex = 0.75, lwd = c(1.5,3,3,3))

dev.off()


# Fig.S3.B
pdf("Fig.S3.B_ZA_ratio_quantile_cutoff_gene_percentage_4in1.pdf", width = 12, height = 9, pointsize = 20)
par(mfrow = c(2, 2), mar=c(1.1,1.1,1.1,1.1), oma=c(4.5,5,1.5,6), mgp = c(3, 1, 0))
# p1
plot(F.quant$RaA, type = "l", col = "darkorange", xaxt = "n", las = 1, lwd = 4,
     xlab = "Quantile cut-off threshold", ylab = "Z:A ratio", ylim = c(0.5, 1.15))
lines(F.quant$RnA, lwd = 4, col = "olivedrab")
axis(1, labels = F)
axis(2, at = c(0.5, 1), las = 1, font = 2, lwd = 0)
grid(col = "grey60", lwd = 1.5)
abline(h = c(0.5, 1), lwd = 1.5, lty = 5)
mtext("Z:A ratio", side = 2, line = 4)
# p2
plot(M.quant$RaA, type = "l", col = "darkorange", xaxt = "n", las = 1, lwd = 4, yaxt = "n",
     xlab = "Quantile cut-off threshold", ylab = "Z:A ratio", ylim = c(0.5, 1.15), )
lines(M.quant$RnA, lwd = 4, col = "olivedrab")
axis(1, labels = F)
axis(2, labels = F) 
grid(col = "grey60", lwd = 1.5)
abline(h = c(0.5, 1), lwd = 1.5, lty = 5)
# p3 gene fraction
plot(F.quant$frac_A, type = "l", col = "grey20", ylim = c(15,100),
     ylab = "percentage genes remaining", las = 1, lwd = 3)
lines(F.quant$frac_Zn, lwd = 3, col = "olivedrab")
lines(F.quant$frac_Za, lwd = 3, col = "darkorange")
grid(col = "grey60", lwd = 1.5)
mtext("percentage genes remaining", side = 2, line = 4)
mtext("quantile cut-off in Female (WZ)", side = 1, line = 3)
# p4
plot(M.quant$frac_A, type = "l", col = "grey20", ylim = c(15,100),
     xlab = "", ylab = "", las = 1, lwd = 3, yaxt = "n")
lines(M.quant$frac_Za, lwd = 3, col = "darkorange")
lines(M.quant$frac_Zn, lwd = 3, col = "olivedrab")
axis(2, labels = F) 
grid(col = "grey60", lwd = 1.5)
mtext("quantile cut-off in Male (ZZ)", side = 1, line = 3)

legend(90, 130, c("aut", "anc-Z", "neo-Z"), col = c("grey20", "darkorange", "olivedrab"), 
       bty = "n", lwd = 4, seg.len = 1, xpd = NA, horiz = F, x.intersp = 0.75)

dev.off()

#title("Z:A ratio under sliding quantile cut-off", cex.main = 1.25, outer = T)

# Fig.S5
pdf("Fig.S5_plotProfile_FM416.420.bs50.chr.pdf", width = 10, height = 8, pointsize = 16)
par(mfrow=c(2,2), mar=c(4.1,2.1,2.1,0), oma = c(2,5,3,3))
# p1
plot(t(FM416.chr[1,3:222]), type = "l", xlim = c(0, 220), ylim = c(-0.4, 1.2), col = "olivedrab", lwd = 3, 
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "H4K16ac")
axis(2, at = seq(-0.4, 1.2, by = 0.4), tck = 0.02, las = 1)
axis(4, at = seq(-0.4, 1.2, by = 0.4), tck = 0.02, labels = F)
for (i in 3:31) {lines(t(FM416.chr[i,3:222]), type = "l", col = "grey60")}
lines(t(FM416.chr[2,3:222]), type = "l", col = "darkorange", lwd = 3)
abline(v = c("60","160"), lty = 5)
mtext(expression("log"[2]*" (ChIP:Input)"), side = 2, line = 4)
# p2
plot(0, type = "l", xlim = c(0, 220), ylim = c(-0.4, 1.2), col = "olivedrab", lwd = 3, 
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
axis(2, at = seq(-0.4, 1.2, by = 0.4), tck = 0.02, labels = F)
axis(4, at = seq(-0.4, 1.2, by = 0.4), tck = 0.02, labels = F)
for (i in 34:62) {lines(t(FM416.chr[i,3:222]), type = "l", col = "grey60")}
lines(t(FM416.chr[32,3:222]), type = "l", col = "olivedrab", lwd = 3)
lines(t(FM416.chr[33,3:222]), type = "l", col = "darkorange", lwd = 3)
abline(v = c("60","160"), lty = 5)
# p3
plot(0, type = "l", xlim = c(0, 220), ylim = c(0, 1.8), 
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", 
     xlab = "gene distance", ylab = "", main = "H4K20me1")
axis(1, at = c("0","60","160","220"), c("-3Kb","TSS","TES","+3Kb"), tck = 0, lwd = 0, lwd.ticks = 1, line = -0.5)
axis(2, at = seq(0, 1.8, by = 0.4), tck = 0.02, las = 1)
axis(4, at = seq(0, 1.8, by = 0.4), tck = 0.02, labels = F)
for (i in 3:31) {lines(t(FM420.chr[i,3:222]), type = "l", col = "grey60")}
lines(t(FM420.chr[1,3:222]), type = "l", col = "olivedrab", lwd = 3)
lines(t(FM420.chr[2,3:222]), type = "l", col = "darkorange", lwd = 3)
abline(v = c("60","160"), lty = 5)
mtext(expression("log"[2]*" (ChIP:Input)"), side = 2, line = 4)
# p4
plot(0, type = "l", xlim = c(0, 220), ylim = c(0, 1.8),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "gene distance",
     ylab = "", main = "")
axis(1, at = c("0","60","160","220"), c("-3Kb","TSS","TES","+3Kb"), tck = 0, lwd = 0, lwd.ticks = 1, line = -0.5)
axis(2, at = seq(0, 1.8, by = 0.4), tck = 0.02, labels = F)
axis(4, at = seq(0, 1.8, by = 0.4), tck = 0.02, labels = F)
for (i in 34:62) {lines(t(FM420.chr[i,3:222]), type = "l", col = "grey60")}
lines(t(FM420.chr[32,3:222]), type = "l", col = "olivedrab", lwd = 3)
lines(t(FM420.chr[33,3:222]), type = "l", col = "darkorange", lwd = 3)
abline(v = c("60","160"), lty = 5)

dev.off()


