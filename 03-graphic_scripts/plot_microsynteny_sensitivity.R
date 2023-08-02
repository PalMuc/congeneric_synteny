# plot length parameter vs synteny blocks
# using microsynteny.py script at:
# https://github.com/wrf/genomeGTFtools
# 2023-07-31 by WRF


################################################################################
# check for intergenic length sensitivity
# with parameter -z in microsynteny.py
z_length_data = read.table("~/git/speciation_synteny/summary_data/microsynteny_z_length_dependence.tab", header = TRUE, sep = "\t")

obi_genome_size = 2342534168
talb_genome_size = 792101331
tmac_genome_size = 782423979
bufo_genome_size = 4545482719
# magic number is 1 / 2e-4, so 1/5000

genus = z_length_data$genus
genome_sizes = rep(c( bufo_genome_size, obi_genome_size, talb_genome_size ), table(genus) )
sp_colors = rep(c( "#ab099daa", "#590aacaa", "#0a4facaa" ), table(genus) )
sp_points = rep(c(15, 16, 17), table(genus) )

# make double plot for supplemental fig
pdf(file = "~/git/speciation_synteny/supplements_for_paper/microsynteny_z_length_dependence.pdf", width = 8, height = 10, paper = "a4")
par(mfrow=c(2,1), mar=c(4.5,4.5,2,1))
plot(z_length_data$z/genome_sizes, z_length_data$s_blocks, 
     xlab = "Length to next gene (fraction of genome length)", ylab = "Number of genes in blocks of 3 or more", log = "x", 
     cex = 2, pch = sp_points, col = sp_colors,
     axes = FALSE, cex.lab=1.4 )
axis(1, cex.axis=0.8 )
axis(2, cex.axis=1.3 )
legend("topleft", legend = unique(genus), cex = 1.2,
       pch = c(15,16,17), col = c( "#ab099daa", "#590aacaa", "#0a4facaa" ))
segments(2e-4,0,2e-4,25000, lwd = 2, lty = 2, col = "#0000007e")
mtext("A",at = 1e-6, cex = 2)
plot(z_length_data$z/genome_sizes, z_length_data$n_blocks, 
     xlab = "Length to next gene (fraction of genome length)", ylab = "Number of blocks", log = "x", 
     cex = 2, pch = sp_points, col = sp_colors,
     axes = FALSE, cex.lab=1.4 )
axis(1, cex.axis=0.8 )
axis(2, cex.axis=1.3 )
legend("topright", legend = unique(genus), cex = 1.2,
       pch = c(15,16,17), col = c( "#ab099daa", "#590aacaa", "#0a4facaa" ))
mtext("B",at = 1e-6, cex = 2)
dev.off()


################################################################################
# check for intervening gene sensitivity
# with parameter -s in microsynteny.py
s_length_data = read.table("~/git/speciation_synteny/summary_data/microsynteny_skip_dependence.tab", header = TRUE, sep = "\t")


pdf(file = "~/git/speciation_synteny/supplements_for_paper/microsynteny_s_skip_dependence.pdf", width = 8, height = 8, paper = "a4")
par(mfrow=c(2,2), mar=c(4.5,4.5,2,1))
plot(s_length_data$s, s_length_data$q_blocks/s_length_data$q_total,
     ylim=c(0,1), axes = FALSE, cex.lab=1.4,
     xlab="Allowed intervening genes (-s)", ylab="Fraction of total genes in blocks",
     pch = ifelse( s_length_data$m==2, 15, 17), cex=ifelse( s_length_data$z==200000, 1.3 , 2 ), 
     col = ifelse( s_length_data$random=="y", ifelse( s_length_data$m==2, "#e2bfafe7", "#ac260aaa" ), ifelse( s_length_data$m==2, "#afc4e2ee", "#0a4facaa" ) ) )
# legend("right", legend = c("m = 2", "m = 3", "random, m = 2", "random, m = 3"), cex = 1.2,
#        pch = c(15,17), col = c( "#afc4e2ee", "#0a4facaa" , "#e2bfafe7", "#ac260aaa"))
axis(1, at = c(1,2,3,4,5,6,10), cex.axis=1.3 )
axis(2, cex.axis=1.3 )
legend(1,0.85, legend = c("m = 2", "m = 3", "random, m = 2", "random, m = 3"), cex = 1.2,
       pch = c(15,17), col = c( "#afc4e2ee", "#0a4facaa" , "#e2bfafe7", "#ac260aaa"))
legend(6.5,0.85, legend = c("z = 200kb", "z = 500kb"), cex = 1.2,
       pch = c(17), pt.cex = c(1.3,2), col = "black")
mtext("A",at = 0, cex = 2)
#
plot(s_length_data$s, s_length_data$n_blocks,
     axes = FALSE, cex.lab=1.4,
     xlab="Allowed intervening genes (-s)", ylab="Total number of blocks (>= m genes)",
     pch = ifelse( s_length_data$m==2, 15, 17), cex=ifelse( s_length_data$z==200000, 1.3 , 2 ), 
     col = ifelse( s_length_data$random=="y", ifelse( s_length_data$m==2, "#e2bfafe7", "#ac260aaa" ), ifelse( s_length_data$m==2, "#afc4e2ee", "#0a4facaa" ) ) )
axis(1, at = c(1,2,3,4,5,6,10), cex.axis=1.3 )
axis(2, cex.axis=1.3 )
mtext("B",at = 0, cex = 2)
#
plot(s_length_data$s, s_length_data$longest,
     axes = FALSE, cex.lab=1.4, ylim=c(0,1300),
     xlab="Allowed intervening genes (-s)", ylab="Length of longest block (n genes)",
     pch = ifelse( s_length_data$m==2, 15, 17), cex=ifelse( s_length_data$z==200000, 1.3 , 2 ), 
     col = ifelse( s_length_data$random=="y", ifelse( s_length_data$m==2, "#e2bfafe7", "#ac260aaa" ), ifelse( s_length_data$m==2, "#afc4e2ee", "#0a4facaa" ) ) )
# legend("right", legend = c("m = 2", "m = 3", "random, m = 2", "random, m = 3"), cex = 1.2,
#        pch = c(15,17), col = c( "#afc4e2ee", "#0a4facaa" , "#e2bfafe7", "#ac260aaa"))
axis(1, at = c(1,2,3,4,5,6,10), cex.axis=1.3 )
axis(2, cex.axis=1.3 )
mtext("C",at = 0, cex = 2)
#
plot(s_length_data$s, s_length_data$mean_block,
     axes = FALSE, cex.lab=1.4, ylim=c(0,180),
     xlab="Allowed intervening genes (-s)", ylab="Average block length (n genes)",
     pch = ifelse( s_length_data$m==2, 15, 17), cex=ifelse( s_length_data$z==200000, 1.3 , 2 ), 
     col = ifelse( s_length_data$random=="y", ifelse( s_length_data$m==2, "#e2bfafe7", "#ac260aaa" ), ifelse( s_length_data$m==2, "#afc4e2ee", "#0a4facaa" ) ) )
# legend("right", legend = c("m = 2", "m = 3", "random, m = 2", "random, m = 3"), cex = 1.2,
#        pch = c(15,17), col = c( "#afc4e2ee", "#0a4facaa" , "#e2bfafe7", "#ac260aaa"))
axis(1, at = c(1,2,3,4,5,6,10), cex.axis=1.3 )
axis(2, cex.axis=1.3 )
mtext("D",at = 0, cex = 2)
dev.off()


################################################################################
# plot intergenic distance

intergenic_data1 = read.table("~/git/speciation_synteny/summary_data/GCF_910596095.1_fThuMac1.1_genomic.intergenic.tab.gz", header = FALSE, sep = "\t")
intergenic_data2 = read.table("~/git/speciation_synteny/summary_data/GCF_914725855.1_fThuAlb1.1_genomic.intergenic.tab.gz", header = FALSE, sep = "\t")
ig_distance1 = sort(intergenic_data1$V5)
ig_distance2 = sort(intergenic_data2$V5)
ig_h1 = hist(ig_distance1, breaks = 2000)
ig_h2 = hist(ig_distance2, breaks = 2000)

#table(ig_distance1 > tmac_genome_size/5000)
# FALSE  TRUE 
# 29621   152

intergenic_data3 = read.table("~/git/speciation_synteny/summary_data/GCF_020631705.1_ASM2063170v1.1_genomic.intergenic.tab.gz", header = FALSE, sep = "\t")
intergenic_data4 = read.table("~/git/speciation_synteny/summary_data/GCF_021134715.1_ASM2113471v1_genomic.intergenic.tab.gz", header = FALSE, sep = "\t")
ig_distance3 = sort(intergenic_data3$V5)
ig_distance4 = sort(intergenic_data4$V5)
ig_h3 = hist(ig_distance3, breaks = 1000)
ig_h4 = hist(ig_distance4, breaks = 1000)
dmagna_genome_size = 161467056
dpulex_genome_size = 133196385

intergenic_data5 = read.table("~/git/speciation_synteny/summary_data/GCF_001194135.2_ASM119413v2_genomic.intergenic.tab.gz", header = FALSE, sep = "\t")
intergenic_data6 = read.table("~/git/speciation_synteny/summary_data/GCF_006345805.1_ASM634580v1_genomic.intergenic.tab.gz", header = FALSE, sep = "\t")
ig_distance5 = sort(intergenic_data5$V5)
ig_distance6 = sort(intergenic_data6$V5)
ig_h5 = hist(ig_distance5, breaks = 10000)
ig_h6 = hist(ig_distance6, breaks = 10000)
osin_genome_size = 2719151902

intergenic_data7 = read.table("~/git/speciation_synteny/summary_data/GCF_025612915.1_ASM2561291v2_genomic.intergenic.tab.gz", header = FALSE, sep = "\t")
intergenic_data8 = read.table("~/git/speciation_synteny/summary_data/GCF_002022765.2_C_virginica-3.0_genomic.intergenic.tab.gz", header = FALSE, sep = "\t")
ig_distance7 = sort(intergenic_data7$V5)
ig_distance8 = sort(intergenic_data8$V5)
ig_h7 = hist(ig_distance7, breaks = 2000)
ig_h8 = hist(ig_distance8, breaks = 2000)
cang_genome_size = 624273838
cvir_genome_size = 684741128

pdf(file = "~/git/speciation_synteny/supplements_for_paper/sp_intergenic_length_histogram.pdf", width = 8, height = 9, paper = "a4")
par(mfrow=c(2,2) )
# daphnia
plot( ig_h4$mids/dmagna_genome_size, ig_h4$counts , type="n", log="x",
      xlab="Intergenic distance / genome size (log)",
      ylab="Number of introns",
      cex.lab=1.4, cex.axis=1.0)
segments(ig_h3$mids/dmagna_genome_size, rep(0,length(ig_h3$counts)), 
         ig_h3$mids/dmagna_genome_size, ig_h3$counts, lwd=3, col="#220a7eaa")
segments(ig_h4$mids/dpulex_genome_size*1.1, rep(0,length(ig_h4$counts)), 
         ig_h4$mids/dpulex_genome_size*1.1, ig_h4$counts, lwd=3, col="#520aaeaa")
legend("topright", legend = c("D. magna", "D. pulex"), lwd = 5, col = c("#220a7eaa","#520aaeaa"), text.font = 3)
points(2e-4,500, pch=6, lwd=4)
text(2e-4,1200, "z-length\ncutoff\n26kb")
points(mean(ig_distance3)/dmagna_genome_size,1000, pch=6, col="#220a7eaa", lwd=4)
points(mean(ig_distance4)/dpulex_genome_size,600, pch=6, col="#520aaeaa", lwd=4)
text(2.3e-5,1700, "mean\n3.2kb")
# crassostrea
plot( ig_h8$mids/cvir_genome_size, ig_h8$counts , type="n", log="x",
      xlab="Intergenic distance / genome size (log)",
      ylab="Number of introns",
      cex.lab=1.4, cex.axis=1.0)
segments(ig_h7$mids/cang_genome_size, rep(0,length(ig_h7$counts)), 
         ig_h7$mids/cang_genome_size, ig_h7$counts, lwd=3, col="#ec60bdaa")
segments(ig_h6$mids/cvir_genome_size*1.1, rep(0,length(ig_h6$counts)), 
         ig_h6$mids/cvir_genome_size*1.1, ig_h6$counts, lwd=3, col="#bc206daa")
legend("topright", legend = c("C. angulata", "C. virginica"), lwd = 5, col = c("#ec60bdaa","#bc206daa"), text.font = 3)
points(2e-4,200, pch=6, lwd=4)
text(2e-4,600, "z-length\ncutoff\n136kb")
points(mean(ig_distance7)/cang_genome_size,600, pch=6, col="#ec60bdaa", lwd=4)
points(mean(ig_distance8)/cvir_genome_size,700, pch=6, col="#bc206daa", lwd=4)
text(1.2e-5,1100, "mean\n8.9kb")
# octopus
plot( ig_h6$mids/osin_genome_size, ig_h6$counts , type="n", log="x",
      xlab="Intergenic distance / genome size (log)",
      ylab="Number of introns",
      cex.lab=1.4, cex.axis=1.0)
segments(ig_h5$mids/obi_genome_size, rep(0,length(ig_h5$counts)), 
         ig_h5$mids/obi_genome_size, ig_h5$counts, lwd=3, col="#590aacaa")
segments(ig_h6$mids/osin_genome_size*1.1, rep(0,length(ig_h6$counts)), 
         ig_h6$mids/osin_genome_size*1.1, ig_h6$counts, lwd=3, col="#a20aeeaa")
legend("topright", legend = c("O. bimaculoides", "O. sinensis"), lwd = 5, col = c("#590aacaa","#a20aeeaa"), text.font = 3)
points(2e-4,200, pch=6, lwd=4)
text(2e-4,500, "z-length\ncutoff\n506kb")
points(mean(ig_distance5)/obi_genome_size,300, pch=6, col="#590aacaa", lwd=4)
points(mean(ig_distance6)/osin_genome_size,400, pch=6, col="#a20aeeaa", lwd=4)
text(3e-5,700, "mean\n67kb")
# tuna
plot( ig_h1$mids/talb_genome_size, ig_h1$counts , type="n", log="x",
      xlab="Intergenic distance / genome size (log)",
      ylab="Number of introns",
      cex.lab=1.4, cex.axis=1.0)
segments(ig_h1$mids/tmac_genome_size, rep(0,length(ig_h1$counts)), 
         ig_h1$mids/tmac_genome_size, ig_h1$counts, lwd=3, col="#0a4facaa")
segments(ig_h2$mids/talb_genome_size*1.1, rep(0,length(ig_h2$counts)), 
         ig_h2$mids/talb_genome_size*1.1, ig_h2$counts, lwd=3, col="#3a4fecaa")
legend("topright", legend = c("T. maccoyii", "T. albacares"), lwd = 5, col = c("#0a4facaa","#3a4fecaa"), text.font = 3)
points(2e-4,400, pch=6, lwd=4)
text(2e-4,800, "z-length\ncutoff\n156kb")
points(mean(ig_distance1)/tmac_genome_size,800, pch=6, col="#0a4facaa", lwd=4)
points(mean(ig_distance2)/talb_genome_size,600, pch=6, col="#3a4fecaa", lwd=4)
text(1.3e-5,1200, "mean\n10.2kb")
dev.off()


#