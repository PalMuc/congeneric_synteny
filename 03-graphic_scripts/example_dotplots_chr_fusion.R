# make dot plot of synteny between two genomes, based on unidirectional blast hits (i.e. not reciprocal)
# last modified 2023-07-27

#pdf(file = "~/git/speciation_synteny/figures_for_paper/combined_dotplots_v1.pdf", width = 12, height = 12)
#png(file = "~/git/speciation_synteny/figures_for_paper/combined_dotplots_v1.png", width = 12, height = 12, units = "in", res = 90)
#par(mfrow = c(2,2))
all2Dfile = "~/git/speciation_synteny/04-macrosynteny_plots/Dmel6_vs_DereRS2.scaffold_synteny.tab.gz"
all2Ddata = read.table(all2Dfile, sep="\t", stringsAsFactors=FALSE)
categories = all2Ddata[,1]
is_scaf1 = which(categories=="s1")
scafdata1 = all2Ddata[is_scaf1,]
longestscaf1 = max(scafdata1[,6])
is_longscafs1 = which(as.numeric(scafdata1[,5]) > 0.02)
longscafs1 = c(0, scafdata1[,6][is_longscafs1] )
longscafs1_names = scafdata1[is_longscafs1,2]
is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
longestscaf2 = max(scafdata2[,6])
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.02)
scafdata2_long = scafdata2[is_longscafs2,]
longscafs2 = c(0, scafdata2[,6][is_longscafs2] )
length(longscafs2)
longscafs2_names = scafdata2_long[,2]
is_points = which(categories=="g")
pointsdata = all2Ddata[is_points,]
is_both_longscaf = !is.na(match(pointsdata[,3], longscafs1_names)) & !is.na(match(pointsdata[,5], longscafs2_names))
pointsdata_long = pointsdata[is_both_longscaf,]

genome_x = pointsdata_long[,6]
genome_y = pointsdata_long[,7]
xmax = tail( pretty(longscafs1), n=1)
ymax = tail( pretty(longscafs2), n=1)
longscafs_x = longscafs1
longscafs_y = longscafs2
nscafs_x = length(longscafs1)
nscafs_y = length(longscafs2)
xlab = genome1_lab
ylab = genome2_lab
xmax_mb = round(xmax / 1000000)
ymax_mb = round(ymax / 1000000)
dotcolor = "#cd9714aa" # color blue

# make PDF
# outputfile = "~/git/speciation_synteny/figures_for_paper/drosophila_dotplot_v1.pdf"
# pdf(file=outputfile, width=6, height=6) # a5 size
# outputfile = "~/git/speciation_synteny/figures_for_paper/drosophila_dotplot_v1.png"
# png(file=outputfile, width=6, height=6, units = "in", res = 90) # a5 size
par( mar=c(4.5,4.5,1,1) )
plot(genome_x, genome_y, pch=16, 
     xlim = c(0,135000000), ylim = c(0,135000000),
     col=dotcolor, cex=0.5, cex.lab=1.4, 
     main="", xlab="D. melanogaster", ylab="D. erecta", font.lab=3,
     axes=FALSE )
mtext("A",at = max(longscafs_x)*-0.08, cex = 2, line = -1)
axis(1, at=longscafs_x, labels=round(longscafs_x/1000000,ifelse(longscafs_x<100000000, 1, 0)), cex.axis=0.6, gap.axis = -1 )
segments(longscafs_x, 0, longscafs_x, longscafs_y[nscafs_y], lwd=0.3, col="#00000088")
axis(2, at=longscafs_y, labels=round(longscafs_y/1000000,ifelse(longscafs_y<100000000, 1, 0)), cex.axis=0.6, line=0, las=1, gap.axis = -1 )
segments( 0, longscafs_y, longscafs_x[nscafs_x], longscafs_y, lwd=0.8, col="#00000088")
rect(longscafs_x[3], longscafs_y[6], longscafs_x[3]+6e6, longscafs_y[7], border = "#0000ea88", lwd = 3)
rect(longscafs_x[5]+16e6, longscafs_y[3], longscafs_x[6], longscafs_y[5], border = "#0000ea88", lwd = 3)
# dev.off()

#

all2Dfile = "~/git/speciation_synteny/04-macrosynteny_plots/Obi_vs_Osi.scaffold_synteny.tab.gz"
all2Ddata = read.table(all2Dfile, sep="\t", stringsAsFactors=FALSE)
categories = all2Ddata[,1]
is_scaf1 = which(categories=="s1")
scafdata1 = all2Ddata[is_scaf1,]
longestscaf1 = max(scafdata1[,6])
is_longscafs1 = which(as.numeric(scafdata1[,5]) > 0.005)
longscafs1 = c(0, scafdata1[,6][is_longscafs1] )
longscafs1_names = scafdata1[is_longscafs1,2]
is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
longestscaf2 = max(scafdata2[,6])
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.005)
scafdata2_long = scafdata2[is_longscafs2,]
longscafs2 = c(0, scafdata2[,6][is_longscafs2] )
length(longscafs2)
longscafs2_names = scafdata2_long[,2]
is_points = which(categories=="g")
pointsdata = all2Ddata[is_points,]
is_both_longscaf = !is.na(match(pointsdata[,3], longscafs1_names)) & !is.na(match(pointsdata[,5], longscafs2_names))
pointsdata_long = pointsdata[is_both_longscaf,]

genome_x = pointsdata_long[,6]
genome_y = pointsdata_long[,7]
xmax = tail( pretty(longscafs1), n=1)
ymax = tail( pretty(longscafs2), n=1)
longscafs_x = longscafs1
longscafs_y = longscafs2
nscafs_x = length(longscafs1)
nscafs_y = length(longscafs2)
xlab = genome1_lab
ylab = genome2_lab
xmax_mb = round(xmax / 1000000)
ymax_mb = round(ymax / 1000000)
dotcolor = "#590aacaa" # color blue

# make PDF
# outputfile = "~/git/speciation_synteny/figures_for_paper/octopus_dotplot_v1.pdf"
# pdf(file=outputfile, width=6, height=6) # a5 size
# outputfile = "~/git/speciation_synteny/figures_for_paper/octopus_dotplot_v1.png"
# png(file=outputfile, width=6, height=6, units = "in", res = 90) # a5 size
par( mar=c(4.5,4.5,1,1) )
plot(genome_x, genome_y, pch=16, 
     xlim = c(0,2350000000), ylim = c(0,2350000000),
     col=dotcolor, cex=0.5, cex.lab=1.4, 
     main="", xlab="O. bimaculoides", ylab="O. sinensis", font.lab=3,
     axes=FALSE )
mtext("B",at = max(longscafs_x)*-0.08, cex = 2, line = -1)
axis(1, at=longscafs_x, labels=round(longscafs_x/1000000,ifelse(longscafs_x<100000000, 1, 0)), cex.axis=0.6, gap.axis = -1, las=2 )
segments(longscafs_x, 0, longscafs_x, longscafs_y[nscafs_y], lwd=0.3, col="#00000088")
axis(2, at=longscafs_y, labels=round(longscafs_y/1000000,ifelse(longscafs_y<100000000, 1, 0)), cex.axis=0.6, line=0, las=1, gap.axis = -1 )
segments( 0, longscafs_y, longscafs_x[nscafs_x], longscafs_y, lwd=0.8, col="#00000088")
#rect(longscafs_x[8], longscafs_y[9], longscafs_x[9], longscafs_y[10], border = "#00000088", lwd = 3)
# dev.off()

#


all2Dfile = "~/git/speciation_synteny/04-macrosynteny_plots/Dpul_vs_Dmag.scaffold_synteny.tab.gz"

all2Ddata = read.table(all2Dfile, sep="\t", stringsAsFactors=FALSE)
categories = all2Ddata[,1]
is_scaf1 = which(categories=="s1")
scafdata1 = all2Ddata[is_scaf1,]
longestscaf1 = max(scafdata1[,6])
is_longscafs1 = which(as.numeric(scafdata1[,5]) > 0.0009)
longscafs1 = c(0, scafdata1[,6][is_longscafs1] )
longscafs1_names = scafdata1[is_longscafs1,2]
is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
longestscaf2 = max(scafdata2[,6])
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.02)
scafdata2_long = scafdata2[is_longscafs2,]
longscafs2 = c(0, scafdata2[,6][is_longscafs2] )
length(longscafs2)
longscafs2_names = scafdata2_long[,2]
is_points = which(categories=="g")
pointsdata = all2Ddata[is_points,]
is_both_longscaf = !is.na(match(pointsdata[,3], longscafs1_names)) & !is.na(match(pointsdata[,5], longscafs2_names))
pointsdata_long = pointsdata[is_both_longscaf,]

genome_x = pointsdata_long[,6]
genome_y = pointsdata_long[,7]
xmax = tail( pretty(longscafs1), n=1)
ymax = tail( pretty(longscafs2), n=1)
longscafs_x = longscafs1
longscafs_y = longscafs2
nscafs_x = length(longscafs1)
nscafs_y = length(longscafs2)
xlab = genome1_lab
ylab = genome2_lab
xmax_mb = round(xmax / 1000000)
ymax_mb = round(ymax / 1000000)
dotcolor = "#220abe88" # color blue

# make PDF
# outputfile = "~/git/speciation_synteny/figures_for_paper/daphnia_dotplot_v1.pdf"
# pdf(file=outputfile, width=6, height=6) # a5 size
# outputfile = "~/git/speciation_synteny/figures_for_paper/daphnia_dotplot_v1.png"
# png(file=outputfile, width=6, height=6, units = "in", res = 90) # a5 size
par( mar=c(4.5,4.5,1,1) )
plot(genome_x, genome_y, pch=16, 
     xlim = c(0,135000000), ylim = c(0,135000000),
     col=dotcolor, cex=0.5, cex.lab=1.4, 
     main="", xlab="D. pulex", ylab="D. magna", font.lab=3 ,
     axes=FALSE )
mtext("C",at = max(longscafs_x)*-0.08, cex = 2, line = -1)
axis(1, at=longscafs_x, labels=round(longscafs_x/1000000,ifelse(longscafs_x<100000000, 1, 0)), cex.axis=0.6, gap.axis = -1 )
segments(longscafs_x, 0, longscafs_x, longscafs_y[nscafs_y], lwd=0.3, col="#00000088")
axis(2, at=longscafs_y, labels=round(longscafs_y/1000000,1), cex.axis=0.6, line=0, las=1)
segments( 0, longscafs_y, longscafs_x[nscafs_x], longscafs_y, lwd=0.8, col="#00000088")
rect(longscafs_x[10], longscafs_y[2], longscafs_x[11], longscafs_y[3], border = "#ea000088", lwd = 3)
rect(longscafs_x[11], longscafs_y[2], longscafs_x[12], longscafs_y[3], border = "#ea000088", lwd = 3)
rect(longscafs_x[9], longscafs_y[1], longscafs_x[10], longscafs_y[2], border = "#ea000088", lwd = 3)
rect(longscafs_x[12], longscafs_y[1], longscafs_x[13], longscafs_y[2], border = "#ea000088", lwd = 3)
# dev.off()

#

all2Dfile = "~/git/speciation_synteny/04-macrosynteny_plots/Cang_vs_Cvir3.scaffold_synteny.tab.gz"
all2Ddata = read.table(all2Dfile, sep="\t", stringsAsFactors=FALSE)
categories = all2Ddata[,1]
is_scaf1 = which(categories=="s1")
scafdata1 = all2Ddata[is_scaf1,]
longestscaf1 = max(scafdata1[,6])
is_longscafs1 = which(as.numeric(scafdata1[,5]) > 0.02)
longscafs1 = c(0, scafdata1[,6][is_longscafs1] )
longscafs1_names = scafdata1[is_longscafs1,2]
is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
longestscaf2 = max(scafdata2[,6])
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.02)
scafdata2_long = scafdata2[is_longscafs2,]
longscafs2 = c(0, scafdata2[,6][is_longscafs2] )
length(longscafs2)
longscafs2_names = scafdata2_long[,2]
is_points = which(categories=="g")
pointsdata = all2Ddata[is_points,]
is_both_longscaf = !is.na(match(pointsdata[,3], longscafs1_names)) & !is.na(match(pointsdata[,5], longscafs2_names))
pointsdata_long = pointsdata[is_both_longscaf,]

genome_x = pointsdata_long[,6]
genome_y = pointsdata_long[,7]
xmax = tail( pretty(longscafs1), n=1)
ymax = tail( pretty(longscafs2), n=1)
longscafs_x = longscafs1
longscafs_y = longscafs2
nscafs_x = length(longscafs1)
nscafs_y = length(longscafs2)
xlab = genome1_lab
ylab = genome2_lab
xmax_mb = round(xmax / 1000000)
ymax_mb = round(ymax / 1000000)
dotcolor = "#fc60bd88" # color blue

# make PDF
# outputfile = "~/git/speciation_synteny/figures_for_paper/crassostrea_dotplot_v1.pdf"
# pdf(file=outputfile, width=6, height=6) # a5 size
# outputfile = "~/git/speciation_synteny/figures_for_paper/crassostrea_dotplot_v1.png"
# png(file=outputfile, width=6, height=6, units = "in", res = 90) # a5 size
par( mar=c(4.5,4.5,1,1) )
plot(genome_x, genome_y, pch=16, 
     xlim = c(0,700000000), ylim = c(0,700000000),
     col=dotcolor, cex=0.5, cex.lab=1.4, 
     main="", xlab="C. angulata", ylab="C. virginica", font.lab=3 ,
     axes=FALSE )
mtext("D",at = max(longscafs_x)*-0.08, cex = 2, line = -1)
axis(1, at=longscafs_x, labels=round(longscafs_x/1000000,ifelse(longscafs_x<100000000, 1, 0)), cex.axis=0.6, gap.axis = -1 )
segments(longscafs_x, 0, longscafs_x, longscafs_y[nscafs_y], lwd=0.3, col="#00000088")
axis(2, at=longscafs_y, labels=round(longscafs_y/1000000,ifelse(longscafs_y<100000000, 1, 0)), cex.axis=0.6, line=0, las=1)
segments( 0, longscafs_y, longscafs_x[nscafs_x], longscafs_y, lwd=0.8, col="#00000088")
rect(longscafs_x[3], longscafs_y[2], longscafs_x[4], longscafs_y[3], border = "#ea000088", lwd = 3)
rect(longscafs_x[9], longscafs_y[2], longscafs_x[10], longscafs_y[3], border = "#ea000088", lwd = 3)
rect(longscafs_x[1], longscafs_y[1], longscafs_x[2], longscafs_y[2], border = "#0000ea88", lwd = 3)
rect(longscafs_x[1], longscafs_y[10], longscafs_x[2], longscafs_y[11], border = "#0000ea88", lwd = 3)
dev.off()


#
