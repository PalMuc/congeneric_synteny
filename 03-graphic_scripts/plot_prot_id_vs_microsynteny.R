# plot pairwise protein identity, and other stats
# by WRF 2023-08-02


id_clust_data = read.table("~/git/speciation_synteny/summary_data/identity_and_microsynteny_pairwise.tab", header =TRUE, sep = "\t")
id_clust_data = id_clust_data[2:dim(id_clust_data)[1],]

pair_color_list = c( #"#a0e499aa", 
                     "#ac0a18aa", "#590aacaa", "#ec60bdaa", "#ec60bdaa", "#220a7eaa",
                     rep("#cd9714aa",8),
                     rep("#ab099daa",3),
                     rep("#098eabaa",5), 
                     rep("#666666aa",4), rep("#9999eeaa",2), 
                     rep("#9fcd1aaf",6))
pair_color_list.no_alpha = substr(pair_color_list,1,7) # trim alpha


# plot all vs all
plot(id_clust_data[,3:23],
     pch = 16, col=pair_color_list, cex=2 )

pair_genus_names = c(#"Tethya", 
                     "Acropora", "Octopus", "Crassostrea", "Crassostrea", "Daphnia", 
                     "Anastrepha", "Culex", "Drosophila", "", "Drosophila", "Vespa", "Bombus", 
                     "Lytechinus", "Mauremys", "Bufo", "Cervus",
                     "Perca", "Epinephelus", "", "Oreochromis", "Thunnus",
                     "Insect", "Insect", "Insect", "",
                     "Placozoa", "Hexactinellida",
                     "","","","","","")
                     #"H-C", "H-G", "H-P", "H-N", "H-M", "H-O" )

pair_genus_names = c(#"Tethya", 
                     "Acropora", "Octopus", "Crassostrea", "", "Daphnia", 
                     "", "", "Drosophila", "", "Drosophila", "", "", 
                     "", "", "", "",
                     "", "", "", "", "Thunnus",
                     "Insect", "Insect", "Insect", "",
                     "Placozoa", "Hexactinellida",
                     "","","","Mammal","Mammal","Mammal")

fraction_macro = id_clust_data$on_main_scaf / (id_clust_data$on_main_scaf + id_clust_data$off_main_scaf)
fraction_micro = id_clust_data$t_gene / id_clust_data$t_total

#
pdf(file = "~/git/speciation_synteny/figures_for_paper/sp_pair_identity_vs_microsyn_v5.pdf", width = 5, height = 5, useDingbats = FALSE)
par(mar=c(4.5,4.5,1,1))
plot(id_clust_data$prot_id, fraction_micro,
     xlim = c(0.45,1), ylim = c(0,1), frame.plot = FALSE,
     xlab = "Average pairwise protein identity",
     ylab = "Fraction of microsyntenic genes",
     #ylab = "Fraction of genes in microsyntenic blocks (3 or more)",
     pch = 16, col=pair_color_list, cex=3,
     cex.axis =1.3, cex.lab=1.3)
segments(0.5,0,0.84,1, col="#00000044", lty=2)
#abline(a = -1.470588 , b = 2.94, col="#00000044", lty=2)
#abline(a = -1.870588 , b = 2.94, col="#00000044", lty=2)
text(id_clust_data$prot_id, fraction_micro,
     pair_genus_names, pos = 2, font = c(rep(3,22),rep(1,12)), col = pair_color_list.no_alpha)
dev.off()

pdf(file = "~/git/speciation_synteny/figures_for_paper/sp_macrosyn_vs_microsyn_v5.pdf", width = 5, height = 5, useDingbats = FALSE)
par(mar=c(4.5,4.5,1,1))
plot( fraction_macro, fraction_micro ,
      xlim = c(0,1), ylim = c(0,1), frame.plot = FALSE,
      xlab = "Fraction of macrosyntenic genes",
      ylab = "Fraction of microsyntenic genes",
      pch = 16, col=pair_color_list, cex=3,
      cex.axis =1.3, cex.lab=1.3)
text(fraction_macro, fraction_micro,
     pair_genus_names, pos = 2, font = c(rep(3,22),rep(1,12)), col = pair_color_list.no_alpha)
dev.off()

pdf(file = "~/git/speciation_synteny/figures_for_paper/sp_pair_identity_vs_macrosyn_v5.pdf", width = 5, height = 5, useDingbats = FALSE)
par(mar=c(4.5,4.5,1,1))
plot(id_clust_data$prot_id, fraction_macro,
     xlim = c(0.45,1), ylim = c(0,1), frame.plot = FALSE,
     xlab = "Average pairwise protein identity",
     ylab = "Fraction of macrosyntenic genes",
     pch = 16, col=pair_color_list, cex=3,
     cex.axis =1.3, cex.lab=1.3)
text(id_clust_data$prot_id, fraction_macro,
     pair_genus_names, pos = 2, font = c(rep(3,22),rep(1,12)), col = pair_color_list.no_alpha)
dev.off()

################################################################################

# supplemental figure showing divergence time
pdf(file = "~/git/speciation_synteny/supplements_for_paper/sp_pair_identity_vs_div_time_v2.pdf", width = 10, height = 5)
par(mfrow=c(1,2), mar=c(4.5,4.5,3,1))
plot(id_clust_data$prot_id, id_clust_data$div_time,
     xlim = c(0.5,1), ylim = c(350,0), frame.plot = FALSE,
     xlab = "Average pairwise protein identity",
     ylab = "Estimated divergence time (Ma)",
     #ylab = "Fraction of genes in microsyntenic blocks (3 or more)",
     pch = 16, col=pair_color_list, cex=3,
     cex.axis =1.3, cex.lab=1.3)
text(id_clust_data$prot_id, id_clust_data$div_time,
     pair_genus_names, pos = 2, font = c(rep(3,22),rep(1,12)), col = pair_color_list.no_alpha)
mtext("A", side = 3, line = 1, at = 0.4, cex = 2)
plot(fraction_micro, id_clust_data$div_time,
     xlim = c(0,1), ylim = c(350,0), frame.plot = FALSE,
     ylab = "Estimated divergence time (Ma)",
     xlab = "Fraction of macrosyntentic genes",
     pch = 16, col=pair_color_list, cex=3,
     cex.axis =1.3, cex.lab=1.3)
text(fraction_micro, id_clust_data$div_time,
     pair_genus_names, pos = c(rep(2,22),rep(4,13)), font = c(rep(3,22),rep(1,12)), col = pair_color_list.no_alpha)
mtext("B", side = 3, line = 1, at = -0.1, cex = 2)
dev.off()


################################################################################

# define line plotting
get_protein_identity = function(filename){
  d = read.table(filename, sep="\t", stringsAsFactors = FALSE )
  id_pct_no_gap = d$V5/(d$V5 + d$V7)
  id_pct_no_gap.s = sort(id_pct_no_gap ,decreasing = TRUE)
  print( paste(basename(filename), mean(id_pct_no_gap.s) ) )
  return(id_pct_no_gap.s)
}

pair_color_list = c("#a0e499aa", "#ac0a18aa", "#590aacaa", "#ec60bdaa", "#ec60bdaa","#ec60bdaa", "#220a7eaa",
                    "#2ecd14aa", "#44cd14aa", "#cd9714aa", "#cd9714aa", "#cd9714aa", 
                    "#ac0a18aa", "#cd1f14aa",
                    "#0da730aa", "#a6ab09aa", "#ab099daa", "#9fcd1aaf",
                    "#098eabaa", "#e38601aa", "#e38601aa", "#e33601aa", "#0a4facaa",
                    "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa",
                    "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa" )
pair_color_list.no_alpha = substr(pair_color_list,1,7) # trim alpha


id_file_list_dir = dir( "~/git/speciation_synteny/06-prot_id_tables/", "*.homologs_identity.tab.gz", recursive = TRUE)
#id_file_list = id_file_list_dir[c(33,1,31,6,8,5,17,2,9,12,14,13,  27,4,28,29,3,7,32,18,19,30,20)]
id_file_list = c("Tmi_V4b17_hintsutr_vs_TwiV4_AUG.homologs_identity.tab.gz", 
                 "AhyaV1_vs_Amilv2.homologs_identity.tab.gz", 
                 "Obi_vs_Osi.homologs_identity.tab.gz", 
                 "Cang_vs_Cvir3.homologs_identity.tab.gz", 
                 "Cgig1_vs_Cvir3.homologs_identity.tab.gz", 
                 "Cang_vs_Cgig1.homologs_identity.tab.gz", 
                 "Dpul_vs_Dmag.homologs_identity.tab.gz", 
                 "AnaObli1_vs_AnaLude1.homologs_identity.tab.gz", 
                 "Cqui1_vs_CPPV2.homologs_identity.tab.gz", 
                 "Dmel6_vs_DereRS2.homologs_identity.tab.gz", 
                 "Dmel6_vs_DpseMV25.homologs_identity.tab.gz", 
                 "Dmel6_vs_Dgrim.homologs_identity.tab.gz", 
                 "iyVesCrab1_vs_iVesVel2.homologs_identity.tab.gz", 
                 "Bpyro_vs_iyBomTerr1.homologs_identity.tab.gz", 
                 "Lvar3_vs_Lpic2.homologs_identity.tab.gz", 
                 "Mreeves_vs_Mmutica.homologs_identity.tab.gz", 
                 "Bgarg_vs_aBufBuf1.homologs_identity.tab.gz", 
                 "Ccan_vs_mCerEla1.homologs_identity.tab.gz", 
                 "PFLA1_vs_Pfluv1.homologs_identity.tab.gz", 
                 "Efusco_vs_Elanceo.homologs_identity.tab.gz", 
                 "Efusco_vs_Emoara.homologs_identity.tab.gz", 
                 "OaurZZ_vs_OnilUMD.homologs_identity.tab.gz", 
                 "fThuAlb1_vs_fThuMac1.homologs_identity.tab.gz")

pair_genus_names_short = c("Tethya", "Acropora", "Octopus", 
                           "Crassostrea", "Crassostrea", "Crassostrea", "Daphnia", 
                           "Anastrepha", "Culex", "D.mel-D.ere", "D.mel-D.pse", "D.mel-D.gri", 
                           "Vespa", "Bombus", 
                           "Lytechinus", "Mauremys", "Bufo", "Cervus",
                           "Perca", "Epinephelus", "Epinephelus", "Oreochromis", "Thunnus" )
# sorted curve
pdf(file = "~/git/speciation_synteny/figures_for_paper/sp_pair_protein_identity_v3.pdf", width = 6, height = 5)
par(mar=c(1,4.5,4.5,1))
plot( 0,0, type = 'n', xlim=c(0,14000), ylim=c(0.2,1), frame.plot=FALSE, axes=FALSE,
      xlab="", ylab="Percent identity between orthologous proteins" )
axis(2, at=seq(0.2,1,0.2), labels = seq(20,100,20))
axis(3)
mtext("Proteins ranked from highest identity to lowest identity",3,line = 2)
for (i in seq(1,length(id_file_list),1) ){
  sp_id = get_protein_identity(paste0("~/git/speciation_synteny/06-prot_id_tables/",id_file_list[i]) )
  lines( 1:length(sp_id), sp_id , lwd=4, col=pair_color_list[i] )
  text( length(sp_id), min(sp_id), pair_genus_names_short[i], col=pair_color_list.no_alpha[i] , 
        srt = 45, pos = 2, font=3, )
}
dev.off()


