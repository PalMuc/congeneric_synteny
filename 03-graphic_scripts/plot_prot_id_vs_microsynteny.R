# plot pairwise protein identity, and other stats
# by WRF 2023-08-02


id_clust_data = read.table("~/git/speciation_synteny/summary_data/identity_and_microsynteny_pairwise.tab", header =TRUE, sep = "\t")

pair_color_list = c("#a0e499aa", "#ac0a18aa", "#590aacaa", "#ec60bdaa", "#ec60bdaa", "#220a7eaa",
                    "#2ecd14aa", "#44cd14aa", "#cd9714aa", "#cd9714aa", "#cd9714aa", 
                    "#ac0a18aa", "#cd1f14aa",
                    "#0da730aa", "#a6ab09aa", "#ab099daa", "#9fcd1aaf",
                    "#098eabaa", "#e38601aa", "#e38601aa", "#e33601aa", "#0a4facaa",
                    "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa",
                    "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa" )
pair_color_list.no_alpha = substr(pair_color_list,1,7) # trim alpha

pair_color_list = c( "#a0e499aa", "#ac0a18aa", "#590aacaa", "#ec60bdaa", "#ec60bdaa", "#220a7eaa",
                     rep("#cd9714aa",8),
                     rep("#ab099daa",3),
                     rep("#098eabaa",5), 
                     rep("#666666aa",6), 
                     rep("#9fcd1aaf",6))
pair_color_list.no_alpha = substr(pair_color_list,1,7) # trim alpha


# plot all vs all
plot(id_clust_data[,3:23],
     pch = 16, col=pair_color_list, cex=2 )

pair_genus_names = c("Tethya", "Acropora", "Octopus", "Crassostrea", "Crassostrea", "Daphnia", 
                     "Anastrepha", "Culex", "Drosophila", "", "Drosophila", "Vespa", "Bombus", 
                     "Lytechinus", "Mauremys", "Bufo", "Cervus",
                     "Perca", "Epinephelus", "", "Oreochromis", "Thunnus",
                     "Insect", "Insect", "Insect", "",
                     "Placozoa", "Hexactinellida",
                     "","","","","","")
                     #"H-C", "H-G", "H-P", "H-N", "H-M", "H-O" )

pair_genus_names = c("Tethya", "Acropora", "Octopus", "Crassostrea", "", "Daphnia", 
                     "", "", "Drosophila", "", "Drosophila", "", "", 
                     "", "", "", "",
                     "", "", "", "", "Thunnus",
                     "Insect", "Insect", "Insect", "",
                     "Placozoa", "Hexactinellida",
                     "","","","Mammal","Mammal","Mammal")

#
pdf(file = "~/git/speciation_synteny/figures_for_paper/sp_pair_identity_vs_microsyn_v4.pdf", width = 5, height = 5)
par(mar=c(4.5,4.5,1,1))
plot(id_clust_data$prot_id, id_clust_data$t_gene / id_clust_data$t_total,
     xlim = c(0.45,1), ylim = c(0,1), frame.plot = FALSE,
     xlab = "Average pairwise protein identity",
     ylab = "Fraction of microsyntenic genes",
     #ylab = "Fraction of genes in microsyntenic blocks (3 or more)",
     pch = 16, col=pair_color_list, cex=3,
     cex.axis =1.3, cex.lab=1.3)
segments(0.5,0,0.84,1, col="#00000044", lty=2)
#abline(a = -1.470588 , b = 2.94, col="#00000044", lty=2)
#abline(a = -1.870588 , b = 2.94, col="#00000044", lty=2)
text(id_clust_data$prot_id, id_clust_data$t_gene / id_clust_data$t_total,
     pair_genus_names, pos = 2, font = c(rep(3,22),rep(1,12)), col = pair_color_list.no_alpha)
dev.off()

pdf(file = "~/git/speciation_synteny/figures_for_paper/sp_macrosyn_vs_microsyn_v4.pdf", width = 5, height = 5)
par(mar=c(4.5,4.5,1,1))
plot( id_clust_data$on_main_scaf / (id_clust_data$on_main_scaf + id_clust_data$off_main_scaf), 
      id_clust_data$t_gene / id_clust_data$t_total ,
      xlim = c(0,1), ylim = c(0,1), frame.plot = FALSE,
      xlab = "Fraction of macrosyntenic genes",
      ylab = "Fraction of microsyntenic genes",
      pch = 16, col=pair_color_list, cex=3,
      cex.axis =1.3, cex.lab=1.3)
text(id_clust_data$on_main_scaf / (id_clust_data$on_main_scaf + id_clust_data$off_main_scaf), 
     id_clust_data$t_gene / id_clust_data$t_total,
     pair_genus_names, pos = 2, font = c(rep(3,22),rep(1,12)), col = pair_color_list.no_alpha)
dev.off()

fraction_macro = id_clust_data$on_main_scaf / (id_clust_data$on_main_scaf + id_clust_data$off_main_scaf)
fraction_micro = id_clust_data$t_gene / id_clust_data$t_total



################################################################################

# define line plotting
get_protein_identity = function(filename){
  d = read.table(filename, sep="\t", stringsAsFactors = FALSE )
  id_pct_no_gap = d$V5/(d$V5 + d$V7)
  id_pct_no_gap.s = sort(id_pct_no_gap ,decreasing = TRUE)
  print( paste(basename(filename), mean(id_pct_no_gap.s) ) )
  return(id_pct_no_gap.s)
}
pair_genus_names_short = c("Tethya", "Acropora", "Octopus", 
                           "Crassostrea", "Crassostrea", "Crassostrea", "Daphnia", 
                           "Anastrepha", "Culex", "D.mel-D.ere", "D.mel-D.pse", "D.mel-D.gri", 
                           "Vespa", "Bombus", 
                           "Lytechinus", "Mauremys", "Bufo", "Cervus",
                           "Perca", "Epinephelus", "Epinephelus", "Oreochromis", "Thunnus" )
# sorted curve
pdf(file = "~/git/speciation_synteny/figures_for_paper/sp_pair_protein_identity_v2.pdf", width = 6, height = 5)
par(mar=c(1,4.5,4.5,1))
plot( 0,0, type = 'n', xlim=c(0,14000), ylim=c(0.2,1), frame.plot=FALSE, axes=FALSE,
      xlab="", ylab="Percent identity between orthologous proteins" )
axis(2)
axis(3)
mtext("Protein ranked from highest identity to lowest",3,line = 2)
for (i in seq(1,length(id_file_list),1) ){
  sp_id = get_protein_identity(paste0("~/git/speciation_synteny/06-prot_id_tables/",id_file_list[i]) )
  lines( 1:length(sp_id), sp_id , lwd=4, col=pair_color_list[i] )
  text( length(sp_id), min(sp_id), pair_genus_names_short[i], col=pair_color_list.no_alpha[i] , 
        srt = 45, pos = 2, font=3, )
}
dev.off()


