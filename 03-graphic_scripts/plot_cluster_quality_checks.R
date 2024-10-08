# plot cluster quality checks
# examining which species contribute more to gaps in the alignment
# WRF 2023-07-26

###
# This script generates only supplemental figures
###

homolog_id_cols = c("target", "taxon", "gaps", "gap_pct", 
                    "identities", "id_pct", "differences", "diff_pct",  
                    "ids_diffs", "q_gaps", "s_gaps")

id_file_list_dir = dir( "~/git/speciation_synteny/06-prot_id_tables/", "*.homologs_identity.tab.gz", recursive = TRUE)
#id_file_list = id_file_list_dir[c(33,1,31,6,8,5,17,2,9,12,14,13,  27,4,28,29,3,7,32,18,19,30,20)]
id_file_list = c(#"Tmi_V4b17_hintsutr_vs_TwiV4_AUG.homologs_identity.tab.gz", 
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
pair_genus_names = c(#"Tethya wilhelma-minuta", 
                     "Acropora hyacinthus-millepora", "Octopus bimaculoides-sinensis", 
                     "Crassostrea angulata-virginica", "Crassostrea gigas-virginica", "Crassostrea angulata-gigas",
                     "Daphnia pulex-magna", "Anastrepha obliqua-ludens", "Culex quinquefasciatus-pipiens", 
                     "Drosophila melanogaster-erecta", "Drosophila melanogaster-pseudoobscura", "Drosophila melanogaster-grimshawi", 
                     "Vespa crabro-velutina", "Bombus pyrosoma-terrestris", "Lytechinus variegatus-picta", 
                     "Mauremys reevesii-mutica", 
                     "Bufo gargarizans-bufo", "Cervus hanglu-elaphus",
                     "Perca flavescens-fluviatilis", "Epinephelus fuscoguttatus-lanceolatum", "Epinephelus fuscoguttatus-moara", 
                     "Oreochromis aureus-niloticus", "Thunnus maccoyii-albacares" )

pair_color_list = c(#"#a0e499aa", 
                    "#ac0a18aa", "#590aacaa", "#ec60bdaa", "#ec60bdaa", "#ec60bdaa", "#220a7eaa",
                    "#2ecd14aa", "#44cd14aa", "#cd9714aa", "#cd9714aa", "#cd9714aa", 
                    "#ac0a18aa", "#cd1f14aa",
                    "#0da730aa", "#a6ab09aa", "#ab099daa", "#9fcd1aaf",
                    "#098eabaa", "#e38601aa", "#e38601aa", "#e33601aa", "#0a4facaa",
                    "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa")
pair_color_list.no_alpha = substr(pair_color_list,1,7) # trim alpha
letter_list = unlist(strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ",""))

pdf(file = "~/git/speciation_synteny/supplements_for_paper/cluster_alignment_lengths_v2.pdf", width = 8, height = 10, paper = "a4")
par(mfrow=c(4,2), mar=c(4.5,4.5,4,1) )
for (i in 1:length(id_file_list) ){
  cluster_id = read.table( paste0("~/git/speciation_synteny/06-prot_id_tables/",id_file_list[i]), 
                           header = FALSE, sep = "\t", col.names = homolog_id_cols, stringsAsFactors = FALSE)

  sp1_header = sapply(strsplit(cluster_id$target[1],"|", fixed = TRUE), getElement, 1)
  sp2_header = sapply(strsplit(cluster_id$taxon[1],"|", fixed = TRUE), getElement, 1)
  is_q = grepl(sp1_header, cluster_id$target)
  is_s = grepl(sp1_header, cluster_id$taxon)
  
  protein_sizes = cluster_id$identities + cluster_id$differences
  total_proteome = sum(protein_sizes)
  total_q_gaps = sum(c(cluster_id$q_gaps[is_q],cluster_id$s_gaps[is_s]))
  total_s_gaps = sum(c(cluster_id$q_gaps[is_s],cluster_id$s_gaps[is_q]))
  all_gaps = sum(cluster_id$gaps)
  #all_gaps
  all_gaps - total_q_gaps - total_s_gaps
  bv = c( total_proteome + all_gaps,  sum( cluster_id$identities ), sum( cluster_id$differences ), 
         total_q_gaps , total_s_gaps )
  
  b = barplot( bv , col = pair_color_list[i], main = pair_genus_names[i],
           ylab = "Total amino acids",
           cex.names = 0.9, cex.axis = 1.2, cex.lab = 1.2, font.main=4,
           names.arg = c("Alignment\nlength", "Identities", "Differences", paste0(sp1_header,"\ngaps"), paste0(sp2_header,"\ngaps") ) )
  text(b[2:5], bv[2:5] + bv[1] * c(-0.05,0.05,0.05,0.05), round(bv[2:5]/bv[1],3) )
  mtext(letter_list[i], side = 3, line = 2, at = -0.5, cex = 2)
}
dev.off()


################################################################################

pdf(file = "~/git/speciation_synteny/supplements_for_paper/sp_percent_id_histograms_v3.pdf", width = 8, height = 10, paper = "a4")
par(mfrow=c(4,2), mar=c(4.5,4.5,4,1) )
for (i in 1:length(id_file_list) ){
  cluster_id = read.table( paste0("~/git/speciation_synteny/06-prot_id_tables/",id_file_list[i]), 
                           header = FALSE, sep = "\t", col.names = homolog_id_cols, stringsAsFactors = FALSE)
  id_pct_no_gap = 100*round(sort( cluster_id$identities/(cluster_id$identities + cluster_id$differences) , decreasing = TRUE), digits = 3)
  protein_id_table = table(id_pct_no_gap)
  protein_id_mean = mean(id_pct_no_gap)
  plot( names(protein_id_table), protein_id_table, type = 'h', 
        xlim = c(100,30), main = pair_genus_names[i], font.main=4,
        xlab = "Percent protein identity", ylab = "N protein pairs",
        lwd=3, col = pair_color_list[i],
        axes = FALSE , cex.lab=1.4  )
  axis(1,cex.axis = 1.3)
  axis(2,cex.axis = 1.3)
  segments(protein_id_mean, 0, protein_id_mean, max(protein_id_table), lwd = 2, lty=2, col = "#00000066")
  text(protein_id_mean, max(protein_id_table) * 0.9, paste("Mean:",round(protein_id_mean, digits = 3)), pos=4 )
  mtext(letter_list[i], side = 3, line = 2, at = 111, cex = 2)
}
dev.off()

cras_dmel_numbers = c(6,4,5,10,11,12)
id_file_subset = id_file_list[cras_dmel_numbers]
pair_genus_subset = c( "Crassostrea angulata-gigas", "Crassostrea angulata-virginica", "Crassostrea gigas-virginica",
                      "Drosophila melanogaster-\nerecta", "Drosophila melanogaster-\npseudoobscura", "Drosophila melanogaster-\ngrimshawi" )

pdf(file = "~/git/speciation_synteny/figures_for_paper/cras_dmel_percent_id_histograms_v2.pdf", width = 8, height = 6 )
#png(file = "~/git/speciation_synteny/figures_for_paper/cras_dmel_percent_id_histograms_v2.png", width = 720, height = 540 , res = 90)
par(mfrow=c(2,3), mar=c(4.5,4.5,4,1) )
for (i in 1:length(id_file_subset) ){
  cluster_id = read.table( paste0("~/git/speciation_synteny/06-prot_id_tables/",id_file_subset[i]), 
                           header = FALSE, sep = "\t", col.names = homolog_id_cols, stringsAsFactors = FALSE)
  id_pct_no_gap = 100*round(sort( cluster_id$identities/(cluster_id$identities + cluster_id$differences) , decreasing = TRUE), digits = 3)
  protein_id_table = table(id_pct_no_gap)
  protein_id_mean = mean(id_pct_no_gap)
  plot( names(protein_id_table), protein_id_table, type = 'h', 
        xlim = c(100,30), main = pair_genus_subset[i], font.main=4,
        xlab = "Percent protein identity", ylab = "N protein pairs",
        lwd=3, col = pair_color_list[cras_dmel_numbers][i],
        axes = FALSE , cex.lab=1.4  )
  axis(1,cex.axis = 1.3)
  axis(2,cex.axis = 1.3)
  segments(protein_id_mean, 0, protein_id_mean, max(protein_id_table), lwd = 2, lty=2, col = "#00000066")
  text(protein_id_mean, max(protein_id_table) * 0.9, paste("Mean:",round(protein_id_mean, digits = 3)), pos=4 )
  mtext(letter_list[i], side = 3, line = 1, at = 116, cex = 2)
}
dev.off()


################################################################################

id_clust_data = read.table("~/git/speciation_synteny/summary_data/identity_and_microsynteny_pairwise.tab", header =TRUE, sep = "\t")
id_clust_data = id_clust_data[2:dim(id_clust_data)[1],]
pdf(file = "~/git/speciation_synteny/supplements_for_paper/clustering_and_blast_overview_v2.pdf" , width = 8, height = 10, paper = "a4")
par(mfrow=c(4,2), mar=c(4.5,4.5,4,1))
for ( i in 1:length(pair_genus_names) ){
  bm = matrix( data = c(id_clust_data$q_total[i], id_clust_data$t_total[i], 
                        id_clust_data$q_blast[i], id_clust_data$t_blast[i], 
                        id_clust_data$q_gene[i] , id_clust_data$t_gene[i],
                        id_clust_data$n_ortho[i], id_clust_data$n_ortho[i] ), ncol = 4)
  b = barplot( bm , col = pair_color_list[i], beside = TRUE, main = pair_genus_names[i],
               font.main = 4,
               names.arg = c("Total\ngenes", "Genes w/\nblast hits", "Syntenic\ngenes", "One-to-one\northologs"))
  text(b[3:8], c(id_clust_data$q_blast[i], id_clust_data$t_blast[i], 
                 id_clust_data$q_gene[i] , id_clust_data$t_gene[i], 
                 id_clust_data$n_ortho[i], id_clust_data$n_ortho[i]) - 
         rep(c(id_clust_data$q_total[i], id_clust_data$t_total[i]),3)*0.1,
       round(c(id_clust_data$q_blast[i], id_clust_data$t_blast[i], 
               id_clust_data$q_gene[i] , id_clust_data$t_gene[i], 
               id_clust_data$n_ortho[i], id_clust_data$n_ortho[i]) / 
               rep(c(id_clust_data$q_total[i], id_clust_data$t_total[i]),3),digits = 2)
  )
  mtext(letter_list[i], side = 3, line = 2, at = -0.5, cex = 2)
}
dev.off()


################################################################################

acropora_id_file.g = "~/project/speciation-synteny/datasets/acropora/AhyaV1_vs_Amilv2.homologs_identity.g.tab"
acropora_id_file.m = "~/project/speciation-synteny/datasets/acropora/AhyaV1_vs_Amilv2.homologs_identity.m.tab"

acropora_id.g = read.table( acropora_id_file.g, header = FALSE, sep = "\t", col.names = homolog_id_cols)
acropora_id.m = read.table( acropora_id_file.m, header = FALSE, sep = "\t", col.names = homolog_id_cols)

head(acrorpora_id.g)

acropora_id_pct.g = round( as.numeric(acropora_id.g$identities / (acropora_id.g$identities + acropora_id.g$differences )),digits = 3)
acropora_id_pct.m = round( as.numeric(acropora_id.m$identities / (acropora_id.m$identities + acropora_id.m$differences )),digits = 3)
acropora_id_table.g = table(acropora_id_pct.g)
acropora_id_table.m = table(acropora_id_pct.m)
mean(acropora_id_pct.g)
mean(acropora_id_pct.m)


daphnia_id_file.g = "~/project/speciation-synteny/datasets/daphnia/Dpul_vs_Dmag.homologs_identity.g.tab"
daphnia_id_file.m = "~/project/speciation-synteny/datasets/daphnia/Dpul_vs_Dmag.homologs_identity.m.tab"
daphnia_id.g = read.table( daphnia_id_file.g, header = FALSE, sep = "\t", col.names = homolog_id_cols)
daphnia_id.m = read.table( daphnia_id_file.m, header = FALSE, sep = "\t", col.names = homolog_id_cols)
daphnia_id_pct.g = round( as.numeric(daphnia_id.g$identities / (daphnia_id.g$identities + daphnia_id.g$differences )),digits = 3)
daphnia_id_pct.m = round( as.numeric(daphnia_id.m$identities / (daphnia_id.m$identities + daphnia_id.m$differences )),digits = 3)
daphnia_id_table.g = table(daphnia_id_pct.g)
daphnia_id_table.m = table(daphnia_id_pct.m)
mean(daphnia_id_pct.g)
mean(daphnia_id_pct.m)


plot( 0,0, type = 'n', cex.lab=1.4,
      xlab = "Percent protein identity", ylab = "N protein pairs",
      xlim = c(1,0.3), ylim = c(0,100), axes = FALSE )
axis(1, cex.lab=1.3, cex.axis=1.3)
axis(2, cex.lab=1.3, cex.axis=1.3)
lines(names(acropora_id_table.g), acropora_id_table.g, lwd=2, col="#ac0a18aa" )
lines(names(acropora_id_table.m), acropora_id_table.m, lwd=2, col="#560a08aa" )
lines(names(daphnia_id_table.g), daphnia_id_table.g, lwd=2, col="#220a7eaa" )
lines(names(daphnia_id_table.m), daphnia_id_table.m, lwd=2, col="#426aceaa" )




#