# plot synteny block lengths across species
# using microsynteny.py script at:
# https://github.com/wrf/genomeGTFtools
# 2023-07-25 by WRF


ms_file_list_dir = dir( "~/git/speciation_synteny/05-microsynteny_plots/", "*.microsynteny.*tab.gz", recursive = TRUE)
ms_file_list = ms_file_list_dir[c(33,1,31,6,8,5,17,2,9,12,14,13,
             27,4,28,29,3,7,32,18,19,30,20)]

pair_color_list = c("#a0e499aa", "#ac0a18aa", "#590aacaa", "#ec60bdaa", "#ec60bdaa", "#ec60bdaa", "#220a7eaa",
                    "#2ecd14aa", "#44cd14aa", "#cd9714aa", "#cd9714aa", "#cd9714aa", 
                    "#ac0a18aa", "#cd1f14aa",
                    "#0da730aa", "#a6ab09aa", "#ab099daa", "#9fcd1aaf",
                    "#098eabaa", "#e38601aa", "#e38601aa", "#e33601aa", "#0a4facaa",
                    "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa", "#666666aa")
pair_color_list.no_alpha = substr(pair_color_list,1,7) # trim alpha

block_headers = c("q_scaffold", "s_scaffold", "block_ID", 
                  "q_gene_ID", "q_gene_pos_start", "q_gene_pos_end", "q_gene_strand", 
                  "s_gene_ID", "s_gene_pos_start", "s_gene_pos_end", "s_gene_strand", 
                  "bitscore" )
pair_genus_names = c("Tethya wilhelma-minuta", "Acropora hyacinthus-millepora", "Octopus bimaculoides-sinensis", 
                     "Crassostrea angulata-virginica", "Crassostrea gigas-virginica", "Crassostrea angulata-gigas",
                     "Daphnia pulex-magna", "Anastrepha obliqua-ludens", "Culex quinquefasciatus-pipiens", 
                     "Drosophila melanogaster-erecta", "Drosophila melanogaster-pseudoobscura", "Drosophila melanogaster-grimshawi", 
                     "Vespa crabro-velutina", "Bombus pyrosoma-terrestris", "Lytechinus variegatus-picta", 
                     "Mauremys reevesii-mutica", 
                     "Bufo gargarizans-bufo", "Cervus hanglu-elaphus",
                     "Perca flavescens-fluviatilis", "Epinephelus fuscoguttatus-lanceolatum", "Epinephelus fuscoguttatus-moara", 
                     "Oreochromis aureus-niloticus", "Thunnus maccoyii-albacares" )

# make PDF
pdf(file = "~/git/speciation_synteny/supplements_for_paper/microsynteny_blocks_overview_v1.pdf" , width = 8, height = 10, paper = "a4")
par(mfrow=c(4,2), mar=c(4.5,4.5,3,1))
# loop through each species, and make 1/8-page-size plot
for (i in 1:(length(ms_file_list))){
    blockdata = read.table(paste0("~/git/speciation_synteny/05-microsynteny_plots/", ms_file_list[i]), 
                           col.names = block_headers, sep="\t", fill = TRUE)
    counttable = table(blockdata[,3])
    totalblk = length(counttable)
    totalgenes = length(unique(blockdata[,4]))
    counttable_hist = table(counttable)
    longest_block = max(as.integer(names(counttable_hist)))
    mostcommon = max(counttable_hist)
    block_len_pos = c(2,5,10,20,50,100,200,500,1000)
    block_len_labs = block_len_pos[c(TRUE,block_len_pos<longest_block)]
    
    plot( as.integer(names(counttable_hist)) , counttable_hist , log="xy",
         xlim = range(block_len_labs),
         main=pair_genus_names[i], font.main = 4, 
         axes = FALSE, cex.lab=1.4,
         xlab="Number of genes in block", ylab="Number of blocks",
         pch=16, col=pair_color_list[i], cex=2)
    axis(1, at = block_len_labs,  cex.axis=1.3 )
    axis(2, cex.axis=1.3 )
    totallab = paste(totalblk, "blocks", sep=" ")
    genelab = paste(totalgenes, "genes in blocks", sep=" ")
    text(longest_block, mostcommon*0.50, totallab, cex=1.3, pos=2)
    text(longest_block, mostcommon*0.85, genelab, cex=1.3, pos=2)
    # calculate regression as log-log, plot as gray line
    bcor = cor( log10(as.integer(names(counttable_hist))) , log10(counttable_hist) )
    bcor2 = bcor^2
    bc_lm = lm( log10(counttable_hist) ~ log10(as.integer(names(counttable_hist))) )
    text(2,2,paste("R2 =",round(bcor2,digits = 3)), pos=4, cex=1.2)
    abline(bc_lm, lty=2, lwd=4, col="#00000066", untf=FALSE)
}
dev.off()



#