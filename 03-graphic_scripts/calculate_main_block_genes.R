# calculate main block genes
# created by WRF 2023-07-28
# determine number of genes that are on the dominant matching chromosome
# print totals
# and then print a table of those genes that are NOT on the main block

args = commandArgs(trailingOnly=TRUE)

all2Dfile = args[1]
outputfile = gsub("([\\w/]+)\\....$","\\1.off_main.tab",gsub(".gz$","",all2Dfile,perl=TRUE),perl=TRUE)
#all2Dfile = "macrosynteny_plots/Cqui1_vs_CPPV2.scaffold_synteny.tab.gz"
all2Ddata = read.table(all2Dfile, sep="\t", stringsAsFactors=FALSE)
categories = all2Ddata[,1]
is_scaf1 = which(categories=="s1")
scafdata1 = all2Ddata[is_scaf1,]
longestscaf1 = max(scafdata1[,6])
is_longscafs1 = which(as.numeric(scafdata1[,5]) > 0.01)
longscafs1 = c(0, scafdata1[,6][is_longscafs1] )
longscafs1_names = scafdata1[is_longscafs1,2]
is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
longestscaf2 = max(scafdata2[,6])
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.01)
scafdata2_long = scafdata2[is_longscafs2,]
longscafs2 = c(0, scafdata2[,6][is_longscafs2] )
length(longscafs2)
longscafs2_names = scafdata2_long[,2]
is_points = which(categories=="g")
pointsdata = all2Ddata[is_points,]
is_both_longscaf = !is.na(match(pointsdata[,3], longscafs1_names)) & !is.na(match(pointsdata[,5], longscafs2_names))
pointsdata_long = pointsdata[is_both_longscaf,]

scaffold_match_data = data.frame( sc1 = pointsdata_long[,3], sc2 = pointsdata_long[,5], stringsAsFactors = FALSE )

match_freq_table = table(scaffold_match_data)
#match_freq_table
scaffold_best_match_xi = apply(match_freq_table, 1, which.max)
#scaffold_best_match_xi
scaffold_best_match_yi = apply(match_freq_table, 2, which.max)
#scaffold_best_match_yi
x_best_of_y_name = names(scaffold_best_match_yi)[scaffold_best_match_xi]
#x_best_of_y_name
point_best_y = x_best_of_y_name[match(pointsdata_long$V3,names(scaffold_best_match_xi))]
is_on_best_match = pointsdata_long$V5==point_best_y
counts_is_best = table(is_on_best_match)
print( paste( basename(all2Dfile), counts_is_best[["TRUE"]], counts_is_best[["FALSE"]], sep = "   ") )

off_main_genes = pointsdata_long[!is_on_best_match,]

write.table(off_main_genes, file = outputfile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


