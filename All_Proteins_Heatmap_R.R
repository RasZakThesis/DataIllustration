#Import of data

uniprot = read.table('./uniprot-reviewed_yes_headers.txt', header =T, sep = "\t")
blast_chromosomes = read.table('genemapping_usearch/chromosomes/blastout_function_selected.txt', header =F, sep = "\t")
blast_plasmids = read.table('genemapping_usearch/plasmids/blastout_function_selected.txt', header =F, sep = "\t")
blast_plasmids_d = read.table('genemapping_diamond_blast/plasmids/blastout_function_selected.txt', header =F, sep = "\t")
taxa = read.table('./taxonomy2.txt', header =T, sep = "\t") #taxa needs to be ajusted to the amount of samples in the data set

###############################################
#Data treatment

pathways1 = as.data.frame(merge(blast_plasmids_d,uniprot))

enzyme.table1 = as.data.frame.matrix(table(pathways1$V1,pathways1$enzyme))

enzyme.table2 = enzyme.table1[,'PI-type proteinase']

enzyme.table2 = as.data.frame(enzyme.table2)

names(enzyme.table2)[1] <- 'PI-type proteinase'

rownames = rownames(enzyme.table1)

rownames(enzyme.table2)<-rownames


###########################################
blast = rbind(blast_chromosomes, blast_plasmids)

pathways = as.data.frame.matrix(merge(blast,uniprot))

enzyme.table = as.data.frame.matrix(table(pathways$V1,pathways$enzyme))

enzyme.table$`PI-type proteinase` <- NULL

enzyme.table = cbind(enzyme.table, enzyme.table2)

rownames(enzyme.table) <- taxa$isolate

trans.selected.enzyme = t(enzyme.table)


###########################################
#Heatmap done with trans.selected.enzyme

trans.selected.enzyme = apply(trans.selected.enzyme, MARGIN = 2, FUN = function(x) {log(x + 20/10)})


#creating heatmap

ht1 = Heatmap(trans.selected.enzyme, name = "Abundance (log)", col = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
              heatmap_legend_param = list(grid_height = unit(10, "mm"), grid_width = unit(10,"mm"), 
                                          labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 20, fontface = "bold")), 
              rect_gp = gpar(col = "black", lwd = 1.2),
              clustering_distance_rows = "canberra", 
              clustering_method_rows = "ward.D2",
              column_km = 3,
              clustering_distance_columns = "canberra", 
              clustering_method_columns = "ward.D2",
              row_names_gp = gpar(fontsize = 12, fontface = "bold"), row_names_side = "left",
              row_title = "Protein Name", row_title_gp = gpar(fontsize = 20, fontface = "bold"), row_title_side = "left",
              row_dend_side = "right",
              row_names_max_width = max_text_width(rownames(trans.selected.enzyme)),
              column_names_gp = gpar(col = "black", fontsize = 14, fontface = "bold"), column_names_rot = 45,
              column_title = "Isolate Name", column_title_gp = gpar(fontsize = 20, fontface = "bold"), column_title_side = "bottom",
              column_dend_height = unit(5, "cm"), row_dend_width = unit(4, "cm"),
              heatmap_width = unit(65, "cm"),
              heatmap_height = unit(600, "cm"),
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("green", "azure4", "cadetblue1")), 
                                                                  labels = c("azure4" = "L. cremoris", 
                                                                             "cadetblue1" = "L. laudensis", 
                                                                             "green" = "Le. mesenteroides"), 
                                                                  labels_gp = gpar(col = "black", fontsize = 11, fontface = "bold")))) 




jpeg(filename = "./heatmapallproteins.jpeg", width = 80, height = 650, units = "cm", res = 150)

draw(ht1)

dev.off()





