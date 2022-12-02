#Import of data

uniprot = read.table('./uniprot-reviewed_yes_headers.txt', header =T, sep = "\t")
blast_chromosomes = read.table('genemapping_usearch/chromosomes/blastout_function_selected.txt', header =F, sep = "\t")
blast_plasmids = read.table('genemapping_usearch/plasmids/blastout_function_selected.txt', header =F, sep = "\t")
blast_plasmids_d = read.table('genemapping_diamond_blast/plasmids/blastout_function_selected.txt', header =F, sep = "\t")
Padloc = read.table("./padloc_ds1.txt", sep = "\t")
proteins_chromosomes = read.table("./chromosomes_protein.txt", sep = "\t")
proteins_plasmids = read.table("./plasmid_protein.txt", sep = "\t")
proteins_duplicates = read.table("./duplicates_proteins.txt", sep = "\t")
protein_class = read.table("./protein_class.txt", sep = "\t")
taxa = read.table('./taxonomy.txt', header =T, sep = "\t") #taxa needs to be ajusted to the amount of samples in the data set

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
names(protein_class)[1] <- 'protein'
names(protein_class)[2] <- 'Protein_Location'

proteins_all = rbind(proteins_chromosomes, proteins_plasmids)

pathways_padloc = as.data.frame.matrix(merge(Padloc,uniprot))
padloc_table = as.data.frame.matrix(table(pathways_padloc$V1,pathways_padloc$enzyme))

padloc_table[padloc_table>0] <-1

blast = rbind(blast_chromosomes, blast_plasmids)

pathways = as.data.frame.matrix(merge(blast,uniprot))

enzyme.table = as.data.frame.matrix(table(pathways$V1,pathways$enzyme))

enzyme.table = cbind(padloc_table, enzyme.table)

enzyme.table$`PI-type proteinase` <- NULL

enzyme.table = cbind(enzyme.table, enzyme.table2)

selected.enzyme = enzyme.table[ ,colnames(enzyme.table) %in% proteins_all$V1]

selected.enzyme = selected.enzyme[-49,]

rownames(selected.enzyme) <- taxa$isolate

trans.selected.enzyme = t(selected.enzyme)


################################################
#PCoA plot done with selected.enszyme

mod3 = capscale(selected.enzyme ~ 1, dist = "canberra")
pcs = (as.data.frame(mod3$CA$u))[,(1:6)]
eig.vals = mod3$CA$eig
ggplot = ggplot(pcs, aes(x = MDS1, y = MDS2, shape = taxa$Species, color = taxa$Species)) + 
  geom_point(size = 3.5) + 
  labs(color = "Species") +
  labs(shape = "Species") +
  scale_color_manual(values=c("green", "azure4", "cadetblue1")) + 
  theme_bw() + geom_vline(xintercept = 0, lty = 2, color = alpha(colour = "orange", alpha = 0.5)) +
  geom_hline(yintercept = 0, lty = 2, color = alpha(colour = "orange", alpha = 0.5))  +
  xlab(glue("MDS1, {round(eig.vals[[1]]/sum(eig.vals) * 100, 1)}%")) + 
  ylab(glue("MDS2, {round(eig.vals[[2]]/sum(eig.vals) * 100, 1)}%")) + ggtitle("PCoA") + coord_fixed()
ggplot

ggsave("./PCoA.jpeg", device = "jpeg", dpi = 400, width = 10, height = 10)
dev.off()


###########################################
#Heatmap done with trans.selected.enzyme

trans.selected.enzyme = apply(trans.selected.enzyme, MARGIN = 2, FUN = function(x) {log(x + 20/10)})

n = 3 
set.seed(1000)
palette = distinctColorPalette(n)
Protein_Location = sort(unique(protein_class$Protein_Location))
color.key.taxonomy = data.frame(Protein_Location,palette)

rownames = (protein_class$protein)
rownames(protein_class)<-rownames

df <- data.frame(protein = rownames(protein_class), 
                 Protein_Location = ifelse(protein_class$Protein_Location == "Chromosome", "Chromosome", 
                                           ifelse(protein_class$Protein_Location == "Plasmid", "Plasmid", "Chromosome&Plasmid")))
df = data.frame(df, palette = ifelse(df$Protein_Location == "Chromosome", "deeppink4", 
                                     ifelse(df$Protein_Location== "Plasmid", "darkgreen", "blue4")))
#df = as.matrix(df) 
rownames(df) <- df$protein
protein_lgd = Legend(df$palette, labels = unique(df$Protein_Location), 
                     legend_gp = gpar(fill = unique(df$palette)), title = "Protein Gene Location", 
                     grid_width = unit(10, "mm"), grid_height = unit(10, "mm"),
                     title_gp = gpar(fontsize = 20, fontface = "bold"),
                     labels_gp = gpar(fontsize = 15))


#creating heatmap

ht1 = Heatmap(trans.selected.enzyme, name = "Abundance (log)", col = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
              heatmap_legend_param = list(grid_height = unit(10, "mm"), grid_width = unit(10,"mm"), 
                                          labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 20, fontface = "bold")), 
              rect_gp = gpar(col = "black", lwd = 1.2),
              clustering_distance_rows = "canberra", 
              clustering_method_rows = "ward.D2",
              show_row_dend = F,
              column_km = 3,
              clustering_distance_columns = "canberra", 
              clustering_method_columns = "ward.D2",
              row_names_gp = gpar(col = df$palette, fontsize = 16, fontface = "bold"), row_names_side = "left",
              row_dend_side = "right",
              row_names_max_width = max_text_width(rownames(trans.selected.enzyme)),
              column_names_gp = gpar(col = "black", fontsize = 14, fontface = "bold"), column_names_rot = 45,
              column_title = "Isolate Name", column_title_gp = gpar(fontsize = 20, fontface = "bold"), column_title_side = "bottom",
              column_dend_height = unit(5, "cm"), row_dend_width = unit(4, "cm"),
              heatmap_width = unit(54, "cm"),
              heatmap_height = unit(77, "cm"),
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("green", "azure4", "cadetblue1")), 
                                                                  labels = c("azure4" = "L. cremoris", 
                                                                             "cadetblue1" = "L. laudensis", 
                                                                             "green" = "Le. mesenteroides"), 
                                                                  labels_gp = gpar(col = "black", fontsize = 13, fontface = "bold")))) 




jpeg(filename = "./heatmap.jpeg", width = 72, height = 83, units = "cm", res = 300)

draw(ht1, heatmap_legend_list = list(protein_lgd))

dev.off()





