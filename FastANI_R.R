ANIchromchrom <- readxl::read_xlsx("./FastANI.xlsx", sheet = 1 )
ANIchromref <- readxl::read_xlsx("./FastANI.xlsx", sheet = 2 )
taxa = read.table('./taxonomy.txt', header =T, sep = "\t")


###########################################
n = 3 
set.seed(1000)
palette = distinctColorPalette(n)
Species = sort(unique(taxa$Species))
color.key.taxonomy = data.frame(Species,palette)
taxonomy.table = join(taxa,color.key.taxonomy)
taxonomy.table <- column_to_rownames(taxonomy.table, "isolate")


###########################################
ANIchromchrom <- ANIchromchrom %>% column_to_rownames("...1")

ANIchromchrom <- as.matrix(ANIchromchrom)

ANIchromchrom = t(ANIchromchrom)


###########################################
r_an1 = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c( "green", "azure4", "cadetblue1")), 
                                           labels = c("green" = "L. cremoris",
                                                      "cadetblue1" = "L. laudensis",
                                                      "azure4" = "Le. mesenteroides"), 
                                           labels_gp = gpar(col = "black", fontsize = 11, fontface = "bold")))


ht1 = Heatmap(ANIchromchrom, name = "Percent ID", top_annotation = r_an1, col = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
              heatmap_legend_param = list(grid_height = unit(10, "mm"), grid_width = unit(10,"mm"), 
                                          labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 20, fontface = "bold")), 
              rect_gp = gpar(col = "black", lwd = 1.2),
              column_km = 3,
              clustering_distance_rows = "canberra", 
              clustering_method_rows = "ward.D2",
              clustering_distance_columns = "canberra", 
              clustering_method_columns = "ward.D2",
              show_row_dend = F,
              show_column_dend = F,
              row_names_gp = gpar(fontsize = 20, fontface = "bold"), row_names_side = "left",
              row_dend_side = "right",
              max_text_width(rownames(df)),
              column_names_gp = gpar(col = "Black", fontsize = 14, fontface = "bold"), column_names_rot = 45,
              column_title = "Isolate Name", column_title_gp = gpar(fontsize = 20, fontface = "bold"), column_title_side = "bottom",
              column_dend_height = unit(10, "cm"), row_dend_width = unit(5, "cm"),
              heatmap_width = unit(54, "cm"),
              heatmap_height = unit(8, "cm"))

jpeg(filename = "./ANI.jpeg", width = 72, height = 10, units = "cm", res = 400)

draw(ht1)

dev.off()


###########################################
ANIchromref <- column_to_rownames(ANIchromref, "...1")

ANIchromref <- as.matrix(ANIchromref)

ANIchromref[upper.tri(ANIchromref, diag = T)] <- -1

ANIchromref[upper.tri(ANIchromref, diag = F)] <- ANIchromref[lower.tri(ANIchromref, diag = F)]

ANIchromref[ANIchromref== -1] <- 100

ANIchromref[upper.tri(ANIchromref, diag = F)]<-NA


###########################################
Species = taxonomy.table %>% dplyr::select(Species)
species = rownames_to_column(Species, "db")
ANIchromref1 <- ANIchromref %>% data.frame
ANIchromref1 <- rownames_to_column(ANIchromref1, "db")

ANIchromref1 <- left_join(ANIchromref1, species, keep = T,by = "db")
ANIchromref2 = ANIchromref1 %>% dplyr::select(1:59)
spec.col = ANIchromref1 %>% dplyr::select(60:dim(ANIchromref1)[2]) 

ANIchromref2 <- column_to_rownames(ANIchromref2, "db.x")
ANIchromref2 <- as.matrix(ANIchromref2)

spec.col <- column_to_rownames(spec.col, "db.y")


###########################################
r_an2 = HeatmapAnnotation(Species = spec.col$Species,
                          which = "row", 
                          annotation_name_rot = 45, 
                          col = list(Species = c("L. cremoris" = "green", 
                                                 "L. laudensis"  = "azure4",
                                                 "Le. mesenteroides" = "cadetblue1")))


jpeg(filename = "./ANI2.jpeg", width = 68, height = 68, units = "cm", res = 300) 

Heatmap(ANIchromref2, border = F, split = spec.col, left_annotation = r_an2 , name = "Percent ID", color_space = "RGB", col = RColorBrewer::brewer.pal(n = 9, name = "Blues"),
        rect_gp = gpar(col = "white", lwd = 1.2),
        row_names_gp = gpar(fontsize = 20, fontface = "bold"), row_names_side = "left",
        cluster_columns = FALSE,
        cluster_rows = FALSE, na_col = "White",
        max_text_width(rownames(ANIchromref2)),
        column_names_gp = gpar(col = "black", fontsize = 20, fontface = "bold"), column_names_rot = 45,
        column_title = "Isolate Name", column_title_gp = gpar(fontsize = 20, fontface = "bold"), column_title_side = "top",
        heatmap_width = unit(54, "cm"),
        heatmap_height = unit(54, "cm"))

dev.off()

