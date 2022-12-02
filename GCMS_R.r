gcms = read.table('./gcms.txt', header =T, sep = "\t", row.names = 1 )

mod = prcomp(gcms, scale = T, center = T)

fviz_eig(mod, title = "GCMS Scree Plot")

ggsave("./GCMS_scree_plot.jpeg", device = "jpeg", dpi = 400, width = 7.5, height = 3.75, bg="white")

dev.off()

fviz_pca_var(mod, repel = TRUE, title = "GCMS Loadings", col.var = "darkgreen", col.circle = "black")


ggsave("./GCMS_rate_loadings.jpeg", device = "jpeg", dpi = 400, width = 7.5, height = 7.5, bg="white")

dev.off()

fviz_pca_ind(mod, repel = TRUE, title = "GCMS Scores", col.ind = "darkgreen")

ggsave("./GCMS_rate_scores.jpeg", device = "jpeg", dpi = 400, width = 7.5, height = 7.5, bg="white")

dev.off()


