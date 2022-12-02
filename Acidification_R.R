#Accidification curve
data <- readxl::read_xlsx("./Samlet.xlsx", sheet = 4, col_names = T )

ggp <- ggplot(data, aes(time, ph, color = Blend)) +
  geom_line() +
  labs(x = "Time (h)", y = "pH") + 
  theme_bw()

ggpcont <- ggp + scale_y_continuous(expand = c(0,0), breaks = c(4,4.5,5,5.5,6,6.5,7), limits = c(4.1,6.8)) + 
  scale_x_continuous(expand = c(0,0), breaks = c(0,2,4,6,8,10,12,14,16,18,20), limits = c(0,21))

ggpcont

ggsave("./acidi_Curve.jpeg", device = "jpeg", dpi = 400, width = 10, height = 5)
dev.off()


#######################################
#Accidification rate
f1 <- readxl::read_xlsx("./Samlet.xlsx", sheet = 8, col_names = T )
f2 <- readxl::read_xlsx("./Samlet.xlsx", sheet = 9, col_names = T )
f3 <- readxl::read_xlsx("./Samlet.xlsx", sheet = 10, col_names = T )
f4 <- readxl::read_xlsx("./Samlet.xlsx", sheet = 11, col_names = T )
f5 <- readxl::read_xlsx("./Samlet.xlsx", sheet = 12, col_names = T )
f6 <- readxl::read_xlsx("./Samlet.xlsx", sheet = 13, col_names = T )
f9 <- readxl::read_xlsx("./Samlet.xlsx", sheet = 14, col_names = T )
f10 <- readxl::read_xlsx("./Samlet.xlsx", sheet = 15, col_names = T )

sg1 <- as.data.frame(savitzkyGolay(f1$deri,0,1,13))
sg2 <- as.data.frame(savitzkyGolay(f2$deri,0,1,13))
sg3 <- as.data.frame(savitzkyGolay(f3$deri,0,1,13))
sg4 <- as.data.frame(savitzkyGolay(f4$deri,0,1,13))
sg5 <- as.data.frame(savitzkyGolay(f5$deri,0,1,13))
sg6 <- as.data.frame(savitzkyGolay(f6$deri,0,1,13))
sg9 <- as.data.frame(savitzkyGolay(f9$deri,0,1,13))
sg10 <- as.data.frame(savitzkyGolay(f10$deri,0,1,13))

write_xlsx(sg1,"./df1.xlsx")
write_xlsx(sg2,"./df2.xlsx")
write_xlsx(sg3,"./df3.xlsx")
write_xlsx(sg4,"./df4.xlsx")
write_xlsx(sg5,"./df5.xlsx")
write_xlsx(sg6,"./df6.xlsx")
write_xlsx(sg9,"./df9.xlsx")
write_xlsx(sg10,"./df10.xlsx")

data2 <- readxl::read_xlsx("./Samlet.xlsx", sheet = 17, col_names = T )


ggp <- ggplot(data2, aes(time, deri, color = Blend)) +
  geom_line() +
  labs(x = "Time (h)", y = "dpH/dt (pH/min)") + 
  theme_bw()
ggpcont <- ggp + scale_y_continuous(expand = c(0,0), breaks = c(-0.007, -0.006, -0.005, -0.004, -0.003, -0.002, -0.001, 0.00), limits = c(-0.007,0.0005)) + 
  scale_x_continuous(expand = c(0,0), breaks = c(0,2,4,6,8,10,12,14,16,18,20), limits = c(0,20.5))

ggpcont

ggsave("./acidi_rate.jpeg", device = "jpeg", dpi = 400, width = 10, height = 5)
dev.off()
