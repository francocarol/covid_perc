#############
### Load packages
#############

require(tidyverse)
if(!require(ggpubr)){install.packages("ggpubr"); library(ggpubr)}

#############
### Load data
#############

c_home <- data.frame(read.csv("DATA/5. matrix_cont_home_PONTO.csv", stringsAsFactors = FALSE))
c_school <- data.frame(read.csv("DATA/6. matrix_cont_school_PONTO.csv", stringsAsFactors = FALSE))
c_work <- data.frame(read.csv("DATA/7. matrix_cont_work_PONTO.csv", stringsAsFactors = FALSE))
c_other <- data.frame(read.csv("DATA/8. matrix_cont_other_loc_PONTO.csv", stringsAsFactors = FALSE))

load("PREM/prem2020_bra.Rdata")

#############
### Plot function script
#############

plot_cmatrix <- function(dat, grid = FALSE, normalise = FALSE, labels = FALSE, grid_color = "grey", 
                         text_size = 2.5, round_number = 2, legend.pos="right") {

if(normalise) {
  e <- as.numeric(eigen(dat)$values[1])
  dat <- dat/e
}

colnames(dat) <- paste0("S",seq(1,nrow(dat)))
rownames(dat) <- paste0("f",seq(1,nrow(dat)))

R <- dat %>% 
  as.data.frame() %>%
  rownames_to_column("f_id") %>%
  pivot_longer(-c(f_id), names_to = "samples", values_to = "counts") %>%
  mutate(samples= fct_relevel(samples,colnames(dat))) %>%
  mutate(f_id= fct_relevel(f_id,rownames(dat)))

if(grid == FALSE) { 
gplot <- R %>% ggplot(aes(x=samples, y=f_id, fill=counts)) + 
  geom_raster() + scale_fill_gradient2("", low = "white", high = "royalblue4") +
  labs(x = "Age of individual", y = "Age of contact") +
  scale_x_discrete(labels = (0:ncol(dat))*5) +
  scale_y_discrete(labels = (0:ncol(dat))*5)
} else {
gplot <- ggplot(R, aes(samples, f_id)) + geom_tile(aes(fill = counts),
                colour = grid_color) + scale_fill_gradient2("", low = "white", high = "royalblue4") +
  labs(x = "", y = "") + 
  scale_x_discrete(labels = (0:ncol(dat))*5) +
  scale_y_discrete(labels = (0:ncol(dat))*5)
}

if (labels) gplot <- gplot + geom_text(aes(label = round(counts,round_number)), size = text_size)

if (legend.pos=="right"){
  gplot<- gplot +
  theme(legend.position = "right")}
if (legend.pos=="bottom") {
  gplot <- gplot +
  theme(legend.position = "bottom", legend.direction = "vertical")
  }

return(gplot)
}

#############
### Plot Contact Matrices
#############

normalise = FALSE
if (normalise) {file.ending = "_norm.png"} else {file.ending = ".png"}

#######
# All
d = 0.2
nh <- plot_cmatrix(new_home, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
ns <- plot_cmatrix(new_school, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
nw <- plot_cmatrix(new_work, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
no <- plot_cmatrix(new_others, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))

prem2020_plots <- ggarrange(nh, ns, nw, no,
                    labels = c("Home", "School", "Work", "Other"),
                    ncol = 2, nrow = 2)

filename = paste0("PREM_2020",file.ending)
ggsave(prem2020_plots, file = filename, dpi = 100, width = 40, height = 30, units = "cm")

# Urban
uh <- plot_cmatrix(u_home, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
us <- plot_cmatrix(u_school, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
uw <- plot_cmatrix(u_work, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
uo <- plot_cmatrix(u_others, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))

prem2020_urban_plots <- ggarrange(uh, us, uw, uo,
                            labels = c("Home", "School", "Work", "Other"),
                            ncol = 2, nrow = 2)

filename = paste0("PREM_2020_urban",file.ending)
ggsave(prem2020_urban_plots, file = filename, dpi = 100, width = 40, height = 30, units = "cm")

# Rural
rh <- plot_cmatrix(r_home, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
rs <- plot_cmatrix(r_school, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
rw <- plot_cmatrix(r_work, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
ro <- plot_cmatrix(r_others, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))

prem2020_rural_plots <- ggarrange(rh, rs, rw, ro,
                                  labels = c("Home", "School", "Work", "Other"),
                                  ncol = 2, nrow = 2)

filename = paste0("PREM_2020_rural",file.ending)
ggsave(prem2020_rural_plots, file = filename, dpi = 100, width = 40, height = 30, units = "cm")

# Prem 2017
ch <- plot_cmatrix(c_home, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
cs <- plot_cmatrix(c_school, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
cw <- plot_cmatrix(c_work, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))
co <- plot_cmatrix(c_other, normalise = normalise) + theme(plot.margin = margin(t = d, r = d, b = d, l = 2.2, unit = "cm"))

prem2017_plots <- ggarrange(ch, cs, cw, co,
                                  labels = c("Home", "School", "Work", "Other"),
                                  ncol = 2, nrow = 2)

filename = paste0("PREM_2017",file.ending)
ggsave(prem2017_plots, file = filename, dpi = 100, width = 40, height = 30, units = "cm")

#############
### Correlation matrix
#############

c1 = cor(c(as.matrix(c_home)), c(new_home))
c2 = cor(c(as.matrix(c_school)), c(new_school))
c3 = cor(c(as.matrix(c_work)), c(new_work))
c4 = cor(c(as.matrix(c_other)), c(new_others))

c5 = cor(c(as.matrix(c_home)), c(u_home))
c6 = cor(c(as.matrix(c_school)), c(u_school))
c7 = cor(c(as.matrix(c_work)), c(u_work))
c8 = cor(c(as.matrix(c_other)), c(u_others))

c9 = cor(c(as.matrix(c_home)), c(r_home))
c10 = cor(c(as.matrix(c_school)), c(r_school))
c11 = cor(c(as.matrix(c_work)), c(r_work))
c12 = cor(c(as.matrix(c_other)), c(r_others))

cormatrix <- matrix(0, 2, 6)
cormatrix[1,] = c(c1,c2,c5,c6,c9,c10)
cormatrix[2,] = c(c3,c4,c7,c8,c11,c12)

colnames(cormatrix) = c("All_2020.1", "All_2020.2", "Urban_2020.1", "Urban_2020.2", "Rural_2020.1", "Rural_2020.2")
rownames(cormatrix) = c("All_2017.1", "All_2017.2")

write.table(round(cormatrix, 2), file = "cor_matrix.txt", quote = FALSE, sep = "\t")

