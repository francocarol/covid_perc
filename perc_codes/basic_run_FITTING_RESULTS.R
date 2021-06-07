################################################
################################################
################################################

# Basic Run to make plots and get AIC 

################################################
################################################
################################################

# All necessary source codes are in this folder:
# REMEMBER TO ERASE IT WHEN ZIPPING IT AND SHARING IT 
setwd("/home/caroline/Documents/paper_covid/codes/")
#IN CLUSTER
#setwd("~/extraspace/Carol/paper_covid/")
### Definir cidade:
estado = "SP"
cidade = "Sao_Paulo"

# Loading data files corresponding to the fitting results for the 3 model versions:
#old fitting results:
file1 <- paste0("fits/new_model1_srag.1.CM.log0-2021_05_17_SP.csv")
file2 <- paste0("fits/model2_scn50_srag_2021_05_24_SP.csv")
file3 <- paste0("fits/model3_scn54_30npi_srag_2021_05_24_SP.csv")
Scn1 <- 50
Scn2 <- 50
Scn3 <- 54

# RUN CALC_AIC FUNCTION:
source("AIC_calc.R")
AICtable <- calc_AIC(file1, file2, file3)
AICtable
write.csv(AICtable, "output/AIC_table_3models.csv", row.names=FALSE)

# RUN PLOTTING FUNCTION
source("plots_fitting_results.R")
combined_plots <- plot_together(file1, file2, file3, Scn1, Scn2, Scn3,
                                y_lim_cases=6000, y_lim_mort=1750)
print(combined_plots)
ggsave(combined_plots, file = "output/fitting_all_versions.png", device = "png", dpi = 100, width = 30, height = 10, units = "cm")

