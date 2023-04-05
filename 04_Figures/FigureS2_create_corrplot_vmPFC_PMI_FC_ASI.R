# Load necessary libraries, install using 'install.packages("LIBRARY")' if needed
library(tidyverse)
library(ggplot2)
library(DescTools)
library(ggpubr)

# Load data. First command will prompt a window to locate the file with subjects demographics "subject_correlations.csv" on your machine
# For the second prompt please locate the file with subjects information "subject_info.csv" on your machine
data.path.corr <- file.choose()
df.correlations <- read_csv(data.path.corr)
data.path.info <- file.choose()
df.info <- read_csv(data.path.info)

# normalize Pearson's r using Fisher r-to-z transformation
df.correlations$cor_z <- FisherZ(df.correlations$cor_r)

# combine both dataframes
df.full <- inner_join(df.correlations, df.info, by = c("ID", "Group"))


# make plot
df.full %>% 
  filter(Group == "MA") %>%
  filter(roi_1 == "vmPFC" & roi_2 == "pm_ins") %>% 
  ggscatter(x = "ASI.Total", y = "cor_z",
            add = "reg.line",
            conf.int = TRUE, 
            add.params = list(color = "black",
                              fill = "lightgray",
                              cor.coef.size = 20),
            xlab = "ASI Total score",
            ylab = "vmPFC-PMI correlation coefficient (z-score)") +
  theme(axis.line = element_line(colour = 'black', size = 1.1),
        text = element_text(size = 14),
        axis.text = element_text(face = "bold", color = "black"),
        axis.ticks.x = element_line(size = 1.3, color = "black"),
        axis.ticks.y = element_line(size = 1.3, color = "black"),
        legend.key=element_blank()) +
  stat_cor(method = "pearson", label.x = 37, label.y = 0.7, cor.coef.name = "r", size = 5) + 
  annotate("text", size = 5, x = 45, y = 0.63, label = "95% CI [-0.69, -0.05]") + 
  annotate("text", size = 5, x = 40.7, y = 0.56, label = "BF = 3.03")

ggsave("C:\\Users\\Jonas\\Documents\\GitHub\\GAD_MBA_FC\\04_Figures\\FigureS2_vmPFC-PMI-FC_ASI_corr.pdf", device = "pdf", dpi = "print")