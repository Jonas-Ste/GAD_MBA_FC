# Load necessary libraries, install using 'install.packages("LIBRARY")' if needed
library(tidyverse)

# Load resting-state data
# First command will prompt a window to locate the file with subjects demographics "subject_correlation.csv" on your machine
corrdata.path <- file.choose()
df.correlations <- read_csv(corrdata.path)

# normalize Pearson's r using Fisher r-to-z transformation
df.correlations$cor_z <- FisherZ(df.correlations$cor_r)

# export resting-state data in correct format as input for AFNIs MBA
df.correlations %>%
  rename(Subj = ID,
         ROI1 = roi_1,
         ROI2 = roi_2, 
         group = Group) %>% 
  select(Subj, ROI1, ROI2, cor_z, group) %>% 
  # adjust path to match data structure on your local machine
  write_delim("MBA_input_full.txt", delim = " ")