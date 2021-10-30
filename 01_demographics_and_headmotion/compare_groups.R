# Load necessary libraries, install using 'install.packages("LIBRARY")' if needed
library(tidyverse)


# set up function to always show mean and standard deviation in demographics table - analogous to the table in the manuscript 
mean_sd_builder <- function(variable){
  paste(round(mean(variable, na.rm = TRUE), 1), round(sd(variable, na.rm = TRUE), 1), sep = " \u00b1 ")
}


# Load data. First command will prompt a window to locate the file with subjects demographics "subject_info.csv" on your machine
data.path <- file.choose()
df.demographics <- read_csv(data.path)

# this builds the descriptive demographic table as in the manuscript: mean, standard deviation, number of participants
df.demographics %>% 
  # this function performs all future calculations on the HC and the GAD (MA) group separately 
  group_by(Group) %>% 
  summarise(mean_sd_age = mean_sd_builder(Age),
            n_age = sum(!is.na(Age)),
            mean_sd_bmi = mean_sd_builder(BMI),
            n_bmi = sum(!is.na(BMI)),
            mean_sd_PHQ = mean_sd_builder(PHQ.9),
            n_PHQ = sum(!is.na(PHQ.9)),
            mean_sd_OASIS = mean_sd_builder(OASIS),
            n_OASIS = sum(!is.na(OASIS)),
            mean_sd_GAD7 = mean_sd_builder(GAD7),
            n_GAD7 = sum(!is.na(GAD7)),
            mean_sd_STAI_S = mean_sd_builder(STAI.State),
            n_STAI_S = sum(!is.na(STAI.State)),
            mean_sd_STAI_T = mean_sd_builder(STAI.Trait),
            n_STAI_T = sum(!is.na(STAI.Trait)),
            mean_sd_ASI_T = mean_sd_builder(ASI.Total),
            n_ASI_T = sum(!is.na(ASI.Total)),
            mean_sd_ASI_P = mean_sd_builder(ASI.Physical),
            n_ASI_P = sum(!is.na(ASI.Physical)),
            mean_sd_ASI_C = mean_sd_builder(ASI.Cognitive),
            n_ASI_C = sum(!is.na(ASI.Cognitive)),
            mean_sd_ASI_S = mean_sd_builder(ASI.Social),
            n_ASI_S = sum(!is.na(ASI.Social))) %>% 
  # convert everything to character, otherwise pivoting will fail
  mutate(across(everything(), as.character)) %>% 
  # these two functions  transform the dataframe from wide format to long and back to wide for easier viewing
  pivot_longer(-Group, names_to = "key", values_to = "val", values_ptypes = list(val = "character")) %>% 
  pivot_wider(names_from=Group, values_from=val) %>%
  select(key, MA, HC) %>% 
  print(n = 22)

# test group differences for demographic variables in the order they appear in the demographics table
t.test(Age ~ Group, data = df.demographics)
t.test(BMI ~ Group, data = df.demographics)
t.test(PHQ.9 ~ Group, data = df.demographics)
t.test(OASIS ~ Group, data = df.demographics)
t.test(GAD7 ~ Group, data = df.demographics)
t.test(STAI.State ~ Group, data = df.demographics)
t.test(STAI.Trait ~ Group, data = df.demographics)
t.test(ASI.Total ~ Group, data = df.demographics)
t.test(ASI.Physical ~ Group, data = df.demographics)
t.test(ASI.Cognitive ~ Group, data = df.demographics)
t.test(ASI.Social ~ Group, data = df.demographics)
  
  
# compare groups based on average head motion during resting state scan
df.demographics %>% 
  group_by(Group) %>% 
  summarize(mean = mean(avg_motion),
            standard_error = sd(avg_motion)/sqrt(length(Group)))

# test if groups differ significantly on average head motion
wilcox.test(avg_motion ~ Group, data = df.demographics)