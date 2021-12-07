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
mean(subset(df.demographics, Group == "HC")$Age) - mean(subset(df.demographics, Group == "MA")$Age)

t.test(BMI ~ Group, data = df.demographics)
mean(subset(df.demographics, Group == "HC")$BMI) - mean(subset(df.demographics, Group == "MA")$BMI)

t.test(PHQ.9 ~ Group, data = df.demographics)
mean(subset(df.demographics, Group == "HC")$PHQ.9) - mean(subset(df.demographics, Group == "MA")$PHQ.9)

t.test(OASIS ~ Group, data = df.demographics)
mean(subset(df.demographics, Group == "HC")$OASIS) - mean(subset(df.demographics, Group == "MA")$OASIS)

t.test(GAD7 ~ Group, data = df.demographics)
mean(subset(df.demographics, Group == "HC")$GAD7) - mean(subset(df.demographics, Group == "MA")$GAD7)

# omitting NA because STAI date one from HC participant is missing
t.test(STAI.State ~ Group, data = df.demographics)
mean(subset(df.demographics, Group == "HC")$STAI.State, na.rm = TRUE) - mean(subset(df.demographics, Group == "MA")$STAI.State)

# omitting NA because STAI date one from HC participant is missing
t.test(STAI.Trait ~ Group, data = df.demographics)
mean(subset(df.demographics, Group == "HC")$STAI.Trait, na.rm = TRUE) - mean(subset(df.demographics, Group == "MA")$STAI.Trait)

t.test(ASI.Total ~ Group, data = df.demographics)
mean(subset(df.demographics, Group == "HC")$ASI.Total) - mean(subset(df.demographics, Group == "MA")$ASI.Total)

t.test(ASI.Physical ~ Group, data = df.demographics)
mean(subset(df.demographics, Group == "HC")$ASI.Physical) - mean(subset(df.demographics, Group == "MA")$ASI.Physical)

t.test(ASI.Cognitive ~ Group, data = df.demographics)
mean(subset(df.demographics, Group == "HC")$ASI.Cognitive) - mean(subset(df.demographics, Group == "MA")$ASI.Cognitive)

t.test(ASI.Social ~ Group, data = df.demographics)
mean(subset(df.demographics, Group == "HC")$ASI.Social) - mean(subset(df.demographics, Group == "MA")$ASI.Social)

  
# compare groups based on average head motion during resting state scan
df.demographics %>% 
  group_by(Group) %>% 
  summarize(mean = mean(avg_motion),
            standard_error = sd(avg_motion)/sqrt(length(Group)))

# test if groups differ significantly on average head motion
# mean and 95% CI of mean average motion in GAD group
df.demographics %>% 
  filter(Group == "MA") %>% 
  {t.test(avg_motion ~ 1, data = .)}

# mean and 95% CI of mean average motion in HC group
df.demographics %>% 
  filter(Group == "HC") %>% 
  {t.test(avg_motion ~ 1, data = .)}

wilcox.test(avg_motion ~ Group, data = df.demographics, conf.int = TRUE, conf.level = 0.95)
paste0("Difference between group means (average motion): ",
       mean(subset(df.demographics, Group == "HC")$avg_motion) - mean(subset(df.demographics, Group == "MA")$avg_motion))