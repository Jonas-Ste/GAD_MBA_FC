# Load necessary libraries, install using 'install.packages("LIBRARY")' if needed
library(tidyverse)
library(DescTools)
library(broom)

# Load data. First command will prompt a window to locate the file with subjects demographics "subject_correlations.csv" on your machine
data.path <- file.choose()
df.correlations <- read_csv(data.path)

# dataframe contains Pearson's r each possible ROI combination for each individual subject
glimpse(df.correlations)
head(df.correlations)

# normalize Pearson's r using Fisher r-to-z transformation
df.correlations$cor_z <- FisherZ(df.correlations$cor_r)


#### Hypothesis 1: PCC - vmPFC - Decreased ####
# mean z-score and standard error of the mean
df.correlations %>% 
  filter(roi_1 == "vmPFC" & roi_2 == "PCC") %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(cor_z),
            SE = sd(cor_z)/sqrt(length(cor_z)))

# maximum number of hypotheses the tested ROIs are included in 
n_hyp <- 4

# t-test for group difference
df.correlations %>% 
  filter(roi_1 == "vmPFC" & roi_2 == "PCC") %>% 
  t.test(cor_z ~ Group, data = .) %>% 
  # display results in tidy format
  tidy() %>% 
  # Bonferroni-correct by number of hypotheses for tested ROIs
  mutate(p_adj = ifelse(p.value * n_hyp > 1, 1, p.value * n_hyp)) %>% 
  # order results and rename estimates for easy interpretation
  rename(df = parameter,
         t_value = statistic,
         mean_HC = estimate1,
         mean_GAD = estimate2,
         mean_difference = estimate) %>% 
  relocate(p_adj, .after = p.value) %>% 
  relocate(mean_difference, conf.low, conf.high, .after = mean_GAD)



#### Hypothesis 2: PCC - dmPFC - Decreased ####
# mean z-score and standard error of the mean
df.correlations %>% 
  filter(roi_1 == "dmPFC" & roi_2 == "PCC") %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(cor_z),
            SE = sd(cor_z)/sqrt(length(cor_z)))

# maximum number of hypotheses the tested ROIs are included in 
n_hyp <- 2

# t-test for group difference
df.correlations %>% 
  filter(roi_1 == "dmPFC" & roi_2 == "PCC") %>% 
  t.test(cor_z ~ Group, data = .) %>% 
  # display results in tidy format
  tidy() %>% 
  # Bonferroni-correct by number of hypotheses for tested ROIs
  mutate(p_adj = ifelse(p.value * n_hyp > 1, 1, p.value * n_hyp)) %>% 
  # order results and rename estimates for easy interpretation
  rename(df = parameter,
         t_value = statistic,
         mean_HC = estimate1,
         mean_GAD = estimate2,
         mean_difference = estimate) %>% 
  relocate(p_adj, .after = p.value) %>% 
  relocate(mean_difference, conf.low, conf.high, .after = mean_GAD)



#### Hypothesis 3: AI - dACC - Decreased ####
# mean z-score and standard error of the mean
df.correlations %>% 
  filter(roi_1 == "a_ins" & roi_2 == "dACC") %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(cor_z),
            SE = sd(cor_z)/sqrt(length(cor_z)))

# maximum number of hypotheses the tested ROIs are included in 
n_hyp <- 3

# t-test for group difference
df.correlations %>% 
  filter(roi_1 == "a_ins" & roi_2 == "dACC") %>% 
  t.test(cor_z ~ Group, data = .) %>% 
  # display results in tidy format
  tidy() %>% 
  # Bonferroni-correct by number of hypotheses for tested ROIs
  mutate(p_adj = ifelse(p.value * n_hyp > 1, 1, p.value * n_hyp)) %>% 
  # order results and rename estimates for easy interpretation
  rename(df = parameter,
         t_value = statistic,
         mean_HC = estimate1,
         mean_GAD = estimate2,
         mean_difference = estimate) %>% 
  relocate(p_adj, .after = p.value) %>% 
  relocate(mean_difference, conf.low, conf.high, .after = mean_GAD)



#### Hypothesis 4: dlPFC - Amygdala - Increased ####
# mean z-score and standard error of the mean
df.correlations %>% 
  filter(roi_1 == "amyg" & roi_2 == "dlPFC") %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(cor_z),
            SE = sd(cor_z)/sqrt(length(cor_z)))

# maximum number of hypotheses the tested ROIs are included in 
n_hyp <- 6

# t-test for group difference
df.correlations %>% 
  filter(roi_1 == "amyg" & roi_2 == "dlPFC") %>% 
  t.test(cor_z ~ Group, data = .) %>% 
  # display results in tidy format
  tidy() %>% 
  # Bonferroni-correct by number of hypotheses for tested ROIs
  mutate(p_adj = ifelse(p.value * n_hyp > 1, 1, p.value * n_hyp)) %>% 
  # order results and rename estimates for easy interpretation
  rename(df = parameter,
         t_value = statistic,
         mean_HC = estimate1,
         mean_GAD = estimate2,
         mean_difference = estimate) %>% 
  relocate(p_adj, .after = p.value) %>% 
  relocate(mean_difference, conf.low, conf.high, .after = mean_GAD)



#### Hypothesis 5: AI - vmPFC - Decreased ####
# mean z-score and standard error of the mean
df.correlations %>% 
  filter(roi_1 == "vmPFC" & roi_2 == "a_ins") %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(cor_z),
            SE = sd(cor_z)/sqrt(length(cor_z)))

# maximum number of hypotheses the tested ROIs are included in
n_hyp <- 4

# t-test for group difference
df.correlations %>% 
  filter(roi_1 == "vmPFC" & roi_2 == "a_ins") %>% 
  t.test(cor_z ~ Group, data = .) %>% 
  # display results in tidy format
  tidy() %>% 
  # Bonferroni-correct by number of hypotheses for tested ROIs
  mutate(p_adj = ifelse(p.value * n_hyp > 1, 1, p.value * n_hyp)) %>% 
  # order results and rename estimates for easy interpretation
  rename(df = parameter,
         t_value = statistic,
         mean_HC = estimate1,
         mean_GAD = estimate2,
         mean_difference = estimate) %>% 
  relocate(p_adj, .after = p.value) %>% 
  relocate(mean_difference, conf.low, conf.high, .after = mean_GAD)



#### Hypothesis 6: PMI - vmPFC - Decreased ####
# mean z-score and standard error of the mean
df.correlations %>% 
  filter(roi_1 == "vmPFC" & roi_2 == "pm_ins") %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(cor_z),
            SE = sd(cor_z)/sqrt(length(cor_z)))

# maximum number of hypotheses the tested ROIs are included in 
n_hyp <- 4

# t-test for group difference
df.correlations %>% 
  filter(roi_1 == "vmPFC" & roi_2 == "pm_ins") %>% 
  t.test(cor_z ~ Group, data = .) %>% 
  # display results in tidy format
  tidy() %>% 
  # Bonferroni-correct by number of hypotheses for tested ROIs
  mutate(p_adj = ifelse(p.value * n_hyp > 1, 1, p.value * n_hyp)) %>% 
  # order results and rename estimates for easy interpretation
  rename(df = parameter,
         t_value = statistic,
         mean_HC = estimate1,
         mean_GAD = estimate2,
         mean_difference = estimate) %>% 
  relocate(p_adj, .after = p.value) %>% 
  relocate(mean_difference, conf.low, conf.high, .after = mean_GAD)



#### Hypothesis 7: Amygdala - dACC - Increased ####
# mean z-score and standard error of the mean
df.correlations %>% 
  filter(roi_1 == "dACC" & roi_2 == "amyg") %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(cor_z),
            SE = sd(cor_z)/sqrt(length(cor_z)))

# maximum number of hypotheses the tested ROIs are included in 
n_hyp <- 6

# t-test for group difference
df.correlations %>% 
  filter(roi_1 == "dACC" & roi_2 == "amyg") %>% 
  t.test(cor_z ~ Group, data = .) %>% 
  # display results in tidy format
  tidy() %>% 
  # Bonferroni-correct by number of hypotheses for tested ROIs
  mutate(p_adj = ifelse(p.value * n_hyp > 1, 1, p.value * n_hyp)) %>% 
  # order results and rename estimates for easy interpretation
  rename(df = parameter,
         t_value = statistic,
         mean_HC = estimate1,
         mean_GAD = estimate2,
         mean_difference = estimate) %>% 
  relocate(p_adj, .after = p.value) %>% 
  relocate(mean_difference, conf.low, conf.high, .after = mean_GAD)



#### Hypothesis 8: Amygdala - TP - Increased ####
# mean z-score and standard error of the mean
df.correlations %>% 
  filter(roi_1 == "amyg" & roi_2 == "TP") %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(cor_z),
            SE = sd(cor_z)/sqrt(length(cor_z)))

# maximum number of hypotheses the tested ROIs are included in 
n_hyp <- 6

# t-test for group difference
df.correlations %>% 
  filter(roi_1 == "amyg" & roi_2 == "TP") %>% 
  t.test(cor_z ~ Group, data = .) %>% 
  # display results in tidy format
  tidy() %>% 
  # Bonferroni-correct by number of hypotheses for tested ROIs
  mutate(p_adj = ifelse(p.value * n_hyp > 1, 1, p.value * n_hyp)) %>% 
  # order results and rename estimates for easy interpretation
  rename(df = parameter,
         t_value = statistic,
         mean_HC = estimate1,
         mean_GAD = estimate2,
         mean_difference = estimate) %>% 
  relocate(p_adj, .after = p.value) %>% 
  relocate(mean_difference, conf.low, conf.high, .after = mean_GAD)



#### Hypothesis 9: Amygdala - vmPFC - Increased ####
# mean z-score and standard error of the mean
df.correlations %>% 
  filter(roi_1 == "vmPFC" & roi_2 == "amyg") %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(cor_z),
            SE = sd(cor_z)/sqrt(length(cor_z)))

# maximum number of hypotheses the tested ROIs are included in 
n_hyp <- 6

# t-test for group difference
df.correlations %>% 
  filter(roi_1 == "vmPFC" & roi_2 == "amyg") %>% 
  t.test(cor_z ~ Group, data = .) %>% 
  # display results in tidy format
  tidy() %>% 
  # Bonferroni-correct by number of hypotheses for tested ROIs
  mutate(p_adj = ifelse(p.value * n_hyp > 1, 1, p.value * n_hyp)) %>% 
  # order results and rename estimates for easy interpretation
  rename(df = parameter,
         t_value = statistic,
         mean_HC = estimate1,
         mean_GAD = estimate2,
         mean_difference = estimate) %>% 
  relocate(p_adj, .after = p.value) %>% 
  relocate(mean_difference, conf.low, conf.high, .after = mean_GAD)



#### Hypothesis 10: Amygdala - AI - Increased ####
# mean z-score and standard error of the mean
df.correlations %>% 
  filter(roi_1 == "a_ins" & roi_2 == "amyg") %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(cor_z),
            SE = sd(cor_z)/sqrt(length(cor_z)))

# maximum number of hypotheses the tested ROIs are included in 
n_hyp <- 6

# t-test for group difference
df.correlations %>% 
  filter(roi_1 == "a_ins" & roi_2 == "amyg") %>% 
  t.test(cor_z ~ Group, data = .) %>% 
  # display results in tidy format
  tidy() %>% 
  # Bonferroni-correct by number of hypotheses for tested ROIs
  mutate(p_adj = ifelse(p.value * n_hyp > 1, 1, p.value * n_hyp)) %>% 
  # order results and rename estimates for easy interpretation
  rename(df = parameter,
         t_value = statistic,
         mean_HC = estimate1,
         mean_GAD = estimate2,
         mean_difference = estimate) %>% 
  relocate(p_adj, .after = p.value) %>% 
  relocate(mean_difference, conf.low, conf.high, .after = mean_GAD)



#### Hypothesis 11: Amygdala - PMI - Increased ####
# mean z-score and standard error of the mean
df.correlations %>% 
  filter(roi_1 == "pm_ins" & roi_2 == "amyg") %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(cor_z),
            SE = sd(cor_z)/sqrt(length(cor_z)))

# maximum number of hypotheses the tested ROIs are included in 
n_hyp <- 6

# t-test for group difference
df.correlations %>% 
  filter(roi_1 == "pm_ins" & roi_2 == "amyg") %>% 
  t.test(cor_z ~ Group, data = .) %>% 
  # display results in tidy format
  tidy() %>% 
  # Bonferroni-correct by number of hypotheses for tested ROIs
  mutate(p_adj = ifelse(p.value * n_hyp > 1, 1, p.value * n_hyp)) %>% 
  # order results and rename estimates for easy interpretation
  rename(df = parameter,
         t_value = statistic,
         mean_HC = estimate1,
         mean_GAD = estimate2,
         mean_difference = estimate) %>% 
  relocate(p_adj, .after = p.value) %>% 
  relocate(mean_difference, conf.low, conf.high, .after = mean_GAD)
