# Load necessary libraries, install using 'install.packages("LIBRARY")' if needed
library(tidyverse)
library(DescTools)
library(broom)
library(BayesFactor)


# Load demographical data 
# First command will prompt a window to locate the file with subjects demographics "subject_info.csv" on your machine
demdata.path <- file.choose()
df.demographics <- read_csv(demdata.path)

# Load resting-state data
# First command will prompt a window to locate the file with subjects demographics "subject_correlation.csv" on your machine
corrdata.path <- file.choose()
df.correlations <- read_csv(corrdata.path)

# normalize Pearson's r using Fisher r-to-z transformation
df.correlations$cor_z <- FisherZ(df.correlations$cor_r)

# build new dataframe with results z-score of vmPFC-PMI connectivity and scores on clinical questionnaires
df.clinical <- df.correlations %>%
  filter(roi_1 == "vmPFC" & roi_2 == "pm_ins") %>% 
  inner_join(df.demographics, by = "ID") %>% 
  rename(Group = Group.x) %>% 
  select(-Group.y)


### correlations between vmPFC-PMI connectivity and clinical questionnaires are exploratory ###
### correlate for all clinical questionnaires for each group individually                   ###
### correct for multiple comparisons using FDR and save results in tidy dataframes          ###


### Correlations for GAD-Group ###

# cor: ASI-Total in GAD
clinical.cor.GAD <- 
  df.clinical %>%
  filter(Group == "MA") %>% 
  cor.test(.$cor_z, .$ASI.Total, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "ASI-Total",
         BF = NA_real_) %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter)
# calculate Bayes Factor
clinical.cor.GAD[clinical.cor.GAD$questionnaire == "ASI-Total", "BF"] <- df.clinical %>%
  filter(Group == "MA") %>% 
  {correlationBF(.$cor_z, .$ASI.Total)} %>% 
  extractBF() %>% 
  {.$bf}

# cor: GAD-7 in GAD
clinical.cor.GAD <- 
  df.clinical %>%
  filter(Group == "MA") %>% 
  cor.test(.$cor_z, .$GAD7, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "GAD-7") %>% 
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.GAD)

# calculate Bayes Factor
clinical.cor.GAD[clinical.cor.GAD$questionnaire == "GAD-7", "BF"] <- df.clinical %>%
  filter(Group == "MA") %>% 
  {correlationBF(.$cor_z, .$GAD7)} %>% 
  extractBF() %>% 
  {.$bf}

# cor: PHQ-9 in GAD
clinical.cor.GAD <- 
  df.clinical %>%
  filter(Group == "MA") %>% 
  cor.test(.$cor_z, .$PHQ.9, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "PHQ-9") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.GAD)

# calculate Bayes Factor
clinical.cor.GAD[clinical.cor.GAD$questionnaire == "PHQ-9", "BF"] <- df.clinical %>%
  filter(Group == "MA") %>% 
  {correlationBF(.$cor_z, .$PHQ.9)} %>% 
  extractBF() %>% 
  {.$bf}

# cor: STAI-State in GAD
clinical.cor.GAD <- 
  df.clinical %>%
  filter(Group == "MA") %>% 
  cor.test(.$cor_z, .$STAI.State, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "STAI-State") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.GAD)

# calculate Bayes Factor
clinical.cor.GAD[clinical.cor.GAD$questionnaire == "STAI-State", "BF"] <- df.clinical %>%
  filter(Group == "MA") %>% 
  {correlationBF(.$cor_z, .$STAI.State)} %>% 
  extractBF() %>% 
  {.$bf}

# cor: STAI-Trait in GAD
clinical.cor.GAD <- 
  df.clinical %>%
  filter(Group == "MA") %>% 
  cor.test(.$cor_z, .$STAI.Trait, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "STAI-Trait") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.GAD)

# calculate Bayes Factor
clinical.cor.GAD[clinical.cor.GAD$questionnaire == "STAI-Trait", "BF"] <- df.clinical %>%
  filter(Group == "MA") %>% 
  {correlationBF(.$cor_z, .$STAI.Trait)} %>% 
  extractBF() %>% 
  {.$bf}

# cor: OASIS in GAD
clinical.cor.GAD <- 
  df.clinical %>%
  filter(Group == "MA") %>% 
  cor.test(.$cor_z, .$OASIS, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "OASIS") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.GAD)
  

# calculate Bayes Factor
clinical.cor.GAD[clinical.cor.GAD$questionnaire == "OASIS", "BF"] <- df.clinical %>%
  filter(Group == "MA") %>% 
  {correlationBF(.$cor_z, .$OASIS)} %>% 
  extractBF() %>% 
  {.$bf}

# correct for multiplicity using the Bonferroni-method
clinical.cor.GAD <- 
  clinical.cor.GAD %>% 
  mutate(p.value.adj = p.adjust(.$p.value, method = "bonferroni"))


# reorder to match order in manuscript
clinical.cor.GAD <- 
  clinical.cor.GAD %>% 
  relocate(BF, .after = p.value) %>% 
  relocate(questionnaire, .before = pearson_r) %>% 
  relocate(p.value.adj, .after = p.value) %>% 
  arrange(-row_number())

clinical.cor.GAD


### Correlations for HC-Group ###

# cor: ASI-Total in HC
clinical.cor.HC <- 
  df.clinical %>%
  filter(Group == "HC") %>% 
  cor.test(.$cor_z, .$ASI.Total, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "ASI-Total",
         BF = NA_real_) %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter)

# calculate Bayes Factor
clinical.cor.HC[clinical.cor.HC$questionnaire == "ASI-Total", "BF"] <- df.clinical %>%
  filter(Group == "HC") %>% 
  {correlationBF(.$cor_z, .$ASI.Total)} %>% 
  extractBF() %>% 
  {.$bf}

# cor: GAD-7 in HC
clinical.cor.HC <- 
  df.clinical %>%
  filter(Group == "HC") %>% 
  cor.test(.$cor_z, .$GAD7, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "GAD-7") %>% 
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.HC)

# calculate Bayes Factor
clinical.cor.HC[clinical.cor.HC$questionnaire == "GAD-7", "BF"] <- df.clinical %>%
  filter(Group == "HC") %>% 
  {correlationBF(.$cor_z, .$GAD7)} %>% 
  extractBF() %>% 
  {.$bf}

# cor: PHQ-9 in HC
clinical.cor.HC <- 
  df.clinical %>%
  filter(Group == "HC") %>% 
  cor.test(.$cor_z, .$PHQ.9, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "PHQ-9") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.HC)

# calculate Bayes Factor
clinical.cor.HC[clinical.cor.HC$questionnaire == "PHQ-9", "BF"] <- df.clinical %>%
  filter(Group == "HC") %>% 
  {correlationBF(.$cor_z, .$PHQ.9)} %>% 
  extractBF() %>% 
  {.$bf}

# cor: STAI-State in HC  - missing data from one HC
clinical.cor.HC <- 
  df.clinical %>%
  filter(Group == "HC") %>% 
  cor.test(.$cor_z, .$STAI.State, method = "pearson", na.action = "na.omit", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "STAI-State") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.HC)

# calculate Bayes Factor
clinical.cor.HC[clinical.cor.HC$questionnaire == "STAI-State", "BF"] <- df.clinical %>%
  filter(Group == "HC") %>% 
  {correlationBF(.$cor_z, .$STAI.State)} %>% 
  extractBF() %>% 
  {.$bf}

# cor: STAI-Trait in HC - missing data from one HC
clinical.cor.HC <- 
  df.clinical %>%
  filter(Group == "HC") %>% 
  cor.test(.$cor_z, .$STAI.Trait, method = "pearson", na.action = "na.omit", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "STAI-Trait") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.HC)

# calculate Bayes Factor
clinical.cor.HC[clinical.cor.HC$questionnaire == "STAI-Trait", "BF"] <- df.clinical %>%
  filter(Group == "HC") %>% 
  {correlationBF(.$cor_z, .$STAI.Trait)} %>% 
  extractBF() %>% 
  {.$bf}

# cor: OASIS in HC
clinical.cor.HC <- 
  df.clinical %>%
  filter(Group == "HC") %>% 
  cor.test(.$cor_z, .$OASIS, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "OASIS") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.HC)

# calculate Bayes Factor
clinical.cor.HC[clinical.cor.HC$questionnaire == "OASIS", "BF"] <- df.clinical %>%
  filter(Group == "HC") %>% 
  {correlationBF(.$cor_z, .$OASIS)} %>% 
  extractBF() %>% 
  {.$bf}

# correct for multiplicity using the Bonferroni-method
clinical.cor.HC <- 
  clinical.cor.HC %>% 
  mutate(p.value.adj = p.adjust(.$p.value, method = "bonferroni"))

# reorder to match order in manuscript
clinical.cor.HC <- 
  clinical.cor.HC %>% 
  relocate(BF, .after = p.value) %>% 
  relocate(questionnaire, .before = pearson_r) %>% 
  relocate(p.value.adj, .after = p.value) %>% 
  arrange(-row_number())

clinical.cor.HC