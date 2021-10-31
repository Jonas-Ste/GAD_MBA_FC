# Load necessary libraries, install using 'install.packages("LIBRARY")' if needed
library(tidyverse)
library(DescTools)
library(broom)


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
  select(-Group) %>% 
  inner_join(df.demographics, by = "ID")


### correlations between vmPFC-PMI connectivity and clinical questionnaires are exploratory ###
### correlate for all clinical questionnaires for each group individually                   ###
### correct for multiple comparisons using FDR and save results in tidy dataframes          ###


### Correlations for GAD-Group ###

# cor: GAD-7 in GAD
clinical.cor.GAD <- 
  df.clinical %>%
  filter(Group == "MA") %>% 
  cor.test(.$cor_z, .$GAD7, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "GAD-7") %>% 
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter)

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

# cor: ASI-Total in GAD
clinical.cor.GAD <- 
  df.clinical %>%
  filter(Group == "MA") %>% 
  cor.test(.$cor_z, .$ASI.Total, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "ASI-Total") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.GAD)

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

# correct for multiple comparisons using FDR, reorder to match order in manuscript

clinical.cor.GAD <- 
  clinical.cor.GAD %>% 
  mutate(p_adj = p.adjust(.$p.value, method = "fdr")) %>%
  relocate(p_adj, .after = p.value) %>% 
  relocate(questionnaire, .before = pearson_r) %>%
  arrange(-row_number())

clinical.cor.GAD


### Correlations for HC-Group ###

# cor: GAD-7 in HC
clinical.cor.HC <- 
  df.clinical %>%
  filter(Group == "HC") %>% 
  cor.test(.$cor_z, .$GAD7, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "GAD-7") %>% 
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter)

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

# cor: ASI-Total in HC
clinical.cor.HC <- 
  df.clinical %>%
  filter(Group == "HC") %>% 
  cor.test(.$cor_z, .$ASI.Total, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "ASI-Total") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.HC)

# cor: STAI-State in HC
clinical.cor.HC <- 
  df.clinical %>%
  filter(Group == "HC") %>% 
  cor.test(.$cor_z, .$STAI.State, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "STAI-State") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.HC)

# cor: STAI-Trait in HC
clinical.cor.HC <- 
  df.clinical %>%
  filter(Group == "HC") %>% 
  cor.test(.$cor_z, .$STAI.Trait, method = "pearson", data = .) %>% 
  tidy() %>% 
  mutate(questionnaire = "STAI-Trait") %>%
  rename(pearson_r = estimate,
         t_value = statistic,
         df = parameter) %>% 
  bind_rows(clinical.cor.HC)

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

# correct for multiple comparisons using FDR, reorder to match order in manuscript

clinical.cor.HC <- 
  clinical.cor.HC %>% 
  mutate(p_adj = p.adjust(.$p.value, method = "fdr")) %>%
  relocate(p_adj, .after = p.value) %>% 
  relocate(questionnaire, .before = pearson_r) %>%
  arrange(-row_number())

clinical.cor.HC