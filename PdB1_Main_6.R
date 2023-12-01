# PdB1 BHT VTS experiment 1 analysis # 
# Bartosz Majchrowicz, majchrowicz.b@gmail.com #


# Load libraries and settings ---------------------------------------------

{
  # load libraries, use install.packages('package_name') if library not yet installed
  library(tidyverse)
  library(easystats)
  library(rstatix)
  library(ggplot2)
  library(ggdist)
  library(cowplot)
  library(readbulk)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(Hmisc)
  library(RColorBrewer) # display.brewer.all()
  library(ggthemes)
  library(grid)
  library(gridExtra)
  library(scales)
  library(knitr)
  library(faux)
  library(sjPlot)
  library(janitor)
  library(beepr)
  library(tictoc)
  library(ggeffects)

  # some settings for convenience
  options(scipen=999, width = 150)
  Sys.setenv(LANG = "en")
  theme_set(theme_bw(base_size = 14) +
              theme(panel.grid.minor = element_blank()))
  '%nin%' <- Negate("%in%")
  
  library(default)
  default(get_summary_stats) <- list(type = 'common') # get_summary_stats <- reset_default(get_summary_stats)
  default(kable) <- list(digits = 3, format = 'simple')
  
  myprint <- function(model, logOdds = F, sig = T) {
    print_myprint <- coef(summary(model)) %>% as_tibble(rownames='Term') %>% 
      left_join(logoddsratio_to_d(fixef(model)) %>% as_tibble(rownames='Term'), by = 'Term') %>% 
      mutate(OddsRatio = exp(Estimate), .after = 'Estimate', .keep = 'unused') %>% 
      relocate(value, .after = 'OddsRatio') %>% 
      rename(SE = 'Std. Error', z = 'z value', p = 'Pr(>|z|)', d = 'value') %>%
      mutate(p = round(p,5))
    
    if (sig == T){
      print_myprint <- print_myprint %>% 
        mutate(Sig = case_when(p <= 0.05 & p > 0.01 ~ '*', p <= 0.01 & p > 0.001 ~ '**',
                               p <= 0.001 ~ '***', TRUE ~ ''), .after = 'p')
    }
    if (logOdds == T){
      print_myprint <- print_myprint %>% 
        left_join(coef(summary(model)) %>% as_tibble(rownames='Term') %>% 
                    select(Term, Estimate)) %>% 
        relocate(Estimate, .before = 'OddsRatio') %>% rename(LogOdds = Estimate)
    }
    return(print_myprint %>% kable(digits = 3))
  }
  
}

# Vars to clear 
# rm(pi1,pi2,all1,all2,p0,p1,p2,p3,p4,p6,p5,b1,i1) # free up memo

# Vars to save
# save(p7,c3,d2,q2,m1, file = 'Rdata_1.RData')
# load('Rdata_1.RData')

# Data wrangling ----------------------------------------------------------

# read Biostat data
p0 <- read_bulk(directory = "C:/Users/barto/Psych/Badania/PdB Exp1/BHTandVTSb/data/data_biostat",
                subdirectories = F, extension = 'csv', verbose = F, fun = read.csv) %>% # load data from multiple files 
  as_tibble()

p1 <- p0 %>% 
  mutate_if(is.character, list(~na_if(.,""))) %>%  # replace blank spaces with NAs
  rename(token = id) %>% 
  rename(id = participant) %>%
  mutate(countBal = case_when(id %% 8 == 1 ~ 1,  # code counterbalance
                              id %% 8 == 2 ~ 2,
                              id %% 8 == 3 ~ 3,
                              id %% 8 == 4 ~ 4,
                              id %% 8 == 5 ~ 5, 
                              id %% 8 == 6 ~ 6,
                              id %% 8 == 7 ~ 7,
                              id %% 8 == 0 ~ 8,
                              TRUE ~ NA_real_), .after = 'id') %>% 
  mutate(group = case_when(id %% 2 == 1 ~ 0,  # code group (1: with control, 0: lack control)
                           id %% 2 == 0 ~ 1,
                           TRUE ~ NA_real_), .after = 'countBal') %>% 
  mutate(groupType = case_when(group == 0 ~ 'lack control',
                               group == 1 ~ 'full control'), .after = 'group') %>% 
  rename(type = Type, stimuli = Stimuli, location = Location, shape = Shape, task = Task, cue = Cue, block = Block, cuefol = cue_follow,
         sex = płeć, age = wiek, hand = dominująca.ręka, file = File,
         MotivSlider.response = WaznoscSlider.response, ProcSlider.response = ProcedureCheckQuestionSlider.response,
         lackControl.keys = lcontrol1.keys, lackControl.rt = lcontrol1.rt, lackControl.corr = lcontrol1.corr,
         lackControlFinal.keys = final_lack_control.keys, lackControlFinal.rt = final_lack_control.rt, lackControlFinal.corr = final_lack_control.corr,
         withControl.keys = key_control1.keys, withControl.rt = key_control1.rt, withControl.corr = key_control1,
         withControlFinal.keys = final_key_control1.keys, withControlFinal.rt = final_key_control1.rt, withControlFinal.corr = final_key_control1.corr) %>% 
  mutate(age = case_when(age == 'odmowa' ~ NA,
                         age == '.' ~ NA,
                         TRUE ~ age)) %>% 
  mutate(age = as.integer(age)) %>%
  mutate(sex = tolower(sex)) %>% 
  mutate(sex = case_when(sex == 'dsadsa' ~ NA,
                         sex == 'false' ~ NA,
                         sex == '.' ~ NA,
                         TRUE ~ sex)) %>% 
  mutate(sex = str_sub(sex, 1, 1)) %>% # extract 1st letter from string
  mutate(sex = case_when(sex == 'ż' ~ 'k',
                         sex == 'd' ~ 'k',
                         TRUE ~ sex)) %>% 
  mutate(source = 'biostat', .before = 'token')

biostat <- read_csv2('biostat_braki_uzupelnione_2.csv') %>%  # read missing demographic data
  select(-token) %>% 
  mutate(sex = str_sub(sex, 1, 1)) %>%
  mutate(age = as.integer(age),
         sex = tolower(sex)) %>% 
  mutate(source = 'biostat')  

# read IPs data
i0 <- read_bulk(directory = "C:/Users/barto/Psych/Badania/PdB Exp1/BHTandVTSb/data/data_ips",
                subdirectories = F, extension = 'csv', verbose = F, fun = read.csv) %>% # load data from multiple files 
  as_tibble() 

i1 <- i0 %>% 
  mutate_if(is.character, list(~na_if(.,""))) %>%  # replace blank spaces with NAs
  rename(id = participant) %>%
  mutate(countBal = case_when(id %% 8 == 1 ~ 1, id %% 8 == 2 ~ 2, id %% 8 == 3 ~ 3, id %% 8 == 4 ~ 4,
                              id %% 8 == 5 ~ 5, id %% 8 == 6 ~ 6, id %% 8 == 7 ~ 7, id %% 8 == 0 ~ 8,
                              TRUE ~ NA_real_), .after = 'id') %>% 
  mutate(group = case_when(id %% 2 == 1 ~ 0,  # code group (1: with control, 0: lack control)
                           id %% 2 == 0 ~ 1,
                           TRUE ~ NA_real_), .after = 'countBal') %>%
  mutate(groupType = case_when(group == 0 ~ 'lack control',
                               group == 1 ~ 'full control'), .after = 'group') %>% 
  rename(type = Type, stimuli = Stimuli, location = Location, shape = Shape, task = Task, cue = Cue, block = Block, cuefol = cue_follow,
         pseud = Pseudonim, sex = Płeć, age = Wiek, hand = Dominująca.ręka, file = File,
         MotivSlider.response = WaznoscSlider.response, ProcSlider.response = ProcedureCheckQuestionSlider.response,
         lackControl.keys = lcontrol1.keys, lackControl.rt = lcontrol1.rt, lackControl.corr = lcontrol1.corr,
         lackControlFinal.keys = final_lack_control.keys, lackControlFinal.rt = final_lack_control.rt, lackControlFinal.corr = final_lack_control.corr,
         withControl.keys = key_control1.keys, withControl.rt = key_control1.rt, withControl.corr = key_control1,
         withControlFinal.keys = final_key_control1.keys, withControlFinal.rt = final_key_control1.rt, withControlFinal.corr = final_key_control1.corr) %>% 
  mutate(age = if_else(age == '20 lat ', '20', age)) %>% 
  mutate(age = as.integer(age),
         sex = tolower(sex)) %>%
  mutate(sex = str_sub(sex, 1, 1)) %>% # extract 1st letter from string
  mutate(source = 'ips', .after = 'groupType')


# combine Biostats and Ips data

compare_df_cols(p1,i1, return = 'mismatch', bind_method = 'rbind')
intersect(unique(i1$id), unique(p1$id)) # mind same ids across ips and biostat datasets

pi1 <- p1 %>% 
  bind_rows(i1 %>% 
              mutate(id = as.integer(id+3000))) %>% # change ips ids to not duplicate ids
  relocate(pseud, .after = 'token')

compare_df_cols(p1,i1, pi1, return = 'mismatch', bind_method = 'rbind')
max(p1$id); max(i1$id); max(pi1$id)
nrow(p1) + nrow(i1) == nrow(pi1) # should be T

pi2 <- pi1 %>% # remove duplicated tokens; see duplicateCases
  filter(!(id == 605  & token == '18367609b03cd27a2acf6a7d855b9b34f7eebf80'),
         !(id == 1621 & token == '7694268b2c29afb39297350fd478fd6ee59c0844'),
         !(id == 847  & token == '8e90d20102aab0b1724980016abbf98ff6a63aa9'),
         !(id == 1805 & token == 'e5ad54d80f0dee6bbc80d7cc2aa61ddeec62d7fe'))

pi1 %>% group_by(id) %>% slice(1) %>% ungroup() %>% count()
pi2 %>% group_by(id) %>% slice(1) %>% ungroup() %>% count() # 327 in biostat


# remove 999s
# 
# remove999s <- pi1 %>% # filter any 999s
#   filter(if_any(ends_with('.corr'),
#                 ~.==999)) 
# 
# pi2 <- pi1 %>% anti_join(remove999s)  # remove any 999s from main df
# 
# nrow(pi1) - nrow(remove999s) == nrow(pi2) # check, T
# pi2 %>% filter(if_any(ends_with('.corr'), ~.==999)) # check, should be 0

# specify columns
bht_cols <- p1 %>% 
  select(c('id','countBal','group','source',  # specify BHT columns
           
           'lackControl.keys':'lackControl.corr',
           'LackControlLoopTrials.thisTrialN':'LackControlLoopTrials.ran', 
           'withControl.keys':'withControl.corr',
           'WithControlLoopTrials.thisTrialN':'WithControlLoopTrials.ran', 
           'CorrAns', 'Feedback', 
           
           'LackControlLoopProblems.thisTrialN':'LackControlLoopProblems.ran', 
           'lackControlFinal.keys':'lackControlFinal.corr',
           'WithControlLoopProblems.thisTrialN':'WithControlLoopProblems.ran',
           'withControlFinal.keys':'withControlFinal.corr',

           'FinalCorrAns', 'FinalFeedback',
           'BHTglobalClockTime'
           )) %>% 
  colnames()

q_cols <-  p1 %>% 
  select(c('id','countBal','group','source',
           contains('Slider.response'))) %>%  # specify Questionnaires columns
  colnames()

vts_cols <- p1 %>% 
  select(c('id','countBal','group','source',
           'type':'shape', 'task', 'cue','cuefol', 'block',  # specify VTS columns
           'correct_loc_left':'correct_shape_right', 
            contains('timebreak'),
           
           'nBlocksMixCb1.thisRepN':'nBlocksMixCb1.thisIndex',
           'nTrialsMixCb1.thisRepN':'nTrialsMixCb1.thisIndex',
           'nBlocksMixCb2.thisRepN':'nBlocksMixCb2.thisIndex',
           'nTrialsMixCb2.thisRepN':'nTrialsMixCb2.thisIndex',
           'nBlocksMixCb3.thisRepN':'nBlocksMixCb3.thisIndex',
           'nTrialsMixCb3.thisRepN':'nTrialsMixCb3.thisIndex',
           'nBlocksMixCb4.thisRepN':'nBlocksMixCb4.thisIndex',
           'nTrialsMixCb4.thisRepN':'nTrialsMixCb4.thisIndex',
           
           'nBlocksFreeCb1.thisRepN':'nBlocksFreeCb1.thisIndex',
           'nTrialsFreeCb1.thisRepN':'nTrialsFreeCb1.thisIndex',
           'nBlocksFreeCb2.thisRepN':'nBlocksFreeCb2.thisIndex',
           'nTrialsFreeCb2.thisRepN':'nTrialsFreeCb2.thisIndex',
           'nBlocksFreeCb3.thisRepN':'nBlocksFreeCb3.thisIndex',
           'nTrialsFreeCb3.thisRepN':'nTrialsFreeCb3.thisIndex',
           'nBlocksFreeCb4.thisRepN':'nBlocksFreeCb4.thisIndex',
           'nTrialsFreeCb4.thisRepN':'nTrialsFreeCb4.thisIndex',
           
           'both_shape_cb1_mix.keys', 'both_shape_cb1_mix.rt', 'both_shape_cb1_mix.corr', 
           'both_loc_cb1_mix.keys',   'both_loc_cb1_mix.rt',   'both_loc_cb1_mix.corr', 
           'both_shape_cb1_free.keys', 'both_shape_cb1_free.rt', 'both_shape_cb1_free.corr', 
           'both_loc_cb1_free.keys',   'both_loc_cb1_free.rt',   'both_loc_cb1_free.corr', 
           
           'both_shape_cb2_mix.keys', 'both_shape_cb2_mix.rt', 'both_shape_cb2_mix.corr',
           'both_loc_cb2_mix.keys',   'both_loc_cb2_mix.rt',   'both_loc_cb2_mix.corr',
           'both_shape_cb2_free.keys', 'both_shape_cb2_free.rt', 'both_shape_cb2_free.corr',
           'both_loc_cb2_free.keys',   'both_loc_cb2_free.rt',   'both_loc_cb2_free.corr',
           
           'both_shape_cb3_mix.keys', 'both_shape_cb3_mix.rt', 'both_shape_cb3_mix.corr',
           'both_loc_cb3_mix.keys',   'both_loc_cb3_mix.rt', 'both_loc_cb3_mix.corr',
           'both_shape_cb3_free.keys', 'both_shape_cb3_free.rt', 'both_shape_cb3_free.corr',
           'both_loc_cb3_free.keys',  'both_loc_cb3_free.rt', 'both_loc_cb3_free.corr',
           
           'both_shape_cb4_mix.keys', 'both_shape_cb4_mix.rt', 'both_shape_cb4_mix.corr',
           'both_loc_cb4_mix.keys',   'both_loc_cb4_mix.rt',   'both_loc_cb4_mix.corr',
           'both_shape_cb4_free.keys', 'both_shape_cb4_free.rt', 'both_shape_cb4_free.corr',
           'both_loc_cb4_free.keys',   'both_loc_cb4_free.rt',   'both_loc_cb4_free.corr',
  )) %>% 
  colnames()
  

# all1 <- pi2 %>% 
#   mutate(type = tolower(type)) %>% 
#   select(all_of(c(bht_cols, q_cols, vts_cols))) # select predefined columns
# 
# colnames(all1)
# unique(all1$id) # check ids

# Demographics ------------------------------------------------------------

d1 <- pi2 %>% group_by(id) %>% slice(1) %>% ungroup()

# d1 %>% select(id, token, sex, age) %>%
#   filter(is.na(age)) %>%
#   full_join(d1 %>%
#               select(id, token, sex, age) %>%
#               filter(is.na(sex))) %>%
#   write_csv('biostat_braki_2.csv')

d2 <- d1 %>% select(id, sex, age, source) %>% # join with biostat-provided missing data
  filter(id %nin% biostat$id) %>% 
  bind_rows(biostat) %>% 
  arrange(id)

d2 %>% get_summary_stats(age) %>% kable() # age
d2 %>% group_by(source) %>% get_summary_stats(age) # age
d2 %>% filter(is.na(age))
ggplot(d2, aes(x = age, fill = source)) + geom_histogram()

ggsave('age_hist.png', path = 'plots_1/', w=5,h=4)

d2 %>% filter(age>40)

d2 %>% filter(is.na(age)) %>% count() %>% rename('age unknown' = 1) %>%
  bind_cols(d2 %>% filter(age <= 40) %>% count() %>% rename('age<=40' = 1)) %>% 
  bind_cols(d2 %>% filter(age>40) %>% count() %>% rename('age>40' = 1)) %>% 
  kable()

d2 %>% count() # n
d2 %>% group_by(source) %>%  count()

d2 %>% freq_table(sex, na.rm=F) %>% kable() # sex
d2 %>% freq_table(source, sex, na.rm=F) %>% kable() # sex
d2 %>% filter(is.na(sex))
ggplot(d2, aes(x = age, fill = sex)) + geom_histogram()

d2 %>% filter(sex=='n')

# group assignment 
d1 %>% freq_table(groupType) %>% kable()  
d1 %>% freq_table(countBal) %>% kable()


# duplicates in biostat data
{p1 %>% group_by(id) %>% slice(1) %>% ungroup() %>% count() # 331 ids
p1 %>% group_by(token) %>% slice(1) %>% ungroup() %>% count() # 327 tokens

p1 %>% filter(is.na(token)) # no missing tokens
p1 %>% filter(is.na(id)) # no missing ids

p1 %>% group_by(id) %>% slice(1) %>% ungroup() %>% select(token) %>% # find duplicates
  group_by_all() %>%
  filter(n()>1) %>%
  ungroup() %>% arrange(token) %>% 
  slice(1,3,5,7) -> dupl

p1 %>% group_by(id) %>% slice(1) %>% ungroup() %>% select(id, token) %>% 
  filter(token %in% dupl$token) -> dupl2

p1 %>% group_by(id, token) %>% slice(1) %>% ungroup() %>% select(id, token, sex, age, file) %>% 
  mutate(
         file_nr = as.integer(str_sub(file, 1, 6)),
         file_date = str_sub(file, -27, -5),
         .keep = 'unused') %>%
  filter(id %in% dupl2$id) %>%
  arrange(token) -> dupl3
  
write_csv2(dupl3, 'duplicates.csv')

duplicateCases <- dupl3 %>% select(id, token) %>% 
  slice(2,4,6,8) # select 2nd case per token, to be removed
}

# BHT ---------------------------------------------------------------------

b1 <-  pi2 %>% 
  select(all_of(bht_cols)) # subset BHT columns

b2 <- pi2 %>% 
  select(all_of(bht_cols)) %>% 
  filter_at(vars(contains('ControlLoop')),
            any_vars(!is.na(.))) # remove rows with NAs (we infer which rows are practice based on presence of NAs in above-specified non-practice columns)

b2 %>% get_summary_stats(contains('rt')) %>% 
  ggplot(., aes(y = mean, x = variable)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.15) +
  # geom_point(size = 3)+
  geom_point(aes(size = n)) +
  # geom_text(aes(label = paste('N:', n)), nudge_x = .2, nudge_y = .2, show.legend = FALSE) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_size(breaks= c(1000,6000,12000)) +
  labs(y = 'RT', x = 'Response type')

b2 %>% # check corr logging
  filter(group == 1) %>% # with control
  select(id, group, CorrAns, withControl.keys, withControl.corr) %>% 
  mutate(checkCorr = case_when(withControl.keys == CorrAns ~ 1,
                               (withControl.keys != CorrAns ~ 0),
                               TRUE ~ NA_real_)) %>% 
  mutate(check = if_else(checkCorr == withControl.corr, 0, 1)) %>% 
  freq_table(check) # should be 0 in 100

table(b3$lackControl.corr); table(b3$withControl.corr); table(b3$lackControlFinal.corr); table(b3$withControlFinal.corr)

b2 %>% group_by(group) %>% get_summary_stats(contains('.corr'))

b2 %>% freq_table(group, lackControl.corr)
b2 %>% freq_table(group, withControl.corr)
b2 %>% freq_table(group, lackControlFinal.corr)
b2 %>% freq_table(group, withControlFinal.corr)

b2 %>% 
  select(contains('corr')) %>% 
  count(any_vars(is.na(.)))

b2 %>% get_summary_stats() %>% print(n=Inf)

# b3 <- b2 %>% # remove 999s
  

# Questionnaires ----------------------------------------------------------

q1 <- pi2 %>% 
  select(all_of(q_cols)) %>% 
  filter(if_any(ends_with('.response'), # exclude 'id':'group' from filtering
                ~ !is.na(.)))  # remove rows with NAs
  
colnames(q1) <- gsub("Slider.response", "", colnames(q1))
colnames(q1) <- tolower(colnames(q1))

q2 <- q1 %>% 
  group_by(id) %>% mutate(indx = row_number()) %>% 
  mutate(task = case_when(indx == 1 ~ 'bht',
                           indx == 2 ~ 'vts'),
         .keep = 'unused', .after = 'group') %>% 
  ungroup() %>% mutate(group = as.integer(group))


q3 <- q2 %>% 
  select(-proc) %>%
  pivot_longer(cols = c('mental':'motiv'), names_to = 'q', values_to = 'ans')

q3 %>% 
  group_by(task, group) %>% 
  get_summary_stats(ans)

colnames(q2)

# procedure check question
q1 %>% freq_table(source, proc) %>% 
  mutate(proc = if_else(proc == 1, 'no', 'yes')) %>% 
  rename(ans = proc, propWithinSource = prop) %>% kable()
(68+2)/358

# NASA
nasa_q <- c('mental', 'physical', 'hurried', 'successful', 'hard', 'insecure')

nasaAll <- q2 %>% 
  select(id, group, task, all_of(nasa_q)) %>% 
  rowwise(id, task) %>% 
  summarise(nasaSum = sum(c(mental, physical, hurried, successful, hard, insecure)),
            nasaMean=mean(c(mental, physical, hurried, successful, hard, insecure))) %>% 
  pivot_wider(names_from = task, values_from = c(nasaSum, nasaMean)) %>% 
  ungroup() %>%
  left_join(q2 %>% select(id,group) %>% group_by(id) %>% slice(1)) %>% relocate(group, .after = id)

nasaAll %>% 
  group_by(group) %>% 
  get_summary_stats(-id)

nasaAll %>% mutate(group = as.factor(group)) %>% 
  group_by(group) %>% 
  get_summary_stats(nasaSum_bht) %>% 
  ggplot(., aes(y = mean, x = group, group = group)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1)  +
  labs(y = 'NASA BHT sum', x = 'Group')

nasaAll %>% mutate(group = as.factor(group)) %>% 
  group_by(group) %>% 
  get_summary_stats(nasaMean_bht) %>% 
  ggplot(., aes(y = mean, x = group, group = group)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1)  +
  labs(y = 'NASA BHT mean', x = 'Group')

q2 %>% mutate(group = as.factor(group)) %>% 
  select(id, group, task, mental) %>% 
  group_by(group) %>% 
  get_summary_stats(mental) %>% 
  ggplot(.,aes(y = mean, x = group, group = group)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1)  +
  labs(y = 'NASA BHT mental item', x = 'Group')

nasaAll %>% mutate(group = as.factor(group)) %>% 
  group_by(group) %>% 
  ggplot(., aes(x= nasaMean_bht)) +
  facet_wrap(~group) +
  geom_histogram()

# NASA factor analysis

nasaQs <- q2 %>% 
  select(all_of(nasa_q))

factanal(nasaQs, factors = 2) # factor 1: item 1,2,3,5
factanal(nasaQs, factors = 3) # factor 1: item 2,3,6, factor 2: item 1,5

factanal(nasaQs, factors = 2, rotation = 'promax') 
factanal(nasaQs, factors = 3, rotation = 'promax') 

nasa3q <- q2 %>% # use 1st out of 2 factors only
  select(id, group, task, mental, physical, hurried, hard) %>% 
  rowwise(id, task) %>% 
  summarise(nasaSum = sum(c(mental, physical, hurried, hard)),
            nasaMean=mean(c(mental, physical, hurried, hard))) %>% 
  pivot_wider(names_from = task, values_from = c(nasaSum, nasaMean)) %>% 
  ungroup() %>%
  left_join(q2 %>% select(id,group) %>% group_by(id) %>% slice(1)) %>% relocate(group, .after = id) %>%
  mutate(nasa_bht = nasaMean_bht, nasa_vts = nasaMean_vts, .after = group, .keep = 'unused')

nasa3q %>% mutate(group = as.factor(group)) %>% 
  group_by(group) %>% 
  get_summary_stats(nasaSum_bht) %>% 
  ggplot(., aes(y = mean, x = group, group = group)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1)  +
  labs(y = 'NASA BHT sum', x = 'Group')

nasa3q %>% mutate(group = as.factor(group)) %>% 
  group_by(group) %>% 
  get_summary_stats(nasaMean_bht) %>% 
  ggplot(., aes(y = mean, x = group, group = group)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1)  +
  labs(y = 'NASA BHT mean', x = 'Group')

# VTS ---------------------------------------------------------------------

p2 <- pi2 %>% # remove practice trials
  select(all_of(vts_cols), groupType) %>% 
  filter_at(vars(
    'both_shape_cb1_mix.corr', 'both_loc_cb1_mix.corr', 'both_shape_cb2_mix.corr', 'both_loc_cb2_mix.corr',
    'both_shape_cb3_mix.corr', 'both_loc_cb3_mix.corr', 'both_shape_cb4_mix.corr', 'both_loc_cb4_mix.corr',
    'both_shape_cb1_free.corr', 'both_loc_cb1_free.corr', 'both_shape_cb2_free.corr', 'both_loc_cb2_free.corr',
    'both_shape_cb3_free.corr', 'both_loc_cb3_free.corr', 'both_shape_cb4_free.corr', 'both_loc_cb4_free.corr',
    'nBlocksMixCb1.thisTrialN', 'nBlocksMixCb2.thisTrialN', 'nBlocksMixCb3.thisTrialN', 'nBlocksMixCb4.thisTrialN',
    'nBlocksFreeCb1.thisTrialN', 'nBlocksFreeCb2.thisTrialN', 'nBlocksFreeCb3.thisTrialN', 'nBlocksFreeCb4.thisTrialN',
  ),
  any_vars(!is.na(.))) # remove rows with NAs (we infer which rows are practice based on presence of NAs in above-specified non-practice columns)

p3 <- p2  %>%  # collapse (unite) block number (so that it's coded in a single column, not split to separate columns based on cb)
  bind_cols(p2 %>% 
              select(starts_with('nBlocks') & ends_with('.thisRepN')) %>% 
              unite('blockNr', everything(), sep='', na.rm=T) %>% # unite columns
              mutate(blockNr = as.integer(blockNr) + 1) %>% 
              fill(blockNr, .direction = "up")) %>% # fill-up block nr from row in which it was specified to all other rows (trials)
  relocate(blockNr, .after = block) # relocate for convenience

p3 %>% group_by(id) %>% count() # check nr of trials per id
p3 %>% group_by(id) %>% count() %>% filter(n != 292) # if 999 not removed!

p3 %>% group_by(id) %>% count(blockNr) # check nr of trials per id and block
p3 %>% group_by(id) %>% count(blockNr) %>% filter(n != 146) # if 999 not removed!

p4 <- p3 %>% 
  filter_at(vars('type':'block'), any_vars(!is.na(.))) # remove remaining rows with NAs (already used to fill-up block nr)

p5 <- p4 %>% 
  select(-c('correct_loc_left':'correct_shape_right', 'stimuli')) # remove useless columns
      
p6 <- p5 %>% # get collapsed task per trial across tasks and cb's (actually performed and cued) - quite slow part of the code
  bind_cols(p5 %>% 
              select(starts_with('both_loc_') & ends_with('.keys')) %>%  # selects loc task across all cbs and free/mixed blocks
              mutate(across(everything(), ~replace(., . == 'None', NA))) %>%  # replace 'None' with NAs (case_when works too but is slower)
              mutate(task1 = case_when(if_all(everything(), is.na) ~ 'shape')) %>% 
              bind_cols(p5 %>% # bind columns processed similarly as those above
                          select(starts_with('both_shape_') & ends_with('.keys')) %>%   
                          mutate(across(everything(), ~replace(., . == 'None', NA))) %>%   
                          mutate(task2 = case_when(if_all(everything(), is.na) ~ 'loc')) %>% 
                          select(task2)) %>% 
              unite(task_perf, c(task1,task2), sep = "", na.rm = TRUE) %>% 
              select(task_perf)) %>% 
  # select(!(starts_with('both_loc_') & ends_with('.keys'))) %>%  # remove columns that have just been united
  mutate(task_cued = case_when(cue == 'JAKI' ~ 'shape', # actually just renamed 'cue' column
                               cue == 'GDZIE'~ 'loc',
                               cue == 'XXXX' ~ 'any',
                               TRUE ~ 'unknownTask')) %>% 
  relocate(c(task_cued, task_perf), .after = 'block') # relocate for convenience

p5a <- p5 %>%  # collapse (unite) all trial and block data across tasks and cb's
  select(starts_with('nBlocks') & ends_with('.thisRepN')) %>% # block columns
  unite('nBlocks.thisRepN', everything(), sep='', na.rm=T) %>% 
  bind_cols(p5 %>% 
              select(starts_with('nBlocks') & ends_with('.thisTrialN')) %>% 
              unite('nBlocks.thisTrialN', everything(), sep='', na.rm=T)) %>% 
  bind_cols(p5 %>% 
              select(starts_with('nBlocks') & ends_with('.thisN')) %>% 
              unite('nBlocks.thisN', everything(), sep='', na.rm=T)) %>% 
  bind_cols(p5 %>% 
              select(starts_with('nBlocks') & ends_with('.thisIndex')) %>% 
              unite('nBlocks.thisIndex', everything(), sep='', na.rm=T)) %>% 
  bind_cols(p5 %>% 
              select(starts_with('nTrials') & ends_with('.thisRepN')) %>% # trial columns
              unite('nTrials.thisRepN', everything(), sep='', na.rm=T) %>% 
              bind_cols(p5 %>% 
                          select(starts_with('nTrials') & ends_with('.thisTrialN')) %>% 
                          unite('nTrials.thisTrialN', everything(), sep='', na.rm=T)) %>% 
              bind_cols(p5 %>% 
                          select(starts_with('nTrials') & ends_with('.thisN')) %>% 
                          unite('nTrials.thisN', everything(), sep='', na.rm=T)) %>% 
              bind_cols(p5 %>% 
                          select(starts_with('nTrials') & ends_with('.thisIndex')) %>% 
                          unite('nTrials.thisIndex', everything(), sep='', na.rm=T))) %>% 
  mutate_if(is.character,as.integer)

p5b <- p5 %>% # collapse response (keys) data
  select(starts_with('both_') & ends_with('.keys')) %>% 
  mutate(across(everything(), ~replace(., . == 'None', NA))) %>% 
  unite(resp, everything(), sep='', na.rm=T) %>% 
  bind_cols(p5 %>% # collapse rt data
              select(starts_with('both_') & ends_with('.rt')) %>% 
              unite(rt, everything(), sep='', na.rm=T) %>% 
              mutate(rt = as.numeric(rt))) %>% 
  bind_cols(p5 %>% # collapse accuracy (corr) data
              select(starts_with('both_') & ends_with('.corr')) %>% 
              mutate(across(everything(), ~replace(., . == '888', NA))) %>% 
              unite(corr, everything(), sep='', na.rm=T) %>% 
              mutate(corr = as.integer(corr)) )

nrow(p6) == nrow(p5a) && nrow(p6) == nrow(p5b) # before merging make sure nr of rows is identical (should be TRUE)

p5b %>% filter(is.na(rt))

p7 <- p6 %>% 
  select(id:blockNr) %>% # select remaining non-redundant columns 
  bind_cols(p5a, p5b) %>%  # bind with collapsed trial, block, resp, rt and corr dataframes 
  group_by(id, blockNr) %>% # grouping needed to code task switch type (to not consider block starts)
  mutate(switch_type = case_when(task_perf == 'shape' & lag(task_perf) == 'shape' ~ 'stay_shape', # code task switch type
                                 task_perf == 'loc' & lag(task_perf) == 'loc' ~ 'stay_loc',
                                 task_perf == 'shape' & lag(task_perf) == 'loc' ~ 'switch_shape',
                                 task_perf == 'loc' & lag(task_perf) == 'shape' ~ 'switch_loc')) %>%  
  mutate(switch = case_when(switch_type == 'stay_shape' | switch_type == 'stay_loc' ~ 0, # code switch logical
                            switch_type == 'switch_shape' | switch_type == 'switch_loc' ~ 1,
                            TRUE ~ NA_real_)) %>% 
  relocate(resp:switch, .after = blockNr) %>% 
  relocate(cuefol, .after = task_perf) %>% 
  ungroup() %>% 
  filter(corr < 999) %>% # remove unknown corr - - - - - - - - - - - - - - - - - <- consider!
  filter(!nchar(resp) > 1) %>% # remove rare cases of two responses registered
  mutate(trialNr = nTrials.thisN + 1, .after = blockNr) %>% # recode trial nr per block
  select(-c('cue', 'nBlocks.thisRepN':'nTrials.thisIndex')) %>% # again remove some columns
  rename(trial_type = task) %>% 
  mutate(cBal = case_when( # simplified counterbalance (blocks order only)
    countBal == 1 | countBal == 2 | countBal == 5 | countBal == 6 ~ 'mix-free', # mix > free 
    countBal == 3 | countBal == 4 | countBal == 7 | countBal == 8 ~ 'free-mix', # free > mix
    TRUE ~ NA), .after = countBal)  %>% 
  mutate(blockNrC = case_when(cBal == 'mix-free' & (block == 'mixed' & blockNr == 1) ~ 1, # combined block nr, accounting for block type
                              cBal == 'mix-free' & (block == 'mixed' & blockNr == 2) ~ 2,
                              cBal == 'mix-free' & (block == 'freeOnly' & blockNr == 1) ~ 3,
                              cBal == 'mix-free' & (block == 'freeOnly' & blockNr == 2) ~ 4,
                              cBal == 'free-mix' & (block == 'freeOnly' & blockNr == 1) ~ 1,
                              cBal == 'free-mix' & (block == 'freeOnly' & blockNr == 2) ~ 2,
                              cBal == 'free-mix' & (block == 'mixed' & blockNr == 1) ~ 3,
                              cBal == 'free-mix' & (block == 'mixed' & blockNr == 2) ~ 4,
                              TRUE ~ NA_real_))

p7 %>% freq_table(blockNrC, na.rm = F)
p7 %>% freq_table(countBal, na.rm = F)

# Post-wrangling checks and descriptives

colnames(p7)
p0 %>% distinct(frameRate, .keep_all = TRUE) %>% select(participant, frameRate) %>% arrange(frameRate) %>% 
  ggplot(., aes(x=frameRate)) + 
  geom_histogram(color = 'grey20') +
  geom_vline(xintercept = c(60,75,144), linetype = 'dashed') # most common refresh rates are 60Hz, 75Hz, 144Hz, and 240Hz

# trials nr, counts
p7 %>% group_by(id) %>% count() # mind if you removed corr 999 cases
p7 %>% group_by(id, block) %>% count() %>% pivot_wider(names_from = block, values_from = n) %>% print(n=Inf)

p7 %>% group_by(id) %>% count(blockNr)
p7 %>% count(block)

p7 %>% freq_table(resp)
p7 %>% freq_table(corr)
p7 %>% freq_table(task_perf)
p7 %>% freq_table(block, task_perf)
p7 %>% freq_table(task_cued)


# Durations ---------------------------------------------------------------


t1 <- pi2 %>% select(id) %>% # time of between-blocks breaks 
  bind_cols(pi2 %>% 
              select(starts_with('timebreak')) %>% 
              unite('timeBreak', everything(), sep='', na.rm=T)) %>% 
  mutate(timeBreak = (as.numeric(timeBreak))/60) %>% 
  filter(!is.na(timeBreak))
get_summary_stats(t1 %>% select(timeBreak))

t2 <- pi2 %>% select(id, globalClockTime) %>% # total time
  filter(!is.na(globalClockTime)) %>% 
  mutate(globalClockTime = globalClockTime/60) 

t3 <- pi2 %>%
  select(id, BHTglobalClockTime ) %>% # 4 times: BHT, VTS, between tasks, total
  filter(!is.na(BHTglobalClockTime)) %>%
  left_join(pi2 %>%
              select(id, BeforeVTSglobalClockTime) %>% 
              filter(!is.na(BeforeVTSglobalClockTime))) %>% 
  left_join(pi2 %>%
              select(id, globalClockTime ) %>% 
              filter(!is.na(globalClockTime))) %>% 
  rename(postBHT = BHTglobalClockTime, preVTS = BeforeVTSglobalClockTime, postVTS = globalClockTime) %>% 
  mutate(betweenTasks = preVTS - postBHT,
         pureVTS = postVTS - preVTS) %>% 
  mutate_at(vars(-id), funs(./60)) %>% # convert to minutes
  left_join(tibble(id = allOutIDs4c$id, # add ID outlier flags (but note that allOutIDs is created later in the script!)
                   out = T, 
                   outType = allOutIDs4c$outType)) %>% 
  mutate(out = if_else(is.na(out), F, T))

tibble(id = c(1002, insOutIDs1$id, insOutIDs2$id, durOutIDs4c$id), out = T)
allOutIDs5c
allOutIDs4c

t3b <- t1 %>%  # add breaks
  bind_cols(t1 %>% arrange(-timeBreak) %>% get_summary_stats(timeBreak) %>% select(mean,sd) %>% 
              mutate(up = mean+(sd*3), lo = mean-(sd*3), .keep = 'unused')) %>% 
  mutate(out_break = if_else(timeBreak > up | timeBreak < lo, T, F)) %>% 
  arrange(-timeBreak)

t4 <- t3 %>% 
  left_join(t3b %>% # out = all IDs to be removed (all durations, instr, age); out_break = break outliers only
              group_by(id) %>% 
              summarise(max_break= max(timeBreak)) %>% 
              left_join(t3b %>% filter(out_break == T) %>% 
                          group_by(id) %>% 
                          summarise(out_break= max(timeBreak)) %>% 
                          mutate(out_break = T)) %>% 
              mutate(out_break = if_else(is.na(out_break), F, out_break))) %>% 
  relocate(out, .before = out_break)
  
# some checks
breaksOutIDs <- t4 %>% filter(out_break==T) %>% distinct(id)
t4 %>% filter(out_break ==T)
t3 %>% filter(out==T) # nrow 41 including, 34 excluding break time outliers
t4 %>% filter(out==T)
t4 %>% filter(out==T & outType == 'dur')
t4 %>% filter(out_break==T & out == F & outType == 'dur')
t4 %>% filter(out_break==T & out == T & outType == 'dur')
t4 %>% filter(out==F & out_break==T) # break time only outliers
allOutIDs5c %>% group_by(id) %>% count() %>% filter(n>1) # no subject is outType of more than 1 type, in both 4c and 5c

# Plots

# break durations
pt1 <- t1 %>% 
  left_join(t2) %>% 
  mutate(id = fct_reorder(as.factor(id), -desc(globalClockTime))) %>% 
  ggplot(., aes(y = timeBreak, x = as.factor(id))) +
  geom_point(alpha = .3) +
  labs(y = 'Break duration (min)', x = 'ID') +
  geom_hline(yintercept = durVals[5], linetype = 'dashed', colour = 'blue1') +
  scale_y_continuous(breaks=seq(0,35,5))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())

pt1o <- t1 %>% 
  left_join(t2) %>% 
  left_join(t3) %>% 
  mutate(id = fct_reorder(as.factor(id), -desc(globalClockTime))) %>% 
  ggplot(., aes(y = timeBreak, x = as.factor(id), colour = out)) +
  geom_point(alpha = .3) +
  labs(y = 'Break duration (min)', x = 'ID') +
  geom_hline(yintercept = durVals[5], linetype = 'dashed', colour = 'blue1') +
  scale_y_continuous(breaks=seq(0,35,5)) +
  scale_colour_manual(values = c('black','red'), name = 'Outlier') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())

# total experiment time (postVTS)
pt2 <- t3 %>% 
  mutate(id = fct_reorder(as.factor(id), -desc(postVTS))) %>% 
  ggplot(., aes(y = postVTS, x = as.factor(id))) +
  geom_point(size = 3, alpha = .3) +
  labs(y = 'Experiment duration (min)', x = 'ID') +
  geom_hline(yintercept = c(durVals[4], durVals2[4]), linetype = 'dashed', colour = c('blue1', 'blue4')) +
  scale_y_continuous(breaks=seq(20,160,20)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())

pt2o <- t3 %>% 
  mutate(id = fct_reorder(as.factor(id), -desc(postVTS))) %>% 
  ggplot(., aes(y = postVTS , x = as.factor(id), colour = out)) +
  geom_point(size = 3, alpha = .3) +
  labs(y = 'Experiment duration (min)', x = 'ID') +
  scale_colour_manual(values = c('black','red'), name = 'Outlier') +
  geom_hline(yintercept = c(durVals[4], durVals2[4]), linetype = 'dashed', colour = c('blue1', 'blue4')) +
  scale_y_continuous(breaks=seq(20,160,20)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())

pta1 <- grid.arrange(pt1, pt2, nrow = 2)

ggsave('dursId.png', pta1, path = 'plots_1/', w=5,h=8)
 
# VTS vs BHT
pt3 <- ggplot(t3, aes(x = postBHT, y = pureVTS)) +
  geom_point(size = 2, alpha = .3) +
  labs(x = 'BHT time (min)', y = 'VTS time (min)') +
  scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks())

pt3o <- ggplot(t3, aes(x = postBHT, y = pureVTS, colour = out)) +
  geom_point(size = 2, alpha = .3) +
  scale_colour_manual(values = c('black','red'), name = 'Outlier') +
  labs(x = 'BHT time (min)', y = 'VTS time (min)') +
  scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks())

# time between the tasks
pt4 <- t3 %>% 
  mutate(id = fct_reorder(as.factor(id), -desc(postVTS))) %>% # order by experiment time
  ggplot(., aes(y = betweenTasks, x = as.factor(id))) +
  geom_point(size = 3, alpha = .3) +
  labs(y = 'Time between the tasks (min)', x = 'ID') +
  geom_hline(yintercept = c(durVals[3], durVals2[3]), linetype = 'dashed', colour = c('blue1', 'blue4')) +
  scale_y_continuous(breaks=seq(0,80,10)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())

pt4o <- t3 %>% 
  mutate(id = fct_reorder(as.factor(id), -desc(postVTS))) %>% # order by experiment time
  ggplot(., aes(y = betweenTasks, x = as.factor(id), colour = out)) +
  geom_point(size = 3, alpha = .3) +
  scale_colour_manual(values = c('black','red'), name = 'Outlier') +
  labs(y = 'Time between the tasks (min)', x = 'ID') +
  geom_hline(yintercept = c(durVals[3], durVals2[3]), linetype = 'dashed', colour = c('blue1', 'blue4')) +
  scale_y_continuous(breaks=seq(0,80,10)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())

# arrange plots
pta2 <- grid.arrange(pt1, pt4, pt2, pt3,   
                     bottom = textGrob("Durations data. Subjects ordered by total experiment time.", gp = gpar(fontsize = 18)),
                     nrow = 2)

ggsave('dursId.png', pta2, path = 'plots_1/', w=16,h=8)

pta2o <- grid.arrange(pt1o + theme(legend.position="none"), pt4o+ theme(legend.position="none"), 
                      pt2o+ theme(legend.position="none"), pt3o+ theme(legend.position=c(0.9, 0.8),
                                                                       legend.background = element_rect(fill = "white", color = "grey80")),   
                     bottom = textGrob("Durations data. Subjects ordered by total experiment time.", gp = gpar(fontsize = 18)),
                     nrow = 2)

ggsave('dursId_out.png', pta2o, path = 'plots_1/', w=16,h=8)

t3 %>% cor_test('postBHT', 'pureVTS')


# histograms, post cleaning

ggplot(t1, aes(x=timeBreak*60)) + 
  geom_histogram(color = 'grey20', binwidth = 2) +
  geom_vline(aes(xintercept=median(timeBreak)), linetype = 'dashed') +
  annotate('text', label = paste('Med = ', t1 %>% get_summary_stats(timeBreak) %>% select(median) %>% round(1), '(s)\n',
                                 'Mean = ', t1 %>% get_summary_stats(timeBreak) %>% select(mean) %>% round(1), '(s)'),
           x=45, y=10, size = 5) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = 'Count', x = 'Break time (s)') 

ggsave('breakDur.png', path = 'plots_1/', w=5,h=4)

ggplot(t1, aes(x = timeBreak)) + 
  geom_histogram(color = "grey20", fill = 'grey50', binwidth = 0.01) +
  scale_x_continuous(breaks = pretty_breaks()) + xlab('timeBreak (min)') +
  geom_vline(xintercept = (durVals[5]), linetype = 'dashed', colour = 'blue1')

ggsave('timeBreak_hist.png', path = 'plots_1/', w=10,h=5)

ggplot(t3, aes(x = postBHT)) + 
  geom_histogram(color = "grey20", fill = 'grey50', binwidth = 0.25) +
  scale_x_continuous(breaks = pretty_breaks()) +
  geom_vline(xintercept = c(durVals[1], durVals2[1]), linetype = 'dashed', colour = c('blue1', 'blue4'))

ggsave('postBHT_hist.png', path = 'plots_1/', w=10,h=5)

ggplot(t3, aes(x = pureVTS)) + 
  geom_histogram(color = "grey20", fill = 'grey50', binwidth = 0.25) +
  scale_x_continuous(breaks = pretty_breaks()) +
  geom_vline(xintercept = c(durVals[2], durVals2[2]), linetype = 'dashed', colour = c('blue1', 'blue4'))

ggsave('pureVTS_hist.png', path = 'plots_1/', w=10,h=5)

ggplot(t3, aes(x = betweenTasks)) + 
  geom_histogram(color = "grey20", fill = 'grey50', binwidth = 0.25) +
  scale_x_continuous(breaks = seq(0,90,10)) +
  geom_vline(xintercept = c(durVals[3], durVals2[3]), linetype = 'dashed', colour = c('blue1', 'blue4'))

ggsave('betweenTasks_hist.png', path = 'plots_1/', w=10,h=5)

ggplot(t3, aes(x = postVTS)) + 
  geom_histogram(color = "grey20", fill = 'grey50', binwidth = 0.25) +
  scale_x_continuous(breaks = seq(0,90,10)) +
  geom_vline(xintercept = c(durVals[4], durVals2[4]), linetype = 'dashed', colour = c('blue1', 'blue4'))

ggsave('postVTS_hist.png', path = 'plots_1/', w=10,h=5)

ggplot(t2, aes(x = globalClockTime)) +
  geom_histogram(color = 'grey20', binwidth = 0.1) +
  geom_vline(aes(xintercept=median(globalClockTime)), linetype = 'dashed') +
  annotate('text', label = paste('Med = ', t2 %>% get_summary_stats(globalClockTime) %>% select(median) %>% round(1), '(min)\n',
                                 'Mean = ', t2 %>% get_summary_stats(globalClockTime) %>% select(mean) %>% round(1), '(min)'),
           x = 27, y= 2.1, size = 5) +
  labs(y = 'Count', x = 'Experiment time (min)') +
  scale_x_continuous(breaks = pretty_breaks())

ggsave('expDur.png', path = 'plots_1/', w=5,h=4)

ggplot(t3, aes(x = postVTS, fill = out)) + 
  geom_histogram(binwidth = 0.25) +
  scale_x_continuous(breaks = seq(0,90,10)) +
  scale_fill_manual(values = c('grey50','red'), name = 'Outlier') +
  geom_vline(xintercept = c(durVals[4], durVals2[4]), linetype = 'dashed', colour = c('blue1', 'blue4'))


# Durations outliers

durOut <- t3 %>% 
  left_join(t3 %>% select(id, postBHT) %>% 
              identify_outliers(postBHT) %>% 
              rename(extreme_postBHT = is.extreme) %>% 
              filter(extreme_postBHT == T) %>% 
              select(1,2,4)) %>% 
  left_join(t3 %>% select(id, betweenTasks) %>% 
              identify_outliers(betweenTasks) %>% 
              rename(extreme_betweenTasks = is.extreme) %>% 
              filter(extreme_betweenTasks == T) %>% 
              select(1,2,4)) %>% 
  left_join(t3 %>% select(id, pureVTS) %>% 
              identify_outliers(pureVTS) %>% 
              rename(extreme_pureVTS = is.extreme) %>% 
              filter(extreme_pureVTS == T) %>% 
              select(1,2,4)) %>% 
  left_join(t3 %>% select(id, postVTS) %>% 
              identify_outliers(postVTS) %>% 
              rename(extreme_postVTS = is.extreme) %>% 
              filter(extreme_postVTS == T) %>% 
              select(1,2,4)) %>% ungroup() %>% 
  mutate(across(contains('extreme'),
                ~ case_when(is.na(.) ~ F, # or: ~ replace(., is.na(.), F)))
                            T ~ .))) 

# durOut %>% 
#   group_by(extreme_postBHT ) %>%
#   get_summary_stats(postBHT) %>% 
#   select(1:n, min, max) %>% 
#   bind_rows()
# durOut %>% 
#   group_by(extreme_betweenTasks) %>%
#   get_summary_stats(betweenTasks) %>% 
#   select(1:n, min, max)
# durOut %>% 
#   group_by(extreme_pureVTS) %>%
#   get_summary_stats(pureVTS) %>% 
#   select(1:n, min, max)
# durOut %>% 
#   group_by(extreme_postVTS) %>%
#   get_summary_stats(postVTS) %>% 
#   select(1:n, min, max)

durOut2 <- durOut %>% 
  filter(if_any(contains('extreme'),
                ~ .==T))

durOut3 <- durOut %>% 
  mutate(across(postBHT:pureVTS, 
                list(mean = mean, sd = sd))) %>%  # named list of funs
  mutate(across(postBHT:pureVTS,  # calculate 3SD threshold value
                ~3*sd(.),
                .names = "{.col}_thr")) %>% 
  mutate(out_postBHT = case_when(postBHT < postBHT_mean - postBHT_thr | postBHT > postBHT_mean + postBHT_thr ~ T, TRUE ~ F), # mark outliers
         out_betweenTasks = case_when(betweenTasks < betweenTasks_mean - betweenTasks_thr | betweenTasks > betweenTasks_mean + betweenTasks_thr ~ T, TRUE ~ F),
         out_pureVTS = case_when(pureVTS < pureVTS_mean - pureVTS_thr | pureVTS > pureVTS_mean + pureVTS_thr ~ T, TRUE ~ F),
         out_postVTS = case_when(postVTS < postVTS_mean - postVTS_thr | postVTS > postVTS_mean + postVTS_thr ~ T, TRUE ~ F))

durOut3 %>% count(out_postBHT == T, out_pureVTS == T, out_betweenTasks == T, out_postVTS == T) %>% kable()
durOut2 %>% count(extreme_postBHT == T, extreme_pureVTS == T, extreme_betweenTasks == T, extreme_postVTS == T) %>% kable()

durOut3 %>% filter(out_postVTS == T) %>% count()
durOut2 %>% filter(extreme_postVTS == T) %>% count()

durOut3 %>% group_by(across(ends_with('_postBHT'))) %>% get_summary_stats(postBHT)
durOut3 %>% group_by(out_postBHT) %>% get_summary_stats(postBHT)
durOut3 %>% group_by(extreme_postBHT) %>% get_summary_stats(postBHT)

durOut3 %>% group_by(across(ends_with('_betweenTasks'))) %>% get_summary_stats(postBHT)


# 3SD min values
durOut4 <- durOut3 %>% 
  filter(out_postBHT == T) %>% 
  select(postBHT, out_postBHT) %>% 
  get_summary_stats(postBHT) %>% select(min, max) %>% 
  mutate(var = 'postBHT') %>% 
  bind_rows(durOut3 %>% 
              filter(out_pureVTS == T) %>% 
              select(pureVTS, out_pureVTS) %>% 
              get_summary_stats(pureVTS) %>% select(min, max) %>% 
              mutate(var = 'pureVTS')) %>% 
  bind_rows(durOut3 %>% 
              filter(out_betweenTasks == T) %>% 
              select(betweenTasks, out_betweenTasks) %>% 
              get_summary_stats(betweenTasks) %>% select(min, max) %>% 
              mutate(var = 'betweenTasks')) %>% 
  bind_rows(durOut3 %>% 
              filter(out_postVTS == T) %>% 
              select(postVTS, out_postVTS) %>% 
              get_summary_stats(postVTS) %>% select(min, max) %>% 
              mutate(var = 'postVTS')) %>% 
  bind_rows(t3b %>% filter(out_break==T) %>% 
              select(timeBreak, out_break) %>% 
              get_summary_stats(timeBreak) %>% select(min, max) %>% 
              mutate(var = 'timeBreak'))

durOut3 %>% # mark based on betweenTasks and postVTS, N=9
  select(id, betweenTasks, postVTS, out_betweenTasks, out_postVTS) %>% 
  filter(out_betweenTasks == T | out_postVTS == T) %>% 
  arrange(postVTS)

durOut8 <- durOut3 %>% # mark based on all 5 durations: BHT, VTS, between tasks, total, breaks
  left_join(t4 %>% select(id, max_break, out_break)) %>% 
  relocate(max_break, .after = postVTS) %>% 
  select(id, postBHT, pureVTS, betweenTasks, postVTS, max_break, 
         out_postBHT, out_pureVTS, out_betweenTasks, out_postVTS, out_break) %>% 
  filter(out_postBHT == T |  out_pureVTS == T |  
           out_betweenTasks == T | out_postVTS == T | out_break == T) %>% 
  arrange(postVTS) 

durOutIDs5c <- durOut8 %>% 
  filter(if_any(starts_with('out'),
                ~.==T)) %>% 
  select(id, (starts_with('out')))

durOutIDs4c <- durOut8 %>% 
  filter(if_any(c(starts_with('out'), -out_break),
                ~.==T)) %>% 
  select(id, (starts_with('out')), -out_break)

durVals <- durOut4 %>% select(min) %>% t() # get min values of outliers per duration type
durOut4 %>% select(var, min) %>% kable()

durOut9 <- durOut8 %>% # create table with counts of outliers per duration type
  count(out_postBHT == T, out_pureVTS == T, out_betweenTasks == T, out_postVTS == T, out_break == T) %>% 
  arrange(desc(n))
names(durOut9) <- gsub(" == T", "", names(durOut9))

durOut9 %>% # print nicely
  mutate(n = as.character(n)) %>%
  mutate_if(is.logical, as.character) %>%
  add_row(out_postBHT= as.character(durVals[1]),
          out_pureVTS= as.character(durVals[2]),
          out_betweenTasks= as.character(durVals[3]),
          out_postVTS= as.character(durVals[4]),
          out_break = as.character(durVals[5]),
          n = 'min value')%>% kable()

durOut8 %>% # remove based on all 5 durations, N=29 
  select(id, postBHT, pureVTS, betweenTasks, postVTS, max_break,
         out_postBHT, out_pureVTS, out_betweenTasks, out_postVTS, out_break) %>% 
  arrange(postVTS) %>% kable()


# extreme min values
# durOut6 <- durOut3 %>% 
#   filter(extreme_postBHT == T) %>% 
#   select(postBHT, extreme_postBHT) %>% 
#   get_summary_stats(postBHT) %>% select(min, max) %>% 
#   mutate(var = 'postBHT') %>% 
#   bind_rows(durOut3 %>% 
#               filter(extreme_pureVTS == T) %>% 
#               select(pureVTS, extreme_pureVTS) %>% 
#               get_summary_stats(pureVTS) %>% select(min, max) %>% 
#               mutate(var = 'pureVTS')) %>% 
#   bind_rows(durOut3 %>% 
#               filter(extreme_betweenTasks == T) %>% 
#               select(betweenTasks, extreme_betweenTasks) %>% 
#               get_summary_stats(betweenTasks) %>% select(min, max) %>% 
#               mutate(var = 'betweenTasks')) %>% 
#   bind_rows(durOut3 %>% 
#               filter(extreme_postVTS == T) %>% 
#               select(postVTS, extreme_postVTS) %>% 
#               get_summary_stats(postVTS) %>% select(min, max) %>% 
#               mutate(var = 'postVTS'))
# 
# 
# durVals2 <- durOut6 %>% select(min) %>% t()
# durOut7 <- durOut3 %>% count(extreme_postBHT == T, extreme_pureVTS == T, extreme_betweenTasks == T, extreme_postVTS == T) %>% 
#   filter(!if_all(starts_with('extreme'), ~.==F))
# names(durOut7) <- gsub(" == T", "", names(durOut7))
# 
# durOut7 %>% 
#   mutate(n = as.character(n)) %>% 
#   mutate_if(is.logical, as.character) %>% 
#   add_row(extreme_postBHT = as.character(durVals2[1]), 
#           extreme_pureVTS = as.character(durVals2[2]),
#           extreme_betweenTasks = as.character(durVals2[3]),
#           extreme_postVTS = as.character(durVals2[4]),
#           n = 'min value') %>% kable()
# 
# durOut3 %>% # remove based on betweenTasks and postVTS: N=24
#   select(id, betweenTasks, postVTS, extreme_betweenTasks, extreme_postVTS) %>% 
#   filter(extreme_betweenTasks == T | extreme_postVTS == T) %>% 
#   arrange(postVTS) %>% kable()


# Outliers ----------------------------------------------------------------

# instruction check
q3 %>% 
  filter(q == 'controlover' & task == 'bht') %>% # no control group but answered full control
  filter(group == 0 & ans == 6) %>% 
  select(id, group, ans) -> insOutIDs1
4/358*100 # %

q3 %>% 
  filter(q == 'controlover' & task == 'bht') %>% # full control group but answered no control
  filter(group == 1 & ans == 1) %>% 
  select(id, group, ans) -> insOutIDs2
3/358*100 # %

# durations
# - see section above

# rt
p7 %>% 
  group_by(id, block, switch_type) %>%
  identify_outliers(rt) %>% 
  filter(is.extreme == T)

rtOut <- p7 %>% 
  group_by(id, block, task_perf) %>% 
  get_summary_stats(rt) %>% 
  complete(id, block, task_perf, variable) # explicit NAs

rtOut2 <- rtOut %>% 
  select(id:task_perf, mean, sd) %>% 
  rename(rtMean = mean, rtSD = sd) %>% 
  mutate(rtUp = rtMean + 3*rtSD,
         rtLo = rtMean - 3*rtSD)


# p8 <- p7 %>% left_join(rtOut2, by = c('id', 'block', 'task_perf')) %>% 
#   filter(!(rt > rtUp | rt < rtLo))
# 
# 100-(nrow(p8)/nrow(p7)*100) # % removed by 3SD


rtOut3 <- p7 %>% 
  group_by(id, block, task_perf) %>% 
  identify_outliers(rt) 

100-((nrow(p7) - rtOut3 %>% filter(is.extreme == T) %>% count()) / p7 %>% filter(!is.na(rt)) %>% count())*100 # % removed by extreme


# accuracy 

accOutAvg <- p7 %>% # calculate 65% per subject only (not per task/block)
  group_by(id) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% 
              group_by(id) %>% 
              count()) %>% 
  mutate(acc = sumCorr/n) %>% ungroup()

ggplot(accOutAvg, aes(x = acc)) + 
  geom_histogram(color = "grey20", fill = 'grey50', binwidth = 0.01) +
  geom_vline(xintercept = 0.65, linetype = 'dashed', color = 'dark red') +
  labs(y='Count', x='Accuracy')

ggsave('accXid.png', path = 'plots_1/', w=4,h=4)

accOutIDs <- accOutAvg %>% # mark <65% subjects
  filter(acc<.65) %>% select(id,acc) %>% arrange(acc)

p10 <- p7 %>% filter(id %nin% accOutIDs$id)

(1-(nrow(p10))/nrow(p7))*100 # % trials 
4/358


accOut %>% 
  ggplot(., aes(x = as.factor(id), y = acc, colour = task_perf)) +
  geom_point()+
  facet_wrap(~block) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())

accOut %>% 
  ggplot(., aes(x = task_perf, y = acc, colour = as.factor(id), group = as.factor(id))) +
  geom_point(aes(size=n),alpha = .5, position = position_dodge(width = .05))+
  geom_line(alpha = .3, position = position_dodge(width = .05))+
  facet_wrap(~block) +
  scale_size(breaks= c(1,50,100), range=c(1,4)) +
  guides(colour="none")

p7 %>% count(is.na(corr))

# Cleaning pipeline -------------------------------------------------------

c1 <- p7 %>%  # final data cleaning
  arrange(id) %>% 
  filter(id %nin% 1002) %>% # age 61
  filter(id %nin% c(insOutIDs1$id, insOutIDs2$id)) %>% # sense of control reversed answer
  filter(id %nin% durOutIDs4c$id) %>% # prolonged durations (mind if using 4 or 5 criteria)
  filter(id %nin% accOutIDs$id) # accuracy <65% per subject only

c1(1-(nrow(c1)/nrow(p7)))*100 # % removed trials
(1-(c1 %>% distinct(id) %>% nrow()/354))*100 # % removed subjects

# allOutIDs5c <- insOutIDs1 %>% select(id) %>% mutate(outType = 'ins') %>% # 5 criteria of dur out
#   bind_rows(insOutIDs2 %>% select(id) %>% mutate(outType = 'ins')) %>% 
#   bind_rows(durOutIDs5c %>% select(id) %>% mutate(outType = 'dur')) %>% 
#   bind_rows(tibble(id = 1002, outType = 'age')) %>% 
#   bind_rows(accOutIDs %>% select(id) %>% mutate(outType = 'acc')) 
# length(allOutIDs5c$id) 

allOutIDs4c <- insOutIDs1 %>% select(id) %>% mutate(outType = 'ins') %>% # 4 criteria of dur out
  bind_rows(insOutIDs2 %>% select(id) %>% mutate(outType = 'ins')) %>% 
  bind_rows(durOutIDs4c %>% select(id) %>% mutate(outType = 'dur')) %>% 
  bind_rows(tibble(id = 1002, outType = 'age')) %>% 
  bind_rows(accOutIDs %>% select(id) %>% mutate(outType = 'acc')) 
length(allOutIDs4c$id) 

allOutIDs4c %>% group_by(id) %>% count() %>% filter(n>1) # no subject is outType of more than 1 type


# Prepare final dataset ---------------------------------------------------

colnames(c1) 

fct2lvl <- c('group', 'block', 'trial_type', 'task_perf', 'switch')

# set contrasts with faux 
c2 <- c1 %>% 
  left_join(nasa3q) %>% 
  # mutate(across(fct2lvl,
  #               ~contr_code_anova(.))) %>%  # change all 2-level factors to anova coding using faux
  filter(id %nin% vsrOutIDs$id) # remove 0 VSR cases

contrasts(c2$block)

# set position dodge params
pd = position_dodge(0.2)
pdtxt = position_dodge(1)
facet_group <- c('0' = 'No control', '1' = 'Full control')


# VSR ---------------------------------------------------------------------


# check VSR outliers

s1a <- c1 %>% 
  group_by(id) %>% 
  filter(!is.na(switch)) %>%
  filter(trial_type == 'free') %>% # discard forced trials
  summarise(switch_sum = sum(switch)) %>% # get nr of all switches per id and block 
  left_join(c1 %>% # join/merge with another df in which we get nr of all observations
              group_by(id) %>% 
              filter(trial_type == 'free') %>% # discard forced trials
              filter(!is.na(switch)) %>%
              count()) %>% 
  mutate(vsr = round(switch_sum/n, 3)) # get VSR (nr of switches / nr of observations) 

vsrOutIDs <- s1a %>% filter(switch_sum == 0) %>% # subjects to be removed due to 0 vsr/switches
  left_join(c1 %>% group_by(id,group) %>% slice(1) %>% select(id,group)) 
11/358

vsrOutIDs %>% freq_table(group)

ggplot(s1a, aes(x=switch_sum)) + 
  geom_histogram(binwidth = 1) +
  scale_x_continuous(breaks=c(seq(0,10,5), seq(20,150,20))) +
  labs(x = 'Switch sum')
ggplot(s1a, aes(x=vsr)) + 
  geom_histogram(binwidth = 0.005) +
  scale_x_continuous(breaks=c(seq(0,.05,0.05), seq(.1,.7,.1))) +
  labs(x = 'VSR')

ggsave('switchSum_hist.png', path = 'plots_1/', w=5,h=4)


# vsrOutIDs2 <- s1a %>% filter(switch_sum <= 2) %>%
#   left_join(c1 %>% group_by(id,group) %>% slice(1) %>% select(id,group))

s1a %>% # cleaning by 3sd doesn't help here
  cross_join(s1a %>% summarise(mean = mean(vsr), thr = 3*(sd(vsr)), up = mean+thr, lo = mean-thr)) %>% 
  mutate(out_vsr = if_else(vsr < lo | vsr > up, T, F)) %>% 
  filter(out_vsr == T)

# consider cuefol == 0

s1b <- c1 %>% 
  group_by(id) %>% 
  filter(trial_type == 'forced', cuefol == 0) %>% 
  count() %>% rename(cuefolNo = n) %>% 
  left_join(c1 %>% group_by(id) %>% 
              filter(trial_type == 'forced') %>% 
              count()) %>% 
  mutate(prop_cuefolNo = cuefolNo/n) %>% 
  arrange(-prop_cuefolNo) %>% 
  left_join(c1 %>% group_by(id,group) %>% slice(1) %>% select(id,group))


s1c <- s1b %>% 
  full_join(c1 %>% group_by(id) %>% slice(1) %>% select(id,group)) %>% 
  mutate(across(everything(),
                ~ if_else(is.na(.), 0, .))) %>% 
  arrange(id) %>% ungroup()
  

ggplot(s1c, aes(x = prop_cuefolNo)) +
  geom_histogram() +
  scale_x_continuous(breaks = pretty_breaks()) +
  geom_vline(xintercept = 0.503, linetype = 'dashed')

s1c %>% summarise(mean = mean(prop_cuefolNo), thr = 3*(sd(prop_cuefolNo, na.rm=T)), up = mean+thr, lo = mean-thr)

c1 %>% filter(id %nin% (s1c %>% filter(cuefolNo==0) %>% select(id)))

c1a <- c1 %>% 
  filter(block == 'mixed') %>% 
  group_by(id) %>% 
  mutate(cueToSwitch = if_else(task_cued == lag(task_perf), F, T)) %>% 
  filter(id==46)
  
# Overall VSR per id 

s1 <- c2 %>% 
  group_by(id, block, group) %>% 
  filter(!is.na(switch)) %>%
  filter(trial_type == 'free') %>% # discard forced trials
  summarise(switch_sum = sum(switch)) %>% # get nr of all switches per id and block 
  left_join(c1 %>% # join/merge with another df in which we get nr of all observations
              group_by(id, block, group) %>% 
              filter(trial_type == 'free') %>% # discard forced trials
              filter(!is.na(switch)) %>%
              count()) %>% 
  mutate(group = as.factor(group)) %>%
  mutate(vsr = round(switch_sum/n, 3)) %>%  # get VSR (nr of switches / nr of observations) 
  filter(id %nin% vsrOutIDs$id)
  
s1 %>% # plot
  group_by(block) %>%
  get_summary_stats(vsr) %>%
  ggplot(., aes(y = mean, x = block)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1) +
  geom_text(aes(label = round(mean,2)), nudge_x = 0.2, show.legend = FALSE) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(y = 'VSR (%)', x = 'Block')

ggsave('vsr_block.png', path = 'plots_1/', w=4,h=4)

s1 %>% # plot
  group_by(block, group) %>%
  get_summary_stats(vsr) %>%
  ggplot(., aes(y = mean, x = block, colour = group)) +
  geom_point(size=4, position = pd) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = pd) +
  geom_text(aes(label = round(mean,2)), position = pdtxt, show.legend = FALSE) +
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  labs(y = 'VSR (%)', x = 'Block', colour = 'Group')

ggsave('vsr_blockXgroup.png', path = 'plots_1/', w=5.5,h=4)


# VSR based on switch difficulty, per id (VSR decomposed into easy and difficult)
s2 <- s1 %>% 
  left_join(c2 %>% 
              group_by(id, block, group) %>% 
              filter(trial_type == 'free') %>% # discard forced trials
              filter(switch_type == 'switch_shape') %>% 
              count(name = 'switch_diff') %>% # count difficult switches (location -> shape)
              left_join(c2 %>% 
                          group_by(id, block, group) %>%  
                          filter(trial_type == 'free') %>%
                          filter(switch_type == 'switch_loc') %>% 
                          count(name = 'switch_easy')) %>% 
              mutate(group = as.factor(group))) %>%  # count easy switches (shape -> location)
  relocate(c('switch_easy', 'switch_diff'), .before = 'n') %>% # rearrange for convenience
  mutate(vsr_easy = switch_easy/n, # get VSR for easy and diff switches
         vsr_diff = switch_diff/n) %>% 
  replace(is.na(.), 0)

s2 %>%  # plot
  group_by(block, group) %>% 
  get_summary_stats(vsr_easy, vsr_diff) %>% 
  ggplot(., aes(y = mean, x = block, fill = variable, colour = variable)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = pd) +
  geom_text(aes(label = round(mean,3)), position = pdtxt, show.legend = FALSE) +
  facet_wrap(~group, labeller = as_labeller(facet_group))+
  labs(y = 'VSR', colour = 'VSR component', x = 'Block') + 
  scale_colour_discrete(labels = c('Easy', 'Difficult')) +
  guides(fill="none")

ggsave('vsr_blockXdiffXgroup.png', path = 'plots_1/', w=8,h=4)

s2 %>%  # plot
  group_by(block, group) %>% 
  get_summary_stats(vsr_easy, vsr_diff) %>% 
  ggplot(., aes(y = mean, x = block, fill = group, colour = group)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = pd) +
  geom_text(aes(label = round(mean,3)), position = pdtxt, show.legend = FALSE) +
  facet_wrap(~variable, labeller = as_labeller(c('vsr_easy' = 'VSR easy', 'vsr_diff' = 'VSR difficult')))+
  labs(y = 'VSR', colour = 'Group', x = 'Block') + 
  scale_colour_discrete(labels = c('no control', 'full control')) +
  guides(fill="none")

ggsave('vsr_blockXdiffXgroup2.png', path = 'plots_1/', w=8,h=4)


# check if VSR can be dissociated in freeOnly and mixed block
s2 %>% ungroup() %>%
  group_by(block) %>% 
  mutate(switch_difference = switch_easy - switch_diff, .before = 'n') %>% 
  get_summary_stats(switch_difference)

s3 <- s2 %>% 
  left_join(s2 %>% # add VSR per id only (averaged across blocks)
              ungroup() %>% 
              group_by(id) %>% 
              summarise(vsr_avg = mean(vsr),
                        vsr_easy_avg = mean(vsr_easy),
                        vsr_diff_avg = mean(vsr_diff))) %>% 
  mutate(block = as.factor(block))

# proportion of VSR easy to VSR difficult (another way to plot vsr_blockXdiff)
s2 %>% # averaged across id
  ungroup() %>% 
  mutate(vsr_easy_prop = vsr_easy/(vsr_diff + vsr_easy),
         vsr_diff_prop = vsr_diff/(vsr_diff + vsr_easy)) %>% 
  pivot_longer(cols = c('vsr_easy_prop', 'vsr_diff_prop'), values_to = 'vsr_prop', names_to = 'vsr_type') %>% 
  group_by(block, vsr_type) %>% 
  get_summary_stats(vsr_prop) %>% 
  ggplot(., aes(y = mean, x = block, colour = vsr_type)) +
  geom_point(size = 3, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.3)) +
  geom_text(aes(label = round(mean,3)), position = position_dodge(1), show.legend = FALSE) +
  labs(y = 'Proportion', colour = 'VSR component prop', x = 'Block') + 
  scale_colour_discrete(labels = c('Difficult', 'Easy'))

s2 %>% # per id
  ungroup() %>% 
  mutate(vsr_easy_prop = vsr_easy/(vsr_diff + vsr_easy),
         vsr_diff_prop = vsr_diff/(vsr_diff + vsr_easy)) %>% 
  pivot_longer(cols = c('vsr_easy_prop', 'vsr_diff_prop'), values_to = 'vsr_prop', names_to = 'vsr_type') %>% 
  group_by(id, block, vsr_type) %>% 
  get_summary_stats(vsr_prop) %>% 
  ggplot(., aes(y = mean, x = block, colour = vsr_type)) +
  geom_point(size = 3, position = position_dodge(0.3)) +
  facet_wrap(~id) +
  geom_text(aes(label = round(mean,3)), position = position_dodge(1), show.legend = FALSE) +
  labs(y = 'Proportion', colour = 'VSR component prop', x = 'Block') + 
  scale_colour_discrete(labels = c('Difficult', 'Easy'))


# VSR by id and diff in mixed block only
s4 <- s3 %>% # reshape and summarise for plot 2
  filter(block == 'mixed') %>% 
  pivot_longer(cols = c('vsr_easy', 'vsr_diff'), names_prefix = 'vsr_', names_to = 'difficulty', # reshape (to long format)
               values_to = 'vsr_difficulty') %>% 
  mutate(difficulty = as.factor(difficulty)) %>% 
  relocate(difficulty, .after = block)

ggplot(s4, aes(y = vsr_difficulty, x = difficulty, colour = difficulty)) + # plot VSR in mixed only
  geom_point(size = 3, position = position_dodge(0.1)) +
  facet_wrap(~id, labeller = labeller(id = label_both)) +
  geom_text(aes(label = round(vsr_difficulty,2)), show.legend = F, nudge_x = 0.3, size=4) +
  scale_x_discrete(labels = c('Easy', 'Difficult')) +
  labs(y = 'VSR (%)', x = 'Difficulty of switch', caption = 'Mixed block only') + 
  guides(colour="none")

ggsave('vsr_diffXid.png', path = 'plots_1/', w=12.5,h=10)


# Counterbalance and VSR

s2 %>% 
  left_join(c2 %>% group_by(id) %>% filter(row_number()==1) %>% select(id, cBal), 
            by = 'id') %>% 
  group_by(cBal, block) %>%
  get_summary_stats(vsr) %>%
  ggplot(., aes(y = mean, x = block)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1) +
  geom_text(aes(label = round(mean,2)), nudge_x = 0.25, show.legend = FALSE) +
  facet_wrap(~cBal) +
  labs(y = 'VSR (%)', x = 'Block')

ggsave('vsr_blockXcb.png', path = 'plots_1/', w=6,h=4)

s2 %>% 
  left_join(c2 %>% group_by(id) %>% filter(row_number()==1) %>% select(id, group, cBal), 
            by = c('id', 'group')) %>% 
  group_by(group, cBal, block) %>%
  get_summary_stats(vsr) %>%
  ggplot(., aes(y = mean, x = block, color = group)) +
  geom_point(size=3, position=pd) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=pd) +
  geom_text(aes(label = round(mean,2)), position=pdtxt, show.legend = FALSE) +
  facet_wrap(~cBal) +
  scale_colour_discrete(labels=c('no control', 'full control')) +
  labs(y = 'VSR (%)', x = 'Block')

ggsave('vsr_blockXcbXgroup.png', path = 'plots_1/', w=8,h=4)

# Task selection ##########

c2 %>% # overall proportion of diff/easy task selection
  filter(task_perf != 'unknownTask') %>% 
  filter(trial_type == 'free') %>% # remove instructed trials
  group_by(group, task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate(easy_prop = loc/(loc+shape),
         diff_prop = shape/(loc+shape))

c2 %>% # plot by block only
  filter(task_perf != 'unknownTask') %>% 
  filter(trial_type == 'free') %>% # remove instructed trials\
  mutate(group=as.factor(group)) %>% 
  group_by(id, group, block, task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  mutate(diff_sel = shape/(loc+shape)) %>% # measure of difficult selection
  group_by(group, block) %>% # average across id
  get_summary_stats(diff_sel) %>% 
  ggplot(., aes(x = block, y = mean, colour = group)) +
  geom_point(size = 3, position=pd) +
  geom_text(aes(label = round(mean,2)), position=pdtxt, show.legend = F) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=pd) +
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  labs(y = 'Selection of difficult task (%)', x = 'Block', colour = 'Group')

ggsave('taskSelXblockXgroup.png', path = 'plots_1/', w=5.5,h=4)


c2 %>% # plot averaged, a bit redundant compared to previous plot
  filter(task_perf != 'unknownTask') %>% 
  filter(trial_type == 'free') %>% # remove instructed trials
  group_by(id, block, task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  mutate(easy_sel = loc/(loc+shape),
         diff_sel = shape/(loc+shape)) %>% 
  pivot_longer(cols = c('easy_sel', 'diff_sel'), names_to = 'task_perf', values_to = 'sel') %>% 
  group_by(block, task_perf) %>% 
  get_summary_stats(sel) %>% 
  ggplot(., aes(x = task_perf, y = mean, colour = block)) +
  geom_point(size = 3, position=position_dodge(0.5)) +
  geom_text(aes(label = round(mean,3)), position=position_dodge(1)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=position_dodge(0.5)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(labels=c('Difficult', 'Easy')) +
  labs(y = 'Proportion (%)', x = 'Task')

ggsave('taskSelXblock2.png', path = 'plots_1/', width = 6, height = 4) # a bit redundant 

c2 %>% # counterbalance effect on overall proportion of diff/easy task selection
  filter(task_perf != 'unknownTask') %>% 
  filter(trial_type == 'free') %>% # remove instructed trials
  group_by(cBal, task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate(easy_prop = loc/(loc+shape),
         diff_prop = shape/(loc+shape))

c2 %>% # counterbalance plot by block only
  filter(task_perf != 'unknownTask') %>% 
  filter(trial_type == 'free') %>% # remove instructed trials
  group_by(id, cBal, block, task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  mutate(diff_sel = shape/(loc+shape)) %>% 
  group_by(cBal, block) %>% # average across id
  get_summary_stats(diff_sel) %>% 
  ggplot(., aes(x = block, y = mean)) +
  geom_point(size = 3, position=position_dodge(0.5)) +
  facet_wrap(~cBal) +
  geom_text(aes(label = round(mean,2)), nudge_x = 0.2) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=position_dodge(0.5)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = 'Selection of difficult task (%)', x = 'Block')

ggsave('taskSelXblockXcb.png', path = 'plots_1/', width = 6, height = 4) 

c2 %>% # counterbalance plot by block only
  filter(task_perf != 'unknownTask') %>% 
  filter(trial_type == 'free') %>% # remove instructed trials
  group_by(id, group, cBal, block, task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  mutate(diff_sel = shape/(loc+shape)) %>% 
  group_by(group, cBal, block) %>% # average across id
  get_summary_stats(diff_sel) %>% 
  ggplot(., aes(x = block, y = mean, color = group)) +
  geom_point(size = 3, position=pd) +
  facet_wrap(~cBal) +
  geom_text(aes(label = round(mean,2)), position=pdtxt, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=pd) +
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_discrete(labels=c('no control', 'full control')) +
  labs(y = 'Selection of difficult task (%)', x = 'Block')

ggsave('taskSelXblockXcbXgroup.png', path = 'plots_1/', width = 8, height = 4) 


# Proportion of switches --------------------------------------------------

y1 <- c2 %>% # overall proportion of switch
  filter(!is.na(switch)) %>% 
  filter(trial_type == 'free') %>%
  group_by(id, group, block, switch) %>%
  count() %>% 
  pivot_wider(names_from = 'switch', values_from = 'n') %>% 
  rename(stay = '0', switch = '1') %>% 
  replace(is.na(.), 0) %>%
  mutate(switch_prop = switch/(switch+stay),
         stay_prop   = stay  /(switch+stay),
         group = as.factor(group)) %>% 
  group_by(group, block) %>% 
  get_summary_stats(switch_prop, stay_prop)

y1 %>% # plot
  rename(switch = variable) %>% 
  mutate(switch = as.factor(switch)) %>%
  ggplot(., aes(y = mean, x = block, colour = switch))+
  geom_point(size = 2, position = position_dodge(.3)) +
  geom_text(aes(label = round(mean,2)), show.legend = F, position = position_dodge(1)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=position_dodge(.3)) +
  facet_wrap(~group, labeller = as_labeller(facet_group))+
  scale_colour_discrete(labels=c('Switch', 'Stay')) +
  labs(y = 'Proportion of stay/switch', x = 'Block', colour = 'Type')

y1 %>% # plot
  rename(switch = variable) %>% 
  mutate(switch = as.factor(switch)) %>%
  ggplot(., aes(y = mean, x = block, colour = group))+
  geom_point(size = 2, position = pd) +
  geom_text(aes(label = round(mean,2)), show.legend = F, position = pdtxt) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=pd) +
  facet_wrap(~switch, labeller = as_labeller(c('switch_prop' = 'Switch', 'stay_prop' = 'Stay')))+
  scale_colour_discrete(labels=c('no control', 'full control')) +
  labs(y = 'Proportion of stay/switch', x = 'Block', colour = 'Group')

ggsave('switchProp_blockXgroup.png', path = 'plots_1/', width = 8, height = 4)

# y2 <- p7 %>% # proportion of switch_type
#   filter(switch_type != 'unknownSwitch') %>% 
#   filter(trial_type == 'free') %>%
#   group_by(id, block, switch_type) %>%
#   count() %>% 
#   pivot_wider(names_from = 'switch_type', values_from = 'n') %>% 
#   replace(is.na(.), 0) %>% 
#   mutate(switch_shape_p = switch_shape/(stay_loc+stay_shape+switch_loc+switch_shape),
#          switch_loc_p = switch_loc/(stay_loc+stay_shape+switch_loc+switch_shape),
#          stay_shape_p = stay_shape/(stay_loc+stay_shape+switch_loc+switch_shape),
#          stay_loc_p = stay_loc/(stay_loc+stay_shape+switch_loc+switch_shape)) %>% 
#   group_by(block) %>% 
#   get_summary_stats(c('switch_shape_p':'stay_loc_p'))
# 
# y2 %>% # plot
#   mutate(switch_type = as.factor(variable), .after = 'block', .keep = 'unused') %>%
#   ggplot(., aes(y = mean, x = block, colour = switch_type))+
#   geom_point(size = 2, position = position_dodge(.3)) +
#   # geom_text(aes(label = round(mean,3)), show.legend = F, nudge_x = 0.2, nudge_y = 0.01) +
#   geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=position_dodge(.3)) +
#   scale_colour_discrete(labels=c('Switch to shape', 'Switch to loc', 'Stay with shape', 'Stay with loc')) +
#   labs(y = 'Proportion of stay/switch', x = 'Block', colour = 'Type')
# 
# ggsave('switchProp_blockXswitchType.png', path = 'plots_1/', width = 6, height = 4)


# Repetition bias #########

z2 <- p7 %>%  # switches combined
  filter(trial_type == 'free') %>% 
  freq_table(id, block, switch_type) %>% select(-n) %>% 
  pivot_wider(names_from = switch_type, values_from = prop) %>% 
  mutate(switch = switch_loc + switch_shape, .keep = 'unused') %>% 
  pivot_longer(names_to = 'type', cols = c('stay_loc':'switch'), values_to = 'prop')

z2 %>% 
  group_by(block, type) %>% 
  get_summary_stats(prop) %>% 
  ggplot(., aes(y=mean, x=block, colour=type)) +
  geom_point(size=3, position=position_dodge(0.5)) +
  scale_colour_discrete(labels = c('Easy repetition', 'Difficult repetition', 'Switches combined')) +
  labs(y='Proporion', x='Block', colour = 'Type') +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=position_dodge(0.5)) +
  # geom_text(aes(label = round(mean,2)), position=position_dodge(1)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6))

ggsave('repBiasXblock.png', path = 'plots_1/', width = 6, height = 4)

  
z3 <- p7 %>% # remove switches
  filter(trial_type == 'free') %>% 
  filter(switch == 0) %>% # reps only
  freq_table(id, block, switch_type) %>% select(-n) %>% 
  pivot_wider(names_from = switch_type, values_from = prop) %>% 
  pivot_longer(names_to = 'type', cols = c('stay_loc','stay_shape'), values_to = 'prop')

z3 %>% 
  group_by(block, type) %>% 
  get_summary_stats(prop) %>% 
  ggplot(., aes(y=mean, x=block, colour=type)) +
  geom_point(size=3, position=position_dodge(0.5)) +
  scale_colour_discrete(labels = c('Easy repetition', 'Difficult repetition')) +
  labs(y='Proportion', x='Block', colour = 'Type') +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=position_dodge(0.5)) +
  geom_text(aes(label = round(mean,1)), position=position_dodge(1.5), show.legend = F) 

ggsave('repBiasXblock2.png', path = 'plots_1/', width = 5.5, height = 4)


# Rigidness ########

r0 <- c2 %>% 
  filter(!is.na(switch)) %>% 
  group_by(id, block, blockNr) %>% 
  mutate(sumSwitch = cumsum(as.double(switch))) %>%  # cumulative sum of switches
  filter(trial_type == 'free')

r0a <- c2 %>% 
  filter(!is.na(switch)) %>% 
  left_join(c2 %>% 
              filter(!is.na(switch)) %>% 
              group_by(id, block, blockNr) %>% 
              mutate(sumSwitch = cumsum(as.double(switch))) %>%  # cumulative sum of switches
              filter(trial_type == 'free')) %>% 
  left_join(c2 %>% 
              filter(!is.na(switch)) %>% 
              filter(trial_type == 'free') %>% 
              group_by(id, block, blockNr) %>% 
              mutate(sumSwitch2 = cumsum(as.double(switch)))) %>% 
  relocate(trial_type, .after=sumSwitch2)


r0 %>% ungroup %>% get_summary_stats(sumSwitch)

r1 <- r0 %>% 
  group_by(id, group, block, blockNr, sumSwitch, task_perf) %>% # count trials per each switch (and other grouping factors)
  count() %>% 
  ungroup() %>% mutate(group = contr_code_anova(group))
 
r1 %>% # single task run length per block
  group_by(group, block, task_perf) %>% 
  get_summary_stats(n) %>% 
  select(group, block, task_perf, n, mean, se) %>% 
  ggplot(., aes(y = mean, x = block, colour = group)) +
  geom_point(size=3, position=pd) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.3, position=pd) +
  geom_text(aes(label = round(mean,1)), show.legend = FALSE, position=pdtxt) +
  facet_grid(block~task_perf)+
  scale_colour_discrete(labels = c('no control', 'full control')) +
  labs(y = 'Average single task run (nr of trials)', x = 'Block', colour = 'Group') 

ggsave('runLengthxDiff.png', path = 'plots_1/', width = 6, height = 5)

r1 %>% # single task run length in freeOnly block
  filter(block == 'freeOnly') %>%
  group_by(group, block, task_perf) %>% 
  get_summary_stats(n) %>% 
  select(group, block, task_perf, n, mean, se) %>% 
  ggplot(., aes(y = mean, x = task_perf, colour = group)) +
  geom_point(size=3, position=pd) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.3, position=pd) +
  geom_text(aes(label = round(mean,1)), show.legend = FALSE, position=pdtxt) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  labs(y = 'Average single task run (nr of trials)', x = 'Task', colour = 'Group') 

r1 %>% # counterbalance single task run length per block
  left_join(c2 %>% group_by(id) %>% filter(row_number()==1) %>% select(id, cBal), 
            by = 'id') %>% 
  ungroup() %>% group_by(cBal, block, task_perf) %>% 
  get_summary_stats(n) %>% 
  # select(block, task_perf, n, mean, se) %>% 
  ggplot(., aes(y = mean, x = block, colour = task_perf)) +
  geom_point(size=3, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=position_dodge(0.5)) +
  geom_text(aes(label = round(mean,1)), show.legend = FALSE, position = position_dodge(1.3)) +
  facet_wrap(~cBal)+
  scale_colour_discrete(labels = c('Easy', 'Difficult')) +
  labs(y = 'Average single task run (nr of trials)', x = 'Block', colour = 'Task') 

ggsave('runLengthXdiffXcb.png', path = 'plots_1/', width = 6.5, height = 4)


# Accuracy ##########

p7 %>% # stay vs switch
  group_by(id, block, switch) %>% 
  filter(!is.na(switch)) %>%
  mutate(corr = as.integer(corr)) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% 
              group_by(id, block, switch) %>% 
              filter(!is.na(switch)) %>% 
              count()) %>% 
  mutate(acc = sumCorr/n) %>% 
  group_by(block, switch) %>% 
  get_summary_stats(acc)

p7 %>% # easy vs difficult task
  group_by(id, block, task_perf) %>% 
  filter(task_perf != 'unknownTask') %>%
  mutate(corr = as.integer(corr)) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% 
              group_by(id, block, task_perf) %>% 
              filter(task_perf != 'unknownTask') %>% 
              count()) %>% 
  mutate(acc = sumCorr/n) %>% 
  group_by(block, task_perf) %>% 
  get_summary_stats(acc)

p7 %>% # plot switch type, averaged across ids
  group_by(id, block, switch_type) %>% 
  filter(switch_type != 'unknownSwitch') %>%
  mutate(corr = as.integer(corr)) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% 
              filter(switch_type != 'unknownSwitch') %>% 
              group_by(id, block, switch_type) %>% 
              count()) %>% 
  mutate(acc = sumCorr/n) %>% 
  filter(!(id == 107 & switch_type == 'stay_loc')) %>% # remove 1 weird case
  group_by(block, switch_type) %>% 
  get_summary_stats(acc) %>% 
  ggplot(., aes(y = mean, x = factor(switch_type, level=c('stay_loc','switch_loc','stay_shape','switch_shape')), colour = block)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.2)) +
  geom_text(aes(label = round(mean,3)), show.legend = FALSE, position = position_dodge(1.5)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6)) +
  # scale_y_continuous(limits = c(NA, 1)) +
  labs(y = 'Accuracy (%)', x = 'Trial type')

ggsave('accuracy.png', path = 'plots_1/', width = 6, height = 5)

p7 %>% # plot switch type, per id
  group_by(id, block, switch_type) %>% 
  filter(switch_type != 'unknownSwitch') %>%
  mutate(corr = as.integer(corr)) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% 
              filter(switch_type != 'unknownSwitch') %>% 
              group_by(id, block, switch_type) %>% 
              count()) %>% 
  mutate(acc = sumCorr/n) %>% 
  group_by(id, block, switch_type) %>% 
  ggplot(., aes(y = acc, x = switch_type, colour = block)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  labs(y = 'Accuracy (%)', x = 'Trial type', colour = 'Block type') +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6)) +
  facet_wrap(~id)

ggsave('accuracyXtrialTypeXid.png', path = 'plots_1/', width = 10, height = 8)


p7 %>% # plot switch type, per id
  group_by(id, blockNrC, block) %>% 
  mutate(corr = as.integer(corr)) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% 
              group_by(id, blockNrC, block) %>% 
              count()) %>% 
  mutate(acc = sumCorr/n) %>% 
  group_by(id, blockNrC) %>% 
  ggplot(., aes(y = acc, x = blockNrC, colour = block)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  facet_wrap(~id) +
  labs(y = 'Accuracy (%)', x = 'Block number', colour = 'Block type') 

ggsave('accuracyXblockNrXid.png', path = 'plots_1/', width = 10, height = 8)


# RTs ############

p7 %>% 
  group_by(block, switch_type) %>% 
  filter(switch_type != 'unknownSwitch') %>% 
  get_summary_stats(rt) %>% 
  ggplot(., aes(y = mean, x = factor(switch_type, level=c('stay_loc','switch_loc','stay_shape','switch_shape')), colour = block)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.2)) +
  geom_text(aes(label = round(mean,3)), show.legend = FALSE, position = position_dodge(1.5)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6)) +
  labs(y = 'RT (s)', x = 'Trial type') 

ggsave('rt.png', path = 'plots_1/', width = 6, height = 5)


ggplot(p7, aes(x = rt)) +
  geom_histogram(color = 'grey20') +
  labs(y = 'Count', x = 'RT (s)')

ggplot(p7, aes(x = rt, colour = trial_type)) +
  geom_histogram(aes(y = ..density..), color = "grey20", fill = 'grey50', binwidth = 0.1) +
  geom_density(color = 'darkgreen', fill = "lightgreen", alpha = 0.3) +
  geom_vline(aes(xintercept=median(rt)), linetype = 'dashed') +
  annotate('text', label = paste('Med = ', p7 %>% get_summary_stats(rt) %>% select(median), '(s)\n',
                                 'Mean = ', p7 %>% get_summary_stats(rt) %>% select(mean), '(s)'),
           x=1.3, y=1.5, size = 4) +
  labs(y = 'Density', x = 'RT (s)')

ggsave('rt_hist.png', path = 'plots_1/', width = 6, height = 5)

p7 %>% # plot switch type, per id
  group_by(id, blockNrC, block) %>% 
  get_summary_stats(rt)  %>% 
  ggplot(., aes(y = mean, x = blockNrC, colour = block)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1) +
  facet_wrap(~id) +
  labs(y = 'RT (s)', x = 'Block number', colour = 'Block type') 

ggsave('rtXblockNrXId.png', path = 'plots_1/', width = 10, height = 8)


# Checks ###########

# check for NAs in all columns
p5 %>% 
  ungroup() %>% 
  select(everything()) %>%  # replace to your needs
  summarise_all(funs(sum(is.na(.))))

# check if VSR for easy/diff equals to overall VSR 
s2 %>% 
  mutate(checkVsr = case_when(round(vsr, 3) == round(vsr_easy + vsr_diff, 3) ~ 1, # should be 1 (unless NA)
                              TRUE ~ 0)) %>% 
  print(n=Inf)


# Prepare to modelling -------------------------------------------------------- 


# save models
save(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15, file = 'models1.RData')

# set contrasts with faux 

fct2lvl <- c('group', 'cBal', 'block', 'trial_type', 'task_perf') # do not include 'switch' as it is measre (~VSR)

c3 <- c2 %>% 
  left_join(s2 %>% select(id, block, contains('vsr'))) %>% # add VSR values, calculated per id and block (free/mixed)
  mutate(across(fct2lvl, 
         ~ as.factor(.))) %>% 
  # mutate(across(fct2lvl,
  #               ~contr_code_anova(.))) %>% # change all 2-level factors to anova coding using faux
  # select(!contains('nasaSum')) %>% 
  mutate(
         switch = as.integer(switch),
         nasa_bhtD = if_else(nasa_bht < 3.75, 0, 1),  # nasa 0 below median, 1 equal or above median
         switch = if_else(trial_type == 'forced', NA, switch), # remove switch info from instructed trials
         logrt = log(rt),
         task_dv = if_else(trial_type == 'forced', NA, task_perf)) %>% 
  mutate(switch_iv = contr_code_anova(switch)) %>% 
  mutate(nasas = standardise(nasa_bht),
         nasasV = standardise(nasa_vts))

contrasts(c3$group)
contrasts(c3$cBal)
contrasts(c3$block)
contrasts(c3$task_perf)
contrasts(c3$switch_iv)

contrasts(c3$group) <- c(-0.5,0.5); contrasts(c3$group)
contrasts(c3$cBal) <- c(-0.5,0.5); contrasts(c3$cBal)
contrasts(c3$block) <- c(-0.5,0.5); contrasts(c3$block)
contrasts(c3$trial_type) <- c(-0.5,0.5); contrasts(c3$trial_type)
contrasts(c3$task_perf)  <- c(-0.5,0.5); contrasts(c3$task_perf)


c3 %>% 
  # group_by(id) %>% 
  get_summary_stats(nasa_bht)

c3 %>% filter(nasa_bht < 3.75) %>% group_by(id) %>% distinct(id) %>% nrow()
c3 %>% filter(nasa_bht >= 3.75) %>% group_by(id) %>% distinct(id) %>% nrow()


# Manipulation checks ------------------------------------------------------

colnames(q2)

mc1 <- c3 %>% 
  left_join(q2 %>% 
              filter(task=='bht') %>% 
              select(id, group, controlover, extentofcontrol, effort1, effort2, motiv) %>% 
              mutate(group = as.factor(group))) %>% 
  group_by(id) %>% slice(1) %>% 
  select(id, group, controlover, extentofcontrol, effort1, effort2, motiv) %>% 
  ungroup 

mc1 %>% t_test(controlover ~ group) %>% 
  left_join(mc1 %>% cohens_d(controlover ~ group))

ggplot(mc1, aes(x = controlover, fill = group)) + geom_histogram()

mc1 %>% t_test(extentofcontrol ~ group) %>% 
  left_join(mc1 %>% cohens_d(extentofcontrol ~ group))

mc1 %>% t_test(effort1 ~ group) %>% 
  left_join(mc1 %>% cohens_d(effort1 ~ group))  

mc1 %>% t_test(effort2 ~ group) %>% 
  left_join(mc1 %>% cohens_d(effort2 ~ group)) 

mc1 %>% t_test(motiv ~ group) %>% 
  left_join(mc1 %>% cohens_d(motiv ~ group)) 


# Modelling VSR ---------------------------------------------------------------


# Models preregistered

m0 <- glmer(switch ~ group * nasas + (1 | id),
            c3, family = binomial)
summary(m0)
tab_model(m0)
tab_model(m0, show.stat = T, show.icc = F)
logoddsratio_to_r(fixef(m0)) %>% as_tibble(rownames='term') %>% rename('d'='value')
myprint(m0)
tab_model(m0, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          file='m0.html',
          pred.labels=c('(Intercept)','Group','NASA','Group x NASA'),
          dv.labels='')

emmip(m0, group ~ nasas, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

# Term            OddsRatio        d      SE         z       p  Sig 
# (Intercept)         0.182   -0.939   0.053   -32.001   0.000  *** 
# group1              1.145    0.075   0.106     1.281   0.200      
# nasas               0.803   -0.121   0.053    -4.094   0.000  *** 
# group1:nasas        0.788   -0.132   0.107    -2.229   0.026  *  

m0a <- glmer(switch ~ group * nasas + (1 | id),
             c3, family = binomial)
summary(m0a)

m0b <- glmer(switch ~ group * nasas + (group | id),
            c3, family = binomial)
summary(m0b)
anova(m0,m0b)

m0c <- glmer(switch ~ group * nasas + (nasas | id),
             c3, family = binomial)
summary(m0c)

anova(m0,m0b,m0c)

m0d <- glmer(switch ~ group * nasas + (block | id),
             c3, family = binomial)
summary(m0d)

# 

m1 <- glmer(switch ~ group * nasa_bht + (group * nasa_bht | id),
           c3, family = binomial) # didn't converge
summary(m1)

m2 <- glmer(switch ~ group * nasas + (group + nasa_bht | id),
            c3, family = binomial,
            glmerControl(optimizer='bobyqa'))
summary(m2)

anova(m0,m2) # m0

m3 <- glmer(switch ~ group * nasas + (group | id),
            c3, family = binomial)
summary(m3)
myprint(m3)
# Term            OddsRatio        d      SE         z       p  Sig 
# (Intercept)         0.182   -0.939   0.054   -31.667   0.000  *** 
# group1              1.153    0.079   0.108     1.324   0.185      
# nasas               0.804   -0.120   0.054    -4.027   0.000  *** 
# group1:nasas        0.788   -0.131   0.108    -2.201   0.028  *  

anova(m0,m3) # m0 > m3, slightly

m4 <- glmer(switch ~ group * nasa_bht + (nasa_bht | id),
            c3, family = binomial)
summary(m4)

anova(m0,m3,m4)

m5 <- glmer(switch ~ group * nasa_bht + (group + nasa_bht || id),
            c3, family = binomial)  # singular
summary(m5)

# # try different optimizers
# ## show available methods
# allFit(show.meth.tab=TRUE) 
# gm_all <- allFit(m2)
# ss <- summary(gm_all)
# ss$msgs                ## convergence messages
# ss$which.OK            ## logical vector: which optimizers worked?
# ## the other components only contain values for the optimizers that worked
# ss$llik                ## vector of log-likelihoods
# ss$fixef               ## table of fixed effects
# ss$sdcor               ## table of random effect SDs and correlations
# ss$theta               ## table of random effects parameters, Cholesky scale


# Complex models

m6 <- glmer(switch ~ group * nasas * task_perf + (1 | id),
            c3, family = binomial) 
summary(m6)
# Estimate Std. Error z value    Pr(>|z|)    
# (Intercept)                            -0.95566    0.18884  -5.061 0.000000418 ***
# group.1-0                               0.93917    0.36774   2.554    0.010653 *  
# nasa_bht                               -0.20049    0.04859  -4.126 0.000036917 ***
# task_perf.shape-loc                    -0.34908    0.09471  -3.686    0.000228 ***
# group.1-0:nasa_bht                     -0.21769    0.09479  -2.296    0.021648 *  
# group.1-0:task_perf.shape-loc           0.36834    0.18709   1.969    0.048977 *  
# nasa_bht:task_perf.shape-loc            0.09583    0.02548   3.761    0.000169 ***
# group.1-0:nasa_bht:task_perf.shape-loc -0.09899    0.05036  -1.966    0.049323 *  

anova(m0,m6)
 
m7 <- glmer(switch ~ group * nasa_bht * task_perf * block + (1 | id),
            c3, family = binomial)  # didn't converge
summary(m7)

m8 <- glmer(switch ~ group * nasa_bht + block + (1 | id),
            c3, family = binomial)  #  
summary(m8)

anova(m6,m8)

m9 <- glmer(switch ~ group * nasa_bht * block + (1 | id),
            c3, family = binomial) # didn't converge
summary(m9)

m10 <- glmer(switch ~ group * nasa_bht + task_perf + block + (1 | id),
            c3, family = binomial)  
summary(m10)

anova(m8,m10)

m11 <- glmer(switch ~ group * nasa_bht * block + (1 | id),
            c3, family = binomial,
            glmerControl(optimizer="bobyqa"))  
summary(m11)
#                                         Estimate    SE        z     p   
# (Intercept)                             -0.88904    0.19696  -4.514 6.37e-06 ***
# group.1-0                                0.92450    0.38034   2.431 0.015069 *  
# nasa_bht                                -0.19851    0.05060  -3.923 8.75e-05 ***
# block.mixed-freeOnly                    -0.23725    0.08405  -2.823 0.004761 ** 
# group.1-0:nasa_bht                      -0.21673    0.09790  -2.214 0.026838 *  
# group.1-0:block.mixed-freeOnly          -0.63962    0.16705  -3.829 0.000129 ***
# nasa_bht:block.mixed-freeOnly            0.28459    0.02229  12.769  < 2e-16 ***
# group.1-0:nasa_bht:block.mixed-freeOnly  0.17064    0.04430   3.851 0.000117 ***

anova(m8,m10,m11,m12)

m12 <- glmer(switch ~ group * nasas * block * task_perf + (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa"))  
summary(m12)
tab_model(m12)
myprint(m12)
tab_model(m12, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          file='m12.html',
          pred.labels=c('(Intercept)','Group','NASA','Block','Task',
                        'Group x NASA','Group x Block','NASA x Block','Group x Task','NASA x Task','Block x Task',
                        'Group x NASA x Block','Group x NASA x Task','Group x Block x Task','NASA x Block x Task',
                        'Group x NASA x Block x Task'),
          dv.labels='')

# Term                     LogOdds   OddsRatio        d      SE         z       p  Sig
# (Intercept)               -1.620       0.198   -0.893   0.055   -29.724   0.000  ***
# group                      0.119       1.126    0.065   0.109     1.091   0.275
# nasas                     -0.211       0.809   -0.117   0.055    -3.856   0.000  ***
# block                      0.814       2.257    0.449   0.023    35.709   0.000  ***
# task                       0.041       1.042    0.023   0.026     1.618   0.106
# group:nasas               -0.231       0.794   -0.127   0.110    -2.110   0.035  *
# group:block                0.001       1.001    0.001   0.046     0.029   0.977
# nasas:block                0.304       1.355    0.167   0.024    12.502   0.000  ***
# group:task                 0.010       1.010    0.005   0.051     0.196   0.845
# nasas:task                 0.099       1.104    0.054   0.028     3.482   0.000  ***
# block:task                -0.059       0.943   -0.032   0.047    -1.258   0.208
# group:nasas:block          0.179       1.196    0.099   0.049     3.688   0.000  ***
# group:nasas:task          -0.096       0.908   -0.053   0.057    -1.701   0.089
# group:block:task           0.153       1.165    0.084   0.093     1.635   0.102
# nasas:block:task          -0.081       0.922   -0.045   0.049    -1.645   0.100
# group:nasas:block:task    -0.060       0.942   -0.033   0.098    -0.606   0.544

m12b <- glmer(switch ~ group * nasas * block * task_perf + (group | id),
               c3, family = binomial,
               glmerControl(optimizer="bobyqa"))  
summary(m12b)
anova(m12,m12b) # m12
anova(m11,m12)

m12c <- glmer(switch ~ group * nasas * block * task_perf + (nasas | id),
              c3, family = binomial,
              glmerControl(optimizer="bobyqa"))  
summary(m12c)
anova(m12,m12c) # m12

m12d <- glmer(switch ~ group * nasas * block * task_perf + (block | id),
              c3, family = binomial,
              glmerControl(optimizer="bobyqa"))  
summary(m12d)
anova(m12,m12d) # m12d

m12e <- glmer(switch ~ group * nasas * block * task_perf + (task_perf | id),
              c3, family = binomial,
              glmerControl(optimizer="bobyqa"))  
summary(m12e)
anova(m12,m12e) 

m12f <- glmer(switch ~ group * nasas * block * task_perf + (group + nasas | id),
              c3, family = binomial,
              glmerControl(optimizer="bobyqa"))  
summary(m12f)
anova(m12,m12d,m12e,m12f) 

m12g <- glmer(switch ~ group * nasas * block * task_perf + (group + nasas + block| id),
              c3, family = binomial,
              glmerControl(optimizer="bobyqa"))  # didnt converge
summary(m12g)

m12h <- glmer(switch ~ group * nasas * block * task_perf + (group + block| id),
              c3, family = binomial,
              glmerControl(optimizer="bobyqa"))  
summary(m12h)
anova(m12d,m12h) #m12d

m13 <- glmer(switch ~ group * nasas * block + (group| id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa"))  # 
summary(m13)

anova(m12,m13)

m14 <- glmer(switch ~ (group * nasa_bht) +
               (group * block) +
               (group * task_perf) + 
               (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa")) 
summary(m14)
anova(m12,m14)

m15 <- glmer(switch ~ (group * nasa_bht * block) +
          (group * block * task_perf) +
          (group * nasa_bht * task_perf) + 
          (1 | id),
        c3, family = binomial,
        glmerControl(optimizer="bobyqa"))  # didnt converge 

summary(m15)

m16 <- glmer(switch ~ (group * nasa_bht * block) +
               (group * block * task_perf) +
               (group * nasa_bht * task_perf) + 
               (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="Nelder_Mead"))  # didnt converge  (tried 3 optimizers)

summary(m16)

# Consider counterbalance
m17 <- glmer(switch ~ group * nasa_bht * block * task_perf * cBal + (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa"))  # didnt converge
summary(m17)


m18 <- glmer(switch ~ group * nasas * block * cBal + (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa"))  #  
summary(m18)
myprint(m18)

emmip(m18, group ~ nasas | block * cBal, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-SE, ymax=yvar+SE, fill = group), alpha = 0.3, colour = NA) 

# Term                         OddsRatio        d      SE         z       p  Sig 
# (Intercept)                      0.199   -0.891   0.055   -29.560   0.000  *** 
# group1                           1.134    0.069   0.109     1.152   0.249      
# nasas                            0.806   -0.119   0.056    -3.890   0.000  *** 
# block1                           2.247    0.446   0.023    35.462   0.000  *** 
# cBal1                            0.996   -0.002   0.109    -0.033   0.974      
# group1:nasas                     0.827   -0.104   0.111    -1.714   0.087      
# group1:block1                    1.006    0.003   0.046     0.137   0.891      
# nasas:block1                     1.385    0.179   0.025    13.273   0.000  *** 
# group1:cBal1                     0.634   -0.252   0.214    -2.131   0.033  *   
# nasas:cBal1                      0.798   -0.125   0.111    -2.043   0.041  *   
# block1:cBal1                     1.242    0.119   0.046     4.755   0.000  *** 
# group1:nasas:block1              1.194    0.098   0.049     3.616   0.000  *** 
# group1:nasas:cBal1               0.775   -0.140   0.217    -1.173   0.241      
# group1:block1:cBal1              1.769    0.314   0.091     6.276   0.000  *** 
# nasas:block1:cBal1               1.285    0.138   0.049     5.114   0.000  *** 
# group1:nasas:block1:cBal1        1.044    0.024   0.098     0.438   0.661    

m19 <- glmer(switch ~ (group * nasa_bht) +
               (group * block * cBal) + (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa"))  #  
summary(m19)



anova(m12,m19,m20)

m20 <- glmer(switch ~ group * block * cBal + (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa"))  #  
summary(m20)
myprint(m20)
anova(m18,m20)
tab_model(m20, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          file='m20.html',
          pred.labels=c('(Intercept)','Group','Block','Cb',
                        'Group x Block','Group x Cb','Block x Cb',
                        'Group x Block x Cb'),
          dv.labels='')
# Term                   OddsRatio        d      SE         z       p  Sig 
# (Intercept)                0.203   -0.880   0.055   -28.861   0.000  *** 
# group1                     1.246    0.121   0.110     1.993   0.046  *   
# block1                     2.144    0.421   0.022    34.975   0.000  *** 
# cBal1                      1.012    0.006   0.110     0.104   0.917      
# group1:block1              0.836   -0.099   0.044    -4.114   0.000  *** 
# group1:cBal1               0.707   -0.191   0.220    -1.575   0.115      
# block1:cBal1               1.216    0.108   0.044     4.488   0.000  *** 
# group1:block1:cBal1        1.430    0.197   0.087     4.115   0.000  *** 

emmip(m20, group ~  block | cBal, cov.reduce = range) +
  geom_errorbar(aes(ymin=yvar-SE, ymax=yvar+SE))   # looks like starting with mix increases effort (instructions increase SoA?)
  
m21 <- glmer(switch ~ (group * nasas * cBal) +
               (group * block * cBal) + (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa"))  #  
summary(m21)

# Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)       -0.84408    0.19861  -4.250           0.00002139 ***
# group              0.78491    0.38890   2.018               0.0436 *  
# nasa              -0.20482    0.05093  -4.021           0.00005784 ***
# cBal               0.73409    0.38837   1.890               0.0587 .  
# block              0.76276    0.02181  34.968 < 0.0000000000000002 ***
# group:nasa        -0.17459    0.09981  -1.749               0.0803 .  
# group:cBal         0.29335    0.72536   0.404               0.6859    
# nasa:cBal         -0.19250    0.09978  -1.929               0.0537 .  
# group:block       -0.17823    0.04360  -4.088           0.00004358 ***
# cBal:block         0.19650    0.04360   4.507           0.00000656 ***
# group:nasa:cBal   -0.19387    0.18703  -1.037               0.2999    
# group:cBal:block   0.35868    0.08701   4.122           0.00003752 ***

anova(m12,m20,m21) # m12 > m21 > m20

m22 <- glmer(switch ~ (group * nasas) + (1 | id),
             c3, family = binomial)  #  
summary(m22)

# Plotting VSR ------------------------------------------------------------

# pre-registered 
mm <- m0 # choose model to plot
c3m <- c3 %>% filter(!is.na(switch))
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(y='Fit', x = 'NASA score', color = 'Group', fill = 'Group')

ggsave('m1.png', path = 'plots_1/', width = 5.5, height = 4)

c3 %>% # raw plot
  group_by(group, nasa_bht) %>% 
  get_summary_stats(switch) %>% 
  ggplot(aes(y = mean, x = nasa_bht, color = group, group = group, fill = group)) +
  geom_ribbon(aes(ymin = mean-2*se, ymax = mean+2*se), colour = NA, alpha = .3) +
  geom_line()

c3 %>% distinct(nasa_bht)

# more complex with nasa
mm <- m12 # choose model to plot
c3m <- c3 %>% anti_join(c3 %>% filter(if_any(c(switch, group, nasa_bht, task_perf, block), ~ is.na(.))))
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  facet_grid(block~task_perf) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(y='Fit', x = 'NASA score', color = 'Group', fill = 'Group') +
  theme(legend.position=c(.85,.15), legend.key.size = unit(.5, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))
 
ggsave('m12.png', path = 'plots_1/', w=5,h=4.5)

# with counterbalance
mm <- m20 # choose model to plot
c3m <- c3 %>% anti_join(c3 %>% filter(if_any(c(switch, group, block, cBal), ~ is.na(.))))
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
ggplot(c3m, aes(y = fit, x = block, color = group, fill = group)) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=plo, ymax=phi), position=pd, width = .4) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  facet_wrap(~cBal) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  labs(y='Fit', x = 'Block type', color = 'Group', fill = 'Group') +
  theme(legend.position=c(.85,.15), legend.key.size = unit(.5, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))

ggsave('m20.png', path = 'plots_1/', h=4, w=5.5)

mm <- m18 # choose model to plot
c3m <- c3 %>% anti_join(c3 %>% filter(if_any(c(switch, group, block, cBal), ~ is.na(.))))
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
  facet_grid(block~cBal) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(y='Fit', x = 'NASA score', color = 'Group', fill = 'Group') +
  theme(legend.position=c(.85,.15), legend.key.size = unit(.5, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))

ggsave('m18.png', path = 'plots_1/', h=4, w=5)


# Contrasts VSR ---------------------------------------------------------------

# Pre-registered models
summary(m0)
tab_model(m0)

contrast(emtrends(m0, ~ group, var = 'nasas'),
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df) %>% kable() # pairwise comparisons

test(emtrends(m0, ~ group, var = 'nasas'), adjust = "fdr") %>% # test slopes against 0
  as_tibble() %>% select(-df) %>% kable()

emmeans(m0, pairwise ~ group | nasas, at = list(nasas = c(-1, 0, 1)))$contrasts %>% 
  as_tibble %>% select(-df) %>% kable  # compare effect of group on 3 points of NASA

emmeans(m0, pairwise ~ group | nasas, at = list(nasas = c(-1, 0, 1)))$contrasts %>% 
  as_tibble %>% select(-df)  %>% mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>% kable

emmeans(m0, pairwise ~ group | nasas, at = list(nasas = c(-2, -1, 0, 1, 2)))$contrasts %>% 
  as_tibble %>% select(-df)  %>% mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>% kable

# Complex models
summary(m12)
emmip(m12, group ~ nasas | block * task_perf, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

contrast(emtrends(m12, ~ group * task_perf * block, var = 'nasas'),
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% # anova/contrast coding (-+0.5), compares btw levels 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE,-z.ratio)

contrast(emmeans(m12, ~ group * task_perf * block), # check 'mini main effect' of group in freeOnly, on averaged nasa
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% # anova/contrast coding (-+0.5), compares btw levels 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE,-z.ratio)

test(emtrends(m12, ~ group * block * task_perf, var = 'nasas'), adjust = "fdr") %>% # test slopes against 0
  as_tibble() %>% select(-df) %>% kable()

emmeans(m12, pairwise ~ group | nasas | task_perf | block, at = list(nasas = c(-1, 0, 1)))$contrasts %>% 
  as_tibble %>% select(-df, -SE, -z.ratio) %>% kable  # compare effect of group on 3 points of NASA

emmeans(m12, pairwise ~ group | nasas | task_perf | block, at = list(nasas = c(-1, 0, 1)))$contrasts %>% 
  as_tibble %>% select(-df, -SE, -z.ratio) %>% mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>% kable # add p adjustment

contrast(emmeans(m12, ~ group * task_perf * block * nasas, at = list(nasas = c(-1, 0, 1))), 
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr")

contrast(emmeans(m12, ~ group | task_perf | block | nasas, at = list(nasas = c(-1, 0, 1))), 
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr")

# Counterbalance models
summary(m20)
emmip(m20, group ~  block | cBal) +
  geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)))   # looks like starting with mix increases effort (instructions increase SoA?)
 
contrast(emmeans(m20, ~ group * block * cBal),
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% # anova/contrast coding (-+0.5), compares btw levels 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE,-z.ratio) %>% kable()

emmip(m21, group ~  nasa_bht  | block * cBal, cov.reduce = range) +
  geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4) 

test(emtrends(m21, ~ group * block * cBal, var = 'nasa_bht'), adjust = "fdr") %>% # test slopes against 0
  as_tibble() %>% select(-df) %>% kable()

# playing with contrasts
emmeans(m12, 'nasa_bht', cov.reduce = range)
emmeans(m12, 'group')
m12tr <- emtrends(m12, var = "nasa_bht", 
         cov.reduce = nasa_bht ~ group * task_perf)
summary(m12tr)  # estimated trends at grid nodes
pairs(emmeans(m12tr, "task_perf", weights = "prop"))
pairs(emmeans(m12tr, "group", weights = "prop"))

emtrends(m12, pairwise ~ group, var = 'nasa_bht')
emtrends(m12, pairwise ~ group * task_perf, var = 'nasa_bht', adjust = 'FDR')
emtrends(m12, pairwise ~ group * task_perf * block, var = 'nasa_bht', adjust = 'FDR')

pairs(emtrends(m12, ~ group, var = 'nasa_bht'))
contrast(emtrends(m12, ~ group, var = 'nasa_bht'), "revpairwise")

emtrends(m12, pairwise ~ group * task_perf * block, var = 'nasa_bht', adjust = 'none')[[2]] %>%
  # summary(infer = TRUE) %>%
  as_tibble() %>%
  select(-df,-SE,-z.ratio) %>%
  adjust_pvalue(p.col = 'p.value', output.col = 'adj.p', method = 'fdr') %>% 
  mutate_if(is.numeric, round, 4) %>% print(n=Inf) 


contrast(emtrends(m12, ~ group * task_perf, var = 'nasa_bht'),
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df)

contrast(emtrends(m12, ~ group * task_perf, var = 'nasa_bht'),
         "consec", simple = 'each', combine = TRUE, adjust = "fdr") %>% 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df)

contrast(emtrends(m12, ~ group * task_perf, var = 'nasa_bht'),
         "eff", simple = 'each', combine = TRUE, adjust = "fdr") %>% 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df)

contrast(emtrends(m12, ~ group * task_perf, var = 'nasa_bht'),
         "trt.vs.ctrl", simple = 'each', combine = TRUE, adjust = "fdr") %>% 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df)
#

contrast(emtrends(m12, ~ group * task_perf * block, var = 'nasa_bht'),
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% # anova/contrast coding (-+0.5), compares btw levels 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE,-z.ratio)

contrast(emtrends(m12, ~ group * task_perf * block, var = 'nasa_bht'),
         "pairwise", simple = 'each', combine = TRUE, adjust = "none") %>% # anova/contrast coding (-+0.5), compares btw levels 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE,-z.ratio)

contrast(emtrends(m12, ~ group * task_perf * block, var = 'nasa_bht'),
         "consec", simple = 'each', combine = TRUE, adjust = "fdr") %>% 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df)

contrast(emtrends(m12, ~ group * task_perf * block, var = 'nasa_bht'),
         "eff", simple = 'each', combine = TRUE, adjust = "fdr") %>% # sum coding (-+1), compares lvl to grand mean, half of anova
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df)

contrast(emtrends(m12, ~ group * task_perf * block, var = 'nasa_bht'),
         "trt.vs.ctrl", simple = 'each', combine = TRUE, adjust = "fdr") %>% 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df)


# Modelling VSR early subset ----------------------------------------------

# Prepare
c3s <- c3 %>% 
  filter(cBal == 'free-mix',
         block == 'freeOnly') %>% 
  mutate(trial = case_when(blockNr == 1 ~ trialNr,
                           blockNr == 2 ~ trialNr + 72),
         trials = standardise(trial))

nrow(c3s)/nrow(c3)

c3s %>% freq_table(blockNr, blockNrC)
c3s %>% get_summary_stats(trial, trials)

ggplot(c3s, aes(x=trialNr)) + geom_histogram(binwidth = 1)
ggplot(c3s, aes(x=trials)) + geom_histogram(binwidth = 0.01)


# Simple
ms0 <- glmer(switch ~ group * nasas + (1 | id),
            c3s, family = binomial)
summary(ms0)
tab_model(ms0, show.stat = T, show.icc = F)
tab_model(ms0, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          # file='m0.html',
          pred.labels=c('(Intercept)','Group','NASA','Group x NASA'),
          dv.labels='')

emmip(ms0, group ~ nasas, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 


# Complex
ms1 <- glmer(switch ~ group * nasas * blockNr + (1 | id),
             c3s, family = binomial,
             glmerControl(optimizer = 'bobyqa'))
summary(ms1)
anova(ms0,ms1)

emmip(ms1, group ~ nasas | blockNr, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

ms2 <- glmer(switch ~ group * nasas * trial + (1 | id),
             c3s, family = binomial,
             glmerControl(optimizer = 'bobyqa')) # trials not scaled
summary(ms2)
anova(ms0,ms1,ms2)

emmip(ms2, group ~  trial, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

ms3 <- glmer(switch ~ group * nasas * trials + (1 | id),
             c3s, family = binomial,
             glmerControl(optimizer = 'bobyqa'))
summary(ms3)
anova(ms1,ms3) # ms3
tab_model(ms3, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          file='ms3.html',
          pred.labels=c('(Intercept)','Group','NASA','Trial', 'Group x NASA', 'Group x Trial', 'NASA x Trial', 'Group x NASA x Trial'),
          dv.labels='')

ms3b <- glmer(switch ~ group * nasas * trials + (group | id),
             c3s, family = binomial,
             glmerControl(optimizer = 'bobyqa'))
summary(ms3b)
anova(ms3, ms3b)

1.14008 + (1.96*4.497); 1.14008 - (1.96*4.497) # heterogeneity interval: fixed effect +/- 1.96 x random effect SD; not sure if it makes sense for btwn-sbj variable though


ms3c <- glmer(switch ~ group * nasas * trials + (nasas | id),
              c3s, family = binomial,
              glmerControl(optimizer = 'bobyqa'))
summary(ms3c)
anova(ms3, ms3c)

ms3d <- glmer(switch ~ group * nasas * trials + (trials | id),
              c3s, family = binomial,
              glmerControl(optimizer = 'bobyqa'))
summary(ms3d)
anova(ms3, ms3d) # ms3d

tab_model(ms3d, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          # file='ms3d.html',
          pred.labels=c('(Intercept)','Group','NASA','Trial', 'Group x NASA', 'Group x Trial', 'NASA x Trial', 'Group x NASA x Trial'),
          dv.labels='')

-0.54849 + (1.96*0.6311); -0.54849 - (1.96*0.6311) # heterogeneity interval: fixed effect +/- 1.96 x random effect SD 

ms3e <- glmer(switch ~ group * nasas * trials + (group * trials | id),
              c3s, family = binomial,
              glmerControl(optimizer = 'bobyqa'))
summary(ms3e)
anova(ms3d, ms3e)

ms3f <- glmer(switch ~ group * nasas * trials + (group * trials + nasas| id),
              c3s, family = binomial,
              glmerControl(optimizer = 'bobyqa')) # not conv
summary(ms3f)
anova(ms3f, ms3d)# ms3d


ms4 <- glmer(switch ~ group * trials + (1 | id),
             c3s, family = binomial,
             glmerControl(optimizer = 'bobyqa'))
summary(ms4)
anova(ms3,ms4) # ms3

ms5 <- glmer(switch ~ group * nasas * trials * task_perf + (1 | id),
             c3s, family = binomial,
             glmerControl(optimizer = 'bobyqa'))
summary(ms5)
anova(ms3,ms5)


# Contrasts & plots VSR early subset --------------------------------------


emmip(ms3, group ~ trials | nasas, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

emmip(ms3, group ~ trials | nasas, cov.reduce = range, at = list(nasas = c(-2, -1, 0, 1, 2, 3))) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

emmip(ms3, trials ~ nasas , cov.reduce = range,         # effect of nasa starts to appear only in later trials
      at = list(nasas = c(-2, -1, 0, 1, 2, 3), trials = c(-1.72,0,1.72))) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = trials), alpha = 0.3, colour = NA) +
  # scale_colour_discrete(labels = c('first', 'middle', 'last')) + 
  guides(fill = 'none', group = 'none')  

facet_group <- c('0' = 'No control', '1' = 'Full control')
facet_trial <- c('2' = 'First trial', '72' = 'Middle trial', '144' = 'Last trial')


emmip(ms3, trials ~ nasas , cov.reduce = range,         # effect of nasa starts to appear only in later trials
      at = list(nasas = c(-2, -1, 0, 1, 2, 3), trials = c(-1.72,0,1.72)), plotit = F)  %>% 
  ggplot(., aes(x = xvar, y = yvar, colour = tvar)) +
  geom_line() +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = tvar), alpha = 0.2, colour = NA) +
  scale_colour_discrete(labels = c('first', 'middle', 'last')) + 
  guides(fill = 'none') +
  labs(y = 'Linear prediction', x = 'NASA', colour = 'Trial')

emmip(ms3, trials ~ nasas | group, cov.reduce = range,         # effect of nasa starts to appear only in later trials
      at = list(nasas = c(-2, -1, 0, 1, 2, 3), trials = c(-1.72,0,1.72)), plotit = F)  %>% 
  ggplot(., aes(x = xvar, y = yvar, colour = tvar)) +
  geom_line() +
  facet_wrap(~group, labeller = as_labeller(facet_group)) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = tvar), alpha = 0.15, colour = NA) +
  # scale_color_viridis(labels = c('first', 'middle', 'last')) +
  scale_colour_brewer(labels = c('first', 'middle', 'last'), palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = 'none') +
  labs(y = 'Linear prediction', x = 'NASA', colour = 'Trial')

# more complex 
mm <- ms3 # choose model to plot
c3m <- c3s %>% anti_join(c3s %>% filter(if_any(c(switch, group, nasa_bht, task_perf, block), ~ is.na(.))))
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
c3m %>% filter(trial %in% c(2,72,144)) %>% 
  mutate(trial = as.factor(trial)) %>% 
  ggplot(., aes(y = fit, x = nasas, color = trial, fill = trial)) +
  geom_line(linewidth = .8) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = trial, fill = trial), alpha = 0.2, colour = NA) +
  facet_wrap(~group, labeller = as_labeller(facet_group)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_colour_brewer(labels = c('first', 'middle', 'last'), palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = 'none') +
  labs(y = 'Linear prediction', x = 'NASA', colour = 'Trial') +
  theme(legend.position=c(.62,.2), legend.key.size = unit(.5, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))

ggsave('ms3v1.png', path = 'plots_1/', width = 5.5, height = 4)

  
c3m %>% filter(trial %in% c(2,72,144)) %>% 
  mutate(trial = as.factor(trial)) %>% 
  ggplot(., aes(y = fit, x = nasas, color = group, fill = group)) +
  geom_line(linewidth = .8) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.2, colour = NA) +
  facet_wrap(~trial, labeller = as_labeller(facet_trial)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_colour_discrete(labels = as_labeller(facet_group)) +
  guides(fill = 'none') +
  labs(y = 'Linear prediction', x = 'NASA', colour = 'Group') +
  theme(legend.position=c(.11,.15), legend.key.size = unit(.5, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))

ggsave('ms3v2.png', path = 'plots_1/', width = 7, height = 4)


emmip(ms2, trial ~ nasas , cov.reduce = range, 
      at = list(nasas = c(-2, -1, 0, 1, 2, 3), trial = c(1,72,144))) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), alpha = 0.3, colour = NA) 

emmip(ms3, trials ~ nasas | group, cov.reduce = range, # no interaction with group but trend of nasa only in full control
      at = list(nasas = c(-2, -1, 0, 1, 2, 3), trials = c(-1.72,0,1.72))) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), alpha = 0.3, colour = NA) 

emmip(ms3, trials ~ nasas , cov.reduce = range, 
      at = list(nasas = c(-2, -1, 0, 1, 2, 3), trials = c(-1,0,1)), plotit=F) 

emmip(ms3, group ~ trials, cov.reduce = range) + # effect of group does not wane due to trial nr
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 


# contrast effect of group on first, middle, last trial
emmeans(ms3, pairwise ~ group | trials, at = list(trials = c(-1.72, 0, 1.72)))$contrasts %>% 
  as_tibble %>% select(-df, -SE, -z.ratio) %>% mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>% kable # add p adjustment

# interactions with task
emmip(ms5, group  ~ nasas | trials * task_perf , cov.reduce = range, 
        at = list(nasas = c(-2, -1, 0, 1, 2, 3), trials = c(-1.72,0,1.72)))

emmip(ms5, group  ~ trials  | nasas * task_perf , cov.reduce = range, 
      at = list(nasas = c(-2, -1, 0, 1, 2), trials = c(-1.72,0,1.72)))


# Modelling VSR with NASA VTS ---------------------------------------------

mv1 <- glmer(switch ~ group * nasas * nasasV + (1 | id),
            c3, family = binomial) 
summary(mv1)

emmip(mv1, group ~ nasasV, cov.reduce = range) + #  
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

mv2 <- glmer(switch ~ group * nasas * nasasV * block + (1 | id),
             c3, family = binomial) # not conv 
summary(mv2)

mv3 <- glmer(switch ~ group * nasasV * block * task_perf + (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa")) 
summary(mv3)
anova(mv1,mv3) # mv3
anova(m12,mv3) # m12

tab_model(mv3, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          # file='mv3.html',
          pred.labels=c('(Intercept)','Group','NASA(VTS)','Block','Task',
                        'Group x NASA(VTS)','Group x Block','NASA(VTS) x Block','Group x Task','NASA(VTS) x Task','Block x Task',
                        'Group x NASA(VTS) x Block','Group x NASA(VTS) x Task','Group x Block x Task','NASA(VTS) x Block x Task',
                        'Group x NASA(VTS) x Block x Task'),
          dv.labels='')
emmip(mv3, group ~ nasasV | block, cov.reduce = range) + #  no control group withdraws effort if NASA VTS is high
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA)  
  

emmip(mv3, group ~ task_perf | block) + #  no control group withdraws effort if NASA VTS is high
  geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4)  

mm <- mv3 # choose model to plot
c3m <- c3 %>% anti_join(c3 %>% filter(if_any(c(switch, group, task_perf, block), ~ is.na(.))))
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
c3m %>% 
  group_by(group, nasasV, block) %>% 
  summarise(fit = mean(fit), plo = mean(plo), phi = mean(phi)) %>% 
  ggplot(., aes(y = fit, x = nasasV, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  facet_grid(~block) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(y='Fit', x = 'NASA score', color = 'Group', fill = 'Group') +
  theme(legend.position=c(.85,.15), legend.key.size = unit(.5, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))

mv3lm  <- lm(fit ~ group * nasasV * block,cc)
summary(mv3lm)


# Modelling VSR early subset with NASA VTS --------------------------------

msv3 <- glmer(switch ~ group * nasasV * trials + (trials | id),
              c3s, family = binomial,
              glmerControl(optimizer = 'bobyqa'))
summary(msv3)
anova(ms3d,msv3) # pretty much equal, perhaps msv3 >= ms3d

emmip(msv3, group ~ nasasV, cov.reduce = range) + #  no control group withdraws effort if NASA VTS is high
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 


# Modelling Task selection ------------------------------------------------

# Pre-registered models
mt1 <- glmer(task_dv ~ group * nasas + (1 | id),
             c3, family = binomial)
summary(mt1)
tab_model(mt1)
tab_model(mt1, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          file='mt1.html',
          pred.labels=c('(Intercept)','Group','NASA','Group x NASA'), 
          dv.labels='')
myprint(mt1)

mt2 <- glmer(task_dv ~ group * nasas + (group | id),
             c3, family = binomial)
summary(mt2)
tab_model(mt2)

anova(mt1,mt2) # mt1

# Complex models

mt3 <- glmer(task_dv ~ group * nasas * block + (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa"))  #  
summary(mt3)
myprint(mt3)
tab_model(mt3, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          # file='mt1.html',
          # pred.labels=c('(Intercept)','Group','NASA','Group x NASA'), 
          dv.labels='')

anova(mt1,mt3) # mt3

mt4 <- glmer(task_dv ~ group * nasas * block * switch_iv + (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa"))  #  
summary(mt4)

anova(mt1,mt3,mt4) #  mt4

mt5 <- glmer(task_dv ~ (group * nasas * block) +
               (group * nasas * switch_iv) +
               (block * nasas * switch_iv) +
               (group * block * switch_iv) +
               (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa"))  #  
summary(mt5)
anova(mt4,mt5) # mt5, slightly
myprint(mt5)
tab_model(mt5, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          file='mt5.html',
          pred.labels=c('(Intercept)','Group','NASA','Block','Switch',
                        'Group x NASA','Group x Block','NASA x Block','Group x Switch','NASA x Switch','Block x Switch',
                        'Group x NASA x Block','Group x NASA x Switch','Group x Block x Switch','NASA x Block x Switch'),
          dv.labels='')

# Term                        OddsRatio    d      SE        z       p  Sig
# (Intercept)                 0.579   -0.302   0.131   -4.177   0.000  ***
# group1                      1.167    0.085   0.258    0.598   0.550
# nasas                       0.979   -0.012   0.131   -0.162   0.871
# block1                      0.891   -0.064   0.025   -4.695   0.000  ***
# switch                      1.013    0.007   0.026    0.484   0.628
# group1:nasas                0.913   -0.050   0.257   -0.354   0.723
# group1:block1               1.096    0.051   0.049    1.866   0.062
# nasas:block1                1.140    0.072   0.027    4.798   0.000  ***
# group1:switch               1.029    0.016   0.052    0.549   0.583
# nasas:switch                1.090    0.048   0.029    3.025   0.002  **
# block1:switch               1.065    0.035   0.049    1.281   0.200
# group1:nasas:block1         1.013    0.007   0.044    0.301   0.763
# group1:nasas:switch         0.899   -0.059   0.057   -1.877   0.061
# nasas:block1:switch         0.923   -0.044   0.055   -1.469   0.142
# group1:block1:switch        1.349    0.165   0.100    3.003   0.003  **

anova(mt4s,mt5)

mt6 <- glmer(task_dv ~ nasa_bht * block * switch_iv +
               (1 | id),
             c3, family = binomial,
             glmerControl(optimizer="bobyqa"))  #  didnt converge
summary(mt6)

anova(mt5,mt6)

mt7 <- glmer(task_dv ~ group * block * switch_iv +
             (1 | id),
             c3, family = binomial)  #  
summary(mt7)    
anova(mt5,mt7)

# Plotting Task selection -------------------------------------------------

# Pre-registered 
mm <- mt1 # choose model to plot
c3m <- c3 %>% filter(!is.na(task_dv))
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(y='Fit (task selection)', x = 'NASA score', color = 'Group', fill = 'Group')

ggsave('mt1.png', path = 'plots_1/', width = 5.5, height = 4)

emmip(mt1, group ~ nasa_bht, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 


# Complex

emmip(mt4, group ~ nasa_bht | block, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

emmip(mt5, group ~ nasa_bht | block * switch_iv, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

emmip(mt5, ~ nasa_bht | block * switch_iv, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), alpha = 0.3, colour = NA) 

mm <- mt5 # choose model to plot
c3m <- c3 %>% filter(!is.na(switch_iv))
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
facet_switch <- c('0' = 'Repeat', '1' = 'Switch')
facet_block_switch <- c('0' = 'Repeat', '1' = 'Switch', 'freeOnly' = 'freeOnly', 'mixed' = 'mixed')


ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
  facet_grid(block~switch_iv, labeller = labeller(switch_iv = as_labeller(facet_switch))) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(y='Fit (task selection)', x = 'NASA score', color = 'Group', fill = 'Group') +
  theme(legend.position='bottom') +
  theme(legend.position=c(.86,.1), legend.key.size = unit(.4, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))

c3m %>% group_by(nasas, switch_iv,block) %>% 
  get_summary_stats(fit,phi,plo) %>% 
  select(block:variable,mean) %>% 
  pivot_wider(names_from = variable, values_from = mean) %>% 
  ggplot(., aes(y = fit, x = nasas)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin=plo, ymax=phi), alpha = 0.3, colour = NA) +
  facet_grid(block~switch_iv, labeller = labeller(switch_iv = as_labeller(facet_switch))) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(y='Fit (task selection)', x = 'NASA score') +
  theme(legend.position='bottom') +
  theme(legend.position=c(.86,.1), legend.key.size = unit(.4, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))

emmip(mt5, ~ block)
emmip(mt5, ~ nasas | block, cov.reduce = range) + geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), alpha = 0.3) 
emmip(mt5, ~ nasas | switch_iv, cov.reduce = range) + geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), alpha = 0.3) 
emmip(mt5, ~ nasas | switch_iv * block, cov.reduce = range) + geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), alpha = 0.3) 
emmip(mt5, ~ nasas | block * switch_iv, cov.reduce = range) + 
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), alpha = 0.3) +
  facet_grid(block~switch_iv, labeller = labeller(switch_iv = as_labeller(facet_switch))) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(y='Fit (task selection)', x = 'NASA score') +
  theme(legend.position='bottom') +
  theme(legend.position=c(.86,.1), legend.key.size = unit(.4, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))

ggsave('mt5.png', path = 'plots_1/', width = 5.5, height = 5.5)

# Contrasts Task selection -------------------------------------------------

# Simple
contrast(emtrends(mt4, ~ group * block, var = 'nasa_bht'),
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% # anova/contrast coding (-+0.5), compares btw levels 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE,-z.ratio)

test(emtrends(mt4, ~ group * block, var = 'nasa_bht'), adjust = "fdr") %>% # test slopes against 0
  as_tibble() %>% select(-df) %>% kable()

# Complex
myprint(mt5)
contrast(emtrends(mt5, ~ group * block * switch_iv, var = 'nasas'),
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% # anova/contrast coding (-+0.5), compares btw levels 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE,-z.ratio) %>% kable

test(emtrends(mt5, ~ group * block * switch_iv, var = 'nasas'), adjust = "fdr") %>% # test slopes against 0
  as_tibble() %>% select(-df) %>% kable()

emmeans(mt5, pairwise ~ group | nasas | switch_iv | block, at = list(nasas = c(-1, 0, 1)))$contrasts %>% 
  as_tibble %>% select(-df, -SE, -z.ratio) %>% mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>% kable 

# contrast(emtrends(m12, ~ group * task_perf * block, var = 'nasas'),
#          "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% # anova/contrast coding (-+0.5), compares btw levels 
#   as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE,-z.ratio)
# 
# test(emtrends(m12, ~ group * block * task_perf, var = 'nasas'), adjust = "fdr") %>% # test slopes against 0
#   as_tibble() %>% select(-df) %>% kable()
# 
# emmeans(m12, pairwise ~ group | nasas | task_perf | block, at = list(nasas = c(-1, 0, 1)))$contrasts %>% 
#   as_tibble %>% select(-df, -SE, -z.ratio) %>% kable  # compare effect of group on 3 points of NASA


# Modelling Task selection early subset -----------------------------------

# Simple
mts0 <- glmer(task_dv ~ group * nasas + (1 | id),
             c3s, family = binomial)
summary(mts0)
tab_model(mts0, show.stat = T, show.icc = F)
tab_model(mts0, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          # file='m0.html',
          pred.labels=c('(Intercept)','Group','NASA','Group x NASA'),
          dv.labels='')

emmip(mts0, group ~ nasas, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 


# Complex
mts1 <- glmer(task_dv ~ group * nasas * blockNr + (1 | id),
             c3s, family = binomial,
             glmerControl(optimizer = 'bobyqa'))
summary(mts1)
anova(mts0,mts1)

emmip(mts1, group ~ nasas | blockNr, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

mts3 <- glmer(task_dv ~ group * nasas * trials + (1 | id),
             c3s, family = binomial,
             glmerControl(optimizer = 'bobyqa'))
summary(mts3)
anova(mts1,mts3) 

mts3b <- glmer(task_dv ~ group * nasas * trials + (group | id),
              c3s, family = binomial,
              glmerControl(optimizer = 'bobyqa'))
summary(mts3b)
anova(mts3, mts3b)

mts3c <- glmer(task_dv ~ group * nasas * trials + (nasas | id),
              c3s, family = binomial,
              glmerControl(optimizer = 'bobyqa'))
summary(mts3c)
anova(mts3, mts3c)

mts3d <- glmer(task_dv ~ group * nasas * trials + (trials | id),
              c3s, family = binomial,
              glmerControl(optimizer = 'bobyqa'))
summary(mts3d)
anova(mts3, mts3d) # mts3d

tab_model(mts3d, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          file='mts3d.html',
          pred.labels=c('(Intercept)','Group','NASA','Trial', 'Group x NASA', 'Group x Trial', 'NASA x Trial', 'Group x NASA x Trial'),
          dv.labels='')

mts3e <- glmer(task_dv ~ group * nasas * trials + (group * trials | id),
              c3s, family = binomial,
              glmerControl(optimizer = 'bobyqa'))
summary(mts3e)
anova(mts3d, mts3e) # mts3e


# Contrasts & plotting Task selection early subset ------------------------

emmip(mts3, trials ~ nasas , cov.reduce = range,         # effect of nasa starts to appear only in later trials
      at = list(nasas = c(-2, -1, 0, 1, 2, 3), trials = c(-1.72,0,1.72)), plotit = F)  %>% 
  ggplot(., aes(x = xvar, y = yvar, colour = tvar)) +
  geom_line() +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = tvar), alpha = 0.2, colour = NA) +
  scale_colour_discrete(labels = c('first', 'middle', 'last')) + 
  guides(fill = 'none') +
  labs(y = 'Linear prediction (task selection)', x = 'NASA', colour = 'Trial')

emmip(mts3, trials ~ nasas | group, cov.reduce = range,         # effect of nasa starts to appear only in later trials
      at = list(nasas = c(-2, -1, 0, 1, 2, 3), trials = c(-1.72,0,1.72)), plotit = F)  %>% 
  ggplot(., aes(x = xvar, y = yvar, colour = tvar)) +
  geom_line() +
  facet_wrap(~group, labeller = as_labeller(facet_group)) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = tvar), alpha = 0.15, colour = NA) +
  scale_colour_brewer(labels = c('first', 'middle', 'last'), palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = 'none') +
  labs(y = 'Linear prediction (task selection)', x = 'NASA', colour = 'Trial')

# more complex 
mm <- mts3 # choose model to plot
c3m <- c3s %>% filter(!is.na(task_dv))
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
c3m %>% filter(trial %in% c(2,72,144)) %>% 
  mutate(trial = as.factor(trial)) %>% 
  ggplot(., aes(y = fit, x = nasas, color = trial, fill = trial)) +
  geom_line(linewidth = .8) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = trial, fill = trial), alpha = 0.2, colour = NA) +
  facet_wrap(~group, labeller = as_labeller(facet_group)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_colour_brewer(labels = c('first', 'middle', 'last'), palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = 'none') +
  labs(y = 'Linear prediction (task selection)', x = 'NASA', colour = 'Trial') +
  theme(legend.position=c(.62,.2), legend.key.size = unit(.5, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))

ggsave('mts3v1.png', path = 'plots_1/', width = 5.5, height = 4)

c3m %>% filter(trial %in% c(2,72,144)) %>% 
  mutate(trial = as.factor(trial)) %>% 
  ggplot(., aes(y = fit, x = nasas, color = group, fill = group)) +
  geom_line(linewidth = .8) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.2, colour = NA) +
  facet_wrap(~trial, labeller = as_labeller(facet_trial)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_colour_discrete(labels = as_labeller(facet_group)) +
  guides(fill = 'none') +
  labs(y = 'Linear prediction (task selection)', x = 'NASA', colour = 'Trial') +
  theme(legend.position=c(.11,.15), legend.key.size = unit(.5, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))

ggsave('mts3v2.png', path = 'plots_1/', width = 7, height = 4)

# Modelling rigidness ------------------------------------------------------

# prepare

r2 <- r1 %>% 
  filter(block == 'freeOnly') %>% 
  left_join(c3 %>% group_by(id) %>% slice(1) %>% select(id, cBal, nasas), by = 'id') %>% 
  rename(ntr = n) %>%
  select(-block) %>% 
  mutate(
         # group = contr_code_anova(group),
         # task_perf = contr_code_anova(task_perf),
         # cBal = contr_code_anova(cBal), 
         task_perf = as.factor(task_perf),
         sumSwitch = as.integer(sumSwitch),
         ntr = as.double(ntr),
         nasas = as.double(nasas))

contrasts(r2$group) <- c(-0.5,0.5); contrasts(r2$group)
contrasts(r2$task_perf) <- c(-0.5,0.5); contrasts(r2$task_perf)
contrasts(r2$cBal) <- c(-0.5,0.5); contrasts(r2$cBal)


r3 <- r2 %>% cross_join(r2 %>% ungroup() %>% # clean by 3sd
                          summarise(mean = mean(ntr), thr = 3*(sd(ntr)), up = mean+thr, lo = mean-thr)) %>% 
  mutate(out = if_else(ntr >= up, T, F)) %>% 
  filter(out==F)

nrow(r3)/nrow(r2)


ggplot(r2 %>% cross_join(r2 %>% ungroup() %>% 
                           summarise(mean = mean(ntr), thr = 3*(sd(ntr)), up = mean+thr, lo = mean-thr)) %>% 
         mutate(out = if_else(ntr >= up, T, F)), 
       aes(x=ntr, fill=out)) + 
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = 46, linetype = 'dashed', color= 'dark red') +
  labs(y='Count', x= 'N trials', fill = 'Outlier')

ggplot(r3, aes(x=ntr, fill=out))  + geom_histogram()

ggsave('nTrials_hist.png', path = 'plots_1/', w=5,h=4)



# model

mr0 <- lmer(ntr ~ group * nasas + (1|id), r3)
summary(mr0)

mr0b <- lmer(ntr ~ group * nasas + (nasas|id), r3)
summary(mr0b)

anova(mr0,mr0b)

mr1 <- lmer(ntr ~ group * task_perf + (1|id), r3)
summary(mr1)

mr2 <- lmer(ntr ~ group * task_perf + (task_perf|id), r3)
summary(mr2)

mr3 <- lmer(ntr ~ group * task_perf + (group|id), r3)  
summary(mr3)

anova(mr0,mr1,mr2,mr3) # mr2

mr4 <- lmer(ntr ~ group * nasas  * task_perf+ (1|id), r3)  #  
summary(mr4)

anova(mr3,mr4) # mr4

mr5 <- lmer(ntr ~ group * nasas * task_perf + (task_perf|id), r3)  #  
summary(mr5)
tab_model(mr5, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          file='mr5.html',
          pred.labels=c('(Intercept)','Group','NASA','Task','Group x NASA','Group x Task','NASA x Task',
                        'Group x NASA x Task'),
          dv.labels='')

#                       Estimate Std. Error      df   t value Pr(>|t|)   
# (Intercept)               6.0436     0.3588 173.2263  16.844 <0.0000000000000002 ***
# group1                   -0.4068     0.7176 173.2263  -0.567              0.5716    
# nasas                     0.6651     0.3617 181.9773   1.839              0.0676 .  
# task_perf1               -0.5346     0.3683 104.3327  -1.451              0.1497    
# group1:nasas              1.2461     0.7234 181.9773   1.723              0.0867 .  
# group1:task_perf1         0.5008     0.7366 104.3327   0.680              0.4981    
# nasas:task_perf1         -0.4734     0.3810 112.1724  -1.243              0.2166    
# group1:nasas:task_perf1  -0.7078     0.7620 112.1724  -0.929              0.3549    

anova(mr2,mr4,mr5) # mr5

emmip(mr5, group ~ nasa | task_perf, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-SE, ymax=yvar+SE, fill = group), alpha = 0.3, colour = NA) 

mr6 <- lmer(ntr ~ group * task_perf * nasas + (group|id), r3)   #  
summary(mr6)
#                       Estimate Std. Error      df   t value Pr(>|t|)   
# (Intercept)           3.56827    1.34292  127.32252   2.657  0.00889 **
# group                -4.87650    2.68583  127.32252  -1.816  0.07178 . 
# task                  0.84871    0.38348 6916.86591   2.213  0.02692 * 
# nasa                  0.68574    0.34630  140.57081   1.980  0.04963 * 
# group:task            0.01610    0.76695 6916.86591   0.021  0.98326   
# group:nasa            1.14653    0.69259  140.57081   1.655  0.10007   
# task:nasa            -0.21380    0.10562 6922.88004  -2.024  0.04298 * 
# group:task:nasa       0.03546    0.21123 6922.88005   0.168  0.86670   


anova(mr4,mr5,mr6) # mr5

mr7 <- lmer(ntr ~ group * nasas * task_perf * cBal + (1 | id), r3)  #  
summary(mr7)
anova(mr5,mr7)

mr8 <- lmer(ntr ~ group * task_perf * nasa_bht + (group + task_perf | id), r3,
            lmerControl(optCtrl=list(maxfun=1e5), optimizer = 'bobyqa'))  #  
summary(mr8)

# Plotting rigidness ------------------------------------------------------

mm <- mr1 # choose model to plot
c3m <- r2  
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
ggplot(c3m, aes(y = fit, x = task_perf, color = group)) +
  geom_point(size=3,position=pd) +
  geom_errorbar(aes(ymin=plo, ymax=phi, color = group), position=pd) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  labs(y='Fit (rigidness)', x = 'Task', color = 'Group', fill = 'Group')

ggsave('mr1.png', path = 'plots_1/')


mm <- mr6 # choose model to plot
c3m <- r3  
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
ggplot(c3m, aes(y = fit, x = task_perf, color = group)) +
  geom_point(size=3, position=pd) +
  geom_errorbar(aes(ymin=plo, ymax=phi, color = group), position=pd) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  labs(y='Fit (rigidness)', x = 'Task', color = 'Group', fill = 'Group')

ggsave('mr6.png', path = 'plots_1/', w=5,h=4)

mm <- mr5 # choose model to plot
c3m <- r3
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  scale_x_continuous(breaks = pretty_breaks()) +
  facet_wrap(~task_perf)+
  labs(y='Fit (rigidness)', x = 'NASA score', color = 'Group', fill = 'Group') +
  theme(legend.position=c(.85,.15), legend.key.size = unit(.5, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))

ggsave('mr5.png', path = 'plots_1/', w=5.5,h=4)


r3 %>% 
  group_by(group, task_perf) %>% 
  get_summary_stats(ntr) %>% 
  ggplot(.,aes(y=mean, x=task_perf, group=group, colour=group)) +
  geom_point(size = 3, position = pd) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), position = pd, width = 0.1) +
  geom_text(aes(label = round(mean,1)), show.legend = FALSE, position=pdtxt) 

ggplot(r3,aes(y=mean, x=task_perf, group=group, colour=group)) +
  geom_point(size = 3, position = pd) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), position = pd, width = 0.1) +
  geom_text(aes(label = round(mean,1)), show.legend = FALSE, position=pdtxt) 


emmip(mr1, group ~  task_perf) +
  geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4) 

emmip(mr5, group ~ nasa_bht | task_perf, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

emmip(mr6, group ~ nasa_bht | task_perf, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

# Contrasts rigidness -----------------------------------------------------

# Simple
contrast(emmeans(mr1, ~ group * task_perf),
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% # anova/contrast coding (-+0.5), compares btw levels 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE,-z.ratio) %>% kable()

# Complex
summary(mr5)
contrast(emtrends(mr5, ~ group * task_perf, var = 'nasas'),
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% # anova/contrast coding (-+0.5), compares btw levels 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE,-z.ratio) %>% kable

test(emtrends(mr5, ~ group * task_perf, var = 'nasas'), adjust = "fdr") %>% # test slopes against 0
  as_tibble() %>% select(-df) %>% kable()

emmeans(mr5, pairwise ~ group | nasas | task_perf, at = list(nasas = c(-1, 0, 1)))$contrasts %>% 
  as_tibble %>% select(-df, -SE, -z.ratio) %>% mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>% kable  # compare effect of group 


# Modelling NASA BHT --------------------------------------------------------

n1 <- c3 %>% 
  group_by(id) %>% slice(1) %>% 
  select(id, group, nasas, nasa_bht) %>% 
  ungroup()

ggplot(n1, aes(x = nasas, fill = group)) +
  geom_histogram(binwidth = 0.25) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  labs(y = 'Count', x ='NASA score', fill = 'Group')

ggplot(n1, aes(y = nasas, x = group)) +
  geom_violin(fill = 'grey90') + 
  geom_boxplot(width=0.1) + 
  stat_summary(fun.y=mean, geom="point", color = 'dark red', size = 3) +
  # geom_jitter(position=position_jitter(0.2), alpha = .2) +
  scale_x_discrete(labels = c('no control', 'full control')) +
  labs(y = 'NASA', x = 'Group')

mn1 <- lm(nasas ~ group, n1)
summary(mn1)

n1 %>% t_test(nasas ~ group) %>% 
  left_join(n1 %>% cohens_d(nasas ~ group))
 

ggplot(mn1, aes(x=group, y=nasas, fill = group))+
  geom_flat_violin(position = position_nudge(x = .15, y = 0),adjust = 1, trim = T, colour = NA)+
  geom_point(position = position_jitter(width = .1), size = .7, alpha = .3)+
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = .1, colour = "black") +
  geom_rangeframe(x=c(10000)) + 
  scale_x_discrete(labels = c('no control', 'full control')) +
  guides(fill = "none", colour = "none")   +
  coord_cartesian(xlim = c(1.2, NA)) +
  coord_flip() 

ggplot(mn1, aes(x = group, y = nasas, fill = group, colour = group)) + 
  stat_halfeye(adjust = .5, width = .6, justification = -.2, .width = 0, point_colour = NA) + # add half-violin from {ggdist} package
  geom_boxplot(width = .13, outlier.color = NA, alpha = .5) +
  stat_dots(side = "left", justification = 1.1, binwidth = .1, colour = NA, 
            dotsize = .7, stackratio = .8, alpha = .4) + # add dot plots from {ggdist} package
  # geom_point(position = position_jitter(width = .1), size = 1, alpha = .2) +
  scale_x_discrete(labels = c('no control', 'full control')) +
  coord_cartesian(xlim = c(1, NA)) + # remove white space on the left
  guides(fill = "none", colour = 'none') + labs(x = 'Group', y = 'NASA')

ggsave('nasaXgroup.png', path = 'plots_1/', w=4.5,h=4)

# Modelling RT -------------------------------------------------------------- 

# Prepare
c3rt <- c3 %>% 
  filter(logrt != -Inf) %>% 
  left_join(c3rt %>%
              group_by(group, task_perf, block) %>%
              get_summary_stats(logrt) %>% 
              mutate(lo = mean-3*sd,
                     hi = mean+3*sd) %>% 
              select(group:task_perf, mean,sd,lo,hi)) %>% 
  filter(logrt > lo & logrt < hi) %>% 
  mutate(logrtC = center(logrt))

1-nrow(c3rt)/nrow(c3) # 1.5% removed

c3rt %>% 
  ggplot(aes(x=logrt)) + geom_histogram()

c3rt %>%
  group_by(group, task_perf, block) %>% 
  get_summary_stats(logrt) %>% 
  mutate(lo = mean-3*sd,
         hi = mean+3*sd) %>% 
  select(group:task_perf, mean,sd,lo,hi)

c3rt %>%
  get_summary_stats(logrt, logrtC)


# Simple
mrt1 <- lmer(logrt ~ group * nasas + (1 | id),
             c3rt)
summary(mrt1)


# Complex
mrt2 <- lmer(logrt ~ group * nasas * block * task_perf + (1 | id),
             c3rt) 
summary(mrt2)
anova(mrt1,mrt2)

tab_model(mrt2, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          # file='mrt2.html',
          pred.labels=c('(Intercept)','Group','NASA','Block','Task',
                        'Group x NASA','Group x Block','NASA x Block','Group x Task','NASA x Task','Block x Task',
                        'Group x NASA x Block','Group x NASA x Task','Group x Block x Task','NASA x Block x Task',
                        'Group x NASA x Block x Task'),
          dv.labels='')

mrt2b <- lmer(logrt ~ group * nasas * block * task_perf + (block | id),
             c3rt) 
summary(mrt2b)
anova(mrt2b,mrt2)

mrt2c <- lmer(logrt ~ group * nasas * block * task_perf + (task_perf | id),
              c3rt) 
summary(mrt2c)
anova(mrt2c,mrt2b)

mrt2d <- lmer(logrt ~ group * nasas * block * task_perf + (task_perf | id),
              c3rt) 
summary(mrt2d)
anova(mrt2d,mrt2b)

# Plotting RT -------------------------------------------------------------


mm <- mrt2 # choose model to plot
c3m <- c3rt
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
ggplot(c3m, aes(y = exp(fit), x = nasas, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin=exp(plo), ymax=exp(phi), color = group, fill = group), alpha = 0.3, colour = NA) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  facet_grid(block~task_perf) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(y='Fit (RT)', x = 'NASA score', color = 'Group', fill = 'Group')

ggsave('mr2.png', path = 'plots_1/')


# Modelling accuracy ------------------------------------------------------


c3 %>% 
  group_by(task_perf, group) %>% 
  get_summary_stats(corr)

ma1 <- glmer(corr ~ group * nasas + (1 | id),
            c3 %>% filter(logrt != -Inf),
            family = binomial)
summary(ma1)


ma2 <- glmer(corr ~ group * nasas * block * task_perf + (1 | id),
             c3 %>% filter(logrt != -Inf),
             family = binomial,
             glmerControl(optimizer = 'bobyqa')) #
summary(ma2)
anova(ma1,ma2)

tab_model(ma2, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, 
          file='ma2.html',
          pred.labels=c('(Intercept)','Group','NASA','Block','Task',
                        'Group x NASA','Group x Block','NASA x Block','Group x Task','NASA x Task','Block x Task',
                        'Group x NASA x Block','Group x NASA x Task','Group x Block x Task','NASA x Block x Task',
                        'Group x NASA x Block x Task'),
          dv.labels='')

ma2b <- glmer(corr ~ group * nasas * block * task_perf + (block | id),
             c3 %>% filter(logrt != -Inf),
             family = binomial,
             glmerControl(optimizer = 'bobyqa')) #
summary(ma2b)
anova(ma2b,ma2)

# Plotting accuracy -------------------------------------------------------

mm <- ma2 # choose model to plot
c3m <- c3 %>% filter(logrt != -Inf)
c3m$fit <- model.matrix(mm) %*% fixef (mm)
pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
c3m$plo <- c3m$fit-1.96*sqrt(pvar) # update data with CI
c3m$phi <- c3m$fit+1.96*sqrt(pvar)
ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
  scale_colour_discrete(labels = c('no control', 'full control')) +
  facet_grid(block~task_perf) +
  scale_fill_discrete(labels = c('no control', 'full control')) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(y='Fit (accuracy)', x = 'NASA score', color = 'Group', fill = 'Group') +
  theme(legend.position=c(.15,.12), legend.key.size = unit(.4, 'cm'), 
        legend.background = element_rect(fill = "white", color = "grey85"),
        legend.title = element_text(size=13))


ggsave('ma2.png', path = 'plots_1/')


# Switching in mixed block ------------------------------------------------

# more strict subset of free trials in mixed block: >n+1 after switch from instructed trials
c4 <- c3 %>% 
  mutate(switchMix = case_when(lag(trial_type) == 'forced' ~ NA_real_,
                               TRUE ~ switch))

distinct(c4, switchMix)  

table(c4$switchMix, useNA = "always")  
table(c4$switch, useNA = "always")  
count(c4, switchMix)

mm12 <- glmer(switchMix ~ group * nasas * block * task_perf + (1 | id),
             c4, family = binomial,
             glmerControl(optimizer="bobyqa"))  
summary(mm12)

mm12d <- glmer(switchMix ~ group * nasas * block * task_perf + (block | id),
              c4, family = binomial,
              glmerControl(optimizer="bobyqa"))  
summary(mm12d)
anova(mm12,mm12d)
