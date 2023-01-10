# PdB1 VTS pilot experiment analysis ###########
# Bartosz Majchrowicz, majchrowicz.b@gmail.com #

# Load libraries and settings  ###########
{
  # load libraries, use install.packages('package_name') if library not yet installed
  library(tidyverse)
  library(rstatix)
  library(ggplot2)
  library(readbulk)
  library(lme4)
  library(lmerTest)
  library(emmeans)

  # some settings for convenience
  options(scipen=999, width = 150)
  Sys.setenv(LANG = "en")
  theme_set(theme_bw(base_size = 14) +
              theme(panel.grid.minor = element_blank()))
}

# Load and transform data ###########
p0 <- read_bulk(directory = "C:/Users/barto/Psych/Badania/PdB Exp1/PilotVTS/Exp analysis pilot/data_pilot_2",
                subdirectories = F, verbose = F, fun = read.csv) %>% as_tibble() # load multiple data files 

p1 <- na_if(p0, "") %>% # replace blank spaces with NAs
  rename(id = participant) %>% 
  mutate(countBal = case_when(id %% 4 == 1 ~ 1,  # code counterbalance
                              id %% 4 == 2 ~ 2,
                              id %% 4 == 3 ~ 3,
                              id %% 4 == 0 ~ 4,
                              TRUE ~ NA_real_)) %>% 
  rename(type = Type, stimuli = Stimuli, location = Location, shape = Shape) %>% 
  mutate(type = tolower(type)) %>% 
  select(c('id','countBal', 'frameRate', 'type':'correct_shape_right',  # select only needed columns, others are removed
           'LoopCb1.thisRepN':'LoopCb1.thisIndex',
           'nMainCb1.thisRepN':'nMainCb1.thisIndex',
           'LoopCb2.thisRepN':'LoopCb2.thisIndex',
           'nMainCb2.thisRepN':'nMainCb2.thisIndex',
           'both_shape_cb1.keys', 'both_shape_cb1.rt', 'both_shape_cb1.corr', 
           'both_loc_cb1.keys', 'both_loc_cb1.rt', 'both_loc_cb1.corr', 
           'both_shape_cb2.keys', 'both_shape_cb2.rt', 'both_shape_cb2.corr', 
           'both_loc_cb2.keys', 'both_loc_cb2.rt', 'both_loc_cb2.corr')) %>% 
  mutate(id = replace(id, id == 70842, 108),
         id = replace(id, id == 133966, 111))

table(p1$id) # check ids
  
p1s <- p1 %>% # remove practice trials
  filter_at(vars('LoopCb1.thisRepN', 'LoopCb1.thisTrialN', 'LoopCb1.thisN', 'LoopCb1.thisIndex',
                 'LoopCb2.thisRepN', 'LoopCb2.thisTrialN', 'LoopCb2.thisN', 'LoopCb2.thisIndex',
                 'both_shape_cb1.corr', 'both_loc_cb1.corr', 'both_shape_cb2.corr', 'both_loc_cb2.corr'),
            any_vars(!is.na(.)))


p2 <- p1s %>% # fill in block number data based on values from specified columns
  group_by(id) %>% 
  fill(c('LoopCb1.thisRepN', 'LoopCb2.thisRepN'), .direction = "up") 

p3 <- p2 %>%  # remove rows with no data in specified columns
  filter_at(vars('type','stimuli','location', 'shape'), all_vars(!is.na(.)))
  

p4 <- p3 %>%  # code/collapse block number (so that it's coded in a single column, not split to separate columns based on cb)
  mutate(LoopCb1.thisRepN = as.numeric(LoopCb1.thisRepN),
         LoopCb2.thisRepN = as.numeric(LoopCb2.thisRepN)) %>% 
  mutate(blockNr = case_when(is.na(LoopCb1.thisRepN) ~ LoopCb2.thisRepN,
                             is.na(LoopCb2.thisRepN) ~ LoopCb1.thisRepN,
                             TRUE ~ NA_real_))

p4 %>% group_by(id) %>% count() %>% print(n=Inf) # check nr of trials per id (63 trial x 4 blocks = 252)
 
p5 <- p4 %>% 
  select(-c('correct_loc_left':'correct_shape_right', 'stimuli')) # remove useless columns

p6 <- p5 %>% # code/collapse data across tasks and cb's
  mutate(task = case_when(is.na(both_loc_cb1.keys)   & is.na(both_loc_cb2.keys) ~ 'shape', # code task in a single column
                          is.na(both_shape_cb1.keys) & is.na(both_shape_cb2.keys) ~ 'loc',
                          TRUE ~  'unknownTask')) %>% 
  mutate(resp = case_when(task == 'shape' & is.na(both_shape_cb1.keys) ~ both_shape_cb2.keys, # code response 
                          task == 'shape' & is.na(both_shape_cb2.keys) ~ both_shape_cb1.keys,
                          task == 'loc'   & is.na(both_loc_cb1.keys) ~ both_loc_cb2.keys,
                          task == 'loc'   & is.na(both_loc_cb2.keys) ~ both_loc_cb1.keys,
                          TRUE ~  'unknownResp')) %>% 
  mutate(corr = case_when(task == 'shape' & is.na(both_shape_cb1.keys) ~ as.character(both_shape_cb2.corr), # code correct
                          task == 'shape' & is.na(both_shape_cb2.keys) ~ as.character(both_shape_cb1.corr),
                          task == 'loc'   & is.na(both_loc_cb1.keys) ~ as.character(both_loc_cb2.corr),
                          task == 'loc'   & is.na(both_loc_cb2.keys) ~ as.character(both_loc_cb1.corr),
                          TRUE ~  'unknownCorr')) %>% 
  mutate(rt   = case_when(task == 'shape' & is.na(both_shape_cb1.keys) ~ as.numeric(both_shape_cb2.rt), # code rt
                          task == 'shape' & is.na(both_shape_cb2.keys) ~ as.numeric(both_shape_cb1.rt),
                          task == 'loc'   & is.na(both_loc_cb1.keys) ~ as.numeric(both_loc_cb2.rt),
                          task == 'loc'   & is.na(both_loc_cb2.keys) ~ as.numeric(both_loc_cb1.rt),
                          TRUE ~  NA_real_))

p7 <- p6 %>% # code task switch vs stay
  mutate(switch_type = case_when(task == 'shape' & lag(task) == 'shape' ~ 'stay_shape', # code switch type in a single column
                                 task == 'loc' & lag(task) == 'loc' ~ 'stay_loc',
                                 task == 'shape' & lag(task) == 'loc' ~ 'switch_shape',
                                 task == 'loc' & lag(task) == 'shape' ~ 'switch_loc',
                                 TRUE ~ 'unknownSwitch')) %>%  
  mutate(switch = case_when(switch_type == 'stay_shape' | switch_type == 'stay_loc' ~ 0, # code switch (logical)
                            switch_type == 'switch_shape' | switch_type == 'switch_loc' ~ 1,
                            TRUE ~ NA_real_)) %>% 
  filter(corr != 999) # renove unknown corr

# Voluntary switch rates  ###########

s1 <- p7 %>% 
  group_by(id, blockNr) %>% 
  filter(!is.na(switch)) %>%
  summarise(switch_sum = sum(switch)) %>% # get nr of all switches per id and block 
  left_join(p7 %>% # join/merge with another df in which we get nr of all observations
              group_by(id, blockNr) %>% 
              count(switch_data = !is.na(switch)) %>% # get nr of all non-NA observations about switch (both switch & stay)
              filter(switch_data == TRUE) %>% 
              select(-switch_data)) %>% 
  mutate(vsr = round(switch_sum/n, 3)) # get VSR (nr of switches / nr of observations) 



s2 <- s1 %>%  # add VSR based on switch difficulty, per id and block
  left_join(p7 %>% 
              group_by(id, blockNr) %>% 
              filter(switch_type == 'switch_shape') %>% 
              count(name = 'switch_diff') %>% # count difficult switches (location -> shape)
              left_join(p7 %>% 
                          group_by(id, blockNr) %>% 
                          filter(switch_type == 'switch_loc') %>% 
                          count(name = 'switch_easy'))) %>%  # count easy switches (shape -> location)
  relocate(c('switch_easy', 'switch_diff'), .before = 'n') %>% # rearrange for convenience
  mutate(vsr_easy = switch_easy/n, # get VSR for easy and diff switches
         vsr_diff = switch_diff/n)

s3 <- s2 %>% 
  left_join(s2 %>% # add VSR per id only (averaged across blocks)
              ungroup() %>% 
              group_by(id) %>% 
              summarise(vsr_id = mean(vsr),
                        vsr_easy_id = mean(vsr_easy),
                        vsr_diff_id = mean(vsr_diff))) 

# VSR by diff only
s4 <- s3 %>% # reshape and summarise for plot 1
  pivot_longer(cols = c('vsr_easy', 'vsr_diff'), names_prefix = 'vsr_', names_to = 'difficulty', # reshape (to long format)
               values_to = 'vsr_difficulty') %>% 
  group_by(difficulty) %>% 
  get_summary_stats(vsr_difficulty) %>% # summarise
  mutate(difficulty = as.factor(difficulty))

ggplot(s4, aes(y = mean, x = difficulty)) + # plot VSR by difficulty
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1) +
  # scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(limits = rev, labels = c('Easy', 'Difficult')) +
  labs(y = 'VSR (%)', x = 'Difficulty of switch')

# VSR by id, diff and block
s5 <- s3 %>% # reshape and summarise for plot 2
  pivot_longer(cols = c('vsr_easy', 'vsr_diff'), names_prefix = 'vsr_', names_to = 'difficulty', # reshape (to long format)
               values_to = 'vsr_difficulty') %>% 
  group_by(id, blockNr, difficulty) %>% 
  get_summary_stats(vsr_difficulty) %>% # summarise
  mutate(blockNr = as.factor(blockNr))

ggplot(s5, aes(y = mean, x = blockNr, group = difficulty, colour = difficulty)) + # plot VSR across blocks
  geom_point(size = 5, position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.1)) +
  facet_wrap(~id, labeller = labeller(id = label_both)) +
  # scale_y_continuous(limits = c(0,1)) +
  scale_colour_discrete(limits = rev, labels = c('Easy', 'Difficult')) +
  labs(y = 'VSR (%)', x = 'Block number', colour = 'Difficulty of switch')

# VSR by diff and block
s6 <- s3 %>% # reshape and summarise for plot 3
  pivot_longer(cols = c('vsr_easy', 'vsr_diff'), names_prefix = 'vsr_', names_to = 'difficulty', # reshape (to long format)
               values_to = 'vsr_difficulty') %>% 
  group_by(blockNr, difficulty) %>% 
  get_summary_stats(vsr_difficulty) %>% # summarise
  mutate(blockNr = as.factor(blockNr),
         difficulty = as.factor(difficulty))

ggplot(s6, aes(y = mean, x = blockNr, group = difficulty, colour = difficulty)) + # plot VSR across blocks
  geom_point(size = 5, position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.1)) +
  # scale_y_continuous(limits = c(0,1)) +
  scale_colour_discrete(limits = rev, labels = c('Easy', 'Difficult')) +
  labs(y = 'VSR (%)', x = 'Block number', colour = 'Difficulty of switch')
 

# Stay vs switch by id

z1 <- p7 %>% 
  filter(!is.na(switch)) %>% 
  group_by(id) %>% 
  count(switch_type) %>% 
  pivot_wider(names_from = switch_type, values_from = n) %>% 
  mutate_at(2:5, ~replace_na(.,0)) %>% 
  pivot_longer(cols = c(2:5), names_to = 'switch_type', values_to = 'n')

ggplot(z1, aes(y = n, x = switch_type, fill = switch_type, colour = switch_type)) +
  geom_point(size= 2) +
  facet_wrap(~id) +
  geom_text(aes(label = n), nudge_y = 30) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
  
  
# Easy vs difficult decision (% of easy task selection)

p7 %>% # overall proportion of easy task selection
  filter(task != 'unknownTask') %>% 
  group_by(task) %>%
  count() %>% 
  pivot_wider(names_from = 'task', values_from = 'n') %>% 
  mutate(switch_prop = loc/(loc+shape))

p7 %>% # plot by block and id
  filter(task != 'unknownTask') %>% 
  group_by(id, blockNr, task) %>%
  # group_by(blockNr, task) %>% 
  count() %>% 
  pivot_wider(names_from = 'task', values_from = 'n') %>% 
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  mutate(easy_sel = loc/(loc+shape)) %>% 
  ggplot(., aes(x = blockNr, y = easy_sel)) +
  geom_point(size = 2) +
  facet_wrap(~id) +
  labs(y = 'Selection of easy task (%)', x = 'Block')

p7 %>% # plot by block only
  filter(task != 'unknownTask') %>% 
  group_by(blockNr, task) %>%
  count() %>% 
  pivot_wider(names_from = 'task', values_from = 'n') %>% 
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  mutate(easy_sel = loc/(loc+shape)) %>% 
  ggplot(., aes(x = blockNr, y = easy_sel)) +
  geom_point(size = 2) +
  labs(y = 'Selection of easy task (%)', x = 'Block')

# Accuracy ##########


p7 %>% # averaged across ids
  group_by(id, switch_type) %>% 
  filter(switch_type != 'unknownSwitch') %>%
  mutate(corr = as.integer(corr)) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% 
              filter(switch_type != 'unknownSwitch') %>% 
              group_by(id, switch_type) %>% 
              count()) %>% 
  mutate(acc = sumCorr/n) %>% 
  group_by(switch_type) %>% 
  get_summary_stats(acc) %>% 
  ggplot(., aes(y = mean, x = switch_type)) +
  geom_point(size = 5, position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.1)) +
  labs(y = 'Accuracy (%)', x = 'Trial type') 

p7 %>% # per id
  group_by(id, switch_type) %>% 
  filter(switch_type != 'unknownSwitch') %>%
  mutate(corr = as.integer(corr)) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% 
              filter(switch_type != 'unknownSwitch') %>% 
              group_by(id, switch_type) %>% 
              count()) %>% 
  mutate(acc = sumCorr/n) %>% 
  group_by(id, switch_type) %>% 
  get_summary_stats(acc) %>% 
  ggplot(., aes(y = mean, x = switch_type)) +
  geom_point(size = 5, position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.1)) +
  labs(y = 'Accuracy (%)', x = 'Trial type') +
  facet_wrap(~id)

# RTs ############

p7 %>% 
  group_by(switch_type) %>% 
  filter(switch_type != 'unknownSwitch') %>% 
  get_summary_stats(rt) %>% 
  ggplot(., aes(y = mean, x = switch_type)) +
  geom_point(size = 5, position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.1)) +
  labs(y = 'RT', x = 'Trial type') 


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

