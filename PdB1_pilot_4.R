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

# Data wrangling ###########
p0 <- read_bulk(directory = "C:/Users/barto/Psych/Badania/PdB Exp1/PilotVTS/ForcedC/Cb1-4",
                subdirectories = F, verbose = F, fun = read.csv) %>% # load data from multiple files 
  as_tibble() %>%
  select(File, everything()) %>% 
  separate(File, into = c('x', 'id', 'z'), sep = c(2,3)) %>% 
  select(-c(x,z)) %>% mutate(id = as.integer(id)) 

p1 <- p0 %>% 
  mutate_if(is.character, list(~na_if(.,""))) %>%  # replace blank spaces with NAs
  # rename(id = participant) %>% 
  mutate(countBal = case_when(id %% 4 == 1 ~ 1,  # code counterbalance
                              id %% 4 == 2 ~ 2,
                              id %% 4 == 3 ~ 3,
                              id %% 4 == 0 ~ 4,
                              TRUE ~ NA_real_)) %>% 
  rename(type = Type, stimuli = Stimuli, location = Location, shape = Shape, task = Task, cue = Cue, block = Block, cuefol = cue_follow) %>% 
  mutate(type = tolower(type)) %>% 
  select(c('id','countBal', 'frameRate',   # select only needed columns, others are removed
           'type':'shape', 'task', 'cue','cuefol', 'block',
           'correct_loc_left':'correct_shape_right', 
           
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
           'both_loc_cb3_mix.keys',   'both_loc_cb3_mix.rt',   'both_loc_cb3_mix.corr', 
           'both_shape_cb3_free.keys', 'both_shape_cb3_free.rt', 'both_shape_cb3_free.corr', 
           'both_loc_cb3_free.keys',   'both_loc_cb3_free.rt',   'both_loc_cb3_free.corr', 
           
           'both_shape_cb4_mix.keys', 'both_shape_cb4_mix.rt', 'both_shape_cb4_mix.corr', 
           'both_loc_cb4_mix.keys',   'both_loc_cb4_mix.rt',   'both_loc_cb4_mix.corr', 
           'both_shape_cb4_free.keys', 'both_shape_cb4_free.rt', 'both_shape_cb4_free.corr', 
           'both_loc_cb4_free.keys',   'both_loc_cb4_free.rt',   'both_loc_cb4_free.corr', 
  )) 

colnames(p1)
table(p1$id) # check ids
  
p2 <- p1 %>% # remove practice trials
  filter_at(vars(
    'both_shape_cb1_mix.corr', 'both_loc_cb1_mix.corr', 'both_shape_cb2_mix.corr', 'both_loc_cb2_mix.corr',
    'both_shape_cb3_mix.corr', 'both_loc_cb3_mix.corr', 'both_shape_cb4_mix.corr', 'both_loc_cb4_mix.corr',
    'both_shape_cb1_free.corr', 'both_loc_cb1_free.corr', 'both_shape_cb2_free.corr', 'both_loc_cb2_free.corr',
    'both_shape_cb3_free.corr', 'both_loc_cb3_free.corr', 'both_shape_cb4_free.corr', 'both_loc_cb4_free.corr',
  ),
  any_vars(!is.na(.))) # remove rows with NAs (we infer which rows are practice based on presence of NAs in above-specified non-practice columns)

p2 %>% 
  group_by(id) %>% count()

p3 <- p2 %>% # just for pre-pilot run
  mutate(block = ifelse(is.na(block), 'mixed', block))

p4 <- p3 %>%  # collapse (unite) block number (so that it's coded in a single column, not split to separate columns based on cb)
  rowwise() %>% 
  mutate(blockNr = sum(c_across(c(
    'nBlocksMixCb1.thisRepN', 'nBlocksMixCb2.thisRepN', 'nBlocksMixCb3.thisRepN', 'nBlocksMixCb4.thisRepN',
    'nBlocksFreeCb1.thisRepN', 'nBlocksFreeCb2.thisRepN', 'nBlocksFreeCb3.thisRepN', 'nBlocksFreeCb4.thisRepN',
  )), na.rm = T)) %>% 
  relocate(blockNr, .after = 'block') # relocate for convenience
  
p4 %>% group_by(id) %>% count() %>% print(n=Inf) # check nr of trials per id
 
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
  mutate(task_cued = case_when(cue == 'JAKI' ~ 'shape', # actually just renamed 'cue' column
                               cue == 'GDZIE'~ 'loc',
                               cue == 'XXXX' ~ 'any',
                               TRUE ~ 'unknownTask')) %>% 
  relocate(c(task_cued, task_perf), .after = 'block') # relocate for convenience

colnames(p6)

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
  ungroup()

# Post-wrangling checks and descriptives

colnames(p7)
p7 %>% group_by(id) %>% count()
p7 %>% group_by(id, block) %>% count() %>% pivot_wider(names_from = block, values_from = n)

p7 %>% group_by(id) %>% get_summary_stats(rt)
p7 %>% freq_table(resp)
p7 %>% freq_table(corr)
p7 %>% freq_table(task_perf)
p7 %>% freq_table(task_cued)


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
  scale_x_discrete(limits = rev, labels = c('Easy', 'Difficult')) +
  labs(y = 'VSR (%)', x = 'Difficulty of switch')

ggsave('vsrDiff.png', path = 'plots_pilot_4/')


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

ggsave('vsrDiffxBlock.png', path = 'plots_pilot_4/')

 
# Switch vs stay (% of switches) ##########

z1 <- p7 %>% # number of switch types
  filter(!is.na(switch)) %>% 
  group_by(id) %>% 
  count(switch_type) %>% 
  pivot_wider(names_from = switch_type, values_from = n) %>% 
  mutate_at(2:5, ~replace_na(.,0)) %>% 
  pivot_longer(cols = c(2:5), names_to = 'switch_type', values_to = 'n')

ggplot(z1, aes(y = n, x = switch_type, fill = switch_type, colour = switch_type)) + # switch types by id
  geom_point(size= 2) +
  facet_wrap(~id) +
  geom_text(aes(label = n), nudge_y = 30) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))

ggsave('switchTypexId.png', path = 'plots_pilot_4/', width = 12, height = 11)

p7 %>% # overall proportion of switches
  filter(!is.na(switch)) %>% 
  group_by(switch) %>%
  count() %>% 
  pivot_wider(names_from = 'switch', values_from = 'n') %>% 
  rename(stay = 1, switch = 2) %>% 
  mutate(switch_prop = switch/(switch+stay),
         stay_prop   = stay  /(switch+stay))

p7 %>% # overall proportion of switch types
  filter(switch_type != 'unknownSwitch') %>% 
  group_by(switch_type) %>%
  count() %>% 
  pivot_wider(names_from = 'switch_type', values_from = 'n') %>% 
  mutate(switch_diff_prop = switch_shape/(stay_loc+stay_shape+switch_loc+switch_shape),
         switch_easy_prop =   switch_loc/(stay_loc+stay_shape+switch_loc+switch_shape),
         stay_diff_prop = stay_shape/(stay_loc+stay_shape+switch_loc+switch_shape),
         stay_easy_prop =   stay_loc/(stay_loc+stay_shape+switch_loc+switch_shape))

p7 %>% # plot proportion of switches per id
  filter(!is.na(switch)) %>% 
  group_by(switch, id) %>%
  count() %>% 
  pivot_wider(names_from = 'switch', values_from = 'n') %>% 
  rename(stay = 2, switch = 3) %>% 
  mutate(switch_prop = switch/(switch+stay),
         stay_prop   = stay  /(switch+stay)) %>% 
  pivot_longer(cols = c('switch_prop', 'stay_prop'), names_to = 'type', values_to = 'prop') %>% 
  ggplot(., aes(y = prop, x = type)) +
  geom_point(size = 2) +
  facet_wrap(~id, ncol = 5) +
  geom_text(aes(label = round(prop,2)), nudge_y = 0.1) +
  scale_x_discrete(labels=c('stay', 'switch')) +
  labs(y = 'Proportion of stay/switch', x = 'Type')
  
ggsave('switchPropxId.png', path = 'plots_pilot_4/', width = 9, height = 8)

p7 %>% # plot proportion of switches (averaged across id's)
  filter(!is.na(switch)) %>% 
  group_by(switch) %>%
  count() %>% 
  pivot_wider(names_from = 'switch', values_from = 'n') %>% 
  rename(stay = 1, switch = 2) %>% 
  mutate(switch_prop = switch/(switch+stay),
         stay_prop   = stay  /(switch+stay)) %>% 
  pivot_longer(cols = c('switch_prop', 'stay_prop'), names_to = 'type', values_to = 'prop') %>% 
  ggplot(., aes(y = prop, x = type)) +
  geom_point(size = 2) +
  geom_text(aes(label = round(prop,2)), nudge_y = 0.03) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(labels=c('stay', 'switch')) +
  labs(y = 'Proportion of stay/switch', x = 'Type')

ggsave('switchProp.png', path = 'plots_pilot_4/', width = 5, height = 5)

# Repetition bias #########
z2 <- p7 %>% 
  filter(switch == 0) %>% 
  group_by(task_perf) %>% 
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  bind_cols(p7 %>% 
              filter(switch == 1) %>% 
              ungroup() %>%
              count() %>% rename(switches=1)) %>% 
  mutate(loc_rep = loc/(loc+shape+switches),
         shape_rep=shape/(loc+shape+switches),
         switch = switches/(loc+shape+switches),
         checkSum=(loc_rep+shape_rep+switch)) %>% 
  pivot_longer(cols = 4:6, names_to = 'type', values_to = 'prop')


ggplot(z2, aes(y=prop, x=type)) +
  geom_point(size=3) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(labels = c('Easy repetition', 'Difficult repetition', 'Switch')) +
  labs(y='Proporion', x='Type') +
  geom_text(aes(label = round(prop,2)), nudge_y = 0.05) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6))

ggsave('repBias.png', path = 'plots_pilot_4/', width = 4, height = 4)

  
  
# Easy vs difficult decision (% of easy task selection) ##########

p7 %>% # overall proportion of easy task selection
  filter(task_perf != 'unknownTask') %>% 
  group_by(task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate(easy_prop = loc/(loc+shape),
         diff_prop = shape/(loc+shape))

p7 %>% # plot by block and id
  filter(task_perf != 'unknownTask') %>% 
  group_by(id, blockNr, task_perf) %>%
  # group_by(blockNr, task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  mutate(easy_sel = loc/(loc+shape)) %>% 
  ggplot(., aes(x = blockNr, y = easy_sel)) +
  geom_point(size = 2) +
  facet_wrap(~id) +
  labs(y = 'Selection of easy task (%)', x = 'Block')

p7 %>% # plot by block only
  filter(task_perf != 'unknownTask') %>% 
  group_by(blockNr, task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  mutate(easy_sel = loc/(loc+shape)) %>% 
  ggplot(., aes(x = blockNr, y = easy_sel)) +
  geom_point(size = 3) +
  geom_text(aes(label = round(easy_sel,2)), nudge_y = 0.04) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = 'Selection of easy task (%)', x = 'Block')

ggsave('easySelxBlock.png', path = 'plots_pilot_4/')

p7 %>% # plot averaged
  filter(task_perf != 'unknownTask') %>% 
  group_by(task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  mutate(easy_sel = loc/(loc+shape),
         diff_sel = shape/(loc+shape)) %>% 
  pivot_longer(cols = c('easy_sel', 'diff_sel'), names_to = 'task_perf', values_to = 'sel') %>% 
  ggplot(., aes(x = task_perf, y = sel)) +
  geom_point(size = 3) +
  geom_text(aes(label = round(sel,3)), nudge_y = 0.04) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(labels=c('Difficult', 'Easy')) +
  labs(y = 'Proportion (%)', x = 'Task')

ggsave('taskSel.png', path = 'plots_pilot_4/', width = 5, height = 5)


# Stiffness ########

p8 <- p7 %>% 
  filter(!is.na(switch)) %>% 
  group_by(id, blockNr) %>% 
  mutate(sumSwitch = cumsum(as.double(switch)))

x2 <- p8 %>% 
  group_by(id, blockNr, sumSwitch, task_perf) %>% 
  count()

x2 %>% 
  ungroup() %>% group_by(task_perf) %>% 
  get_summary_stats(n) %>% 
  select(task_perf, n, mean, se) %>% 
  ggplot(., aes(y = mean, x= task_perf)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1) +
  scale_x_discrete(labels = c('Easy', 'Difficult')) +
  labs(y = 'Average single task run (nr of trials)', x = 'Task') 
 
ggsave('runLengthxDiff.png', path = 'plots_pilot_4/', width = 5, height = 5)

x2 %>% 
  ungroup() %>% group_by(id, task_perf) %>% 
  get_summary_stats(n) %>% 
  select(id, task_perf, n, mean, se) %>% 
  ggplot(., aes(y = mean, x= task_perf)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1) +
  facet_wrap(~id) +
  scale_x_discrete(labels = c('Easy', 'Difficult')) +
  labs(y = 'Average single task run (nr of trials)', x = 'Task') 


# Accuracy ##########

p7 %>% # stay vs switch
  group_by(id, switch) %>% 
  filter(!is.na(switch)) %>%
  mutate(corr = as.integer(corr)) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% 
              group_by(id, switch) %>% 
              filter(!is.na(switch)) %>% 
              count()) %>% 
  mutate(acc = sumCorr/n) %>% 
  group_by(switch) %>% 
  get_summary_stats(acc)

p7 %>% # easy vs difficult task
  group_by(id, task_perf) %>% 
  filter(task_perf != 'unknownTask') %>%
  mutate(corr = as.integer(corr)) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% 
              group_by(id, task_perf) %>% 
              filter(task_perf != 'unknownTask') %>% 
              count()) %>% 
  mutate(acc = sumCorr/n) %>% 
  group_by(task_perf) %>% 
  get_summary_stats(acc)

p7 %>% # plot switch type, averaged across ids
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
  ggplot(., aes(y = mean, x = factor(switch_type, level=c('stay_loc','switch_loc','stay_shape','switch_shape')))) +
  geom_point(size = 5, position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.1)) +
  geom_text(aes(label = round(mean,3)), nudge_y = 0.003, nudge_x = 0.3) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6)) +
  labs(y = 'Accuracy (%)', x = 'Trial type')

ggsave('accuracy.png', path = 'plots_pilot_4/', width = 5, height = 5)

p7 %>% # plot switch type, per id
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
  ggplot(., aes(y = mean, x = factor(switch_type, level=c('stay_loc','switch_loc','stay_shape','switch_shape')))) +
  geom_point(size = 5, position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.1)) +
  geom_text(aes(label = round(mean,3)), nudge_y = 0.01, nudge_x = 0.3) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6)) +
  labs(y = 'RT (s)', x = 'Trial type') 

ggsave('rt.png', path = 'plots_pilot_4/', width = 5, height = 5)


# Checks ########### %>% 

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

