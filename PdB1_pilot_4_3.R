# PdB1 VTS pilot experiment analysis # 
# Bartosz Majchrowicz, majchrowicz.b@gmail.com #


# Load libraries and settings ---------------------------------------------

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
  '%nin%' <- Negate("%in%")
  
  library(default)
  default(get_summary_stats) <- list(type = 'common') # get_summary_stats <- reset_default(get_summary_stats)

}

# Data wrangling ###########
p0 <- read_bulk(directory = "C:/Users/barto/Psych/Badania/PdB Exp1/PilotVTS/Exp analysis pilot/data_pilot_4/",
                subdirectories = F, extension = 'csv', verbose = F, fun = read.csv) %>% # load data from multiple files 
  as_tibble()

p1 <- p0 %>% 
  mutate_if(is.character, list(~na_if(.,""))) %>%  # replace blank spaces with NAs
  rename(id = participant) %>%
  mutate(countBal = case_when(id %% 4 == 1 ~ 1,  # code counterbalance
                              id %% 4 == 2 ~ 2,
                              id %% 4 == 3 ~ 3,
                              id %% 4 == 0 ~ 4,
                              TRUE ~ NA_real_)) %>% 
  rename(type = Type, stimuli = Stimuli, location = Location, shape = Shape, task = Task, cue = Cue, block = Block, cuefol = cue_follow) %>% 
  mutate(type = tolower(type)) %>% 
  select(c('id','countBal',   # select only needed columns, others are removed
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
  relocate(blockNr, .after = block)  # relocate for convenience
  
p3 %>% group_by(id) %>% count() %>% print(n=Inf) # check nr of trials per id
p3 %>% group_by(id) %>% count(blockNr)  %>% print(n=Inf) # check nr of trials per id and block

p3 %>% # check nr of counterbalances
  # filter(id >1000) %>%
  group_by(id) %>% filter(row_number()==1) %>% ungroup() %>% freq_table(countBal)

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

nrow(p6) == nrow(p5a) && nrow(p6) == nrow(p5b) # before merging make sure nr or rows is identical (should be TRUE)

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
  filter(corr < 999) %>% # remove unknown corr
  filter(!nchar(resp) > 1) %>% # remove rare cases of two responses registered
  mutate(trialNr = nTrials.thisN + 1, .after = blockNr) %>% # recode trial nr per block
  select(-c('cue', 'nBlocks.thisRepN':'nTrials.thisIndex')) %>% # again remove some columns
  rename(trial_type = task) %>% 
  mutate(cBal = case_when( # simplified counterbalance (blocks order only)
    countBal == 1 | countBal == 3 ~ 'mix-free', # mix > free 
    countBal == 2 | countBal == 4 ~ 'free-mix', # free > mix
    TRUE ~ NA), .after = countBal)

# Post-wrangling checks and descriptives

colnames(p7)
p0 %>% distinct(frameRate, .keep_all = TRUE) %>% select(participant, frameRate) %>% arrange(frameRate) %>% 
  ggplot(., aes(x=frameRate)) + 
  geom_histogram(color = 'grey20') +
  geom_vline(xintercept = c(60,75,144), linetype = 'dashed') # most common refresh rates are 60Hz, 75Hz, 144Hz, and 240Hz

p7 %>% group_by(id) %>% count()
p7 %>% group_by(id, block) %>% count() %>% pivot_wider(names_from = block, values_from = n) %>% print(n=Inf)

p7 %>% group_by(id) %>% count(blockNr)
p7 %>% count(block)

p7 %>% group_by(id) %>% group_by(id) %>% group_by(id) %>% get_summary_stats(rt)
p7 %>% freq_table(resp)
p7 %>% freq_table(corr)
p7 %>% freq_table(task_perf)
p7 %>% freq_table(block, task_perf)
p7 %>% freq_table(task_cued)


# Durations ---------------------------------------------------------------


b1 <- p0 %>% # breaks time
  # rename(id = participant) %>% 
  select(starts_with('timebreak')) %>% 
  filter_all(any_vars(!is.na(.))) %>% 
  unite('timeBreak', everything(), sep='', na.rm=T) %>% 
  mutate(timeBreak = as.numeric(timeBreak))
  
b1 %>% get_summary_stats() 

ggplot(b1, aes(x=timeBreak)) + 
  geom_histogram(color = 'grey20', binwidth = 2) +
  geom_vline(aes(xintercept=median(timeBreak)), linetype = 'dashed') +
  annotate('text', label = paste('Med = ', b1 %>% get_summary_stats() %>% select(median) %>% round(1), '(s)\n',
                                 'Mean = ', b1 %>% get_summary_stats() %>% select(mean) %>% round(1), '(s)'), 
           x=18, y=14.5, size = 5) +
  scale_x_continuous(breaks=seq(0,55,5)) +
  labs(y = 'Count', x = 'Break time (s)')

ggsave('breakDur.png', path = 'plots_pilot_4/', w=5,h=4)


b2 <- p0 %>% select(globalClockTime) %>% # overall time
  filter(!is.na(globalClockTime)) %>% 
  mutate(globalClockTime = globalClockTime/60) 

b2 %>% get_summary_stats()

ggplot(b2, aes(x = globalClockTime)) +
  geom_histogram(color = 'grey20', binwidth = 0.1) +
  geom_vline(aes(xintercept=median(globalClockTime)), linetype = 'dashed') +
  annotate('text', label = paste('Med = ', b2 %>% get_summary_stats() %>% select(median) %>% round(1), '(min)\n',
                                 'Mean = ', b2 %>% get_summary_stats() %>% select(mean) %>% round(1), '(min)'),
           x=21.7, y=1.9, size = 5) +
  # geom_text(aes(x = 22, y = 2.2, label = paste('M = ', b2 %>% get_summary_stats() %>% select(mean) %>% round(1), '(min)')), size = 5) +
  labs(y = 'Count', x = 'Experiment time (min)')

ggsave('expDur.png', path = 'plots_pilot_4/', w=5,h=4)



# Outliers ----------------------------------------------------------------

p7 %>% 
  group_by(id, block, switch_type) %>%
  identify_outliers(rt) %>% 
  filter(is.extreme == T)

p7 %>% 
  group_by(block) %>%
  mahalanobis_distance() 

p7 %>% 
  doo(~mahalanobis_distance(.)) %>% 
  filter(is.outlier == T)

p7 %>% # accuracy  
  group_by(id, block, switch_type) %>% 
  filter(task_perf != 'unknownTask') %>%
  mutate(corr = as.integer(corr)) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% 
              group_by(id, block, task_perf) %>% 
              filter(task_perf != 'unknownTask') %>% 
              count()) %>% 
  mutate(acc = sumCorr/n) %>% 
  # group_by(id, block, switch_type) %>% 
  identify_outliers(acc) %>% 
  filter(is.outlier == T)


# VSR ---------------------------------------------------------------------

# Overall VSR per id

# # new version which drops 0's
# s1 <- p7 %>% 
#   filter(trial_type == 'free') %>% # discard forced trials
#   filter(!is.na(switch)) %>% 
#   freq_table(id, block, switch) %>%
#   filter(switch == 1) %>% # get prop = VSR (nr of switches / nr of observations) 
#   rename(vsr = prop) %>% select(-n)
# 
# s1 %>% 
#   group_by(block) %>% 
#   get_summary_stats(vsr) %>% 
#   ggplot(., aes(y = mean, x = block)) +
#   geom_point(size=5) +
#   geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1) +
#   geom_text(aes(label = round(mean,2)), nudge_x = 0.2, show.legend = FALSE) +
#   labs(y = 'VSR (%)', x = 'Block')
# 
# ggsave('vsr_block.png', path = 'plots_pilot_4/', w=4,h=4)

# p7 %>% # proportion of VSR easy to VSR difficult ||| 1
#   filter(trial_type == 'free') %>% # discard forced trials
#   filter(switch == 1) %>%
#   freq_table(id, block, switch_type)

# p7 %>% # vsr per difficulty per id | proportion of switch to stay per block and task ||| 2
#   filter(trial_type == 'free') %>% # discard forced trials
#   filter(!is.na(switch)) %>% 
#   freq_table(id, block, task_perf, switch)
# 
# s2 <- p7 %>% # vsr per difficulty per id | proportion of switch to stay per block and task ||| 2
#   filter(trial_type == 'free') %>% # discard forced trials
#   filter(!is.na(switch)) %>% 
#   freq_table(id, block, task_perf, switch) %>% 
#   filter(switch == 1) %>% # get prop = VSR (nr of switches / nr of observations)  | VSR per block and task
#   rename(vsr = prop) %>% select(!c('n', 'switch'))
# 
# s2 %>% 
#   group_by(block, task_perf) %>% 
#   get_summary_stats(vsr) %>% 
#   ggplot(., aes(y = mean, x = block, fill = task_perf, colour = task_perf)) +
#   geom_point(position = position_dodge(0.1)) +
#   geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.1)) +
#   geom_text(aes(label = round(mean,1)), nudge_y = 0.001, nudge_x = 0.15, show.legend = FALSE) +
#   labs(y = 'VSR', colour = 'VSR type', x = 'Block') + guides(fill="none")
# 
# ggplot(s2, aes(y = vsr, x = block, fill = task_perf, colour = task_perf)) +
#   geom_point(position = position_dodge(0.1)) +
#   facet_wrap(~id) +
#   geom_text(aes(label = vsr), nudge_y = 0.001, nudge_x = 0.15, show.legend = FALSE) +
#   labs(y = 'VSR', colour = 'VSR type', x = 'Block') + guides(fill="none")

# Old version
# Overall VSR per id
s1 <- p7 %>% 
  group_by(id, block) %>% 
  filter(!is.na(switch)) %>%
  filter(trial_type == 'free') %>% # discard forced trials
  summarise(switch_sum = sum(switch)) %>% # get nr of all switches per id and block 
  left_join(p7 %>% # join/merge with another df in which we get nr of all observations
              group_by(id, block) %>% 
              filter(trial_type == 'free') %>% # discard forced trials
              filter(!is.na(switch)) %>%
              count()) %>% 
  mutate(vsr = round(switch_sum/n, 3)) # get VSR (nr of switches / nr of observations) 

s1 %>% # plot
  group_by(block) %>%
  get_summary_stats(vsr) %>%
  ggplot(., aes(y = mean, x = block)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1) +
  geom_text(aes(label = round(mean,2)), nudge_x = 0.2, show.legend = FALSE) +
  labs(y = 'VSR (%)', x = 'Block')

ggsave('vsr_block.png', path = 'plots_pilot_4/', w=4,h=4)

# add VSR based on switch difficulty, per id (VSR decomposed into easy and difficult)
s2 <- s1 %>% 
  left_join(p7 %>% 
              group_by(id, block) %>% 
              filter(trial_type == 'free') %>% # discard forced trials
              filter(switch_type == 'switch_shape') %>% 
              count(name = 'switch_diff') %>% # count difficult switches (location -> shape)
              left_join(p7 %>% 
                          group_by(id, block) %>%  
                          filter(trial_type == 'free') %>%
                          filter(switch_type == 'switch_loc') %>% 
                          count(name = 'switch_easy'))) %>%  # count easy switches (shape -> location)
  relocate(c('switch_easy', 'switch_diff'), .before = 'n') %>% # rearrange for convenience
  mutate(vsr_easy = switch_easy/n, # get VSR for easy and diff switches
         vsr_diff = switch_diff/n) %>% 
  replace(is.na(.), 0)

s2 %>%  # plot
  group_by(block) %>% 
  get_summary_stats(vsr_easy, vsr_diff) %>% 
  ggplot(., aes(y = mean, x = block, fill = variable, colour = variable)) +
  geom_point(position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position = position_dodge(0.3)) +
  geom_text(aes(label = round(mean,3)), position = position_dodge(1), show.legend = FALSE) +
  labs(y = 'VSR', colour = 'VSR component', x = 'Block') + 
  scale_colour_discrete(labels = c('Easy', 'Difficult')) +
  guides(fill="none")

ggsave('vsr_blockXdiff.png', path = 'plots_pilot_4/', w=6,h=4)


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

ggsave('vsr_diffXid.png', path = 'plots_pilot_4/', w=12.5,h=10)


# Counterbalance and VSR

s1 %>% 
  left_join(p7 %>% group_by(id) %>% filter(row_number()==1) %>% select(id, cBal), 
            by = 'id') %>% 
  group_by(cBal, block) %>%
  get_summary_stats(vsr) %>%
  ggplot(., aes(y = mean, x = block)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1) +
  geom_text(aes(label = round(mean,2)), nudge_x = 0.25, show.legend = FALSE) +
  facet_wrap(~cBal) +
  labs(y = 'VSR (%)', x = 'Block')

ggsave('vsr_blockXcb.png', path = 'plots_pilot_4/', w=6,h=4)

 
# Proportion of switches --------------------------------------------------

y1 <- p7 %>% # overall proportion of switch
  filter(!is.na(switch)) %>% 
  filter(trial_type == 'free') %>%
  group_by(id, block, switch) %>%
  count() %>% 
  pivot_wider(names_from = 'switch', values_from = 'n') %>% 
  rename(stay = '0', switch = '1') %>% 
  replace(is.na(.), 0) %>%
  mutate(switch_prop = switch/(switch+stay),
         stay_prop   = stay  /(switch+stay)) %>% 
  group_by(block) %>% 
  get_summary_stats(switch_prop, stay_prop)

y1 %>% # plot
  rename(switch = variable) %>% 
  mutate(switch = as.factor(switch)) %>%
  ggplot(., aes(y = mean, x = block, colour = switch))+
  geom_point(size = 2, position = position_dodge(.3)) +
  geom_text(aes(label = round(mean,2)), show.legend = F, position = position_dodge(1)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=position_dodge(.3)) +
  scale_colour_discrete(labels=c('Switch', 'Stay')) +
  labs(y = 'Proportion of stay/switch', x = 'Block', colour = 'Type')

ggsave('switchProp_block.png', path = 'plots_pilot_4/', width = 5, height = 4)

y2 <- p7 %>% # proportion of switch_type
  filter(switch_type != 'unknownSwitch') %>% 
  filter(trial_type == 'free') %>%
  group_by(id, block, switch_type) %>%
  count() %>% 
  pivot_wider(names_from = 'switch_type', values_from = 'n') %>% 
  replace(is.na(.), 0) %>% 
  mutate(switch_shape_p = switch_shape/(stay_loc+stay_shape+switch_loc+switch_shape),
         switch_loc_p = switch_loc/(stay_loc+stay_shape+switch_loc+switch_shape),
         stay_shape_p = stay_shape/(stay_loc+stay_shape+switch_loc+switch_shape),
         stay_loc_p = stay_loc/(stay_loc+stay_shape+switch_loc+switch_shape)) %>% 
  group_by(block) %>% 
  get_summary_stats(c('switch_shape_p':'stay_loc_p'))

y2 %>% # plot
  mutate(switch_type = as.factor(variable), .after = 'block', .keep = 'unused') %>%
  ggplot(., aes(y = mean, x = block, colour = switch_type))+
  geom_point(size = 2, position = position_dodge(.3)) +
  # geom_text(aes(label = round(mean,3)), show.legend = F, nudge_x = 0.2, nudge_y = 0.01) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=position_dodge(.3)) +
  scale_colour_discrete(labels=c('Switch to shape', 'Switch to loc', 'Stay with shape', 'Stay with loc')) +
  labs(y = 'Proportion of stay/switch', x = 'Block', colour = 'Type')


ggsave('switchProp_blockXswitchType.png', path = 'plots_pilot_4/', width = 6, height = 4)


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

ggsave('repBiasXblock.png', path = 'plots_pilot_4/', width = 6, height = 4)

  
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

ggsave('repBiasXblock2.png', path = 'plots_pilot_4/', width = 5.5, height = 4)


# Rigidness ########

p8 <- p7 %>% 
  filter(!is.na(switch)) %>% 
  filter(trial_type == 'free') %>% 
  group_by(id, block, blockNr) %>% 
  mutate(sumSwitch = cumsum(as.double(switch))) # cumulative sum of switches

x2 <- p8 %>% 
  group_by(id, block, blockNr, sumSwitch, task_perf) %>% # count trials per each switch (and other grouping factors)
  count() %>% 
  ungroup()

x2 %>% # single task run length per block
  # filter(id %nin% c(5,9,100,106,498686,941944)) %>%
  group_by(block, task_perf) %>% 
  get_summary_stats(n) %>% 
  select(block, task_perf, n, mean, se) %>% 
  ggplot(., aes(y = mean, x = block, colour = task_perf)) +
  geom_point(size=3, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=position_dodge(0.5)) +
  geom_text(aes(label = round(mean,1)), show.legend = FALSE, position = position_dodge(1.3)) +
  scale_colour_discrete(labels = c('Easy', 'Difficult')) +
  labs(y = 'Average single task run (nr of trials)', x = 'Block', colour = 'Task') 

ggsave('runLengthxDiff.png', path = 'plots_pilot_4/', width = 5, height = 4)

x2 %>%  # single task run length per block, per id | exlucde 0's
  filter(block == 'freeOnly') %>% 
  group_by(id, task_perf) %>% 
  get_summary_stats(n) %>% 
  mutate(across(everything(), ~replace(., . == 0, NA))) %>% # replace 0 with NA's
  select(id, task_perf, n, mean, se) %>% 
  ggplot(., aes(y = mean, x= task_perf)) +
  geom_point(size=3) +
  # geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.2) +
  facet_wrap(~id) +
  geom_text(aes(label = round(mean,2)), show.legend = FALSE, nudge_x = 0.3) +
  scale_x_discrete(labels = c('Easy', 'Difficult')) +
  labs(y = 'Average single task run (nr of trials)', x = 'Task') 


x2 %>%  # single task run length per block, per id | include 0's
  group_by(id, task_perf) %>% filter(row_number()==1) %>% select(id, task_perf) %>% # create frame 
  left_join(x2 %>%
              filter(block == 'freeOnly') %>% 
              group_by(id, task_perf) %>% 
              get_summary_stats(n)) %>% 
  mutate(n = case_when(is.na(n) ~ 0, TRUE ~ n),
         mean = case_when(is.na(mean) ~ 0, TRUE ~ mean)) %>% 
  ggplot(., aes(y = mean, x= task_perf)) +
  geom_point(size=3) +
  facet_wrap(~id) +
  geom_text(aes(label = round(mean,2)), show.legend = FALSE, nudge_x = 0.3) +
  scale_x_discrete(labels = c('Easy', 'Difficult')) +
  labs(y = 'Average single task run (nr of trials)', x = 'Task') 

ggsave('runLengthXdiffXid.png', path = 'plots_pilot_4/', width = 10, height = 8)



x2 %>% # counterbalance single task run length per block
  left_join(p7 %>% group_by(id) %>% filter(row_number()==1) %>% select(id, cBal), 
            by = 'id') %>% 
  # filter(id %nin% c(5,9,100,106,498686,941944)) %>%
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

ggsave('runLengthXdiffXcb.png', path = 'plots_pilot_4/', width = 6.5, height = 4)


# Task selection ##########

p7 %>% # overall proportion of diff/easy task selection
  filter(task_perf != 'unknownTask') %>% 
  filter(trial_type == 'free') %>% # remove instructed trials
  group_by(task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate(easy_prop = loc/(loc+shape),
         diff_prop = shape/(loc+shape))

p7 %>% # plot by block only
  filter(task_perf != 'unknownTask') %>% 
  filter(trial_type == 'free') %>% # remove instructed trials
  group_by(id, block, task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  mutate(diff_sel = shape/(loc+shape)) %>% # measure of difficult selection
  group_by(block) %>% # average across id
  get_summary_stats(diff_sel) %>% 
  ggplot(., aes(x = block, y = mean)) +
  geom_point(size = 3, position=position_dodge(0.5)) +
  geom_text(aes(label = round(mean,2)), nudge_x = 0.2) +
  geom_errorbar(aes(ymin = mean-2*se, ymax = mean+2*se), width = 0.1, position=position_dodge(0.5)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = 'Selection of difficult task (%)', x = 'Block')

ggsave('taskSelXblock.png', path = 'plots_pilot_4/', w=4,h=4)

p7 %>% # plot by block and id
  filter(task_perf != 'unknownTask') %>% 
  filter(trial_type == 'free') %>% # remove instructed trials
  group_by(id, block, task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  mutate(diff_sel = shape/(loc+shape)) %>% 
  ggplot(., aes(x = block, y = diff_sel)) +
  geom_point(size = 2) +
  facet_wrap(~id) +
  labs(y = 'Selection of difficult task (%)', x = 'Block')


p7 %>% # plot averaged, a bit redundant compared to previous plot
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

ggsave('taskSelXblock2.png', path = 'plots_pilot_4/', width = 6, height = 4) # a bit redundant 

p7 %>% # counterbalance effect on overall proportion of diff/easy task selection
  filter(task_perf != 'unknownTask') %>% 
  filter(trial_type == 'free') %>% # remove instructed trials
  group_by(cBal, task_perf) %>%
  count() %>% 
  pivot_wider(names_from = 'task_perf', values_from = 'n') %>% 
  mutate(easy_prop = loc/(loc+shape),
         diff_prop = shape/(loc+shape))

p7 %>% # counterbalance plot by block only
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

ggsave('taskSelXblockXcb.png', path = 'plots_pilot_4/', width = 6, height = 4)  


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
  labs(y = 'Accuracy (%)', x = 'Trial type')

ggsave('accuracy.png', path = 'plots_pilot_4/', width = 6, height = 5)

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
  labs(y = 'Accuracy (%)', x = 'Trial type') +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6)) +
  facet_wrap(~id)

ggsave('accuracyXid.png', path = 'plots_pilot_4/', width = 13, height = 11)


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

ggsave('rt.png', path = 'plots_pilot_4/', width = 6, height = 5)


ggplot(p7, aes(x = rt)) +
  geom_histogram(color = 'grey20') +
  labs(y = 'Count', x = 'RT (s)')

ggplot(p7, aes(rt)) +
  geom_histogram(aes(y = ..density..), color = "grey20", fill = 'grey50', binwidth = 0.1) +
  geom_density(color = 'black', fill = "lightgreen", alpha = 0.3) +
  # geom_density(color = NA, fill = "grey50") +
  # geom_histogram(aes(y = ..density..), fill = '#F85700', binwidth = 0.1, alpha = 0.4) +
  geom_vline(aes(xintercept=median(rt)), linetype = 'dashed') +
  annotate('text', label = paste('Med = ', p7 %>% get_summary_stats(rt) %>% select(median), '(s)\n',
                                 'Mean = ', p7 %>% get_summary_stats(rt) %>% select(mean), '(s)'),
           x=1.5, y=1.5, size = 4) +
  labs(y = 'Density', x = 'RT (s)')

ggsave('rt_hist.png', path = 'plots_pilot_4/', width = 6, height = 5)


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

