# PdB1 analysis: simplified pipeline  -------------------------------------
# Bartosz Majchrowicz, majchrowicz.b@gmail.com

# Libraries, settings, functions ------------------------------------------

{
  # libraries 
  library(tidyverse)
  library(easystats)
  library(rstatix)
  library(ggplot2)
  library(ggdist)
  library(cowplot)
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
  library(afex)
  library(sjPlot)
  library(janitor)
  library(beepr)
  library(default)
  library(ggeffects)

  # settings
  options(scipen=999, width = 150)
  Sys.setenv(LANG = "en")
  theme_set(theme_bw(base_size = 14) +
              theme(panel.grid.minor = element_blank()))
  default(get_summary_stats) <- list(type = 'common') # get_summary_stats <- reset_default(get_summary_stats)
  default(kable) <- list(digits = 3, format = 'simple')
  pd = position_dodge(0.2)
  pdtxt = position_dodge(1)
  
  # functions
  isConv <- function (model) { # helper function to check convergence/singular fit
    if (!inherits(model, "merMod")) stop("Error: must pass a lmerMod object")
    msg <- NULL
    if(is.null(unlist(model@optinfo$conv$lme4))) {
      msg = 'Converged normally'
    } else {
      if (isSingular(model)) {
        msg = 'Singular fit'
      } else {
        msg = 'Not converged'
      }
    }
    return(msg)
  }
  format_labels <- function(model){ # formatted labels for tab_model
    names(fixef(model)) %>% str_replace_all(c("[:digit:]"="",':'=' x ')) %>% 
      str_to_title %>% str_replace_all(c('Nasas'='NASA', 'Task_perf'='Task', 'X'='x')) %>% str_remove_all('_iv.-')
  }
  format_formula <- function(model){ # formatted formula for tab_model
    deparse(formula(model), width.cutoff = 200L) %>% str_to_title %>% 
      str_replace_all(c('Id'='ID', '_perf'='', '_dv'='', '\\*'='x', 'Nasas'='NASA', 'Logrt'='LogRT')) %>% str_remove_all('_iv.-')
  }
}  

# Data load -------------------------------------------------------------- 

save(c3, p7, q2, t3, # main datasets
     c3s, c3rt, c3srt, r3, # data subsets/summaries 
     pl_m1, pl_m12j, pl_m22a, pl_ms3d_v1, pl_ms3d_v2, pl_mt3a, pl_ma2c, mrt2c, pl_mr5, # plots
     file='PdB1_Data.RData') 
load('PdB1_Data.RData') # includes datasets and plots

glimpse(p7) # full long-format data, before cleaning and transformations
glimpse(c3) # full long-format data, after cleaning and transformations, used for main modelling and subsetting
glimpse(q2) # questionnaire data
glimpse(t3) # durations data

{
  contrasts(c3$group) # 0 no control, 1 full control
  contrasts(c3$cBal)
  contrasts(c3$block)
  contrasts(c3$task_perf)
  contrasts(c3$switch_iv) # 0 repeat, 1 switch
} # factors are anova/contrast coded (-0.5,0.5) (ref. https://debruine.github.io/faux/articles/contrasts)


# Data cleaning  ---------------------------------------------------------- 

# Control item reversed answers (7 participants)
q3 <- q2 %>% 
  select(-proc) %>%
  pivot_longer(cols = c('mental':'motiv'), names_to = 'q', values_to = 'ans')

insOutIDs1 <- q3 %>% 
  filter(q == 'controlover' & task == 'bht') %>% # no control group but answered full control
  filter(group == 0 & ans == 6) %>% 
  select(id, group, ans)  
insOutIDs2 <- q3 %>% 
  filter(q == 'controlover' & task == 'bht') %>% # full control group but answered no control
  filter(group == 1 & ans == 1) %>% 
  select(id, group, ans)  

# Accuracy <65% (4 participants)
accOutIDs <- p7 %>%  
  group_by(id) %>% 
  summarise(sumCorr = sum(corr, na.rm = T)) %>% 
  left_join(p7 %>% group_by(id) %>% count) %>% 
  mutate(acc = sumCorr/n) %>% ungroup() %>% 
  filter(acc<.65) %>% select(id,acc) 

# Duration at & between the tasks (22 participants)
t3out <- t3 %>% 
  mutate(across(postBHT:pureVTS, 
                list(mean = mean, sd = sd))) %>%  # named list of funs
  mutate(across(postBHT:pureVTS,  # calculate 3SD threshold value
                ~3*sd(.),
                .names = "{.col}_thr")) %>% 
  mutate(out_postBHT      = case_when(postBHT < postBHT_mean - postBHT_thr | postBHT > postBHT_mean + postBHT_thr ~ T, TRUE ~ F),      # time at BHT 
         out_betweenTasks = case_when(betweenTasks < betweenTasks_mean - betweenTasks_thr | betweenTasks > betweenTasks_mean + betweenTasks_thr ~ T, TRUE ~ F), # BHT-VTS
         out_pureVTS      = case_when(pureVTS < pureVTS_mean - pureVTS_thr | pureVTS > pureVTS_mean + pureVTS_thr ~ T, TRUE ~ F),      # time at VTS
         out_postVTS      = case_when(postVTS < postVTS_mean - postVTS_thr | postVTS > postVTS_mean + postVTS_thr ~ T, TRUE ~ F)) %>%  # total experiment time
  select(id:pureVTS, starts_with('out'))

durOutIDs4c <- t3out %>% 
  filter(if_any(c(starts_with('out')), ~.==T)) %>% # select outliers based on 4 duration criteria (BHT, VTS, between the tasks, and total (but not VTS break times))
  select(id, (starts_with('out')))

# Cleaning before looking at switches
c1 <- p7 %>% 
  arrange(id) %>% 
  filter(id %nin% 1002) %>% # age 61
  filter(id %nin% c(insOutIDs1$id, insOutIDs2$id)) %>% # sense of control reversed answer
  filter(id %nin% durOutIDs4c$id) %>% # prolonged durations  
  filter(id %nin% accOutIDs$id) # accuracy <65% per subject 

# No switches (additional 11 participants)
vsrOutIDs <- c1 %>% 
  group_by(id) %>% 
  filter(!is.na(switch)) %>%
  filter(trial_type == 'free') %>% # discard forced trials
  summarise(switch_sum = sum(switch)) %>% # get nr of all switches per id and block 
  left_join(c1 %>% # join/merge with another df in which we get nr of all observations
              group_by(id) %>% 
              filter(trial_type == 'free') %>% # discard forced trials
              filter(!is.na(switch)) %>%
              count()) %>% 
  mutate(vsr = round(switch_sum/n, 3)) %>% # get VSR (nr of switches / nr of observations) 
  filter(switch_sum == 0) # pick those with 0 

# Cleaning by 0 switches
c2 <- c1 %>% 
  filter(id %nin% vsrOutIDs$id)  
  
nrow(c2) == nrow(c3) # c3 has some additional data transformations/wrangling, but dataset size is the same


# Demographics ------------------------------------------------------------

c3 %>% get_summary_stats(age) %>% kable() # age



# Manipulation checks  ----------------------------------------------------

mc1 <- c3 %>%  # BHT questionnaire data for manipulation checks
  group_by(id) %>% slice(1) %>% select(id,group) %>% 
  left_join(q2 %>% 
              filter(task=='bht') %>% 
              select(id, group, controlover, extentofcontrol, effort1, effort2, motiv, mental) %>% 
              mutate(group = as.factor(group))) %>% ungroup 


mc1 %>% t_test(controlover ~ group) %>% 
  left_join(mc1 %>% cohens_d(controlover ~ group))

mc1 %>% t_test(extentofcontrol ~ group) %>% 
  left_join(mc1 %>% cohens_d(extentofcontrol ~ group))

mc1 %>% t_test(effort1 ~ group) %>% 
  left_join(mc1 %>% cohens_d(effort1 ~ group))  

mc1 %>% t_test(effort2 ~ group) %>% 
  left_join(mc1 %>% cohens_d(effort2 ~ group)) 

mc1 %>% t_test(motiv ~ group) %>% 
  left_join(mc1 %>% cohens_d(motiv ~ group)) 


# NASA ~ Group ------------------------------------------------------------

n1 <- c3 %>% 
  group_by(id) %>% slice(1) %>% 
  select(id, group, nasas, nasa_bht) %>% 
  ungroup()

n1 %>% t_test(nasas ~ group) %>% 
  left_join(n1 %>% cohens_d(nasas ~ group))

# Plot 
{ggplot(n1, aes(x = group, y = nasas, fill = group, colour = group)) + 
    stat_halfeye(adjust = .5, width = .6, justification = -.2, .width = 0, point_colour = NA) + # half-violin from ggdist package
    geom_boxplot(width = .13, outlier.color = NA, alpha = .5) +
    stat_dots(side = "left", justification = 1.1, binwidth = .1, colour = NA, 
              dotsize = .7, stackratio = .8, alpha = .4) + # dot plots from ggdist package
    scale_x_discrete(labels = c('no control', 'full control')) +
    coord_cartesian(xlim = c(1, NA)) + # remove white space on the left
    guides(fill = "none", colour = 'none') + labs(x = 'Group', y = 'NASA')}


# Experiences of control and difficulty -----------------------------------

mc2 <- c3 %>% 
  group_by(id) %>% slice(1) %>% 
  select(id, group) %>% ungroup %>% 
  left_join(q2 %>% 
              # filter(task=='bht') %>% 
              select(id, task, group, controlover, extentofcontrol, effort1, effort2, motiv, mental) %>% 
              mutate(group = as.factor(group))) %>% ungroup 

mc2 %>% group_by(group, task) %>%
  get_summary_stats(controlover) 

mc2 %>% group_by(group, task) %>%
  get_summary_stats(mental) 

# Switching: simple model -------------------------------------------------

m1 <- glmer(switch ~ group * nasas + (1 | id),
            data = c3, 
            family = binomial)
summary(m1)
tab_model(m1, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='m1.html',
          pred.labels = format_labels(m1), title = format_formula(m1))
# model's fit assessed with anova() not improved with group, nasas, group+nasas random term

# Contrasts
contrast(emtrends(m1, ~ group, var = 'nasas'),  # comparison of NASA slopes between the levels
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>% 
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df) %>% kable() 

emmeans(m1, pairwise ~ group | nasas, at = list(nasas = c(-1, 0, 1)))$contrasts %>%   # compare effect of group on 3 points of NASA
  as_tibble %>% select(-df)  %>% mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>% kable

test(emtrends(m1, ~ group, var = 'nasas'), adjust = "fdr") %>% # test NASA slopes against 0
  as_tibble() %>% select(-df) %>% kable()

# Plots
{
  # whole model (ref. https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)
  mm <- m1 # choose model to plot
  pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only; requires lots of memory! 
  c3m <- c3 %>% anti_join(c3 %>% filter(if_any(c(switch, group, nasas, task_perf, block), ~ is.na(.)))) %>% 
    mutate(fit = model.matrix(mm) %*% fixef (mm), 
           plo = fit-1.96*sqrt(pvar), phi = fit+1.96*sqrt(pvar))
  pl_m1 <- ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
    scale_colour_discrete(labels = c('no control', 'full control')) +
    scale_fill_discrete(labels = c('no control', 'full control')) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(y='Fit (switching)', x = 'NASA', color = 'Group', fill = 'Group') +
    theme(legend.position=c(.85,.87), legend.key.size = unit(.5, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('m1.png', pl_m1, path='plots_2/',w=5,h=5)
}

ggpredict(model = m1, terms = c("nasas","group")) %>% 
  plot(alpha=0.3, line.size=1)

ggpredict(model = m1, terms = c("nasas","group")) %>% 
  # as_tibble() %>% 
  mutate(predicted = log(predicted)) 

emmip(m1, group ~ nasas, cov.reduce = range, plotit=F) 

emmip(m1, group ~ nasas, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

# Switching: complex model ------------------------------------------------

# Compare models
m12d <- glmer(switch ~ group * nasas * block * task_perf + (block | id),
              data = c3, family = binomial,
              glmerControl(optimizer="bobyqa")) # bobyqa optimizer converges most often; can also try optCtrl=list(maxfun=1e5) 

m12i <- glmer(switch ~ group * nasas * block * task_perf + (block + task_perf | id),
              data = c3, family = binomial,
              glmerControl(optimizer="bobyqa"))  

m12j <- glmer(switch ~ group * nasas * block * task_perf + (block * task_perf | id),
              data = c3, family = binomial,
              glmerControl(optimizer="bobyqa"))  

anova(m12d,m12i,m12j) # m12j > m12i > m12d

anova(m12i, m12j)

# Winning model
summary(m12j)
tab_model(m12j, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='m12j.html',
          pred.labels = format_labels(m12j), title = format_formula(m12j))

# Contrasts

emmeans(m12j, pairwise ~ group | block * task_perf)$contrasts %>%  # effect of group in interaction
  as_tibble %>% select(-df, -SE) %>% adjust_pvalue(p.col = 'p.value', output.col = 'adj.p', method = 'fdr') %>% kable  

contrast(emtrends(m12j, ~ group * task_perf * block, var = 'nasas'), # pairwise comparison of NASA slopes between the levels, all factors
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>%  
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE)  %>% kable  

contrast(emtrends(m12j, ~ group, var = 'nasas'), # pairwise comparison of NASA slopes between the levels, average factors
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>%  
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE)  %>% kable  

emmeans(m12j, pairwise ~ group | nasas | task_perf | block, at = list(nasas = c(-1, 0, 1)))$contrasts %>%  # effect of group on 3 points of NASA, all factors
  as_tibble %>% select(-df, -SE) %>% adjust_pvalue(p.col = 'p.value', output.col = 'adj.p', method = 'fdr') %>% kable  

emmeans(m12j, pairwise ~ group | nasas, at = list(nasas = c(-1, 0, 1)))$contrasts %>%  # effect of group on 3 points of NASA, average factors
  as_tibble %>% select(-df, -SE) %>% adjust_pvalue(p.col = 'p.value', output.col = 'adj.p', method = 'fdr') %>% kable  

test(emtrends(m12j, ~ group * block * task_perf, var = 'nasas'), adjust = "fdr") %>% # test NASA slopes against 0, all factors
  as_tibble() %>% select(-df) %>% kable()

test(emtrends(m12j, ~ group, var = 'nasas'), adjust = "fdr") %>% # test NASA slopes against 0, average factors
  as_tibble() %>% select(-df) %>% kable()

# Plots  
{
  # whole model (ref. https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)
  mm <- m12j # choose model to plot
  pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only; requires lots of memory! 
  c3m <- c3 %>% anti_join(c3 %>% filter(if_any(c(switch, group, nasas, task_perf, block), ~ is.na(.)))) %>% 
    mutate(fit = model.matrix(mm) %*% fixef (mm), 
           plo = fit-1.96*sqrt(pvar), phi = fit+1.96*sqrt(pvar))
  facet_task_block <- c('loc' = 'Location task', 'shape' = 'Shape task', 'freeOnly' = 'Free block', 'mixed' = 'Mixed block')
  pl_m12j <- ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
    scale_colour_discrete(labels = c('no control', 'full control')) +
    facet_grid(block~task_perf, labeller = as_labeller(facet_task_block)) +
    scale_fill_discrete(labels = c('no control', 'full control')) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(y='Fit (switching)', x = 'NASA', color = 'Group', fill = 'Group') +
    theme(legend.position=c(.85,.15), legend.key.size = unit(.5, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('m12j.png', pl_m12j, path = 'plots_2/', w=5,h=4.5)
  
  # quick interactions 
  emmip(m12j, block ~ task_perf) +
    geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4) 
  emmip(m12j, block ~ nasas, cov.reduce = range) +
    geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = block), alpha = 0.3, colour = NA) 
  emmip(m12j, group ~ task_perf | block) +
    geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4) 
  
  emmip(m12j, group ~ task_perf  | block, plotit=F)  %>% 
    ggplot(., aes(x = xvar, y = yvar, colour = tvar)) +
    geom_point(size = 3, position=pd) +
    geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4, position=pd)  +
    facet_grid(~block, labeller = as_labeller(facet_block)) +
    scale_colour_discrete(labels = c('no control', 'full control')) +
    scale_x_discrete(labels = c('Location', 'Shape')) +
    guides(fill = 'none') +
    labs(y = 'Linear prediction (switching)', x = 'Task', colour = 'Group') +
    theme(legend.position=c(.85,.18), legend.key.size = unit(.5, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('m12j_x.png', path = 'plots_2/', width = 5, height = 3.5)
}

ggpredict(m12j, c("nasas","group",'task_perf')) %>% 
  plot(alpha=0.3, line.size=1)

ggpredict(m12j, c("nasas", "group")) %>% 
  plot(alpha=0.3, line.size=1)

emmip(m12j, group ~ nasas, cov.reduce = range) +
  geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = group), alpha = 0.3, colour = NA) 

# Switching: counterbalance -----------------------------------------------

# Compare models
m20a <- glmer(switch ~ group * block * cBal + (block | id),
              c3, family = binomial,
              glmerControl(optimizer="bobyqa"))

m22a <- glmer(switch ~ group * block * task_perf * cBal + (block * task_perf| id),
              c3, family = binomial,
              glmerControl(optimizer="bobyqa")) 

anova(m20a,m22a) # m22a

# Winning model
summary(m22a)
tab_model(m22a, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='m22a.html',
          pred.labels = format_labels(m22a), title = format_formula(m22a))
# interactions with Cbal are trend level only now, but let's still take a look at them

# Contrasts
contrast(emmeans(m22a, ~ group * block * task_perf * cBal), # pairwise comparison between the levels
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>%  
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE)  %>% kable  

# Plots
{
  # whole model (ref. https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)
  mm <- m22a # choose model to plot
  pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
  c3m <- c3 %>% anti_join(c3 %>% filter(if_any(c(switch, group, nasas, task_perf, block), ~ is.na(.)))) %>% 
    mutate(fit = model.matrix(mm) %*% fixef (mm), 
           plo = fit-1.96*sqrt(pvar), phi = fit+1.96*sqrt(pvar))
  facet_cBal_block <- c('free-mix' = 'Free -> Mixed', 'mix-free' = 'Mixed -> Free', 'freeOnly' = 'Free block', 'mixed' = 'Mixed block')
  pl_m22a <- ggplot(c3m, aes(y = fit, x = task_perf, color = group, fill = group)) +
    geom_point(size = 2, position=pd) +
    geom_errorbar(aes(ymin=plo, ymax=phi, color = group), width = 0.4, position=pd) +
    scale_colour_discrete(labels = c('no control', 'full control')) +
    facet_grid(block~cBal, labeller = as_labeller(facet_cBal_block)) +
    scale_fill_discrete(labels = c('no control', 'full control')) +
    labs(y='Fit (switching)', x = 'Task', color = 'Group', fill = 'Group') +
    theme(legend.position=c(.85,.15), legend.key.size = unit(.5, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('pl_m22a.png', pl_m22a, path = 'plots_2/', w=5,h=5)
  
  # quick interactions 
  emmip(m22a, group ~ block | cBal) +
    geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4) 
}


# Switching: early subset ------------------------------------------------- 

# Prepare data subset
c3s <- c3 %>% 
  filter(cBal == 'free-mix', block == 'freeOnly') %>% 
  mutate(trial = case_when(blockNr == 1 ~ trialNr, blockNr == 2 ~ trialNr + 72),
         trials = standardise(trial))

glimpse(c3s)
nrow(c3s)/nrow(c3)
c3s %>% get_summary_stats(trial, trials)

# Model
ms3d <- glmer(switch ~ group * nasas * trials + (trials | id),
              c3s, family = binomial,
              glmerControl(optimizer = 'bobyqa'))
summary(ms3d)
tab_model(ms3d, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='ms3d.html',
          pred.labels = format_labels(ms3d), title = format_formula(ms3d))

# Contrasts
emmeans(ms3d, pairwise ~ group | trials, at = list(trials = c(-1.72, 0, 1.72)))$contrasts %>% # effect of group across trials
  as_tibble %>% select(-df, -SE, -z.ratio) %>% mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>% kable 

emmeans(ms3d, pairwise ~ group | nasas * trials, at = list(nasas = c(-1,0,1), trials = c(-1.72, 0, 1.72)))$contrasts %>% # effect of group across trials and NASA
  as_tibble %>% select(-df, -SE, -z.ratio) %>% mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>% kable 

contrast(emtrends(ms3d, ~ trials, var = 'nasas', at = list(trials = c(-1.72, 0, 1.72))), # pairwise slopes of NASA between trials
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>%  
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE,-z.ratio)

test(emtrends(ms3d, ~ trials, var = 'nasas', at = list(trials = c(-1.72, 0, 1.72))), adjust = "fdr") %>% # slopes of NASA across trials against 0
  as_tibble() %>% select(-df) %>% kable()

# Plots
{
  mm <- ms3d # choose model to plot
  pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
  c3m <- c3s %>% anti_join(c3s %>% filter(if_any(c(switch, group, nasa_bht, task_perf, block), ~ is.na(.)))) %>% 
    mutate(fit = model.matrix(mm) %*% fixef (mm), 
           plo = fit-1.96*sqrt(pvar), phi = fit+1.96*sqrt(pvar))
  facet_trial <- c('2' = 'First trial', '72' = 'Middle trial', '144' = 'Last trial')
  # v1
  pl_ms3d_v1 <- c3m %>% filter(trial %in% c(2,72,144)) %>% # 1st trial in which switch can be assessed; last trial of 1st block; last trial of 2nd block
    mutate(trial = as.factor(trial)) %>% 
    ggplot(., aes(y = fit, x = nasas, color = trial, fill = trial)) +
    geom_line(linewidth = .8) +
    geom_ribbon(aes(ymin=plo, ymax=phi, color = trial, fill = trial), alpha = 0.2, colour = NA) +
    facet_wrap(~group, labeller = as_labeller(facet_group)) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_colour_brewer(labels = c('first', 'middle', 'last'), palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    guides(fill = 'none') +
    labs(y = 'Linear prediction (switching)', x = 'NASA', colour = 'Trial') +
    theme(legend.position=c(.62,.2), legend.key.size = unit(.5, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('ms3v1.png', pl_ms3d_v1, path = 'plots_2/', width = 5.5, height = 4)
  
  # v2
  pl_ms3d_v2 <- c3m %>% filter(trial %in% c(2,72,144)) %>% # 1st trial in which switch can be assessed; last trial of 1st block; last trial of 2nd block
    mutate(trial = as.factor(trial)) %>% 
    ggplot(., aes(y = fit, x = nasas, color = group, fill = group)) +
    geom_line(linewidth = .8) +
    geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.2, colour = NA) +
    facet_wrap(~trial, labeller = as_labeller(facet_trial)) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_colour_discrete(labels = c('no control', 'full control')) +
    guides(fill = 'none') +
    labs(y = 'Linear prediction (switching)', x = 'NASA', colour = 'Group') +
    theme(legend.position=c(.11,.15), legend.key.size = unit(.5, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('ms3v2.png', pl_ms3d_v2, path = 'plots_2/', width = 7, height = 4)
  
  # interactions 
  emmip(ms3, trials ~ nasas , cov.reduce = range, at = list(nasas = c(-2, -1, 0, 1, 2, 3), trials = c(-1.72,0,1.72)), 
        plotit = F)  %>% 
    ggplot(., aes(x = xvar, y = yvar, colour = tvar)) +
    geom_line() +
    geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = tvar), alpha = 0.2, colour = NA) +
    scale_colour_discrete(labels = c('first', 'middle', 'last')) + 
    guides(fill = 'none') +
    labs(y = 'Linear prediction (switching)', x = 'NASA', colour = 'Trial')
}

ggpredict(model = ms3, terms = c("nasas","group")) %>% 
  plot(alpha=0.3, line.size=1)

# Task selection: simple model --------------------------------------------

mt1 <- glmer(task_dv ~ group * nasas + (1 | id),
             c3, family = binomial)
summary(mt1)
tab_model(mt1, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='mt1.html',
          pred.labels = format_labels(mt1), title = format_formula(mt1))
# nothing

# Task selection: complex model -------------------------------------------

# Compare models
mt3d <- glmer(task_dv ~ group * nasas * block + (block | id),
              c3 %>% filter(!is.na(switch_iv)), # remove switch NAs to compare to the model with switch as a factor
              family = binomial,
              glmerControl(optimizer="bobyqa"))  

mt4c <- glmer(task_dv ~ group * nasas * block * switch_iv + (block| id),
              c3 %>% filter(!is.na(switch_iv)),  
              family = binomial,
              control = glmerControl(optCtrl=list(maxfun=1e5), optimizer = 'bobyqa'))

mt4d <- glmer(task_dv ~ group * nasas * block * switch_iv + (switch_iv | id),
              c3 %>% filter(!is.na(switch_iv)),  
              family = binomial,
              control = glmerControl(optCtrl=list(maxfun=1e5), optimizer = 'bobyqa'))

mt4e <- glmer(task_dv ~ group * nasas * block * switch_iv + (block * switch_iv| id),
              c3, family = binomial,
              control = glmerControl(optCtrl=list(maxfun=1e5), optimizer = 'bobyqa'))
isConv(mt4e) # singular

anova (mt3d, mt4c, mt4d) # mt4c

# Winning model

summary(mt4c)
tab_model(mt4c, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='mt4c.html',
          pred.labels = format_labels(mt4c), title = format_formula(mt4c))

# Contrasts
contrast(emtrends(mt4c, ~ switch_iv, var = 'nasas'), # comparison of NASA slopes between the levels of switch
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>%  
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE)  %>% kable  

contrast(emmeans(mt4c, ~ group * block * switch_iv), # pairwise comparisons of factor levels
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>%  
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE)  %>% kable  

emmeans(mt4c, pairwise ~ group | nasas, at = list(nasas = c(-1, 0, 1)))$contrasts %>%  # effect of group on 3 points of NASA, average factors
  as_tibble %>% select(-df, -SE) %>% adjust_pvalue(p.col = 'p.value', output.col = 'adj.p', method = 'fdr') %>% kable  

# Plots  
{
  # whole model (ref. https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)
  mm <- mt4c # choose model to plot
  pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only; requires lots of memory! 
  c3m <- c3 %>% filter(!is.na(task_dv)) %>% filter(!is.na(switch_iv)) %>% 
    mutate(fit = model.matrix(mm) %*% fixef (mm), 
           plo = fit-1.96*sqrt(pvar), phi = fit+1.96*sqrt(pvar))
  facet_block_switch <- c('0' = 'Repeat', '1' = 'Switch', 'freeOnly' = 'Free block', 'mixed' = 'Mixed block')
  pl_mt4c <- ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
    scale_colour_discrete(labels = c('no control', 'full control')) +
    facet_grid(block~switch_iv, labeller = as_labeller(facet_block_switch)) +
    scale_fill_discrete(labels = c('no control', 'full control')) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(y='Fit (task selection)', x = 'NASA', color = 'Group', fill = 'Group') +
    theme(legend.position=c(.85,.14), legend.key.size = unit(.5, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('mt4c.png', pl_mt4c, path = 'plots_2/', width = 5.5, height = 5)
  
  # quick interactions 
  emmip(mt4c, switch_iv ~ nasas, cov.reduce = range) +
    geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = switch_iv), alpha = 0.3, colour = NA) 
  emmip(mt4c, block ~ switch_iv) +
    geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4) 

  emmip(mt4c, group ~ switch_iv  | block, plotit=F)  %>% 
    ggplot(., aes(x = xvar, y = yvar, colour = tvar)) +
    geom_point(size = 3, position=pd) +
    geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4, position=pd)  +
    facet_grid(~block, labeller = as_labeller(facet_block)) +
    scale_colour_discrete(labels = c('no control', 'full control')) +
    scale_x_discrete(labels = c('Repeat', 'Switch')) +
    guides(fill = 'none') +
    labs(y = 'Linear prediction (task selection)', x = 'Switch', colour = 'Group') +
    theme(legend.position=c(.85,.18), legend.key.size = unit(.5, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('mt4c_x.png', path = 'plots_2/', width = 5, height = 3.5)
  
}


# Task selection: early subset --------------------------------------------- 

mts3d <- glmer(task_dv ~ group * nasas * trials + (trials | id),
               c3s, family = binomial,
               glmerControl(optimizer = 'bobyqa'))
summary(mts3d)
tab_model(mts3d, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='mts3d.html',
          pred.labels = format_labels(mts3d), title = format_formula(mts3d))
# nothing


# Accuracy: complex model -------------------------------------------------
# (nothing in simple model)

ma2c <- glmer(corr ~ group * nasas * block * task_perf + (block * task_perf | id),
              c3 %>% filter(logrt != -Inf),
              family = binomial,
              glmerControl(optimizer = 'bobyqa')) #
summary(ma2c)
tab_model(ma2c, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='ma2c.html',
          pred.labels = format_labels(ma2c), title = format_formula(ma2c))

# Plots
{
  # whole model
  mm <- ma2c # choose model to plot
  pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
  c3m <- c3 %>% filter(logrt != -Inf) %>% 
    mutate(fit = model.matrix(mm) %*% fixef (mm), 
           plo = fit-1.96*sqrt(pvar), phi = fit+1.96*sqrt(pvar))
  pl_ma2c <- ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
    scale_colour_discrete(labels = c('no control', 'full control')) +
    facet_grid(block~task_perf, labeller = as_labeller(facet_task_block)) +
    scale_fill_discrete(labels = c('no control', 'full control')) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(y='Fit (accuracy)', x = 'NASA', color = 'Group', fill = 'Group') +
    theme(legend.position=c(.15,.12), legend.key.size = unit(.4, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('ma2c.png', pl_ma2c, path = 'plots_2/', w=5,h=4.5)
  
  
  # interactions
  emmip(ma2c, block ~ task_perf) +
    geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4) 
  emmip(ma2c, block ~ nasas, cov.reduce = range) +
    geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = block), alpha = 0.3, colour = NA) 
}


# Accuracy: early subset -------------------------------------------------- 

mas2a <- glmer(corr ~ group * nasas * trials + (trials | id),
               c3s,
               family = binomial,
               glmerControl(optimizer = 'bobyqa')) #
summary(mas2a)
tab_model(mas2a, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='mas2a.html',
          pred.labels = format_labels(mas2a), title = format_formula(mas2a))
# nothing

# RTs: complex model -------------------------------------------------------
# (nothing in simple model)

# Prepare
c3rt <- c3 %>% 
  filter(logrt != -Inf) %>% 
  left_join(c3 %>%
              filter(logrt != -Inf) %>% 
              group_by(group, task_perf, block) %>%
              get_summary_stats(logrt) %>% 
              mutate(lo = mean-3*sd,  # trials cleaning by +-3SD
                     hi = mean+3*sd) %>% 
              select(group:task_perf, mean,sd,lo,hi)) %>% 
  filter(logrt > lo & logrt < hi) %>% 
  mutate(logrtC = center(logrt))
1-nrow(c3rt)/nrow(c3) # 0.75% trials removed
glimpse(c3rt)

# Model
mrt2c <- lmer(logrt ~ group * nasas * block * task_perf + (block * task_perf | id),
              c3rt,
              control = lmerControl(optCtrl=list(maxfun=1e5), optimizer = 'bobyqa'))  
summary(mrt2c)
tab_model(mrt2c, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='mrt2c.html',
          pred.labels = format_labels(mrt2c), title = format_formula(mrt2c))

# Plots
{
  # whole model
  mm <- mrt2c # choose model to plot
  pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
  c3m <- c3rt %>% 
    mutate(fit = model.matrix(mm) %*% fixef (mm), 
           plo = fit-1.96*sqrt(pvar), phi = fit+1.96*sqrt(pvar))
  pl_mrt2c <- ggplot(c3m, aes(y = exp(fit), x = nasas, color = group, fill = group)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin=exp(plo), ymax=exp(phi), color = group, fill = group), alpha = 0.3, colour = NA) +
    scale_colour_discrete(labels = c('no control', 'full control')) +
    facet_grid(block~task_perf, labeller = as_labeller(facet_task_block)) +
    scale_fill_discrete(labels = c('no control', 'full control')) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(y='Fit (RT)', x = 'NASA', color = 'Group', fill = 'Group') +
    theme(legend.position=c(.85,.12), legend.key.size = unit(.4, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('mrt2c.png', pl_mrt2c, path = 'plots_2/', w=5,h=4.5)
  
  # interactions
  emmip(mrt2c, block ~ task_perf) +
    geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4) 
  emmip(mrt2c, block ~ nasas | task_perf, cov.reduce = range) +
    geom_ribbon(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE), fill = block), alpha = 0.3, colour = NA) 
}


# RTs: early subset ------------------------------------------------------- 

# Prepare
c3srt <- c3rt %>% 
  filter(cBal == 'free-mix', block == 'freeOnly') %>% 
  mutate(trial = case_when(blockNr == 1 ~ trialNr, blockNr == 2 ~ trialNr + 72),
         trials = standardise(trial))
glimpse(c3srt)

# Model
mrts1a <- lmer(logrt ~ group * nasas * trials + (trials | id),
               c3srt)
summary(mrts1a)
tab_model(mrts1a, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='mrts1a.html',
          pred.labels = format_labels(mrts1a), title = format_formula(mrts1a))
# just effect of trials


# Rigidness ---------------------------------------------------------------

r3 # this df has 1 row per single run of a task before switching to another task, subset free block only; ntr - nr of trials per such single run (before switch)
ggplot(r3,aes(x=ntr)) + geom_histogram(binwidth = 1) # post-cleaning by 3SD (ntr>46)


mr5 <- lmer(ntr ~ group * nasas * task_perf + (task_perf|id), r3)  #  
summary(mr5)
tab_model(mr5, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='mr5.html',
          pred.labels = format_labels(mr5), title = format_formula(mr5))

# Plots
{
  mm <- mr5 # choose model to plot
  pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only
  c3m <- r3 %>% 
    mutate(fit = model.matrix(mm) %*% fixef (mm), 
           plo = fit-1.96*sqrt(pvar), phi = fit+1.96*sqrt(pvar))
  facet_task <- c('loc' = 'Location task', 'shape' = 'Shape task')
  pl_mr5 <- ggplot(c3m, aes(y = fit, x = nasas, color = group, fill = group)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
    scale_colour_discrete(labels = c('no control', 'full control')) +
    scale_fill_discrete(labels = c('no control', 'full control')) +
    scale_x_continuous(breaks = pretty_breaks()) +
    facet_wrap(~task_perf, labeller = as_labeller(facet_task))+
    labs(y='Fit (rigidness)', x = 'NASA', color = 'Group', fill = 'Group') +
    theme(legend.position=c(.85,.15), legend.key.size = unit(.5, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('mr5.png', pl_mr5, path = 'plots_2/', w=5,h=4)
  
  emmip(mr5, ~ task_perf) +
    geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4) 
}

# Subset
r3s <- r3 %>% 
  filter(cBal=='free-mix')

mrs5 <- lmer(ntr ~ group * nasas * task_perf + (task_perf|id), r3s)  #  
summary(mrs5)
tab_model(mrs5, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='',
          # file='mr5.html',
          pred.labels = format_labels(mrs5), title = format_formula(mrs5))


# Switching: NASA VTS -----------------------------------------------------

mv3a <- glmer(switch ~ group * nasasV * block * task_perf + (block | id),
              c3, family = binomial,
              glmerControl(optimizer="bobyqa")) 
summary(mv3a)
tab_model(mv3a, show.re.var=F, show.icc=F, show.obs=F, show.r2=F, show.ngroups=F, dv.labels='', 
          # file='mv3a.html',
          pred.labels = format_labels(mv3a), title = format_formula(mv3a))

# Contrasts
emmeans(mv3a, pairwise ~ group | block * task_perf)$contrasts %>%  # effect of group in interaction
  as_tibble %>% select(-df, -SE) %>% adjust_pvalue(p.col = 'p.value', output.col = 'adj.p', method = 'fdr') %>% kable  

contrast(emtrends(mv3a, ~ block * task_perf, var = 'nasasV'), # pairwise comparison of NASA slopes between the levels 
         "pairwise", simple = 'each', combine = TRUE, adjust = "fdr") %>%  
  as_tibble() %>% mutate_at("p.value", list(~round(.,4))) %>% select(-df,-SE)  %>% kable  
 
emmeans(mv3a, pairwise ~ group | nasasV, at = list(nasasV = c(-1, 0, 1)))$contrasts %>%  # effect of group on 3 points of NASA, average factors
  as_tibble %>% select(-df, -SE) %>% adjust_pvalue(p.col = 'p.value', output.col = 'adj.p', method = 'fdr') %>% kable  
 

# Plots  
{
  # whole model (ref. https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)
  mm <- mv3a # choose model to plot
  pvar <- diag(model.matrix(mm) %*% tcrossprod(vcov(mm), model.matrix(mm))) # fixed-effects uncertainty only; requires lots of memory! 
  c3m <- c3 %>% anti_join(c3 %>% filter(if_any(c(switch, group, nasas, task_perf, block), ~ is.na(.)))) %>% 
    mutate(fit = model.matrix(mm) %*% fixef (mm), 
           plo = fit-1.96*sqrt(pvar), phi = fit+1.96*sqrt(pvar))
  facet_task_block <- c('loc' = 'Location task', 'shape' = 'Shape task', 'freeOnly' = 'Free block', 'mixed' = 'Mixed block')
  pl_mv3a <- ggplot(c3m, aes(y = fit, x = nasasV, color = group, fill = group)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin=plo, ymax=phi, color = group, fill = group), alpha = 0.3, colour = NA) +
    scale_colour_discrete(labels = c('no control', 'full control')) +
    facet_grid(block~task_perf, labeller = as_labeller(facet_task_block)) +
    scale_fill_discrete(labels = c('no control', 'full control')) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(y='Fit (switching)', x = 'NASA (VTS)', color = 'Group', fill = 'Group') +
    theme(legend.position=c(.85,.15), legend.key.size = unit(.5, 'cm'), 
          legend.background = element_rect(fill = "white", color = "grey85"),
          legend.title = element_text(size=13))
  ggsave('mv3a.png', pl_mv3a, path = 'plots_2/', w=5,h=4.5)
  
  # interactions
  emmip(mv3a, group ~ task_perf | block) + #  
    geom_errorbar(aes(ymin=yvar-(2*SE), ymax=yvar+(2*SE)), width= .4) 
}