
# topic models
# https://www.nicholas-ollberding.com/post/applying-topic-models-to-microbiome-data-in-r/

library(topicmodels)
library(ldatuning)

# unrarefied
# taxa detected in 2 or fewer samples removed
# O

PHYLO <- read_rds('phyloseq_obj.rds')
F_phylo <- PHYLO %>% prune_samples(samples=PHYLO@sam_data$tissue == 'F')
Q_phylo <- PHYLO %>% prune_samples(samples=PHYLO@sam_data$tissue == 'Q')

F_count_matrix <-
  otu_table(F_phylo) %>% 
  t() %>%
  as(Class='matrix')


F_result <- 
  FindTopicsNumber(
  F_count_matrix,
  topics = seq(from = 3, to = 50, by = 1),
  metrics = c("Griffiths2004","CaoJuan2009", "Arun2010", "Deveaud2014"),
  method = "Gibbs",
  control = list(seed = 243),
  mc.cores = 48,
  verbose = TRUE
)

F_plot <- F_result %>% ldatuning::FindTopicsNumber_plot()

# LDA(F_count_matrix, k = 15, method = 'Gibbs' )

#15 topics for F



Q_count_matrix <- otu_table(Q_phylo) %>% t() %>% as(Class='matrix')
Q_result <- 
  FindTopicsNumber(
    Q_count_matrix,
    topics = seq(from = 3, to = 50, by = 1),
    metrics = c("CaoJuan2009", "Arun2010", "Deveaud2014"),
    method = "Gibbs",
    control = list(seed = 243),
    mc.cores = 48,
    verbose = TRUE
  )

Q_plot <- Q_result %>% ldatuning::FindTopicsNumber_plot()

(Q_result$Arun2010 - min(Q_result$Arun2010)) / max(Q_result$Arun2010 - min(Q_result$Arun2010))
Q_result_dat <- 
  Q_result %>%
  pivot_longer(cols = -topics,
               names_to = 'metric',
               values_to = 'value') %>% 
  group_by(metric) %>% 
  mutate(value=(value-min(value))/max(value-min(value)), 
         objective=case_when(
           metric == 'Deveaud2014' ~ 'maximize', 
           TRUE ~ 'minimize'), 
         objective=factor(objective, levels = c('minimize', 'maximize')))
Q_result_dat %>% 
  ungroup() %>% 
  group_by(metric) %>% 
  mutate(
    dif_from_best=value - value[which.min(value)]
    
  )

Q_plot2 <-
  Q_result %>%
  pivot_longer(cols = -topics,
               names_to = 'metric',
               values_to = 'value') %>% 
  group_by(metric) %>% 
  mutate(value=(value-min(value))/max(value-min(value)), 
         objective=case_when(
           metric == 'Deveaud2014' ~ 'maximize', 
           TRUE ~ 'minimize'), 
         objective=factor(objective, levels = c('minimize', 'maximize'))) %>% 
  ggplot(aes(x=topics, y=value, color=metric, group=metric)) +
  geom_point(size=3) + 
  geom_path() + 
  scale_x_continuous(breaks = seq(min(Q_result$topics), max(Q_result$topics), by = 1) ) +
  facet_wrap(~objective, ncol=1, strip.position = 'right') +
  theme(panel.grid.major.x = element_line(colour = 'grey'), 
        panel.background = element_blank())

Q_plot2


# 18 topics for Q
