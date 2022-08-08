
# topic models
# https://www.nicholas-ollberding.com/post/applying-topic-models-to-microbiome-data-in-r/

library(topicmodels)
library(ldatuning)
library(tidyverse)

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

F_LDA <- LDA(F_count_matrix, k = 15, method = 'Gibbs' )

#15 topics for F



Q_count_matrix <- otu_table(Q_phylo) %>% t() %>% as(Class='matrix')
Q_result <- 
  FindTopicsNumber(
    Q_count_matrix,
    topics = seq(from = 3, to = 50, by = 1),
    metrics = c("Griffiths2004","CaoJuan2009", "Arun2010", "Deveaud2014"),
    method = "Gibbs",
    control = list(seed = 243),
    mc.cores = 48,
    verbose = TRUE
  )

Q_plot <- Q_result %>% ldatuning::FindTopicsNumber_plot()


# 18 topics for Q

Q_LDA <- LDA(Q_count_matrix, k = 18, method = 'Gibbs' )






# 
# (Q_result$Arun2010 - min(Q_result$Arun2010)) / max(Q_result$Arun2010 - min(Q_result$Arun2010))
# Q_result_dat <- 
#   Q_result %>%
#   pivot_longer(cols = -topics,
#                names_to = 'metric',
#                values_to = 'value') %>% 
#   group_by(metric) %>% 
#   mutate(value=(value-min(value))/max(value-min(value)), 
#          objective=case_when(
#            metric == 'Deveaud2014' ~ 'maximize', 
#            TRUE ~ 'minimize'), 
#          objective=factor(objective, levels = c('minimize', 'maximize')))
# Q_result_dat %>% 
#   ungroup() %>% 
#   group_by(metric) %>% 
#   mutate(
#     dif_from_best=value - value[which.min(value)]
#     
#   )
# 
# Q_plot2 <-
#   Q_result %>%
#   pivot_longer(cols = -topics,
#                names_to = 'metric',
#                values_to = 'value') %>% 
#   group_by(metric) %>% 
#   mutate(value=(value-min(value))/max(value-min(value)), 
#          objective=case_when(
#            metric == 'Deveaud2014' ~ 'maximize', 
#            TRUE ~ 'minimize'), 
#          objective=factor(objective, levels = c('minimize', 'maximize'))) %>% 
#   ggplot(aes(x=topics, y=value, color=metric, group=metric)) +
#   geom_point(size=3) + 
#   geom_path() + 
#   scale_x_continuous(breaks = seq(min(Q_result$topics), max(Q_result$topics), by = 1) ) +
#   facet_wrap(~objective, ncol=1, strip.position = 'right') +
#   theme(panel.grid.major.x = element_line(colour = 'grey'), 
#         panel.background = element_blank())
# 
# Q_plot2
# 
# 
# # 18 topics for Q


###  extract useful things from topic models
# per-topic-per-word probabilities via the matrix = “beta”
# and the per-document-per-topic probabilities via the matrix = “gamma”

# lda_k36 <- LDA(count_matrix, k = 36, method = "VEM", control = list(seed = 243))
# 
# b_df <- data.frame(tidy(lda_k36, matrix = "beta"))
# 
# g_df <- data.frame(tidy(lda_k36, matrix = "gamma")) %>%
#   arrange(document, topic)
# 
# head(b_df)


install.packages('tidytext')


library(tidytext)
F_b_df <- tidy(F_LDA, matrix = "beta")


F_b_df %>% 
  filter(beta > 0.01) %>% 
  ggplot(aes(y=term, x=beta)) +
  geom_point() + 
  facet_wrap(~topic, scales='free_y')


F_g_df <- tidy(Q_LDA, matrix = "gamma") %>%
  arrange(document, topic)
F_g_df

F_topics <- 
  F_g_df %>%
  mutate(topic=paste0('F', topic)) %>% 
  pivot_wider(names_from = topic, values_from = gamma)

F_topics

#### Q now ##
Q_b_df <- tidy(Q_LDA, matrix = "beta")


Q_b_df %>% 
  filter(beta > 0.01) %>% 
  ggplot(aes(y=term, x=beta)) +
  geom_point() + 
  facet_wrap(~topic, scales='free_y')


Q_g_df <- tidy(Q_LDA, matrix = "gamma") %>%
  arrange(document, topic)

Q_topics <- 
  Q_g_df %>%
  mutate(topic=paste0('Q', topic)) %>% 
  pivot_wider(names_from = topic, values_from = gamma)

Q_topics


all_topics <- F_topics %>% inner_join(Q_topics)


######
# Building the topic model phyloseq object
# Now we just need to multiply the per-document-per-topic probabilities by the read count for each sample.


lib_size_df <- data.frame(sample_sums(ps_g)) %>%
  dplyr::rename("read_count" = "sample_sums.ps_g.") %>%
  rownames_to_column(var = "document")

tm_df <- left_join(lib_size_df, g_df) %>%
  mutate(topic_count = read_count * gamma,
         topic_count = round(topic_count, 0)) %>%
  dplyr::select(-read_count, -gamma) %>%
  pivot_wider(names_from = topic, values_from = topic_count) %>%
  dplyr::rename_with(~ paste0("Topic_", .), -document) %>%
  column_to_rownames(var = "document") %>%
  t(.) %>%
  data.frame(.)

(ps_topic_g <- phyloseq(
  sample_data(ps_g),
  otu_table(tm_df, taxa_are_rows = TRUE)))


