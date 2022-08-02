library(tidyverse)
library(phyloseq)
library(microbiome)



# TODO
# look for differences in sequencing depth between F and Q 
# subset to only samples with both F and Q
# clr for each matrix independently?
# other normalization?

# 
# 
# meta <- 
#   read_csv('data/FS12_16S_data/FS12_final_meta.csv') %>% sample_data()
# 
# sample_names(meta) <- meta$sample_ID
# OTU <- phyloseq::import_mothur(mothur_shared_file = 'data/FS12_16S_data/FS12.shared')
# 
# TAX <- import_mothur(mothur_constaxonomy_file = 'data/FS12_16S_data/FS12.taxonomy')

# PHYLO <- phyloseq(meta, OTU, TAX)

# write_rds(x = PHYLO, file = '~/Documents/enrich_predict/phyloseq_obj.rds')
# 
# 
# PHYLO <- 
#   PHYLO %>%
#   prune_samples(samples = 
#                   PHYLO@sam_data$day == 'D0' &
#                   !is.na(PHYLO@sam_data$AULC) &
#                   PHYLO@sam_data$experiment == 'X12b') 
# 
# PHYLO <- prune_taxa(taxa = taxa_sums(PHYLO) > 1, x = PHYLO) # remove singleton taxa
# PHYLO <- prune_taxa(PHYLO, taxa = otu_table(PHYLO) %>% rowSums() > 3) # remove taxa found in fewer that 3 samples
PHYLO <- read_rds('phyloseq_obj.rds')


# are all the F and Q samples matched?
SAMPS <- 
  tibble(ID=names(sample_sums(PHYLO)),
       ssums=sample_sums(PHYLO)) %>% 
  mutate(TYPE=sub('.*([QF])$','\\1',ID), 
         SAMP=sub('(.*)([QF])$','\\1',ID))

BAD_SAMPS <- 
  SAMPS %>% 
  group_by(SAMP) %>%
  tally() %>% 
  arrange((n)) %>% 
  filter(n ==1) %>% 
  pull(SAMP)

# These are samples without a matching Q or F pair
REMOVE_THESE <- 
  SAMPS %>%
  filter(SAMP %in% BAD_SAMPS) %>% 
  pull(ID)



PHYLO <- 
  PHYLO %>%
  prune_samples(samples = !(PHYLO@sam_data$sample_ID %in% REMOVE_THESE))


sample_sums(PHYLO) %>% hist()

# rarefy to even depth for each sample
# maybe try some methods without this...
PHYLO <- PHYLO %>% rarefy_even_depth()

#

tibble(ID=names(sample_sums(PHYLO)),
       ssums=sample_sums(PHYLO)) %>% 
  mutate(TYPE=sub('.*([QF])$','\\1',ID)) %>% 
  ggplot(aes(x=TYPE, y=ssums)) + geom_boxplot()

PHYLO_clr <- microbiome::transform(PHYLO, transform = 'clr')

PHYLO_rel_abund <- PHYLO %>% transform_sample_counts(fun = function(x) x / sum(x))


clr_melt <- psmelt2(PHYLO_clr)





PHYLO_rel_abund
# PHYLO@sam_data$experiment %>% table()


tst <- PHYLO_rel_abund %>% psmelt() 
tst %>% select(-ends_with('ate'), -log_sal, )
# BiocManager::install("microbiome")

# library(microbiome)
# centered log transform
# PHYLO_clr <- microbiome::transform(PHYLO, transform = 'clr')


# divergence()

rel_abund_mat <- abundances(PHYLO_rel_abund) %>% t()

rel_abund_dist <- vegan::vegdist(rel_abund_mat)

# distance between Q and F within each animal
intra_animal_distances <- 
  rel_abund_dist %>% 
  as(Class = 'matrix') %>%
  as.data.frame() %>% 
  rownames_to_column(var='from') %>% 
  pivot_longer(cols = -from, names_to = 'to', values_to = 'bray') %>% 
  mutate(from_ID=sub('(X12bP[0-9]+D0)[QF]','\\1',from), 
         to_ID = sub('(X12bP[0-9]+D0)[QF]','\\1',to), 
         from_type=sub('(X12bP[0-9]+D0)([QF])','\\2',from), 
         to_type = sub('(X12bP[0-9]+D0)([QF])','\\2',to)) %>% 
  filter(from_ID == to_ID) %>% 
  filter(from_type != to_type) %>% 
  filter(from_type == 'F') %>% 
  transmute(pignum=sub('X12b(P[0-9]+)D0','\\1',from_ID), 
            ID=from_ID, 
            Q_F_dist=bray)
  
# some animals Q is very similar to their F
# for others Q is very different than Q



#######


# Q diversity and F diversity
# done on counts
QF_diversity <- 
  diversity(x = PHYLO, index = 'all', zeroes = TRUE) %>% 
  rownames_to_column(var = 'ID') %>% 
  mutate(type = sub('X12b(P[0-9]+)D0([QF])','\\2',ID), 
         pignum=sub('X12b(P[0-9]+)D0([QF])','\\1',ID)) %>% 
  mutate(ID=sub('[QF]','',ID)) %>% 
  pivot_longer(cols = -c(ID, pignum, type), names_to = 'div_type')


QF_diversity %>% 
  ggplot(aes(x=type, y=value, group=pignum)) + 
  geom_point() + 
  geom_path()+ 
  facet_wrap(~div_type, scales = 'free')


intra_animal_distances %>% left_join(QF_diversity)
#
# within each fecal sample, 
# is each taxa more abundant in the Q or F condition?
# hist(log(tst$Abundance), breaks=50)

tst <- 
  psmelt2(PHYLO_clr, feature.column = 'OTU') %>% 
  mutate(TYPE=sub('.*([QF])$','\\1', sample_ID), 
         ID=sub('(.*)([QF])$','\\1', sample_ID)) %>%
  select(-sample_ID, -.SampleID) %>% 
  select(-ends_with('ate')) 
# tst %>% g

tst <- 
  tst %>%
  group_by(OTU, ID) %>% 
  summarise(Q_m_F=value[tissue=='Q'] - value[tissue == 'F'], .groups = 'drop') %>% 
  arrange(desc(Q_m_F)) #%>% pull(OTU) %>% unique()
  # ggplot(aes(x=Q_d_F, y=Q_m_F)) + geom_point()

tax_table(PHYLO_clr)

# negative values mean OTU was more abundant in F types
# positive values mean OTU was more abundant in Q types (tet enrichments)
hist(tst$Q_m_F, breaks = 1000)

OTU_summary <- 
  tst %>%
  group_by(OTU) %>%
  summarise(av=mean(Q_m_F), 
            stdev=sd(Q_m_F), 
            stderr=stdev/n()) %>% 
  arrange(desc(av))

hist(OTU_summary$av, breaks = 100)

MODELS <- 
  tst %>% 
  group_by(OTU) %>%
  nest() %>% 
  mutate(mod=map(.x=data, .f=~lm(data=.x, Q_m_F ~ 1)), 
         SUM=map(.x=mod, .f=~summary(.x)), 
         CONFINT=map(.x=mod, .f=~confint(object = .x)),
         c.low=map_dbl(.x=CONFINT, .f=~pluck(.x, 1)),
         c.high=map_dbl(.x=CONFINT, .f=~pluck(.x, 2)),
         Est=map_dbl(.x=SUM, .f=~.x[['coefficients']][1]), 
         PV=map_dbl(.x=SUM, .f=~.x[['coefficients']][4])) %>% 
  ungroup() %>% 
  mutate(FDR=p.adjust(PV, method = 'fdr')) %>% 
  arrange((PV)) %>% 
  left_join(
    PHYLO@tax_table %>%
      as(Class = 'matrix') %>%
      as.data.frame() %>% 
      rownames_to_column(var = 'OTU')
  ) %>% 
  mutate(phylum=fct_lump_n(f=Rank2, n=7), 
         pos=ifelse(Est >0, 'pos', 'neg'))


MODELS$Rank2 %>% table()

MODELS %>%
  filter(FDR < 0.05) %>% 
  ggplot(aes(x=Est, y=OTU, color=phylum)) + 
  geom_point() + 
  geom_segment(aes(x=c.low, xend=c.high, yend=OTU, color=phylum)) +
  geom_vline(xintercept = 0, color='red') + 
  coord_flip()+
  facet_wrap(~phylum, nrow=1)

# TAX %>%
#   as(Class = 'matrix') %>%
#   as.data.frame() %>% 
#   rownames_to_column(var = 'OTU')

# MODELS$CONFINT[[1]][2]

# TMPPPP <- MODELS$SUM[[1]]
# 
# TMPPPP[['coefficients']][1]
# TMPPPP$coefficients[3]



# LM <- lm(formula = Q_m_F ~ OTU + ID, data = tst)


## create new features for each sample from the 
# interesting that actinobacteria are exclusively higher in the Q condition.
# one firmicutes OTU is very high 
#several proteobacteria very high but many quite low too
# bacteroides more low than high 


# should highlight the allstar OTUs from the FS12b publication...

# how to use these data to try and predict shedding level...
# i feel like I need to engineer some features
# condense the Q communities into one or two features

# 
# classify D0F OTUs by up or down in Q

### OOOH a feature could be magnitide of difference between F and Q communities


MODELS %>% filter(FDR <=0.1)
#



# what proportion of OTUs in Q were detected in F?


DETECTED_IN_F <- 
  PHYLO %>% 
  prune_samples(samples=grepl('.*F$',PHYLO@sam_data$sample_ID)) %>% 
  prune_taxa(taxa=taxa_sums(.) > 0) %>%
  taxa_names()

DETECTED_IN_Q <- 
  PHYLO %>% 
  prune_samples(samples=grepl('.*Q$',PHYLO@sam_data$sample_ID)) %>% 
  prune_taxa(taxa=taxa_sums(.) > 0) %>%
  taxa_names()


### THIS IS WRONG ###
ONLY_F <- DETECTED_IN_F[!(DETECTED_IN_F %in% DETECTED_IN_Q)]
ONLY_Q <- DETECTED_IN_Q[!(DETECTED_IN_Q %in% DETECTED_IN_F)]

length(ONLY_F)
length(ONLY_Q)


MODELS %>%
  mutate(OTU_PRES=case_when(
    OTU %in% ONLY_F  ~ 'only_F',
    OTU %in% ONLY_Q  ~ 'only_Q')) %>%
  # filter(FDR < 0.2) %>%
  # pull(OTU_PRES) %>%
  ggplot(aes(x=Est, y=OTU, color=OTU_PRES)) + 
  geom_point() + 
  geom_segment(aes(x=c.low, xend=c.high, yend=OTU)) +
  geom_vline(xintercept = 0, color='red') + 
  coord_flip()+
  facet_wrap(~phylum, nrow=1)

  

MODELS %>% filter(FDR < 0.1)






