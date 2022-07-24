library(tidyverse)
library(phyloseq)



# TODO
# look for differences in sequencing depth between F and Q 
# subset to only samples with both F and Q
# clr for each matrix independently?
# other normalization?



meta <- 
  read_csv('data/FS12_16S_data/FS12_final_meta.csv') %>% sample_data()

sample_names(meta) <- meta$sample_ID
OTU <- phyloseq::import_mothur(mothur_shared_file = 'data/FS12_16S_data/FS12.shared')

TAX <- import_mothur(mothur_constaxonomy_file = 'data/FS12_16S_data/FS12.taxonomy')

PHYLO <- phyloseq(meta, OTU, TAX)

write_rds(x = PHYLO, file = '~/Documents/enrich_predict/phyloseq_obj.rds')


PHYLO <- 
  PHYLO %>%
  prune_samples(samples = 
                  PHYLO@sam_data$day == 'D0' &
                  !is.na(PHYLO@sam_data$AULC) &
                  PHYLO@sam_data$experiment == 'X12b') 

PHYLO <- prune_taxa(taxa = taxa_sums(PHYLO) > 1, x = PHYLO) # remove singleton taxa
PHYLO <- prune_taxa(PHYLO, taxa = otu_table(PHYLO) %>% rowSums() > 3) # remove taxa found in fewer that 3 samples


PHYLO_clr <- microbiome::transform(PHYLO, transform = 'clr')

PHYLO_rel_abund <- PHYLO %>% transform_sample_counts(fun = function(x) x / sum(x))


clr_melt <- psmelt2(PHYLO_clr)





PHYLO_rel_abund
# PHYLO@sam_data$experiment %>% table()


tst <- PHYLO_rel_abund %>% psmelt() 
tst%>% select(-ends_with('ate'), -log_sal, )
BiocManager::install("microbiome")

library(microbiome)
# centered log transform
PHYLO_clr <- microbiome::transform(PHYLO, transform = 'clr')




# within each fecal sample, 
# is each taxa more abundant in the Q or F condition?
 hist(log(tst$Abundance), breaks=1000)

tst <- 
  psmelt2(PHYLO_clr, feature.column = 'OTU') %>% 
  mutate(TYPE=sub('.*([QF])$','\\1', sample_ID), 
         ID=sub('(.*)([QF])$','\\1', sample_ID)) %>%
  select(-sample_ID, -.SampleID) %>% 
  select(-ends_with('ate')) %>%
  group_by(OTU, ID) %>% 
  summarise(Q_m_F=value[tissue=='Q'] - value[tissue == 'F'], .groups = 'drop') %>% 
  arrange(desc(Q_m_F))# %>% 
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
  tst %>% group_by(OTU) %>% nest() %>% 
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
    TAX %>%
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
  facet_wrap(~phylum, nrow=1, scales = 'free')

TAX %>%
  as(Class = 'matrix') %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'OTU')

MODELS$CONFINT[[1]][2]

confint()
TMPPPP <- MODELS$SUM[[1]]

TMPPPP[['coefficients']][1]
TMPPPP$coefficients[3]



LM <- lm(formula = Q_m_F ~ OTU + ID, data = tst)
