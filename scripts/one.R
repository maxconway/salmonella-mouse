library(dplyr)
library(readr)
library(FluxBalanceAnalyzeR)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(magrittr)

setwd('~/git/salmonella-mouse/')

optimize_and_minimize <- function(mod){
  res1 <- gurobi::gurobi(mod, params = list(OutputFlag=0))
  
  mod1 <- mod
  flux <- res1$x
  mod1$lb[flux >= 0] <- 0
  mod1$ub[flux <= 0] <- 0
  mod1$lb[mod$obj !=0] <- flux[mod$obj !=0]
  mod1$ub[mod$obj !=0] <- flux[mod$obj !=0]
  mod1$obj <- -sign(flux)
  
  res2 <- gurobi::gurobi(mod1, params = list(OutputFlag=0))
  res1$x <- res2$x
  return(res1)
}

# helper function
find_fluxes <- function(reaction_table){
  optimresult <- parse_reaction_table(reaction_table) %>%
    optimize_and_minimize
  if('x' %in% names(optimresult)){
    return(optimresult$x)
  }else{
    return(0)
  }
}

# Read in reaction table
rxns <- read_tsv('data/LT2_react.tsv') %>%
  select(abbreviation=`Rxn name`,
         equation=Formula,
         uppbnd=UB,
         lowbnd=LB,
         geneAssociation=`Gene-reaction association`,
         Subsystem,
         name = `Rxn description`
         ) %>%
  mutate(obj_coef=1*(abbreviation=='biomass_iRR1083')) %>%
  mutate(geneAssociation = str_replace_all(geneAssociation, 'and', '&'),
         geneAssociation = str_replace_all(geneAssociation, 'or', '|'))

# From reaction table, list all genes
genes_in_model <- rxns$geneAssociation %>%
  str_split('[()|& ]+') %>%
  flatten_chr() %>%
  discard(is.na) %>%
  discard(~ str_length(.x)==0)

# map between STM style codes and gene names
genes_mapping <- read_tsv('data/LT2_genes.tsv')

gene_names_in_model <- genes_mapping %>%
  filter(Gene %in% genes_in_model) %>%
  getElement('Alias')

# select which group to use as the control. This will produce a wild type output.
control <- 'Input'

# read in rna counts and preprocess them
rnas <- inner_join(
  read_tsv('data/seq_counts.tsv') %>%
    rename(Alias = Gene),
  genes_mapping %>%
    select(Alias,
           Gene)
) %>%
  filter(Gene %in% genes_in_model) %>%
  gather(measurement, count, -Locus, -Gene, -Product, -Alias) %>%
  separate(measurement, c('group', 'replicate'), '_') %>%
  mutate(replicate = str_c(group, '_', replicate)) %>%
  # next, remove experimental variance between replicates
  group_by(replicate) %>%
  mutate(count = count/mean(count)) %>%
  ungroup %>%
  # now, average and normalize groups against control
  group_by(Gene, group) %>%
  summarise(meancount = mean(count)) %>%
  spread(group, meancount) %>%
  ungroup %>%
  mutate(Group1 = Group1/.[[control]],
         Group2 = Group2/.[[control]],
         Group3 = Group3/.[[control]],
         Input = Input/.[[control]]) %>%
  gather('group', 'nexpression', -Gene) %>%
  mutate(nexpression = pmin(nexpression, qexp(0.99, 1/mean(nexpression)))) %>%
  mutate(presence = qgamma(pexp(nexpression, 1), 20,20)) %T>%
  (function(x){print(x %>%
      ggplot(aes(y=presence, x=nexpression)) +
      geom_point()
  )})

# project the gene expressions onto the fluxes
models_with_expression <- rnas %>%
  group_by(group) %>%
  do(bind_rows(., data_frame(Gene=setdiff(genes_in_model, .$Gene), presence=NA, group=.$group[1]))) %>%
  by_slice(function(x){
    rxns %>%
      mutate(activation = FluxBalanceAnalyzeR::gene_eval(geneAssociation, x$Gene, x$presence))
  }) %>%
  unnest(.out)


# adjust flux bounds in line with activation
models_with_targets <- models_with_expression %>%
  group_by(group) %>%
  by_slice(function(x){
    x %>% mutate(baseflux = find_fluxes(x)) # this is useful if we want the flux bounds to be a function of the initial fluxes
  }) %>%
  unnest(.out) %>%
  #mutate(activation=1) %>% # for testing
  mutate(uppbnd = ifelse(is.na(activation), uppbnd, uppbnd * activation),
         lowbnd = ifelse(is.na(activation), lowbnd, lowbnd * activation)) %T>%
  (function(x){print(x %>%
                       filter(!(is.na(activation))) %>%
                       #filter(abs(baseflux)<900) %>%
                       ggplot(aes(x = baseflux)) + 
                       geom_ribbon(aes(ymin=lowbnd, ymax=uppbnd)) + 
                       geom_point(aes(y=baseflux), colour='blue')
  )})

# calculate final fluxes
results <- models_with_targets %>%
  group_by(group) %>%
  by_slice(function(x){
    x %>% 
      mutate(flux = find_fluxes(x))
  }) %>%
  unnest(.out) %T>%
  (. %>%
  filter(str_detect(abbreviation, 'iomass')) %>%
  select(group, abbreviation, flux) %>% 
  print
  )

# analyse
## filter
results_analysed_2 <- results %>%
  group_by(abbreviation) %>%
  mutate(variety = sd(flux)/mean(abs(flux))) %>%
  ungroup %>%
  filter(variety>1)

## plot
results_analysed_2 %>%  
  ggplot(aes(
    x=name, 
    y=flux, 
    fill = group
  )) + 
  geom_bar(stat='identity', position='dodge') + 
  scale_y_continuous(trans = scales::trans_new('cbrt', function(x){sign(x)*abs(x)^(1/3)}, function(x){x^3}, domain = c(-Inf, Inf))) +
  coord_flip()

results %>%
  select(group, abbreviation, equation, flux) %>%
  inner_join(.,.,by=c('abbreviation'='abbreviation','equation'='equation')) %>%
  mutate(contrast=str_c(group.x, ' - ', group.y), diff = flux.x-flux.y) %>%
  filter(contrast %in% c('Input - Group1','Input - Group3','Group1 - Group2')) %>%
  filter(abs(diff)>0.1*pmax(abs(flux.x), abs(flux.y))) %>%
  ggplot(aes(x=abbreviation, y=diff)) + geom_bar(stat='identity') + facet_grid(~contrast) + coord_flip()
  
# looking by metabolites:
full_list_results <- models_with_targets %>%
  nest(.key='fba_mod', -group) %>%
  mutate(gurobi_mod = map(fba_mod, parse_reaction_table),
         gurobi_result = map(gurobi_mod, optimize_and_minimize),
         fluxmat = map2(gurobi_mod, gurobi_result, function(mod, res){
           mod$A * matrix(res$x, nrow=nrow(mod$A), ncol=ncol(mod$A), byrow=TRUE)
         }))

reaction_metab_results <- full_list_results %>%
  select(group, fluxmat) %>%
  mutate(fluxdf = map(fluxmat, function(x){
    mat <- as(x, 'TsparseMatrix')
    data_frame(reaction = colnames(x)[mat@j+1], metabolite = rownames(x)[mat@i+1], flux = mat@x)
  })) %>%
  select(-fluxmat) %>%
  unnest(fluxdf)



# reaction_metab_results %>%
#   inner_join(.,.,by=c('reaction','metabolite')) %>%
#   mutate(contrast=str_c(group.x, ' - ', group.y), diff = flux.x-flux.y) %>%
#   filter(abs(diff)>(0.1*pmax(abs(flux.x), abs(flux.y)))) %>%
#   filter(contrast %in% c('Input - Group1','Input - Group3','Group1 - Group2', 'Group2 - Group3')) %>%
#   select(-group.x, -group.y, -flux.x, -flux.y) %>%
#   spread(metabolite, diff, fill=0) %>%
#   group_by(contrast) %>%
#   by_slice(function(x){
#     res <- as.matrix(x %>% select(-reaction))
#     rownames(res) <- x$reaction
#     return(res)
#   }) %>%
#   pwalk(function(contrast, .out){
#     heatmap.plus::heatmap.plus(.out, 
#                                distfun = (function(x){dist(abs(x))}), 
#                                col=RColorBrewer::brewer.pal(11,'RdYlGn'), 
#                                main=contrast
#     )
#   })


# # shows that at least some groups show effect
#  read_tsv('data/seq_counts.tsv') %>%
#   gather(measurement, count, -Locus, -Gene, -Product) %>%
#   separate(measurement, c('group', 'replicate'), '_') %>%
#   mutate(replicate = str_c(group, '_', replicate)) %>%
#   # first we'll do a round of outlider detection
#   group_by(replicate) %>%
#   filter(count < 100*mad(count)+median(count)) %>%
#   ungroup %>%
#   # next, remove experimental variance between replicates
#   group_by(replicate) %>%
#   mutate(count = count/mean(count)) %>%
#   ungroup %>%
#   # now, average and normalize groups against control
#   group_by(Gene, group) %>%
#   summarise(meancount = mean(count)) %>%
#   spread(group, meancount) %>%
#   ungroup %>%
#   mutate(Group1 = Group1/.[[control]],
#          Group2 = Group2/.[[control]],
#          Group3 = Group3/.[[control]],
#          Input = Input/.[[control]]) %>%
#   gather('group', 'nexpression', -Gene) %>%
#   mutate(normalized = qgamma(plnorm(nexpression),shape = 10,scale = 0.1)) %>%
#   rename(presence = normalized)
#   gather(measurement, count, -Locus, -Gene, -Product) %>%
#   separate(measurement, c('group', 'replicate'), '_') %>%
#   filter(group %in% c('Group1', 'Group2', 'Group3')) %>%
#   group_by(Gene, group) %>%
#   mutate(intra_group_sd=sd(count)) %>%
#   group_by(Gene) %>%
#   mutate(inter_goup_sd = sd(count)) %>%
#   mutate(intra_inter_ratio = inter_goup_sd/intra_group_sd) %>%
#   ggplot(aes(x=intra_inter_ratio)) + geom_density() + scale_x_log10()