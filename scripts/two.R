inner_join(
  read_tsv('data/seq_counts.tsv') %>%
    gather(measurement, count, -Locus, -Gene, -Product) %>%
    separate(measurement, c('group', 'replicate'), '_') %>%
    mutate(replicate = str_c(group, '_', replicate)) %>%
    select(-Locus, -Product) %>%
    filter(group!='Input') %>%
    nest(-Gene, -group, .key='expressiondata'),
  read_tsv('data/protein_weights.tsv') %>%
    select(-1, -3, -4, -9:-11) %>%
    gather(group, weight, `uninfected mouse`:`Group 3`) %>%
    mutate(group = str_replace(group, ' ', ''))
) %>%
  unnest(expressiondata) %>%
  group_by(replicate) %>%
  mutate(count = count/mean(count)) %>%
  ungroup %>%
  group_by(group) %>%
  mutate(weight = weight/mean(weight),
         count = count) %>%
  ungroup %>%
  group_by(Gene, group) %>%
  summarise(count = mean(count),
            weight = mean(weight)) %>%
  ungroup %>%
  mutate(count = qgamma(pexp(count, 1/mean(count)), 20,20)) %>%
  ggplot(aes(x=count, y=weight, colour=group)) + 
  geom_jitter(alpha=0.2) + 
  scale_y_continuous(trans = scales::trans_new('cbrt', function(x){sign(x)*abs(x)^(1/5)}, function(x){x^5}, domain = c(-Inf, Inf))) +
  scale_x_continuous(trans = scales::trans_new('cbrt', function(x){sign(x)*abs(x)^(1/5)}, function(x){x^5}, domain = c(-Inf, Inf)))
