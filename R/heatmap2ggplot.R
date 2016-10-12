heatmap2ggplot <- function(x){
  rowlevels <- x %>% 
    dist() %>%
    hclust() %>%
    with(labels[order])
  
  collevels <- x %>% 
    t() %>%
    dist() %>%
    hclust() %>%
    with(labels[order])
    
  x %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    gather(key = 'colname', value = 'value', -rowname) %>%
    mutate(colname = factor(colname, collevels),
           rowname = factor(rowname, rowlevels)) %>%
    group_by(colname) %>%
    filter(!all(between(value,-0.01,0.01))) %>%
    ungroup %>%
    group_by(rowname) %>%
    filter(!all(between(value,-0.01,0.01))) %>%
    ungroup %>%
    ggplot(aes(x=colname, y=rowname, fill = value)) +
    geom_tile() +
    scale_fill_gradient2()
}