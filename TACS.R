library(dplyr)
TACS = function( dge, gene1, gene2, cutoffs = NULL, 
                 return_val = "plot", density = F, 
                 facet_by = NULL, col = col, include_panel_with_all = F,
				 facet_level = FetchData(dge, facet_by)[[1]] %>% factor %>% levels %>% c(rep("all", include_panel_with_all), .),
                num_genes_add = 100, genesets_predetermined = F, dge_reference = dge, ... ){
  
  
  
get_similar_genes = function(dge, markers, n, anticorr = F ){
  data.use = dge@scale.data
  if(!all(markers %in% rownames(data.use))){ 
    warning("Some of your markers have no data available. Trying Various CASE Changes.")
    markers = unique( c( markers, toupper(markers), Capitalize( markers ) ) )
  }
  markers = intersect(markers, rownames( data.use) ) 
  correlation = rowSums( data.use %*% t( data.use[markers, , drop = F]) ) 
  correlation = correlation[ setdiff( names( correlation ), markers ) ]
  if( anticorr ){
    similar_genes = names( sort( abs( correlation ), decreasing = T )[ 1:n ] )
  } else {
    similar_genes = names( sort( correlation, decreasing = T )[ 1:n ] )
  }
  return( similar_genes )  
}
  
  
  # Get gene sets to average 
  library(magrittr)
  if(genesets_predetermined){
    g1_similar = gene1
    g2_similar = gene2
  } else {
    g1_similar = get_similar_genes(dge_reference, gene1, num_genes_add) %>% c( gene1, . )
    g2_similar = get_similar_genes(dge_reference, gene2, num_genes_add) %>% c( gene2, . ) 
    shared = intersect(g1_similar, g2_similar)
    g1_similar %<>% setdiff(shared)
    g2_similar %<>% setdiff(shared)
  }
  
  
  
  
  # Average gene sets to get scores
  get_score <- function(dge, gene){
    score_data <- dge@scale.data
    gene <- unique(gene)
    score <- colMeans(score_data[gene, ])
    return(score)
  }
  
  g1_score = get_score(dge, g1_similar)
  g2_score = get_score(dge, g2_similar)
  
  g1_score_name = paste0(gene1[1], "_score")
  g2_score_name = paste0(gene2[1], "_score")
  
  
  #Add scores as metadata. Extract with faceting var into plotting data.
  
  dge %<>% AddMetaData(g1_score, col.name = g1_score_name)
  dge %<>% AddMetaData(g2_score, col.name = g2_score_name)
  plot_df = FetchData(dge, c(g1_score_name, g2_score_name, facet_by))
  
  
  
  # Form plot
  p = ggplot(data = as.data.frame(plot_df)) 
  if(density){ 
    p = p + stat_density2d(aes_string( x = g1_score_name, y = g2_score_name, color = facet_by), bins = 50 ) + 
	scale_alpha_continuous(range = c(0.4, 1)) + scale_color_manual(values = col)
  } else {
    p = p + geom_point(aes_string( x=g1_score_name, y=g2_score_name))      
  }
  p = p + expand_limits(y=0, x=0)
  # Facet if desired
  if(!is.null(facet_by)) {
    p = p + facet_wrap(as.formula(paste0("~", facet_by)), ncol = 2)
  }
    
  # Add quadrants and percentages
  add_quadrants = function(p, g1_score_name, g2_score_name, cutoffs, facet_by = NULL){
  
  # Calculate percentages
  div_by_sum = function(x){if( sum(x) == 0) 0*x else x/sum(x)} 
  percentify = function(x){return(100*round(div_by_sum(x),3))}
 
  p = p + geom_vline(data = data.frame(xint = cutoffs[1]), aes(xintercept = xint))
  p = p + geom_hline(data = data.frame(yint = cutoffs[2]), aes(yintercept = yint))
  percentages = p$data[c(g1_score_name, g2_score_name, facet_by)]
  percentages[[g1_score_name]] %<>% is_greater_than(cutoffs[1])
  percentages[[g2_score_name]] %<>% is_greater_than(cutoffs[2])
  percentages %<>% table %>% (reshape2::melt)
  if(!is.null(facet_by)) {
    percentages = percentages[order(percentages[[facet_by]]), ]
    for( facet_level in unique(p$data[[facet_by]])){
      percentages[percentages[[facet_by]] == facet_level, "value"] %<>% percentify()
    } 
  } else {
    percentages$value %<>% percentify()
  }

  
  # Form annotation DF with correct facet and attempt sensible placement of percentages
  for( i in seq_along(percentages$value)){
    if(percentages$value[i]==0){next}
    annot_df = data.frame(
      x = ifelse( percentages[i, g1_score_name], 0.5, -0.2) ,
      y = ifelse( percentages[i, g2_score_name], 1, -0.8) ,
      label = paste0( round(percentages$value[i], 1), "%") )
    if(!is.null(facet_by)) {
      annot_df[[facet_by]] = percentages[i, facet_by]
    }
    p = p + geom_text( data = annot_df, aes(x=x,y=y,label=label), size = 4 )                
  }
  
  return(p)
}

  
  if(!is.null(cutoffs)){
    p %<>% add_quadrants(g1_score_name = g1_score_name,
                         g2_score_name = g2_score_name, 
                         cutoffs = cutoffs,
                         facet_by = facet_by)
  } 
  
  
  # Return everything or just a plot or just a seurat object
  if( return_val == "all" ){
    return( list( plot = p, 
                  dge = dge,             
                  genes = list( gene1 = gene1, gene2 = gene2 ),
                  score_names = c( g1_score_name, g2_score_name ), 
                  genesets = list( g1_similar, g2_similar ),
                  cutoffs = cutoffs,
                  plot_df = plot_df ) )
  } else if( return_val == "seurat" ){
    return(dge)
  } else if( return_val == "plot" ){
    return( p )
  } else {
    stop(" return_val should be 'all', 'seurat', or 'plot'. ")
  }
}
