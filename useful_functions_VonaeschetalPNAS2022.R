lunique <- function(x){
  #' length of unique observations
  #' @param x vector of values
  #' @return length of unique observations
  y = length(unique(x))
  return( y )
}

fix_donor_id <- function( id ){
  #' matching donorID with id
  #' @param id vector of ids
  id = gsub( "CPB", "ACPB", id )
  id = paste0("1429", id )

  return(id)

}

setup_pcoa <- function( Yd, both ){
  #' sets up the PCoA inputs
  #' @param Yd distance matrix (setup_adonis)
  #' @param both dataset from setup_adonis
  #' @return arr (dataset of arrows), vectorPCOA (plot of 2 components), cent (centroids)

out        = ape::pcoa( Yd )
vectorPCOA = cbind.data.frame( out$vectors[, 1:2 ] )
env.fit    = vegan::envfit( vectorPCOA, both, na.rm = TRUE)
arr        = data.frame(env.fit$vectors$arrows)
cent       = data.frame(env.fit$factors$centroids)
cent$amp   = rownames( cent )
arr$amp    = rownames( arr )

return( list( vectorPCOA = vectorPCOA, arr = arr, cent = cent ))

}

setup_adonis <- function( small, normer = "protein"){
  #' sets up the data structure for PERMANOVA analysis
  #'
  #' centers and scales variables
  #' @param small dataset
  #' @param normer normalization variable
  #' @return Yd (matrix of distances); tiny (AMP variables)

normer_var = paste0("_", normer, "_log")

tiny = small %>%
  dplyr::select( c(  ends_with(  normer_var ) ) ) %>%
  as.matrix() %>%
  scale( center = TRUE, scale = TRUE )

colnames(tiny) = gsub( pattern = normer_var, replacement = "", colnames(tiny) )

#-- linear combinations of the columns
both    = cbind.data.frame( small, tiny )
both$id = as.factor( both$id )

#-- adonis PERMANOVA model
Yd   = vegan::vegdist( x = tiny, method = "euclidian") # other choices available

return( list(Yd = Yd, both = both, tiny = tiny ))

}

unmatched_ids <- function( x, y ){
  #' compares two sets of vectors
  #' @param x vector
  #' @param y second vector
  #' @return vector of entries in one but not the other

  # x = ogens$id
  # y = AMPdata$id

  ##-- which ones aren't matched
  x_in_y = is.element( el = x, set = y )
  y_in_x = is.element( el = y, set = x )

  ##-- unique ids
  xx = unique(x)
  yy = unique(y)

  length(x)  ; length(y)
  length(xx) ; length(yy)

  # elements in x, but not in y
  only_x = unique( x[ which( !x_in_y) ])

  # elements in y, but not in x
  only_y = unique( y[ which( !y_in_x) ])

  sum(x_in_y); sum(y_in_x)
  length(x_in_y); length(y_in_x)
  sum(!x_in_y)
  sum(!y_in_x)

  return(list( only_x = only_x, only_y = only_y))

}


many_interactions <- function( mod, elas_data
                               , pred_x = "age"
                               , cat_z = "stunted"
                               , country = 'pays'){
  #' plots many interaction plots
  #' gives 95% prediction intervals and plots the points
  #' @param mod model fit
  #' @param elas_data dataset for model
  #' @param pred_x vector of 3 predictors (x axes)
  #' @param cat_z categorical variable (in addition to country)
  #' @param country default variable called 'pays'

for (j in seq_len( length(pred_x ))){
g1 = interactions::interact_plot( data = elas_data, model = mod
                             , pred = !!( pred_x[j] ), modx = !!(cat_z)
                             , plot.points = TRUE
                             , mod2 = !!( country ) , point.alpha = 0.2
                             , interval = TRUE
                             , int.type = "prediction"
                             , int.width = .95)
print(g1)
}
  }

excel_table <- function( qs, qs_header ){
  #' wrapper for excelTable
  #' @param qs data.frame
  #' @param qs_header columns of dataframe

rownames(qs) = NULL
excelR::excelTable( qs , colHeaders = qs_header
                    , defaultColWidth = 120, rowDrag = FALSE
                    , wordWrap = TRUE, allowDeleteRow = FALSE, allowDeleteColumn = FALSE
                    , editable = FALSE, allowInsertRow = FALSE, allowInsertColumn = FALSE
                    , autoColTypes = FALSE)
}

correlation_plot <- function( dat, var1, var2, title ){
  #' creates a correlation plot for two variables in a dataset
  #' @param dat dataset with var1 and var2
  #' @param var1 column of variable 1
  #' @param var2 column of variable 2
  #' @param title name of pdf
  #' @param savePDF if true saves PDF
  #' @return returns correlation coefficient
  #' @importFrom ggpubr ggscatter

  corr = cor( as.numeric(dat[, var1])
            , as.numeric(dat[, var2]), method = "spearman", use = "complete.obs")
  out  = ggpubr::ggscatter( dat, x = var1, y = var2, add = "reg.line"
                     , conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman"
                     , xlab = paste(var1, "levels")
                     , ylab = paste(var2, "levels"))

  print( out )
  return( corr )
}

# simple imputation step
crp_cleanup = function( x ){
  #' transforming crp values from 001,08 format to 1.08 format
  #' and replacing '<06,0' with smaller values
  #' @param x vector of crp character values
  # x = AMPdata$crp

  x = gsub(pattern = "<06,0", replacement = "6.0006", x = x)
  x = gsub(pattern =  ",", replacement = ".", x = x)
  x = as.numeric( as.character( x ))
  imputes = which( x == 6.0006 )
  x[ imputes ] = runif( length( imputes ), min = 3, max = 6 )
  return( x )
}


cleanup = function(x, under_value = 9e-11, over_value = 9999){
  #' imputes NA values if over or under a threshold
  #' according to diferent rules
  #' @param x vector of values
  #' @param under_value if smaller than this, imputes with half of minimum
  #' @param over_value if larger than this, imputes by twice the maximum
  x = as.numeric( as.character( x ) )

  # converts NAs to median value
  x[ is.na(x) ] = median( x, na.rm = TRUE )

  #if below the threshold, replace with half of minimum
  x[ x == under_value ] = NA
  x[ is.na(x) ] = 0.5 * min(x , na.rm = TRUE )

  #if above the threshold, replace with twice the maximum
  x[ x == over_value ] = NA
  x[ is.na(x) ] = 2 * max(x, na.rm = TRUE )

  return( x )
}

quantcut <- function(x, q = 4 ){
  #' cutting a variable into quantiles
  #' @param x values
  #' @param q number of quantiles
    qs     = seq(0,1, length.out = q + 1)
    ns     = 1 : q
    quant  = quantile(x, qs, na.rm = TRUE)
    retval = cut( x, breaks = quant, include.lowest = TRUE, labels = ns )
    return(retval)
}


return_ids <- function( AMPdata, compartments = c("duodenal","gastric"), count = FALSE ){
  #' obtain subject IDs with values in both compartments
  #' @param AMPdata dataframe
  #' @param compartments names of compartment (2)
  #' @param count if TRUE, gives the number of IDS, else subsets the data
  #' @importFrom stringr str_detect
  #' @importFrom dplyr filter select
  #--- gastric/duodenal

  if( length(compartments) == 1) compartments = c(compartments, compartments)

ids = AMPdata %>%
  dplyr::filter(
    stringr::str_detect( pattern, compartments[1]) &
    stringr::str_detect( pattern, compartments[2])
    ) %>%
  dplyr::select( id ) %>%
  unique() %>%
  unlist()

sub = AMPdata %>% dplyr::filter( id %in% ids & type %in% c(compartments[1], compartments[2]))

if( count) return( length(ids)) else
  return(sub)
}

normalize_amp <- function( AMPdata, amp_vars, divide_by = "GIF", name = "log"){
  #' does a log transformation after dividing by a variable
  #' @param AMPdata dataset
  #' @param divide_by "GIF" or "TCN1"
  #' @param name name suffix of new variables
  #' @importFrom dplyr mutate ends_with select all_of

  out = AMPdata %>%
    #-- removing missing
    subset(., !is.na(AMPdata[, divide_by]))  %>%
    #-- renaming
    dplyr::mutate( dplyr::across(
      .cols = dplyr::all_of( amp_vars )
    , .fns  = ~ . / !! as.symbol( divide_by )
    , .names = paste0("{.col}_", divide_by)
      )) %>%
  dplyr::mutate( dplyr::across(
     .cols = dplyr::ends_with( divide_by )
    , .fns = log10
    , .names = paste0("{.col}_", name)
  )) %>%
  dplyr::select( -ends_with( divide_by ))

    return( out )
}


filter_transform <- function( dataset, reverse = FALSE ){
  #' does filter_taxa and transform_sample_counts
  #' @param dataset phyloseq dataset
  #' @param reverse if TRUE does transform, then filter
  #' otherwise filters then transforms
  #' @return dataset

  if ( reverse ){
    dataset = phyloseq::transform_sample_counts( dataset, fun = function(x) x * 100 / sum(x)   )
    dataset = phyloseq::filter_taxa( dataset, flist = function(x) mean(x) > 0.1, prune = TRUE )
  } else {
    dataset = phyloseq::filter_taxa( dataset, flist = function(x) var(x) > 0, prune = TRUE )
    dataset = phyloseq::transform_sample_counts( dataset, fun = function(x) x * 100 / sum(x)   )
  }

  return( dataset )
}


correlation_setup <- function( dataset, variables, pattern = NULL
                               , level = "Species"){
#' subsets a dataset and uses microbiome::associate for correlation check
#' @param dataset phyloseq dataset
#' @param variables variables of interest that don't follow a naming pattern
#' @param pattern variables of interest that all end with a naming pattern, NULL otherwise
#' @param level level for tax_glom ("Species","Order","Genus", etc )
#' @param filter_then_transform if FALSE, transforms then filters (see filter_transform function)

  # Species, Genus,
  if (level %in% c("Species","Genus")) filter_then_transform = TRUE
  if (level %in% c("Order","Class")  ) filter_then_transform = FALSE

  message( paste("Level:", level, "corresponds to"
                 , ifelse( filter_then_transform, "filter( var > 0) then transform"
                           , "transform then filter (mean > .1)")))

        # dataset = dffilteredduodenal; variables = c("CALPROTECTINEggdePS", "AATmggdePS" ); pattern = "conc"
        samp   = data.frame( phyloseq::sample_data( dataset ))

        #-- selecting specific variables
        if (!is.null( pattern )){
          cyto = samp %>% dplyr::select( c(dplyr::all_of( variables), dplyr::ends_with( pattern )))
        } else {
          cyto = samp %>% dplyr::select(   dplyr::all_of( variables)                              )
        }

        #-- selecting complete cases in sample data
        cyto = data.frame( cyto[ complete.cases(cyto), ] )

        #- tax_glom then transforming then filtering
        # filter_then_transform = FALSE
        trim_rel = dataset %>%
          phyloseq::tax_glom( taxrank = level ) %>%
          filter_transform( reverse = filter_then_transform )

        #- taxonomy dataset transposed
        bact_unrestrict = trim_rel %>%
          phyloseq::otu_table() %>%
          t() %>%
          as.data.frame()

        #-- which samples are in both datasets
        bact = bact_unrestrict %>%
          dplyr::filter(
            row.names( bact_unrestrict ) %in% row.names( cyto ))

        #-- adding Genus, Species (or Class, Order) to column names
        if( level %in% c("Genus","Species")
            ) order_names = c("Genus", "Species") else order_names = c("Class", "Order"  )

        colnames( bact ) =
          c( apply( phyloseq::tax_table( trim_rel )[, order_names ], MARGIN = 1, paste0, collapse = "_")    )

        #-- outputting two datasets
        return(list( bact = bact, cyto = cyto , bact_unrestrict = bact_unrestrict ))

      }# end function

associate_plot <- function( x, y = NULL, mode = "table", pdfname = "test.pdf", save_pdf = FALSE, height = 10, width = 10 ){
    #' does microbiome::associate and makes a plot if requested
    #' @param mode options are 'table' or 'matrix' for microbiome::associate function
    #' @param x vector
    #' @param y optional vector
    #' @param save_pdf if TRUE saves using 'pdfname'
    #' @param height height of PDF
    #' @param width width of PDF

    correlation.table = microbiome::associate(x, y, method = "spearman", mode = mode, p.adj.threshold = 0.05, n.signif = 1)

      if( mode == "table"){
        if (save_pdf) pdf( pdfname , width = width, height = height )
        p = microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)
        if (save_pdf) dev.off()
      }# end table mode
      return( correlation.table )
    }# end function

##### ordination functions
setup_ordination <- function( dataset
                 , level = "Order"
                 , meta_variables = c("pays","stunted","SampleType","sexe","sibo","anemie")
                 ){
  #' creates dataset for plotting ordination
  #' @param dataset phyloseq dataset
  #' @param level Kingdom, Order, Phylum, Genus, Species, etc.
  #' @param meta_variables variables to include for plotting purposes
  #' @return returns 'dist' of distances and both (values of PCOA variables and meta-variables)

  #dataset = dffiltered5_Mada_rar; level = "Order"

  df   = phyloseq::tax_glom( dataset , level )
  taxa = as.data.frame(phyloseq::tax_table( df ))
  meta = phyloseq::sample_data( df )[, meta_variables] %>%
    as.data.frame()

#-- selecting relevant columns and transposing
test   =
  as.data.frame(phyloseq::otu_table( df )) %>%
  dplyr::select( dplyr::starts_with("Sample_")   ) %>%
  t() %>%
  as.data.frame()

#-- assigning names
colnames(test)  = make.unique( as.character( taxa[, level] ) )

#-- removing NAs and 0 values
test [ is.na(test) ] = 0 #convert NAs to zeros

#-- will later add 'sampleid' to be rownames
test4 = dplyr::filter( test, rowSums( test ) != 0 )

# Here we use Bray-Curtis distance metric
dist  = test4 %>%
  microbiome::transform( transform = "log10") %>%
  vegan::vegdist( method = "bray")

return( list( dist = dist, test4 = test4, meta = meta  ))
}# end function

setup_pcoa_scores <- function( dist, test4, meta, sig_level = 0.001, perms = 499 ){
  #' sets up a PCOA plot
  #'
  #' (Pascale's original version, with scores)
  #' @param dist distance matrix
  #' @param test4 dataset with values used in PCOA
  #' @param meta dataset with values for metadata
  #' @param sig_level significance level for plotting
  #' @param perms number of permutations in vegan::envfit
  #' @return returns vectorPCOA and df_envfit (arrows) dataset

#### PCoA
PCOA       = ape::pcoa(dist)
vectorPCOA = cbind.data.frame( PCOA$vectors[, 1:2], meta, sampleid = rownames( test4 ) ) %>% dplyr::distinct()
env.fit    = vegan::envfit( PCOA$vectors[, 1:2], test4, perm = perms, na.rm = TRUE)
df_envfit  = as.data.frame( vegan::scores(env.fit, display = "vectors"))

# plot the eigenvalues and interpret
# barplot( PCOA$ values$ Relative_eig[ 1:10 ] )
# ape::biplot.pcoa(PCOA, test4, scale = FALSE, center = TRUE)
# print( out )

#### Extracting significant pvalues from envfit taxa
A      = as.list(env.fit$vectors)
pvals  = as.data.frame( A$pvals )
arrows = as.data.frame( A$arrows * sqrt(A$r) )
C      = cbind( arrows, pvals, Order = rownames( pvals) )
Cred   = subset( C, pvals <= sig_level )

#### Format taxa scores for plotting
df_envfit       = vegan::scores(env.fit, display = "vectors" )
df_envfit       = ((df_envfit) * vegan:::ordiArrowMul( df_envfit)) %>% as.data.frame()
df_envfit$Order = rownames(df_envfit)

### what if not-enough significant p-values ?
temp  = subset(df_envfit, df_envfit$Order %in% Cred$Order )

##-- otherwise including them all
if (nrow(temp) > 0 ) df_envfit = temp else warning("no arrows are significant; empty dataset")
rownames( df_envfit ) = NULL

return(list( df_envfit = df_envfit, vectorPCOA = vectorPCOA ))
}# end function

make_biplot <- function( vectorPCOA
                         , arrow_dataset
                         , color_vector
                         , shape_vector
                         , nice_colors = c("cyan3","mediumblue","darkorange3","goldenrod1")
                         , normer = "protein"
                         , mult = 5
                         , label_var = "amp"
                         , title = "PCoA bi-plot of AMP variables"){
  #' creates a biplot with arrows for variables
  #' @param arrow_dataset dataset of Axis.1, Axis.2 of variables (amp)
  #' @param color_vector vector with color grouping (e.g. small$ pays )
  #' @param shape_vector vector with shape grouping (small$ stunted )
  #' @param vectorPCOA dataset of Axis.1, Axis.2 of all subjects
  #' @param nice_colors values of colors to plot
  #' @param normer normalization variable used in PCOA
  #' @param mult multiplication factor for length of arrow segments


#-- renaming the grouping column

gg = ggplot( ) +
  geom_point(data = vectorPCOA
        , aes(y = Axis.2, x = Axis.1, col = color_vector
              , shape = shape_vector ) ) +
  scale_color_manual( values = nice_colors ) +
  geom_segment( data = arrow_dataset
        , aes( x = 0, y = 0
               , xend = mult * Axis.1
               , yend = mult * Axis.2 )
            , color = "#808080", alpha = 0.5
        , arrow = arrow(length = unit(0.2, "cm")) ) +
  geom_text( data = arrow_dataset
        , aes_string(y = "mult * Axis.2", x = "mult * Axis.1", label = label_var )
        , alpha = 0.5, vjust = 0, hjust = 0, cex = 2.5) +
  scale_shape_manual(values= c(16, 4, 1, 3) ) +
  labs( y = "Axis 2", x = "Axis 1", col = "Compartment:"
        , title = title
        , subtitle = paste( nrow(vectorPCOA), "obs normalized to:", normer)
        , shape = "Stunted Status") +
  theme_classic() +
  theme(legend.position = "bottom")

  return(gg)
}

make_pca_vars <- function( AMPdata, amp_vars, normer, number_pcs = 2
                           , logistic = FALSE ){
  #' creates PCA variables (2 of them)
  #' @param AMPdata original dataset
  #' @param amp_vars amp variables
  #' @param normer normalizer (GIF, TCN1)
  #' @param number_pcs number of principle components
  #' @param logistic if TRUE, uses logisticPCA package
  #' @return dataset with two new columns (PC1, PC2 )

  small = AMPdata

  if (logistic ) {

    tiny = small %>%
      dplyr::select( all_of( amp_vars )) %>%
      as.matrix()

    pc_cv = logisticPCA::cv.lpca( tiny, number_pcs )

    logpca_model = logisticPCA::logisticPCA(
      tiny, k = number_pcs, m = which.min(pc_cv))

    # principle components
    vals = logpca_model$PCs
    colnames(vals) = paste0("PC", 1:number_pcs )

  } else
if (!logistic){


  if( !is.null( normer )){
  small  = AMPdata %>%
    normalize_amp(amp_vars = amp_vars, divide_by = normer )

  tiny = small %>%
    dplyr::select( c(  ends_with( paste0("_", normer, "_log")) ) ) %>%
    as.matrix() %>%
    scale( center = TRUE, scale = TRUE )
  }

  if( is.null(normer)){

    tiny = small %>%
      dplyr::select( all_of( amp_vars )) %>%
      as.matrix() %>%
      scale( center = TRUE, scale = TRUE )

  }

  #-- PCA
  pc   = stats::prcomp( tiny, center = TRUE, scale = TRUE)
  #-- dimensions: 12 AMP variables x 2 columns  # number_pcs = 2
  rr   =  pc$rotation[, paste0("PC", 1 : number_pcs ) ]
  ##-- these variables could be used in a separate regression
  vals = tiny %*% rr

}

duo = cbind( small, vals )

return(duo)
}

make_haz_cats <- function( AMPdata ){
  #' creates categorical version of haz_cont
  #' @param AMPdata dataset
  #' @return full dataset with new variable 'haz_cats'

  out = AMPdata %>% dplyr::mutate(
  haz_cats = dplyr::case_when( haz_cont >=  0 ~ "[0, 2)"
                            ,  haz_cont > -1 ~ "(-1, 0]"
                            ,  haz_cont > -2 ~ "(-2, -1]"
                            ,  haz_cont > -3 ~ "(-3, -2]"
                            ,  haz_cont > -7 ~ "(-7, -3]"))

  return(out)
}


plot_predictions <- function( mod, small, response = "y"){
  #' plots model predictions, stratified by country
  #'@param mod model
  #'@param small dataset used in model
  #'@param response variable with response
  library(ggplot2)

small$pred = predict(mod)
lmz        = range( c(small[, response], small$pred))
gg = ggplot( small ) +
  geom_point(aes_string(  y = "pred"
                 , x = response
                 , col = "sex" ))+
  geom_abline( slope = 1, intercept = 0, lty = 3, col = 'red')+
  facet_wrap(. ~ pays ) +
  scale_x_continuous( limits = lmz )+
  scale_y_continuous( limits = lmz )+
  theme_minimal()

return(gg)
}

 all.subsets <- function(set) {
   #' generates all possible combinations of a set
   #' @return all possible combinations of a set
   #' @param set vector
   #' @source http://rsnippets.blogspot.com/2012/04/generating-all-subsets-of-set.html
    n   = length(set)
    bin = expand.grid(  plyr::rlply(n, c(FALSE, TRUE)))
    colnames(bin) = set
      return( bin[-1,] )

}

missing_size <- function( AMPdata, controls ){
  #' calculates the number of missing observations
  #' for each combination of controls
  #' reporting the unique values of sample sizes in a table
  #' @param AMPdata dataset
  #' @param controls vector of controls
  #' @return table of sample size and combination of controls

  subsets = all.subsets( controls)

  # loop over controls
  out = data.frame( n = NA, set = NA )
  for (i in seq_len( nrow(subsets) )){
    # i = 1
    ccc = controls[  unlist( subsets[i,]  ) ]
    tmp = AMPdata %>%
      non_missing( ccc ) %>%
      dplyr::select(id) %>%
      unlist() %>%
      lunique()
  out[i, ] = c(tmp, i)
  } #end loop

  # only the unique ones
  # uni = out[ !duplicated( out$n ), ]

  # matching with the controls of interest
  control_name =  apply( subsets[ out$set, ]
         , MARGIN = 1
         , function(x) paste0( sort( controls[x] ), collapse = ", ")
         )

  # reporting a table
  uni2 = out %>%
    dplyr::mutate( name = control_name ) %>%
    dplyr::arrange( -n ) %>%
    dplyr::select( - set )

  # unique ones
  uni = uni2 %>% dplyr::group_by( n ) %>% dplyr::mutate(
    len  = nchar( name ), last = max( len )) %>%
    dplyr::filter( last == len )

      return( uni[, c("n", "name")] )
}


non_missing <- function( AMPdata, controls){
  #' creates a filter for non-missing data, looping through controls
  #' @param AMPdata dataset
  #' @param controls vector of controls
  #' @return filtered dataset

  AMPdata %>% dplyr::filter(
      dplyr::across(
         .cols = dplyr::all_of( controls )
        , .fns = function(x) !is.na(x) ) )

}

permute_cols <- function( tiny ){
  #' permutes the columns of a dataset
  #'@param tiny dataset
  permm = tiny
  for(j in seq_len(ncol( permm ) ) ){
    permm[, j] = permm[ permute::shuffle( nrow(permm) ), j ]
  } # end loop
  return( permm )
}



lasso_fit <- function(y, x, k = 10){
  #' choosing a lambda parameter from LASSO
  #' also including hypothesis testing inference
  #' @param y observed values
  #' @param x design matrix


  #-- creates design matrix
  model_term = paste("~", paste( colnames(x), collapse = "+"))
  xx         = model.matrix( as.formula( model_term ), data = x )

  #-- centers and scales non-categorical ones
  cat_vars = which( apply( xx
                    , MARGIN = 2
                    , FUN = function(z) sum( range(z) == 0 | range(z) == 1 ) > 0
         ))
  con_vars = xx %>%
    as.data.frame() %>%
    dplyr::select( - all_of(cat_vars) ) %>%
    scale() %>%
    as.matrix()

  #-- joining the scaled and centered variables
  x2 = cbind( con_vars, xx[ , cat_vars])

  mod = glmnet::cv.glmnet(x = x2, y = y, nfolds = k )
  coe = coef( mod, s = "lambda.1se")

  las = islasso::islasso( y ~ xx, lambda = mod$lambda.1se )
  return(list( las = las, glm_mod = mod , lasso_coef = coe ) )
}

lars_best_fit <- function( y, x, k = 10 ){
  #' uses cross validation to find the best fit
  #' @source http://www.biostat.umn.edu/~weip/course/dm/examples/exampleforhighd1.R
  #' @return list of

  lars.fit    = lars::lars(x = x, y = y, type = "lasso")
  cv.lars.fit = lars::cv.lars(x = x, y = y, K = k, type = "lasso")

 # choose fraction based min cv error rule
 min.indx  = which.min( cv.lars.fit$ cv )
 s.cvmin   = cv.lars.fit$ index[ min.indx ]
 yhat.lars = predict(lars.fit, newx = x, s = s.cvmin, type = "fit", mode = "fraction")

# choose fraction based 1-se cv error rule
# largest value of lambda such that
# error is within 1 standard error of the minimum:
  cv1se = min(cv.lars.fit$ cv) + cv.lars.fit$cv.error[min.indx]
  indx2 = cv.lars.fit$ cv < cv1se
s.cvmin = max(cv.lars.fit$ index[indx2])

yhat.lars = predict(lars.fit, newx = x, s = s.cvmin,
                   type = "fit", mode = "fraction")

return( list( pred = yhat.lars
              , min_lambda = s.cvmin
              , fit = lars.fit ))
}


cmh_procedure <- function( small
                           , pathogens = 1:9
                           , outvar = "low_Elas"
                           , covar = "pays"
                           , value = "Madagascar"
                           , folder = NULL
                           , removefiles = TRUE
                           ){

#-- binary outcome, sorted by covariates
# small = dplyr::inner_join( elas2, ogens, by = "id") %>%
#   non_missing(controls = controls) %>%
#   dplyr::mutate( low_Elas = elas < 100) %>%
#   dplyr::arrange( pays )
if (is.null( folder)) folder = getwd()

rr = range(which( small[, covar] == value ))

#-- writing the covariate file
write.table( x = c(rr[2], nrow(small) - rr[2])
             , file = "cov.txt"
             , row.names = FALSE, col.names = FALSE, sep = " ")

#-- writing the data file
temp =  unlist(t(small[, pathogens]))
write.table( x = temp
             , file = "data.txt"
             , sep = " ", row.names = FALSE , col.names = FALSE )

#-- writing the label file
write.table( x = as.numeric(small[, outvar ])
            , file = "label.txt"
             , sep = " ", row.names = FALSE , col.names = FALSE )

out = fastcmh::runfastcmh( folder = folder
          , alpha = 0.05, showProcessing = FALSE
          , saveAllPvals = TRUE
          )
print( cat( out$summary, sep = "\t") )

if (removefiles){
file.remove( "label.txt")
file.remove( "cov.txt")
file.remove( "data.txt")
}
return(out)
}


### function version
richness_wilcox <- function( dataset, variable){
  #' wilcox test comparing up to 3 values of a variable
  #' reports Shannon, Chao1, Observed, and InvSimpson p-values
  #' @param dataset phyloseq dataset
  #' @param variable variable of interest
  #' @return p-values from Wilcoxon rank sum tests with continuity corrections

  # store the result
  results    = phyloseq::estimate_richness( dataset )
  d          = phyloseq::sample_data( dataset )

  # automatically detect the values of the variable
  values     = as.character( unlist( unique( d[, variable ]) ) )
  message( paste( paste(values, collapse = ", "), "are the values of:", variable))

  if (length(values) > 3 ) warning("Code is written only for 3 unique values of a variable")
  if (length(values) < 2 )
  {
    warning(paste("Wilcox test has only one unique value:", values))
    return( NULL)

  } else
    {
  val_one    = results[d[, variable] == values[1],]
  val_two    = results[d[, variable] == values[2],]

  # do the wilcox tests
  s12 = wilcox.test(val_one$ Shannon   , val_two$ Shannon)$p.value
  c12 = wilcox.test(val_one$ Chao1     , val_two$ Chao1)  $p.value
  o12 = wilcox.test(val_one$ Observed  , val_two$ Observed) $p.value
  i12 = wilcox.test(val_one$ InvSimpson, val_two$ InvSimpson)$p.value

  # dataset of results
  out = data.frame(
      dataset    = deparse(substitute( dataset))
    , variable   = variable
    , comparison = paste(values[1], "and", values[2])
    , Shannon   = signif( s12, 3)
    , Chao1     = signif( c12, 3)
    , Observed  = signif( o12, 3)
    , InvSimpson= signif( i12, 3)
  )


  if( length(values) == 3 ){
  val_tre    = results[d[, variable] == values[3],]

  # do the wilcox tests
  s13 = wilcox.test(val_one$ Shannon   , val_tre$ Shannon) $p.value
  c13 = wilcox.test(val_one$ Chao1     , val_tre$ Chao1)  $p.value
  o13 = wilcox.test(val_one$ Observed  , val_tre$ Observed) $p.value
  i13 = wilcox.test(val_one$ InvSimpson, val_tre$ InvSimpson)$p.value

  # do the wilcox tests
  s23 = wilcox.test(val_two$ Shannon   , val_tre$ Shannon) $p.value
  c23 = wilcox.test(val_two$ Chao1     , val_tre$ Chao1)  $p.value
  o23 = wilcox.test(val_two$ Observed  , val_tre$ Observed) $p.value
  i23 = wilcox.test(val_two$ InvSimpson, val_tre$ InvSimpson)$p.value

  out = data.frame(
      dataset    = deparse(substitute( dataset))
    , variable   = variable
    , comparison = c(
       paste(values[1], "and", values[2])
      ,paste(values[1], "and", values[3])
      ,paste(values[2], "and", values[3])
    )
    , Shannon   = signif( c( s12, s13, s23 ), 3)
    , Chao1     = signif( c( c12, c13, c23 ), 3)
    , Observed  = signif( c( o12, o13, o23 ), 3)
    , InvSimpson= signif( c( i12, i13, i23 ), 3)
  )


  }# end three values

  my_data = data.frame( Shannon = results$Shannon, var = d[, variable] )
  res.aov = aov( Shannon ~ ., data = my_data )
  print( summary(res.aov) )

  return( out  )
  }# end more than one value to compare
}


pdf_richness <- function( dataset
                          , title = "Alpha Diversity according to Sample Type"
                          , xlab = "Sample Type"
                          , pdfname = "test"
                          , save_pdf = FALSE, width = 17, height = 8){
  #' plots the plot_richness function with boxplots
  #' @param dataset dataset
  #' @param title name on plot
  #' @param xlab xlab on plot
  #' @param pdfname name of pdf
  #' @param save_pdf if TRUE, saves PDF
  #' @param width width of image
  #' @param height height of image

  if( save_pdf ) pdf( paste0( pdfname, ".pdf"), width = width, height = height)
  p = phyloseq::plot_richness( dataset
                              , x = "SampleType"
                              , measures = c("Observed", "Chao1", "Shannon", "InvSimpson")
                              , title = title )
  p = p +
    geom_boxplot(data = p$data, aes(x = SampleType, y = value, color = NULL), alpha = 0.1) +
    theme_minimal() +
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18, face = "bold")
       , plot.title = element_text(size = 22, face = "bold", margin = margin(t = 0, r = 0, b = 20, l = 0))
       , strip.text.x = element_text(size = 18, face = "bold" )) +
    xlab( xlab )

 if( save_pdf ){
   print(p)
   dev.off()
 }# if pdf
  return( p )
}# end


table_richness <- function( dataset, filename = "DataAlphasampletype.txt", write_table = FALSE ){
  #' writes a table of the richness measures
  #' @param dataset phyloseq dataset
  #' @param filename name of text file to save
  #' @param write_table if TRUE, writes a table

  resultsalpha = phyloseq::estimate_richness( dataset
                      , measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
  DataAlpha    = phyloseq::sample_data(resultsalpha)
  DataAlpha$SampleID = rownames(DataAlpha)

  if( write_table) write.table(DataAlpha, file = filename ,sep = "\t", row.names = FALSE)
  return( DataAlpha )
}


exclude_missing <- function( dataset, variable ){
  #' removes NA, "" values and <NA> values from a phyloseq dataset
  #' @param dataset phyloseq dataset
  #' @param variable variable name
  #' @return subsetted dataset
  values  = c(unlist(phyloseq::sample_data( dataset)[ , variable ]))
  bads    = is.na(values)
  phyloseq::sample_data(dataset) = subset( phyloseq::sample_data( dataset), !bads )

  values = c(unlist(phyloseq::sample_data( dataset)[ , variable ]))
  bads   = (values == "")
  phyloseq::sample_data(dataset) = subset( phyloseq::sample_data( dataset), !bads )
  # print( unique( phyloseq::sample_data(dataset)[, variable ]) )

  return(dataset)
}

richness_functions <- function( dataset, variable
                                , write_table = FALSE
                                , save_pdf = FALSE
                                , pdfname = "test.pdf"
                                , filename = "test.txt"
                                , width = 17, height = 8
                                , title = "Alpha Diversity"){
  #' wrapper for three richness functions
  #' @param dataset dataset
  #' @param variable variable of interest
  #' @param pdfname name of PDF
  #' @param title title of graph
  #' @param write_table if TRUE, writes a table
  #' @param save_pdf if TRUE, saves the pdf
  #' @param filename name of txt file (write_table)
# dataset = df_red_rar_duodenal; variable = "sibo";
# filtering out missing and undesired quantities
  dataset = exclude_missing( dataset, variable )

out_graph = pdf_richness( dataset = dataset
              , title = title
              , xlab = variable
              , save_pdf = save_pdf
              , pdfname = pdfname
              , width = width
              , height = height )

out_table = table_richness( dataset = dataset
                , filename = filename
                , write_table = write_table )

out_wilcox = richness_wilcox( dataset = dataset , variable = variable )

return( list(
  out_graph = out_graph, out_wilcox = out_wilcox)
)
}


richness_loop <- function( dataset, variables, dataset_name
                           , save_pdf = TRUE, write_table = TRUE ){
  #' wrapper function for richness_functions to enable loops
  #' @param dataset phyloseq dataset
  #' @param variables vector of variables to loop through
  #' @param dataset_name relevant name of dataset (e.g. fecal_Mada)
  #' @param save_pdf if TRUE, saves PDF
  #' @param write_table if TRUE, saves table


  outs = list()
  wilcox = data.frame()

  for( j in seq_len( length(variables )) ){
  outs[[j]] = richness_functions( dataset
      , variable = variables[j]
      , pdfname = paste0("alpha_diversity_", variables[j],"_", dataset_name, ".pdf")
      , filename = paste0("alpha_diversity_", variables[j],"_", dataset_name, ".csv")
      , title = paste0("Alpha Diversity according to", variables[j], dataset_name)
      , save_pdf = save_pdf, write_table = write_table
     )

  #-- storing to wilcox
  wilcox    = rbind.data.frame( wilcox, outs[[j]]$out_wilcox)
  outs[[j]]$out_wilcox = NULL
  }# end loop

  wilcox$dataset = dataset_name

    return( list(outs = outs, wilcox = wilcox ))
}



beta_diversity_pdf <- function( res.adonis
                                , save_pdf = FALSE
                                , pdfname = "Contributionfulldataset.pdf"
                                , write_csv = FALSE
                                , csvname = "Dispersion_samples.csv"){
  #' takes results from adonis2 and creates a plot of beta diversity
  #' @param results results from vegan::adonis2( )
  #' @param save_pdf if TRUE saves pdf
  #' @param pdfname name of pdf
  #' @param write_csv if TRUE writes csv
  #' @param csvname name of csv

  results = data.frame( res.adonis ) %>%
    dplyr::mutate( Contribution = 100 * R2 )

  results$ variable   = rownames( results )

  pp = ggplot(
    # plots only the significant ones
    data = results %>% dplyr::filter( Pr..F. < 0.05 )
      , aes(  x = reorder(variable, Contribution)
            , y = Contribution ) ) +
    geom_bar(position = "dodge", stat = "identity", color = "black") +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12))+
    theme(axis.text.y = element_text(size = 12))+
    theme(axis.text.y = element_text(size = 12))+
    labs(   x = ""
          , y = "% Contribution"
          , title = "Sig. contribution to beta-dispersion")+
    theme(axis.title.y = element_text(size = 14))+
    theme(title = element_text(size = 16, face = "bold"))

  if (save_pdf ){
    pdf(pdfname , width = 6, height = 8)
    print(pp)
    dev.off()
  }# end pdf

  if (write_csv)  write.csv(results, csvname )
  return(pp)
}



#### prune singlets ####
prune_singlets = function(x) {
  x[x <= 1] = 0
  return(x)
}


# Tell the program how to interpret the names
parse_taxonomy_simple = function(char.vec) {
 ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
      "Species", "Accession", "ASV_sequence")
 rv = strsplit(char.vec, ";")[[1]]
 rv = c(rv, rep(NA, length(ranks) - length(rv)))
 names(rv) = ranks
 return(rv)
}



prevalence = function(x){
  #' the proportion of samples in which a species shows up (compared to the total number of samples).
  #' @param x numeric value
  x[x >= 2] = 1
  return(x)
}

allones = function(x){
  #' creating a matrix identical to the prevalence counts matrix
  #' @param x numeric value
  x[ x >= 0 ] = 1
 return(x)
}

merge_otu_frame <- function( dataset , name = "prevalence", group = "SampleType"){
  #' common transformation of phyloseq dataset
  #' @param group variable to merge with
  #' @param name prevalence or abundance in column names
  #' @param dataset phyloseq dataset
dat =  dataset %>%
 phyloseq::merge_samples( group = group ) %>%
 phyloseq::t() %>%
 phyloseq::otu_table() %>%
 as.data.frame()

# add something to column headers to distinguish between
# relative abundance and prevalence tables
colnames(dat) = paste(
  name, colnames(dat), sep = ".")

return(dat)
}


calculate_prevalence <- function( dataset
                                  , group = "SampleType"
                                  , csvname = "merge.prev_species_sampletype.csv"
                                  , write_csv = FALSE ){
  #' calculates prevalence for ____ for all samples
  #' where ___ is Species, Phylum, Genus, etc.
  #' @param dataset phyloseq tax_table
  #' @param csvname name of csv file
  #' @param write_csv if TRUE, writes a CSV
  #' @param group merging variable in merge_otu_frame

# this produces prevalence "counts" for each species, but not percentages
prev_counts.sampletype = dataset %>%
 phyloseq::transform_sample_counts(fun = prevalence) %>%
 merge_otu_frame( group = group)

# this produces a maximum possible prevalence count per species
prev_possible.sampletype = dataset %>%
 phyloseq::transform_sample_counts(fun = allones) %>%
 merge_otu_frame( group = group )

#dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-species basis
test.prev = ( prev_counts.sampletype / prev_possible.sampletype ) * 100

tax_table.sampletype = as.data.frame( phyloseq::tax_table( dataset ))

merge.prev_species = merge(tax_table.sampletype, test.prev, by = "row.names")

if( write_csv ) write.csv(merge.prev_species, csvname )

  return( merge.prev_species)
}



fuse_abundance <- function( dataset, prevdata
                            , csvname = "merge.abunanceprev.species.sampletype.csv"
                            , write_csv = FALSE
                            , group = "SampleType"){
  #' export abundance tables for species and sampletype and fuse them with prevalence tables to get a single one
  #' @param dataset phyloseq dataset
  #' @param prevdata prevalence dataset (calculate_prevalence)
  #' @param csvname name of CSV
  #' @param write_csv if TRUE, writes a CSV

dat = dataset %>%
  phyloseq::transform_sample_counts( function(x) { x/sum(x)} ) %>%
  phyloseq::merge_samples( group = group ) %>%
  phyloseq::t() %>%
  phyloseq::transform_sample_counts( fun = function(OTU) OTU * 100 / sum(OTU) ) %>%
  phyloseq::otu_table() %>%
  as.data.frame()

# add something to distinguish between relative abundance and prevalence
colnames(dat) = paste(
  "abundance", colnames(dat), sep = ".")

row.names( prevdata ) = prevdata $ Row.names

merged = merge( prevdata, dat , by = "row.names")

test = merged[,1] == merged[,2] # to see if the two colums are actually similar
# they are actually similar, so we can the two first rows once we actually set the rownames


merged = merged[,-1]
if (write_csv) write.csv(merged, csvname)

return( merged)
}

prevalence_abundance <- function( dataset
                                  , level = "Genus"
                                  , name
                                  , write_csv = TRUE, group = "SampleType" ){
  #' creates CSV files for prevalence and abundance
  #' by sorting at different levels (Phylum, Genus, Species)
  #' @param dataset phyloseq dataset
  #' @param level Kingdom, Phylum, Genus, Species, etc
  #' @param name dataset name
  #' @param group merging variable (e.g. "SampleType", "pays")
  #' @param write_csv if TRUE, writes csv files

  glom       = phyloseq::tax_glom( dataset , level )
  prev_level = calculate_prevalence(
                     dataset   = glom
                   , csvname   = paste0("merge_prev_", level, "_", name ,".csv")
                   , write_csv = write_csv
                   , group     = group )
  abundance  = fuse_abundance(
                    dataset   = glom
                  , prevdata  = prev_level
                  , csvname   = paste0("merge_prev_abund_", level, "_", name, ".csv")
                  , write_csv = write_csv
                  , group     = group )

return( abundance )
}


subset_cleanup <- function( dataset, controls){
  #' ensures no controls are missing and removes rowSums, colSums of otu zero
  #' @param dataset phyloseq dataset
  #' @param controls controls to be non-missing

#- removing missing values
for ( j in seq_len( length(controls))){
  dataset = exclude_missing( dataset, variable = controls[ j ])
}

  bads  = (rowSums( phyloseq::otu_table( dataset )) == 0)
  phyloseq::sample_data(dataset) = subset( phyloseq::sample_data( dataset), !bads )

  bads  = (colSums( phyloseq::otu_table( dataset )) == 0)
  phyloseq::sample_data(dataset) = subset( phyloseq::sample_data( dataset), !bads )

return( dataset )
}# end


beta_diversity_setup <- function( dataset, controls, pdfname = "test.pdf"
                                  , save_pdf = FALSE, write_csv = FALSE
                                  , csvname = "test.csv", perms = 999 ){
  #' cleans up the dataset and does the adonis2 analysis
  #' @param dataset phyloseq dataset
  #' @param controls controls for use in the model
  #' @param save_pdf if TRUe saves pdf
  #' @param write_csv if TRUE saves csv
  #' @param pdfname pdfname
  #' @param csvname name of csv file of adonis2 results
  #' @param perms permutations for adonis2

# are the required columns in here?
goods = !(controls %in% colnames( phyloseq::sample_data( dataset )))
if( sum(goods) > 0 ) stop("not all controls are defined in this dataset")

df           = subset_cleanup(dataset, controls = controls )
project_bray = phyloseq::distance( df,  method = "bray", type = "samples")
sample_df    = data.frame( phyloseq::sample_data(df))
form         = formula( c("project_bray ~", paste(controls, collapse = " + " )) )

print("beginning results")
# now the adonis test to see if there is a signficant difference according to different variables
res.adonis = vegan::adonis2( form
                , data = sample_df, method = "bray", permutations = perms  )

out = beta_diversity_pdf( res.adonis
                          , pdfname = pdfname, csvname = csvname
                          , save_pdf = save_pdf, write_csv = write_csv )

return( out )
}


chi_square_table <- function( dataset, variable1, variable2 ){
  #' outputs a table and chi-squared test
  #' @param dataset phyloseq dataset
  #' @param variable1 name of variable
  #' @param variable2 name of variable


  tab = table(
       unlist(phyloseq::sample_data(dataset)[, variable1])
      ,unlist(phyloseq::sample_data(dataset)[, variable2])       )
  print( tab )

  message(paste("Chi-square test of", variable1, "and", variable2))
  tab = chisq.test(
        unlist(phyloseq::sample_data(dataset)[, variable1])
      , unlist(phyloseq::sample_data(dataset)[, variable2])   )

  return(tab)
}


filter_prune <- function( dataset ){
  #' filters and prunes tax samples
  #' @param dataset phyloseq::tax_glom object
df = phyloseq::filter_taxa( physeq = dataset, flist = function(x) var(x) > 0, prune = TRUE)
df = phyloseq::prune_samples( samples = phyloseq::sample_sums( df ) > 0, df )
return(df)
}

gm_mean <- function(x, na.rm=TRUE){
  #' geometric mean
  #' @param x vector of numbers
  y = exp( sum( log( x[ x > 0] ), na.rm = na.rm ) / length(x) )
  return(y)
}


max_reorder_results <- function( df ){
  #' reorders results based on maximum log2FoldChange
  #' in terms of Order, Class, Family, Genus
  #' @param df results from deseq_names function
  #' @return payscorr in different order

  # Order order
x = tapply(df$log2FoldChange, df$Order, function(x) max(x))
x = sort(x, decreasing = TRUE)

# Class order
x = tapply(df$log2FoldChange, df$Class, function(x) max(x))
x = sort(x, decreasing = TRUE)
df$Class = factor(as.character(df$Class), levels = names(x))

# Family order
x = tapply(df$log2FoldChange, df$Family, function(x) max(x))
x = sort(x, decreasing = TRUE)
df$Family = factor(as.character(df$Family), levels = names(x))

# Genus order
x = tapply(df$log2FoldChange, df$Genus, function(x) max(x))
x = sort(x, decreasing = TRUE)
df$Genus = factor(as.character(df$Genus), levels = names(x))

df$ taxonomy = paste0( df$ Genus , "/", df$ Family)

return( df )

}

deseq_names <- function( res, merged
                         , fold_change = 1, sig_level = 0.05
                         , write_csv = TRUE
                         , name = "duodenal"){
  #' joins res (from deseq_cooks) and merged (from prevalence_abundance)
  #' into one object
  #' @param res output from deseq_cooks
  #' @param merged output from prevalence_abundance
  #' @param fold_change log2FoldChange >= this level
  #' @param sig_level only outputs adjusted p-values <= sig_level
  #' @param name name of CSV file
  #' @param write_csv if TRUE, writes csv

  results = as( res$results_ordered, "data.frame")
  results$ Row.names = rownames(results)

  payscorr = dplyr::inner_join( results, merged, by = "Row.names") %>%
    dplyr::filter( abs(log2FoldChange) >= fold_change & padj <= sig_level )

  if (write_csv) write.csv(x = payscorr, file = paste("payscorr_reduced_", name, ".csv"))
  return(payscorr)

}

deseq_cooks <- function( dataset, level = "Genus"
                         , full_model = c("age","sexe","Run","pays")
                         , variable_of_interest = "pays"
                         , sig_level = 0.05
                         , fold_change = 1){
  #' doing the DESeq2 analysis
  #'
  #' @param dataset phyloseq dataset
  #' @param level Genus, Species
  #' @param full_model terms on right hand side
  #' @param variable_of_interest term to remove from RHS in reduced model
  #' @param sig_level alpha level (0.05) for significance
  #' @param fold_change minimum fold change to display

  # level = "Genus"; dataset = dffilteredgastric
  ### formulas
  form     = as.formula( paste("~", paste( full_model , collapse = " + ")))
  red_form = as.formula( paste("~", paste( full_model[ !(full_model == variable_of_interest)]
                                       , collapse = " + ")))

  ### deseq commands
    df = phyloseq::tax_glom( dataset, level )
    dd = phyloseq::phyloseq_to_deseq2( df, design = form)
    gm = apply( DESeq2::counts( dd ), MARGIN = 1, gm_mean )

  ### analysis
    effect_size = DESeq2::estimateSizeFactors( dd , geoMeans = gm, type = "poscount")
    reduced     = DESeq2::DESeq( effect_size
                                 , test = "LRT"
                                 , fitType = "parametric"
                                 , reduced = red_form)

  ### testing an effect
    values = DESeq2::resultsNames( reduced )
    res    = DESeq2::results( reduced
                , cooksCutoff = FALSE
                , alpha = sig_level
                , pAdjustMethod = "BH"
                , name = values[ grepl( pattern = variable_of_interest
                                , x = values, ignore.case = TRUE) ]
                )

  # view a summary of the results table with a padj value < 0.01
  print( summary(res, alpha = sig_level))

  # filtering the results
  # reorder the results table by adjusted p-value and remove any "NA" entries
  ordered  = res[order(res$ padj, na.last = NA), ]
  filtered = ordered[ ordered$padj <= sig_level & abs(ordered$log2FoldChange) >= fold_change , ]

  return( list( results_ordered = ordered, results_filtered = filtered, glom = df  ))

}# end function

wilcox_strata_group <- function( metadata, sampledata, bacteria, strata, group ){
        #' does wilcox tests by first stratifying by 'strata' and then comparing levels of variable 'group'
        #' for the variable 'bacteria' in 'sampledata'
        #' @param metadata dataset with strata and group variables
        #' @param sampledata dataset with bacteria
        #' @param strata splits the dataset into strata
        #' @param group compares across levels of this variable
        #' @return dataset of strata, comparison, bacteria, and p-value

        # strata = 'pays'; group = 'sibo';
        test = data.frame(   strata   = metadata[, strata]
                           , g = metadata[, group]
                           , y = as.numeric( sampledata[, bacteria]) )

        outs  = data.frame()
        vals = unique( test$strata )
        for (j in seq_len( length(vals) )){# ids = "RCA"

          out = data.frame(
               strat = vals[j]
           ,  comp   = group
           ,  bacteria = bacteria
           ,   p.val = wilcox.test( formula = y ~ g,  data = test[ test$strata == vals[j] , ] )$p.value
        )

        outs = rbind.data.frame( outs, out )
        } # stratified loop

        return(outs)
        }

make_ggscatter_plots <- function( sampledata, metadata, bacteria, controls ){
        #' creates ggscatter plots for bacteria in sampledata versus controls in metadata
        #' @param sampledata dataset with 'bacteria'
        #' @param metadata dataset with 'controls'
        #' @param bacteria variable from sampledata
        #' @param controls variables from metadata
        #' @return list of ggscatter objects

        # bacteria = "Streptococcus"; controls = c("haz_cont", "UFCml")
        temp = cbind.data.frame( sampledata[, bacteria], metadata[, controls ])
        colnames( temp) = c( bacteria, controls )

        gg = list()
        for ( j in seq_len( length(controls))) {
          gg[[j]] = ggpubr::ggscatter( temp
                     , x =  bacteria
                     , y =  controls[j]
                     , add = "reg.line"
                     , conf.int = TRUE
                     , cor.coef = TRUE
                     , cor.method = "spearman"
                     , xlab = paste("Rel. abundance of", bacteria, "spp."
                                  , ylab = controls[j]) )
        }# end loop
        return(gg)
      }

differential_plot <- function( dataset, title = "Significantly different species by SIBO status in gastric samples, Bangui" ){
  #' creates a ggplot of log2FoldChange by taxonomy
  #' @param dataset output from deseq_names() function
  #' @param title ggplot title
  #' @return ggplot object

gg = ggplot(data =  dataset , aes(x = reorder(taxonomy, log2FoldChange), y = log2FoldChange)) +
 geom_bar(position = "dodge",stat = "identity", color = "black") +
 coord_flip() +
 theme(axis.text.x  = element_text(size = 12 ))+
 theme(axis.text.y  = element_text(size = 12 ))+
 theme(axis.text.y  = element_text(size = 12 ))+
 theme(axis.title.y = element_text(size = 14 ))+
 theme(title = element_text(size = 16, face = "bold"))+
    labs( title = title, xlab = "", ylab = "log2 fold change")

return(gg)
}
