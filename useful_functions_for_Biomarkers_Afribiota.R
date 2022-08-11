# USEFUL FUNCTIONS

#### de-identifying ####

deidentify_matching <- function( seed = 2812 ){
  #' NOTE: requires the 'full_ids' object used in 'initial_merging.R'
  set.seed( seed )

  new_ids = data.frame( id = sample( unique(full_ids$id), replace = FALSE ) )
  new_ids$new = paste0("Sample", 1:nrow(new_ids) )
  new_ids = new_ids %>%
    dplyr::full_join(full_ids, by = "id") %>%
    dplyr::mutate( new_afri = ifelse( !is.na(id_afri), paste0( new, "_", type ), NA)
                 , run  = as.numeric( as.factor(run) ) )

    #-- this one is duplicated
    #

  return( new_ids )
}

change_ids <- function( seq_dat, sample_dat, new_ids ){
  #' uses the 'new' and 'new_afri' columns of new_ids to replace the 'id' and 'id_afri'
  #' columns (seq_data) and rownames (sample_data)
  #' @param seq_dat sequencing data.frame (column names as 'id_afri')
  #' @param sample_dat sample data.frame( rownames as 'id_afri')
  #' @param new_ids dataset created by deidentify_matching() function

  # seq_dat = seqs; sample_dat = temp_meta
  goods   = colnames( seq_dat ) %in% new_ids$id_afri %>% which()
  seq_dat = seq_dat[ , goods ]

  #-- creating a matching between the sequencing names and the new names
  s_cols = data.frame( id_afri = colnames( seq_dat ) )
  n_cols = dplyr::inner_join( s_cols, new_ids )
  #-- directly applying the new names (they are in the same order as s_cols)
  colnames( seq_dat ) = n_cols$new_afri

  #-- all the ones are in the dataset
  goods      = rownames( sample_dat ) %in% new_ids$id_afri %>% which()
  sample_dat = sample_dat[ goods, ]

  #-- creating a matching between the sequencing names and the new names
  s_rows = data.frame( id_afri = rownames( sample_dat ) )
  n_rows = dplyr::left_join( s_rows, new_ids )



  new_ids[ new_ids$id_afri == "S0079.0457_S273", ]
 ( new_ids$id_afri == n_rows[879,] %$% id_afri ) %>% which()
  new_ids[787,]
  n_rows$new_afri %>% duplicated() %>% which()
  dim(n_rows)
  sample_dat$id_afri %>% unique() %>% length()
  #-- directly applying the new names (they are in the same order as s_cols)
  rownames(sample_dat)   = n_rows$ new_afri
  sample_dat $ id_afri   = n_rows$ new_afri
  sample_dat $ id        = n_rows$ new

  return( list(   sample_dat = sample_dat, seq_dat = seq_dat
  ))

}

##### data cleaning ####
easier_column_names <- function( x , compartment = "duo", ignore = "id"){
  #' renames the column names to get rid
  #'@param x vector of column names
  #'@param compartment suffix for location of variables
  #'@param ignore column names to ignore in the renaming
  # x     = c("id","il6cv","il1ObsConc","il2conc","il2CV","il2cv")

  x = tolower(x)

  #-- removing dashes
  x[ !x %in% ignore ]  = gsub( pattern = "-"
                               , replacement = ""
                               , x =  x[ !x %in% ignore ]  )

  #-- replacing alpha, beta with a, b
  x[ !x %in% ignore ]  = gsub( pattern = "interferon"
                             , replacement = "ifn"
                             , x =  x[ !x %in% ignore ]  )
  x[ !x %in% ignore ]  = gsub( pattern = "alpha"
                             , replacement = "a"
                             , x =  x[ !x %in% ignore ]  )
  x[ !x %in% ignore ]  = gsub( pattern = "beta"
                             , replacement = "b"
                             , x =  x[ !x %in% ignore ]  )
  #-- concentration
  x[ !x %in% ignore ]  = gsub( pattern = "obs_conc_"
                               , replacement = ""
                               , x = x[ !x %in% ignore ]  )
  x[ !x %in% ignore ]  = gsub( pattern = "obsconc"
                               , replacement = ""
                               , x = x[ !x %in% ignore ]  )
  x[ !x %in% ignore ]  = gsub( pattern = "conc"
                               , replacement = ""
                               , x =  x[ !x %in% ignore ]  )
  #-- cv
  x[ !x %in% ignore ]  = paste0( x[ !x %in% ignore ], "_", compartment )
  x[ !x %in% ignore ]  = gsub( pattern = paste0("cv_", compartment)
                               , replacement = paste0("_", compartment, "_cv")
                               , x =   x[ !x %in% ignore ]  )

  return(x)
  }

find_duplicates <- function( df , id = "id"){
  #' finds the duplicated ids and variables in a dataset
  #' @param df dataset

  #df = immu
  dup_ids = df[ df %>% dplyr::select( id ) %>% duplicated(), id] %>% unlist() %>% c()
  dup_df  = df[ df$id %in% dup_ids, ] %>%
    dplyr::group_by( id ) %>%
    #-- gives the values that are different
    dplyr::summarize(
      dplyr::across( .fns = function(x) ifelse( x == x[1], NA, x) )
    )

  dup_df2  = df[ df$id %in% dup_ids, ] %>%
    dplyr::group_by( id ) %>%
    #-- gives the values that are different
    dplyr::summarize(
      dplyr::across( .fns = function(x) ifelse( x == x[ dplyr::n()], NA, x) )
    )

  both = dplyr::full_join( dup_df2, dup_df)

  #-- filtering on the problematic ones
  filt = both[rowSums( !is.na( both )) > 1, ] %>% dplyr::arrange( id )

  if( nrow(filt) == 0 ) filt = NULL
  return(filt)
}

impute_double_half <- function(x){
  #' replaces values over the detection limit with twice the max
  #' replaces values under the detection limit with half the minimum
  #'@param vector of values with numbers and "< OOD" or "> OOD"

  # x = c(3,5,11,100,0.3,"> OOD", "<OOD")
  min_val = min( as.numeric( x), na.rm = TRUE )
  max_val = max( as.numeric( x), na.rm = TRUE )

  x[ grepl(pattern = "<", x = x ) ] = min_val / 2
  x[ grepl(pattern = ">", x = x ) ] = max_val * 2

  x = as.numeric(x)

  return(x)
}

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
  id = trimws( id )
  goods = grepl( "ACPB", x = id )
  id[ !goods ] = gsub( "CPB", "ACPB", id[ ! goods ] )

  #-- missing the 1429 prefix
  bads     = which( !grepl( "1429", id  ) )
  id[bads] = paste0("1429", id[bads] )

  return(id)
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
crp_cleanup = function( x ){
  #' transforming crp values from 001,08 format to 1.08 format
  #' and replacing '<06,0' with smaller values
  #' @param x vector of crp character values
  # x = AMPdata$crp

  x[ grepl(pattern = "<", x) ] = "6.0006"
  x = gsub(pattern =  ",", replacement = ".", x = x)
  x = as.numeric( as.character( x ))
  imputes = which( x == 6.0006 )
  x[ imputes ] = 3
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

collapse_na <- function( dat ){
  #' collapsing multiple rows (with NAs) into one row
  #' @param dat dataset
  dat = dat %>%
    dplyr::summarize(
      dplyr::across( .cols = dplyr::everything(), .fns = function(x)
        dplyr::first( stats::na.omit( x ) ) ) )

  return(dat)
}

check_if_complete <- function( dataset, controls){
  #' determines if a row is a complete case or has missing values
  #' @param dataset dataset
  #' @param controls vector of controls
  #' @returns dataset with additional column, 'any_missing' as TRUE or FALSE

  tt = dataset %>%
    dplyr::mutate(
      any_missing =  (
        rowSums(
          dplyr::across(
            .cols = dplyr::all_of( controls )
           , .fns = function(x) is.na(x) )
          ) != 0
        )
  )

  return(tt)
}# end function


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

ph_and_log_dataset <- function( small, controls, normer = "protein" ){
 #' creates dataset of controls and AMPs for linear regression
 #' @param small dataset
 #' @param controls controls
 #' @param normer normalization variable

 x = small %>% dplyr::select( c(
     dplyr::all_of( controls)
   , ends_with(paste0("_", normer,"_log"))
   ))

colnames(x)[ colnames(x) == "ph_estomac"] = "ph"
colnames(x)[ colnames(x) == "ph_intestin"] = "ph"

xtemp = x
colnames(xtemp)[ colnames(xtemp) == "cal_duo"] = "log10(cal_duo)"
colnames(xtemp)[ colnames(xtemp) == "AAT_duo"] = "log10(AAT_duo)"

return( list( x = x, xtemp = xtemp ))
}

##### pcoa analysis  #####
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

if ( is.null( normer ) ){

  tiny = small %>%
    base::Filter( is.numeric , x = .) %>%
    as.matrix() %>%
    scale( center = TRUE, scale = TRUE )

} else {

normer_var = paste0("_", normer, "_log")

tiny = small %>%
    dplyr::select( c(  ends_with(  normer_var ) ) ) %>%
    as.matrix() %>%
    scale( center = TRUE, scale = TRUE )

colnames(tiny) = gsub( pattern = normer_var, replacement = "", colnames(tiny) )
}

#-- linear combinations of the columns
both    = cbind.data.frame( small, tiny )
both$id = as.factor( both$id )

#-- adonis PERMANOVA model
Yd   = vegan::vegdist( x = tiny, method = "euclidian") # other choices available

return( list(Yd = Yd, both = both, tiny = tiny ))

}

#### visualization ####

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

  gg = list()
for (j in seq_len( length(pred_x ))){
gg[[j]] = interactions::interact_plot( data = elas_data, model = mod
                             , pred = !!( pred_x[j] ), modx = !!(cat_z)
                             , plot.points = TRUE
                             , mod2 = !!( country ) , point.alpha = 0.2
                             , interval = TRUE
                             , int.type = "prediction"
                             , int.width = .95
                             , robust = FALSE
                             )
}
  return(gg)
  }


diagonal_line <- function(){
  #' creates a diagonal line for ggplot2
  ggplot2::geom_abline( slope = 1, intercept = 0, lty = 3, col = 'red')

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
          speedyseq::tax_glom( taxrank = level ) %>%
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

associate_plot <- function( x, y = NULL, mode = "table"
                            , pdfname = "test.pdf"
                            , save_pdf = TRUE
                            , height = 10, width = 10 ){
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
        print(p)
        if (save_pdf) dev.off()
      }# end table mode
      return( correlation.table )
    }# end function

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

#### ordination functions (and plotting) ####

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

  df   = speedyseq::tax_glom( dataset , level )
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
  #' @note setup_adonis() then setup_pcoa(), then make_biplot()

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

original_biplot <- function(
  vectorPCOA
  , arrow_dataset
  , mult = 0.0003
  , color_vector
  , label_var = "Order"
  , label_text = "Sample Type"
){
  #' original biplot code in paper
  #' @param arrow_dataset dataset of Axis.1, Axis.2 of variables (amp)
  #' @param color_vector vector with color grouping (e.g. small$ pays )
  #' @param vectorPCOA dataset of Axis.1, Axis.2 of all subjects
  #' @param mult multiplication factor for length of arrow segments (DEPRECATED)
  #' @param label_text text of color category
  #' @param label_var grouping variable used in setup_ordination()

#### automatic scaling of the arrows to the plot
mm      = max( abs(  range(     vectorPCOA[, c("Axis.1", "Axis.2")] ) )  )
nn      = max( abs( range(   arrow_dataset[, c("Axis.1", "Axis.2")] ) )  )
scaling = mm / nn * 0.9

#### Make ordination plot with Families as explanatory variables
gg = ggplot() +
  geom_point(data = vectorPCOA,
  aes( x = Axis.1
      ,y = Axis.2
      , colour = color_vector
      )) +
 scale_color_manual(values = c("green","blue","red") ) +
 scale_shape_manual(values =  c(0,1,2,5,16) ) +
 geom_segment( data = arrow_dataset,
   aes(  x = 0
       , y = 0
       , xend = Axis.1 * scaling
       , yend = Axis.2 * scaling)
   , arrow = arrow(length = unit(0.2, "cm"))
   , color = "#808080", alpha = 0.5 )+
 geom_text( data = arrow_dataset
   , aes_string(y = "scaling * Axis.2", x = "scaling * Axis.1", label = label_var )
   , color = "#808080", alpha = 0.5
   , vjust = 2.5, hjust = .5) +
 theme_classic() +
  labs( color = label_text )

return(gg)
}

make_biplot <- function( vectorPCOA
       , arrow_dataset
       , color_vector
       , color_text = "Compartment:"
       , shape_vector
       , label_text = "Country of Origin"
       , nice_colors = c("cyan3","mediumblue","darkorange3","goldenrod1")
       , normer = NULL
       , label_var = "amp"
       , title = "PCoA bi-plot of AMP variables"
       , ... ){
  #' creates a biplot with arrows for variables
  #' @param arrow_dataset dataset of Axis.1, Axis.2 of variables (amp)
  #' @param color_vector vector with color grouping (e.g. small$ pays )
  #' @param shape_vector vector with shape grouping (small$ stunted )
  #' @param vectorPCOA dataset of Axis.1, Axis.2 of all subjects
  #' @param nice_colors values of colors to plot
  #' @param normer normalization variable used in PCOA
  #' @param mult multiplication factor for length of arrow segments (DEPRECATED)
  #' @param label_text text of symbol category
  #' @param color_Text text of color cateogry

#### automatic scaling of the arrows to the plot
mm      = max( abs(  range(     vectorPCOA[, c("Axis.1", "Axis.2")] ) )  )
nn      = max( abs( range(   arrow_dataset[, c("Axis.1", "Axis.2")] ) )  )
scaling = mm / nn * 0.9

#-- renaming the grouping column
gg = ggplot( ) +
  geom_point(data = vectorPCOA
        , aes(  y = Axis.2
              , x = Axis.1
              , col = color_vector
              , shape = shape_vector ) ) +
  scale_color_manual( values = nice_colors ) +
  geom_segment( data = arrow_dataset
        , aes( x = 0
               , y = 0
               , xend = scaling * Axis.1
               , yend = scaling * Axis.2 )
        , color = "#808080"
        , alpha = 0.5
        , arrow = arrow(length = unit(0.2, "cm")) ) +
  geom_text( data = arrow_dataset
        , aes_string(  y = "scaling * Axis.2"
                     , x = "scaling * Axis.1"
                     , label = label_var )
        , alpha = 0.5, vjust = 0, hjust = 0, cex = 2.5) +
  scale_shape_manual(values = c(16, 4, 1, 3) ) +
  labs(   y = "Axis 2"
        , x = "Axis 1"
        , col = color_text
        , title = title
        , subtitle = paste( nrow(vectorPCOA), "observations ", normer)
        , shape = label_text ) +
  theme_classic() +
  theme(legend.position = "bottom")

  return(gg)
}

#### bootstrapping ####
bootstrap_coef <- function( small, form , times = 2000, logistic = FALSE){
  #' wrapper function for getting the coefficients from a glm
  #' using a bootstrap approach
  #' @param small dataset
  #' @param form linear model formula
  #' @param times number of bootstrap samples
  #' @note https://www.tidymodels.org/learn/statistics/bootstrap/

family = ifelse( logistic == TRUE, "binomial", "gaussian")

fit_m = function( splits ){
  glm( formula = form , data  = splits, family = family )
}

boots = rsample::bootstraps( data = small, times = times, apparent = TRUE )

boot_models = boots %>%
  dplyr::mutate(
        model = purrr::map(splits, fit_m )
  , coef_info = purrr::map(model, rsample::tidy ) )
##- percentile intervals


pint = rsample::int_pctl(boot_models, statistics = coef_info) %>%
  dplyr::mutate( term = gsub("_protein_log", "", term )) %>%
  dplyr::select( - c(.alpha, .method))

pint$sig = ifelse( (pint$.lower > 0 & pint$.upper > 0) |
                   (pint$.upper < 0 & pint$.lower < 0), TRUE, FALSE )

return(pint)
}

plot_bootstrap <- function( pint ){
  #' plots the dataset created from bootstrap_coef()
  #' @param pint dataframe with columns, .estimate, term, .upper, .lower
  #' @return ggplot object
  library(ggplot2)

  gg = ggplot(pint %>% dplyr::filter( term != "(Intercept)")) +
    geom_point( aes(x = .estimate, y = term, col = sig ) ) +
    geom_errorbar( aes(x = .estimate, xmin = .lower, xmax = .upper
                       , y = term, col = sig ) ) +
    geom_vline( xintercept = 0 , lty = 3, col = "red" ) +
    theme_minimal() +
    labs( x = "Estimate", y = "Variable", col = "Statistically\nSignificant")

  return(gg)
}

unregister_dopar <- function() {
  #' removing parallel connections
  #' @source https://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster
  env <- foreach:::.foreachGlobals
  rm(list = ls(name = env), pos = env )
}

##### imputation functions ####
imputation_approach <- function( dat, n.chains = 5, n.iter = 50, formula, digits = 5 ){
  #' imputation with 'mi' package
  #'@param dat dataset
  #'@param n.chains number of imputation chains
  #'@param n.iter number of imputation iterations
  #'@param formula formula
  #'@param digits number of digits to display

mdf = mi::missing_data.frame(  y = data.frame( dat ) )
imp = mi::mi( mdf, n.iter = n.iter, n.chains = n.chains, max.minutes = 5 )
mod = mi::pool( formula = formula, data = imp   )

arm::display( mod , digits = digits, detail = TRUE)

return( mod )
}

plot_imputed_model <- function( dd ){
#' plots the dd dataset, text
#' @param dd dataset with columns coef.est, variable, sig, upper, lower
gg = ggplot( dd) +
  geom_point(aes(y = coef.est, x = variable, col = sig )) +
  geom_errorbar( aes(ymax = upper, ymin = lower, x = variable, col = sig)
                 , width = 0.5, alpha = 0.5 ) +
  theme_minimal() +
  coord_flip() +
  labs( x = "Variable", y = "Estimate", col = "Statistically\nSignificant"
        , title = "Results from Imputed Data"
        , subtitle = "(Pooled Estimates)")
return(gg)

}

#### rmarkdown helpful functions #####
excel_table <- function( qs, qs_header = NULL, include.rows = FALSE, width = 120  ){
  #' wrapper for excelTable
  #' @param qs data.frame
  #' @param qs_header columns of dataframe
  #' @param include.rows if TRUE, then include rownames

if ( !include.rows) rownames(qs) = NULL

excelR::excelTable( qs , colHeaders = qs_header
                    , defaultColWidth = width, rowDrag = FALSE
                    , wordWrap = TRUE, allowDeleteRow = FALSE, allowDeleteColumn = FALSE
                    , editable = FALSE, allowInsertRow = FALSE, allowInsertColumn = FALSE
                    , autoColTypes = FALSE
                    , autoFill = TRUE
                    , autoWidth = FALSE
                    )
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

#### table1 helper functions  ####

render_prevalence <- function(x, ... ){
  #' borrowed from table1::render.continuous.default function
  with(
    table1::stats.apply.rounding( table1::stats.default(x, ...), ...)
    , c( "" , `Prevalence (SD)` = sprintf("%s (%s)", MEAN, SD)
      ))
}

pvalue_table1 <- function(x, ...) {
#' Construct vectors of data y, and groups (strata) g
#'@note https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html

    y <- unlist(x)
    g <- factor(rep(1:length(x), times=sapply(x, length)))
    if (is.numeric(y)) {
        # For continuous variables, perform a Wilcoxon test
        p <- wilcox.test(y ~ g)$p.value
    } else {
        # For categorical variables, perform a chi-squared test of independence
        p <- chisq.test(table(y, g))$p.value
    }
    # Format the p-value, using an HTML entity for the less-than sign.
    # The initial empty string places the output on the line below the variable label.
    c("", sub("<", "&lt;", format.pval(p, digits = 3, eps = 0.001)))
}

label_table1 <- function( AMPdata ){
  #' adds labels to table1 output
  #' @param AMPdata dataset
  #' @return returns dataset with labels
table1::units( AMPdata$ age) = "months"
table1::label( AMPdata$ age) = "Age"
table1::units( AMPdata$ ageyears) = "years"
table1::label( AMPdata$ ageyears) = "Age"
table1::label( AMPdata$ haz) = "Stunting Status"
table1::label( AMPdata$ sibo) = "SIBO Status"
table1::label( AMPdata$ sex) = "Sex"
table1::label( AMPdata$ anemia) = "Anemia"
colnames(Metadata)

if( "type" %in% colnames(AMPdata)) table1::label( AMPdata$ type) = "Compartment"
if( "elas" %in% colnames(AMPdata)) table1::label( AMPdata$ elas) = "Elastase"

return( AMPdata)
}

nice_print <- function( mod ){
  #' prints a summary of the model (e.g. mod = lm( ... ) )
  #'@param mod model
  out = sjPlot::tab_model( mod )
  cat(out$knitr,"\n--------\n")
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

#### lasso and lars functions ####
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

#### unclassified function #####
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

LDM_results <- function( res.ldm ){
    #' displays significance and summary of LDM results
    #' @param res.ldm output from LDM::ldm( ) call
    #' @source code from LDM vignette

    #--- p-value of global variables
    global.pval = res.ldm $p.global.omni

    print("Global P-Values")
    print( signif( global.pval, 3) )

    sigs = length( res.ldm$detected.otu.omni )
    summary.tab = NULL
    if( sigs > 0 ){

      # detect significant OTUs
      w1 = match(res.ldm$detected.otu.omni[[1]], colnames(res.ldm$q.otu.omni))
      o  = w1[ order(res.ldm$p.otu.omni[1, w1]) ]

    # summary
    summary.tab = data.frame( raw.pvalue = signif( res.ldm$p.otu.omni[1,o],3)
                          , adj.pvalue = signif(res.ldm$q.otu.omni[1,o],3)
                          , mean.freq = signif(res.ldm$mean.freq[o],3)
                          , direction = t(ifelse(res.ldm$beta[1,]>0, "+", "-"))[o,]
                          , otu.name = colnames(res.ldm$q.otu.omni)[o]
                          , row.names = NULL )
    colnames(summary.tab)[4] = paste("direction.", rownames( res.ldm$beta[1,]), sep = "")

    }# end something significant

    return( summary.tab )
  }

#### richness functions ####
richness_wilcox <- function( dataset, variable, csv_name = NULL, dataset_name = NULL ){
  #' wilcox test comparing up to 3 values of a variable
  #' reports Shannon, Chao1, Observed, and InvSimpson p-values
  #' @param dataset phyloseq dataset
  #' @param variable variable of interest
  #' @param csv_name name of csv of Wilcox p-values
  #' @return p-values from Wilcoxon rank sum tests with continuity corrections
  #' @param dataset_name name of dataset

  # dataset = df_red_rar_fecesaat;  variable = "aatstunted"
  # store the result
  results    = phyloseq::estimate_richness( dataset )
  d          = phyloseq::sample_data( dataset )

  # automatically detect the values of the variable
  values     = as.character( unlist( unique( d[, variable ]) ) )
  message( paste( paste(values, collapse = ", "), "are the", length(values), "values of:", variable))

  if (length(values) > 4 ) warning("code only works for up to 4 comparisons")
  if (length(values) < 2 )
  {
    warning(paste("Wilcox test has only one unique value:", values))

    out = data.frame(
      dataset    = deparse(substitute( dataset))
    , variable   = variable
    , comparison = "Only 1 value of variable."
    , Shannon   = NA
    , Chao1     = NA
    , Observed  = NA
    , InvSimpson= NA
    , anova = NA )

    return( out )

  }

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

if ( length(values) == 3 ){
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
if( length(values) == 4 ){
  val_tre    = results[d[, variable] == values[3],]
  val_for    = results[d[, variable] == values[4],]

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

  # do the wilcox tests
  s24 = wilcox.test(val_two$ Shannon   , val_for$ Shannon) $p.value
  c24 = wilcox.test(val_two$ Chao1     , val_for$ Chao1)  $p.value
  o24 = wilcox.test(val_two$ Observed  , val_for$ Observed) $p.value
  i24 = wilcox.test(val_two$ InvSimpson, val_for$ InvSimpson)$p.value

  # do the wilcox tests
  s14 = wilcox.test(val_one$ Shannon   , val_for$ Shannon) $p.value
  c14 = wilcox.test(val_one$ Chao1     , val_for$ Chao1)  $p.value
  o14 = wilcox.test(val_one$ Observed  , val_for$ Observed) $p.value
  i14 = wilcox.test(val_one$ InvSimpson, val_for$ InvSimpson)$p.value

  # do the wilcox tests
  s12 = wilcox.test(val_one$ Shannon   , val_two$ Shannon) $p.value
  c12 = wilcox.test(val_one$ Chao1     , val_two$ Chao1)  $p.value
  o12 = wilcox.test(val_one$ Observed  , val_two$ Observed) $p.value
  i12 = wilcox.test(val_one$ InvSimpson, val_two$ InvSimpson)$p.value

  # do the wilcox tests
  s34 = wilcox.test(val_tre$ Shannon   , val_for$ Shannon) $p.value
  c34 = wilcox.test(val_tre$ Chao1     , val_for$ Chao1)  $p.value
  o34 = wilcox.test(val_tre$ Observed  , val_for$ Observed) $p.value
  i34 = wilcox.test(val_tre$ InvSimpson, val_for$ InvSimpson)$p.value

  out = data.frame(
      dataset    = deparse(substitute( dataset))
    , variable   = variable
    , comparison = c(
       paste(values[1], "and", values[2])
      ,paste(values[1], "and", values[3])
      ,paste(values[2], "and", values[3])
      ,paste(values[2], "and", values[4])
      ,paste(values[3], "and", values[4])
      ,paste(values[1], "and", values[4])
    )
    , Shannon   = signif( c( s12, s13, s23, s24, s34, s14 ), 3)
    , Chao1     = signif( c( c12, c13, c23, c24, c34, c14 ), 3)
    , Observed  = signif( c( o12, o13, o23, o24, o34, o14 ), 3)
    , InvSimpson= signif( c( i12, i13, i23, i24, i34, i14 ), 3)
  )

  }# end four values

  if (length(values) > 1){ #

  my_data = data.frame( Shannon = results$Shannon, var = d[, variable] )
  res.aov = aov( Shannon ~ ., data = my_data )
    p.val = signif( summary(res.aov)[[1]][1, 'Pr(>F)'], 3)
  print( summary(res.aov) )

  #-- updating output object
  out$ anova   = p.val
  if (!is.null( dataset_name) ) out$ dataset = dataset_name
  #-- saves a CSV file
  if (!is.null( csv_name ) ) write.csv( x = out, file = csv_name )

  return( out  )
  }
}

richness_pdf <- function( dataset
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

  if( save_pdf ) pdf( paste0( pdfname ), width = width, height = height)
  p = phyloseq::plot_richness( dataset
            , x = xlab
            , measures = c("Observed", "Chao1", "Shannon", "InvSimpson")
            , title = title )
  p = p +
    geom_boxplot(data = p$data
                 , aes_string(x = xlab, y = "value", color = "NULL"), alpha = 0.1) +
    # theme_minimal() +
    theme(axis.text = element_text(size = 18)
       , axis.title = element_text(size = 18, face = "bold")
       , plot.title = element_text(size = 22, face = "bold"
          , margin = margin(t = 0, r = 0, b = 20, l = 0))
          , strip.text.x = element_text(size = 18, face = "bold" )) +
    xlab( xlab )

 if( save_pdf ){
   print(p)
   dev.off()
 }# if pdf
  return( p )
}# end


richness_table <- function( dataset, filename = "DataAlphasampletype.txt", write_table = FALSE ){
  #' writes a table of the richness measures
  #' @param dataset phyloseq dataset
  #' @param filename name of text file to save
  #' @param write_table if TRUE, writes a table

  resultsalpha = phyloseq::estimate_richness( dataset
                      , measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
  DataAlpha    = phyloseq::sample_data(resultsalpha)
  DataAlpha$SampleID = rownames(DataAlpha)

  if( write_table) write.csv(DataAlpha, file = filename, row.names = FALSE)
  return( DataAlpha )
}


exclude_missing <- function( dataset, variable ){
  #' removes NA, "" values and <NA> values from a phyloseq dataset
  #' @param dataset phyloseq dataset
  #' @param variable variable name
  #' @return subsetted dataset # dataset = df; variable = "pays"
  values  = c(unlist(phyloseq::sample_data( dataset)[ , variable ]))
  bads    = is.na(values)
  phyloseq::sample_data(dataset) = subset( phyloseq::sample_data( dataset), !bads )

  values = c(unlist(phyloseq::sample_data( dataset)[ , variable ]))
  bads   = (values == "")
  phyloseq::sample_data(dataset) = subset( phyloseq::sample_data( dataset), !bads )
  # print( unique( phyloseq::sample_data(dataset)[, variable ]) )

  return(dataset)
}

factor_no_yes <- function( variable ){
  #' turns French two levels into English two levels
  #' @param variable variable of interest
  variable = factor( variable, levels = c("Non", "Oui"), labels = c("No", "Yes") )
  return( variable )
}

identify_taxa <- function( ASVs, taxa ){
  #' easily selects taxa from vector of ASVs
  #' @param ASVs vector of ASVs
  #' @param taxa database with ASV, Kingdom, Class, Genus, etc
# ASVs = important_asvs
  these = which( rownames(taxa) %in% ASVs )
  tt = data.frame(ASV = NA, Phylum = NA, Class = NA
              , Order = NA, Family = NA, Genus = NA, Species = NA)

  if ( length( these) > 0 ){
  tt = taxa[these, c("Phylum","Class","Order","Family","Genus","Species")] %>%
    data.frame() %>%
    tibble::rownames_to_column( var = "ASV") %>%
    data.frame()
  }

  return(tt)
}

richness_functions <- function( dataset, variable
                                , write_table = FALSE
                                , save_pdf = FALSE
                                , pdfname = "test.pdf"
                                , filename = "test.txt"
                                , csv_name = NULL
                                , width = 17, height = 8
                                , title = "Alpha Diversity"
                                , dataset_name = NULL ){
  #' wrapper for three richness functions
  #' @param dataset dataset
  #' @param variable variable of interest
  #' @param pdfname name of PDF
  #' @param title title of graph
  #' @param write_table if TRUE, writes a table
  #' @param save_pdf if TRUE, saves the pdf
  #' @param filename name of txt file (write_table)
  #' @param csv_name if not NULL, then writes csv_name for Wilcox p-values
  #' @param dataset_name name for dataset storage
# dataset = df_red_rar_duodenal; variable = "sibo";
# filtering out missing and undesired quantities
  if ( nrow( phyloseq::sample_data( dataset )) < 1) return(FALSE) else {


  dataset = exclude_missing( dataset, variable )

out_graph = richness_pdf( dataset = dataset
              , title = title
              , xlab = variable
              , save_pdf = save_pdf
              , pdfname = pdfname
              , width = width
              , height = height )

out_table = richness_table( dataset = dataset
                , filename = filename
                , write_table = write_table )

out_wilcox = richness_wilcox( dataset = dataset, variable = variable, csv_name = csv_name
                              , dataset_name = dataset_name )

return( list(
  out_graph = out_graph, out_wilcox = out_wilcox)
)
}# end non-zero dataset
}


richness_loop <- function( dataset, variables, dataset_name
                           , save_pdf = TRUE, write_table = TRUE ){
  #' wrapper function for richness_functions to enable loops
  #' @param dataset phyloseq dataset
  #' @param variables vector of variables to loop through
  #' @param dataset_name relevant name of dataset (e.g. fecal_Mada)
  #' @param save_pdf if TRUE, saves PDF
  #' @param write_table if TRUE, saves table

  outs   = list()
  wilcox = data.frame()

  for( j in seq_len( length(variables )) ){
  outs[[j]] = richness_functions( dataset
      , variable = variables[j]
      , pdfname = paste0("alpha_diversity_", variables[j],"_", dataset_name, ".pdf")
      , filename = paste0("alpha_diversity_", variables[j],"_", dataset_name, ".csv")
      , title = paste0("Alpha Diversity according to ", variables[j], ", (", dataset_name, ")" )
      , csv_name = paste0("alpha_diversity_wilcox_", variables[j], "_", dataset_name, ".csv")
      , save_pdf = save_pdf, write_table = write_table
      , dataset_name = dataset_name
     )

  #-- storing to wilcox
  wilcox    = rbind.data.frame( wilcox, outs[[j]]$out_wilcox)
  outs[[j]]$out_wilcox = NULL
  }# end loop

  wilcox$dataset = dataset_name

    return( list(outs = outs, wilcox = wilcox ))
}


#### phyloseq functions ####

filter_transform <- function( dataset, reverse = FALSE ){
  #' does filter_taxa and transform_sample_counts
  #' @param dataset phyloseq dataset
  #' @param reverse if TRUE does transform, then filter
  #' otherwise filters then transforms
  #' @return dataset

  if ( reverse ){
    dataset = phyloseq::transform_sample_counts( dataset, fun = function(x) x * 100 / sum(x, na.rm = TRUE)   )
    dataset = phyloseq::filter_taxa( dataset, flist = function(x) mean(x, na.rm = TRUE) > 0.1, prune = TRUE )
  } else {
    dataset = phyloseq::filter_taxa( dataset, flist = function(x) var(x, na.rm = TRUE) > 0, prune = TRUE )
    dataset = phyloseq::transform_sample_counts( dataset, fun = function(x) x * 100 / sum(x, na.rm = TRUE)   )
  }

  return( dataset )
}

initial_filter <- function( df , depth = 5000 ){
  #' initial filtering criteria for taxa
  #'@param df phyloseq object
  #'@param depth minimum sequencing depth per sample
  #'@return database without Mitochondria, Chloroplasts, empty taxa

  df2 = phyloseq::prune_samples( phyloseq::sample_sums( phyloseq::sample_data(df) >= depth ), x =
                                   phyloseq::sample_data(df) )
  df2 = df2 %>% phyloseq::subset_taxa(   Rank5 != "Mitochondria"
                           & Rank4 != "Chloroplast"
                           & Rank1 != "Unassigned") %>%
    phyloseq::filter_taxa( flist =  function(x) mean(x) > 0, prune = TRUE)

  return( df2 )
}

change_factor_level <- function( df, variable, values, replacements, newvariable = NULL ){
      #' changes the values of a variable within a phyloseq object
      #' @param df phyloseq dataset
      #' @param variable variable of interest
      #' @param values original values
      #' @param replacements replacement values
      #' @param newvariable name of new variable
      #' @return phyloseq object

    if( is.null( newvariable )) newvariable = variable
    stopifnot( length(values) == length(replacements) )

    for (j in seq_len( length(values )) ){  # j = 1
    phyloseq::sample_data( df )[  , newvariable][  phyloseq::sample_data( df )[, variable] == values[j] ] = replacements[j]
    }

return(df)
}# end function

prune_singlets = function(x) {
  x[x <= 1] = 0
  return(x)
}

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

  glom       = speedyseq::tax_glom( dataset , level )
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
                                  , csvname = "test.csv", perms = 999
                                  , strata = NULL ){
  #' cleans up the dataset and does the adonis2 analysis
  #' @param dataset phyloseq dataset
  #' @param controls controls for use in the model
  #' @param save_pdf if TRUe saves pdf
  #' @param write_csv if TRUE saves csv
  #' @param pdfname pdfname
  #' @param csvname name of csv file of adonis2 results
  #' @param perms permutations for adonis2
# dataset = df_feces
# origdf = df; df <- origdf
  # are the required columns in here?
goods = !(controls %in% colnames( phyloseq::sample_data( dataset )))
if( sum(goods) > 0 ) stop("not all controls are defined in this dataset")

df           = subset_cleanup(dataset, controls = controls )
project_bray = phyloseq::distance( df,  method = "bray", type = "samples")
sample_df    = data.frame( phyloseq::sample_data(df) )
form         = formula( c("project_bray ~", paste(controls, collapse = " + " )) )
print(form)

##-- adding strata option
if (!is.null(strata)) { # strata = "pays"; perms = 199
  strat = sample_df %>% dplyr::pull( strata )
  perms = permute::how(nperm = perms, blocks = strat )

}

print("beginning results")
# now the adonis test to see if there is a signficant difference according to different variables
res.adonis = vegan::adonis2( form
                , data = sample_df, method = "bray", permutations = perms  )

out = beta_diversity_pdf( res.adonis
                          , pdfname = pdfname, csvname = csvname
                          , save_pdf = save_pdf, write_csv = write_csv )

return( out )
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

  print( results )
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


filter_prune <- function( dataset, min_samples = 0 ){
  #' filters and prunes tax samples
  #' @param dataset speedyseq::tax_glom object
  #' @param min_samples minimum number of samples (integer)
df = phyloseq::filter_taxa( physeq = dataset, flist = function(x) var(x) > 0, prune = TRUE)
df = phyloseq::prune_samples( samples = phyloseq::sample_sums( df ) > min_samples, df )
return(df)
}

gm_mean <- function(x, na.rm = TRUE){
  #' geometric mean
  #' @param x vector of numbers
  y = x[ !is.na( x ) ]
  y = exp( sum( log( y[ y > 0] ), na.rm = na.rm ) / length(y) )
  return(y)
}

max_reorder_results <- function( df, level = "Genus" ){
  #' reorders results based on maximum log2FoldChange
  #' in terms of Order, Class, Family, Genus
  #' @param df results from deseq_names function
  #' @return payscorr in different order

  if (nrow( df ) > 0 ){
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

# Genus order
x = tapply(df$log2FoldChange, df$Species, function(x) max(x))
x = sort(x, decreasing = TRUE)
df$Species = factor(as.character(df$Species), levels = names(x))


df$ taxonomy = paste0( df[, level ] , "/\n", df$ Family)
}
return( df )

}

sample_data_factor <- function( dataset, variable, levels = NULL ){
  #' turns a phyloseq::sample_data() variable into a factor
  #' @param dataset phyloseq object
  #' @param variable variable of interest
  #' @return phyloseq object

  # dataset = dffiltered4 ; variable = "SampleType"
  if ( is.null(levels) ){

  temp = phyloseq::sample_data( dataset)[, variable] %>%
    unlist() %>%
    as.factor()

  } else {

  temp = phyloseq::sample_data( dataset)[, variable] %>%
    unlist() %>%
    factor(levels = levels)
  }

  phyloseq::sample_data( dataset)[, variable] = temp

  return( dataset)
}

glom_transform_melt <- function( dataset , level = "Phylum", multiplier = 100 ){
  #' tax_glom, transform_sample_counts, psmelt a dataset
  #' @param dataset phyloseq dataset
  #' @param level Phylum, Genus, Class, etc
  #' @return dataset

  df = dataset %>%
    speedyseq::tax_glom( taxrank = level ) %>%
    phyloseq::transform_sample_counts( fun = function(x) multiplier * x / sum( x ) ) %>%
    phyloseq::psmelt() %>%
    dplyr::arrange( level )

  return(df)
}

transform_core_tax_table <- function( dataset
                                      , detection = .01/100
                                      , prevalence = 90/100
                                      , write_csv = FALSE
                                      , csv_name = "test.csv"){
  #' uses phyloseq::transform, microbiome::core, and outputs the tax_table
  #' @param dataset phyloseq dataset
  #' @param detection detection in microbiome::core()
  #' @param prevalence prevalence in microbiome::core()
  #' @param write_csv if TRUE, writes csv with 'csv_name'
  #' @param csv_name name of csv

  df = dataset %>%
    phyloseq::transform_sample_counts( fun = function(x) x / sum(x) ) %>%
    microbiome::core( detection = detection
                      , prevalence = prevalence
                      , include.lowest = TRUE ) %>%
    phyloseq::tax_table() %>%
    as.data.frame()

  print( paste(nrow(df),"unique taxa at detection"
               , 100* detection, "(%) and prevalence"
               , 100 * prevalence, "(%)"))

  if (write_csv) write.csv( df, csv_name )
  return(df)
}

#### deseq2 functions ####
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
    dplyr::filter( abs(log2FoldChange) >= fold_change
                                & padj <= sig_level )

  if (write_csv) write.csv(x = payscorr, file = paste("payscorr_reduced_", name, ".csv"))
  return(payscorr)

}

deseq_fold_change <- function( dataset, level = "Genus"
                               , factor_variables = c("sexe","Run","pays")
                               , continuous_variables = c("age","haz")
                               , variable_of_interest = "stunted"
                               , sig_level = 0.05
                               , fold_change = 1 ){

  #' doing the DESeq2 analysis
  #'
  #' @param dataset phyloseq dataset
  #' @param level Genus, Species
  #' @param factor_variables factor variables (categorical)
  #' @param continuous_variables continuous variables (age, haz)
  #' @param variable_of_interest term of interest (reduced model without)
  #' @param sig_level alpha level (0.05) for significance
  #' @param fold_change minimum log2 fold change to display
  # level = "Genus"; dataset = dffiltered4; variable_of_interest = "pays"
  ### formulas
  # factor_variables = c("pays","sexe"); continuous_variables = c("haz","age");
  # variable_of_interest = "sibo"
  dataset  = subset_cleanup( dataset = dataset
       , controls = c(factor_variables, continuous_variables, variable_of_interest) )

  fact_var = paste( "as.factor(", factor_variables, ")", collapse = "+")
  cont_var = paste( continuous_variables , collapse = "+")
  form     = as.formula( paste( "~", fact_var, "+", cont_var, "+ as.factor(", variable_of_interest, ")" ))
  red_form = as.formula( paste( "~", fact_var, "+", cont_var ))

  ### deseq commands
  df = dataset %>% speedyseq::tax_glom( taxrank = level )
  dd = df %>% phyloseq::phyloseq_to_deseq2( design = form)

  gm = apply( DESeq2::counts( dd ), MARGIN = 1, gm_mean )

  ### analysis
  effect_size = DESeq2::estimateSizeFactors(dd
                                 , geoMeans = gm, type = "poscount")
  reduced     = DESeq2::DESeq( effect_size
                                 , test = "LRT"
                                 , fitType = "parametric"
                                 , reduced = red_form)

  ### testing an effect; sig_value = 0.05
  values = DESeq2::resultsNames( reduced )
  namez = values[ grepl( pattern = variable_of_interest
                    , x = values, ignore.case = TRUE) ]

  ordered = filtered = NULL
  for( j in seq_len( length( namez )) ){
    # j = 1
  res    = DESeq2::results( reduced
                , cooksCutoff = FALSE
                , alpha = sig_level
                , pAdjustMethod = "BH"
                , name = namez[j]
                )
  res

  # view a summary of the results table with a padj value < 0.01
  print( DESeq2::summary(res, alpha = sig_level))

  # filtering the results
  # reorder the results table by adjusted p-value and remove any "NA" entries
  ordered  = rbind.data.frame( ordered,
    res[order(res$ padj, na.last = NA), ] %>% data.frame() %>% dplyr::mutate( value = namez[j] )
    )
  # fold_change = 1
  filtered = rbind.data.frame( filtered,
    ordered[ ordered$padj <= sig_level & abs(ordered$log2FoldChange) >= fold_change , ])

  }# end loop

  return( list( results_ordered = ordered, results_filtered = filtered, glom = df  ))
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

  # level = "Genus"; dataset = dffilteredgastric; variable_of_interest = "pays"
  ### formulas
  form     = as.formula( paste("~", paste( full_model,  collapse = " + ")))
  red_form = as.formula( paste("~", paste( full_model[ !(full_model == variable_of_interest)]
                                      , collapse = " + ")))

  ### deseq commands
  df = dataset %>% speedyseq::tax_glom( taxrank = level )
  dd = df %>% phyloseq::phyloseq_to_deseq2( design = form)

  gm = apply( DESeq2::counts( dd ), MARGIN = 1, gm_mean )

  ### analysis
    effect_size = DESeq2::estimateSizeFactors(dd
                                 , geoMeans = gm, type = "poscount")
    reduced     = DESeq2::DESeq( effect_size
                                 , test = "LRT"
                                 , fitType = "parametric"
                                 , reduced = red_form)

  ### testing an effect; sig_value = 0.05
    values = DESeq2::resultsNames( reduced )
    res    = DESeq2::results( reduced
                , cooksCutoff = FALSE
                , alpha = sig_level
                , pAdjustMethod = "BH"
                , name = values[ grepl( pattern = variable_of_interest
                                , x = values, ignore.case = TRUE) ]
                )

  # view a summary of the results table with a padj value < 0.01
  print( DESeq2::summary(res, alpha = sig_level))

  # filtering the results
  # reorder the results table by adjusted p-value and remove any "NA" entries
  ordered  = res[order(res$ padj, na.last = NA), ]
  filtered = ordered[ ordered$padj <= sig_level & abs(ordered$log2FoldChange) >= fold_change , ]

  return( list( results_ordered = ordered, results_filtered = filtered, glom = df  ))

}# end function


#### additional functions ####
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
               strat   = vals[j]
           ,  comp     = group
           ,  bacteria = bacteria
           ,  p.val    = wilcox.test( formula = y ~ g
               ,  data = test[ test$strata == vals[j] , ] )$p.value
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
                     , ylab = controls[j] )
        }# end loop
        return(gg)
      }

differential_plot <- function( dataset
          , title = "Significantly different"
          , subtitle
          , level = "Genus" ){
  #' creates a ggplot of log2FoldChange by taxonomy
  #' @param dataset output from deseq_names() function
  #' @param title ggplot title
  #' @return ggplot object

  level_plural = dplyr::case_when(  tolower( level ) == "genus" ~ "genera"
                                   ,tolower( level ) == "phylum" ~ "phyla"
                                   ,tolower( level ) == "species" ~ "species"
                                   ,tolower( level ) == "order" ~ "orders"
                                   ,tolower( level ) == "class" ~ "classes"
  )

gg = ggplot(data = dataset
        , aes(x = reorder(taxonomy, log2FoldChange)
            , y = log2FoldChange)) +
 geom_bar(position = "dodge", stat = "identity", color = "black"
          , width = 0.5 ) +
 coord_flip() +
 theme(axis.text.x  = element_text(size = 12 ))+
 theme(axis.text.y  = element_text(size = 12 ))+
 theme(axis.text.y  = element_text(size = 12 ))+
 theme(axis.title.y = element_text(size = 14 ))+
 theme(title = element_text(size = 16, face = "bold"))+
    labs( title = paste( title, level_plural)
          , x = "taxonomy"
          , y = "log2 fold change"
          , subtitle = subtitle
          )+
  geom_hline( yintercept = c(-1,1), lty = 3, col = 'red') +
  theme_minimal()

return(gg)
}


#### rarefy functions ####
rarefy_loop <- function( dataset, seed = 3, sample.size = 5000, loops = 50){
  #' loops over the rarefy_even_depth function and stores the results
  #' @param dataset phyloseq object
  #' @param loops number of loops
  #' @param seed random seed to set
  #' @param sample.size number of samples (with replacement)

  otu = phyloseq::otu_table( dataset, taxa_are_rows = TRUE )

  #- the first instance of rarefy (n = 1)
  set.seed( seed )
  avg = phyloseq::rarefy_even_depth( otu
          , sample.size = sample.size
          , verbose = FALSE, replace = TRUE, trimOTUs = FALSE)

  #- looping over subsequent instances of rarefy (n = 2, ..., loops)
  for ( n in 2:(loops - 1)){
    out = phyloseq::rarefy_even_depth( otu
          , sample.size = sample.size
          , verbose = FALSE, replace = TRUE
          , trimOTUs = FALSE)

    #-- updating average one observation at a time
    avg = (1 - 1/n) * avg + (1/n) * out

  }# end loop

  #- returning the entire phyloseq object
  df = phyloseq::merge_phyloseq(
      phyloseq::sample_data( dataset )
    , phyloseq::tax_table(   dataset )
    , phyloseq::otu_table(   avg, taxa_are_rows = TRUE ))

  return( df )
}

rarefy_average <- function( dataset, seed = 3, sample.size = 5000){
  #' uses the mean of the prevalence and rarefy_even_depth function
  #' @param dataset phyloseq object
  #' @param loops number of loops
  #' @param seed random seed to set
  #' @param sample.size number of samples (with replacement)

  otu  = phyloseq::otu_table( dataset, taxa_are_rows = TRUE )

  #-- matrix of column sums for easy division
  csum = matrix(
    rep( colSums( otu ), each  = nrow(otu))
         , ncol = ncol(otu), byrow = FALSE )

  #-- so that all samples have the same sequencing depth
  test  = otu / csum * (sample.size * 1.10 )

  #-- so that all samples have the required sample.size
  #-- while also being close to the original distribution
  set.seed( seed )

  test2 = phyloseq::rarefy_even_depth(
    phyloseq::otu_table(test, taxa_are_rows = TRUE)
          , sample.size = sample.size
          , verbose = FALSE, replace = TRUE)

  df = phyloseq::merge_phyloseq(
      phyloseq::sample_data( dataset )
    , phyloseq::tax_table( dataset )
    , phyloseq::otu_table(test2, taxa_are_rows = TRUE ))

  return( df )
}


#### xgboost functions #####
prepare_xgboost <- function( tiny, y = "haz_cont", proportion = .8){
#' turning categorical variables to 'one-hot' encoding
#' @param tiny initial dataset
#' @param y outcome of interest (continuous)
#' @param proportion proportion of dataset for training purpose

sparse_ = Matrix::sparse.model.matrix( ~ . -1, data = tiny )
labels  = sparse_[, y]
sparse_ = sparse_[, -which(colnames(sparse_) == y)]
samps   = sample( 1:nrow( sparse_), floor(nrow( sparse_) * proportion ) )

#-- here we choose the classification of interest
train   = xgboost::xgb.DMatrix( data = sparse_[ samps,], label = labels[ samps])
test    = xgboost::xgb.DMatrix( data = sparse_[-samps,], label = labels[-samps])

return( list(train = train, test = test, spar = sparse_, samps = samps, labels = labels ))

}# end


predict_xgboost <- function( model, oo, title = NULL ){
  # predicts and plots xgboost
  #'@param model model from xgboost
  #'@param oo output from prepare_xgboost (list with: samps, test, train, labels )
  #'@param title title for ggplot

#-- prediction using XGBOOST
pred = predict( model, oo$test)
out  = data.frame( as.matrix( oo$spar[- oo$samps, ])
                   , label = oo$labels[- oo$samps]
                   , pred = pred )

gg = ggplot(out ) +
  geom_point(aes(y = pred, x = label )) +
  theme_minimal()+
  scale_x_continuous( limits = range( c(out$pred, out$label)))+
  scale_y_continuous( limits = range( c(out$pred, out$label)))+
  labs( x = "True Values", y = "Predicted Values"
        , subtitle = "Using XGBOOST Algorithm", title = title )+
  geom_abline( slope = 1, intercept = 0, col = 'red', lty = 3)

return(gg)

}

#### RJAGS functions #####
write_script <- function( csv_out
                          , script_out
                          , sample_data_out
                          , count_table_out
                          , reps = 5
                          , base_taxa = 0
                          , mc_iter = 500
                          , em_iter = 6 ){
  #' writes the script for divnet-rs
  #' @param csv_out path for output
  #' @param script_out path for this script
  #' @param sample_data_out path of sample data
  #' @param count_table_out path of count data
  #' @param reps replicates
  #' @param base_taxa base taxa (ASV number)

  #-- ensuring an even number
  mc_iter = 2 * ceiling(mc_iter/2)
  em_iter = 2 * ceiling(em_iter/2)

script = cat(file = script_out,
'[model]
em_iter = ', em_iter,'
em_burn = ', em_iter / 2,'
mc_iter = ', mc_iter,'
mc_burn = ', mc_iter / 2,'
stepsize = 0.01
perturbation = 0.05
replicates = ', reps, '
base_taxa = ', base_taxa, '

[io]
count_table = "', paste0(count_table_out),'"
sample_data =  "', paste0(sample_data_out),'"
output = "', paste0(csv_out),'"

[misc]
random_seed = 0
  ', sep = "")
  return(TRUE)
}


prepare_data <- function( df, tax_level = "Species"
                          , variable = "char"
                          , count_table_out
                          , sample_data_out ){
  #' prepares the CSV datasets for divnet-rs
  #' @param df phyloseq object
  #' @param tax_level level for tax_glom
  #' @param variable variable of interest
  #' @param count_table_out file name
  #' @param sample_data_out file name
  lee = speedyseq::tax_glom( df , taxrank =  tax_level)
  d   = phyloseq::sample_data( df ) %>% data.frame()

  # Write the count table.
  lee %>%
    phyloseq::otu_table() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("taxa") %>%
    write.table(count_table_out,
                quote = FALSE,
                sep = ",",
                row.names = FALSE)

  # Write the sample data
  model.matrix( formula(paste0("~ ", variable )) , data = d)[, -1] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    write.table(sample_data_out,
                quote = FALSE,
                sep = ",",
                row.names = FALSE)

  return( TRUE )
}


write_jags_model <- function( vars ){
  #' writes a linear model with measurement error
  #' @param vars number of variables (integer)

vars = 1:max(vars)
yvar = "y[i] ~ dnorm( mu[i], tau)"
mvar = paste0("mu[i] <- beta0 + "
              , paste0( paste0( "beta", vars, " * x", vars,"[i]" )
                        , collapse = " + "))
xvar = as.matrix( paste0( "x", vars
                   , "[i] ~ dnorm( obs", vars
                   , "[i], pow( sd", vars, "[i], -2) )" ) )
bvar = as.matrix( paste0("beta", c(0, vars), " ~ dnorm( 0.0, 0.001) "))
svar = as.matrix( c( "tau ~ dgamma( 0.01, 0.01)", "sigma <- pow( tau, -1/2)"))

cat(
    'model \n{ '
  , '\nfor (i in 1:n){ \n'
  , yvar
  , xvar
  , mvar
  , '\n}# end loop'
  , bvar
  , svar
  , '\n}# end model'
  , sep = "\n"
  , file = "model_JAGS.bug"
)

return( TRUE )
}# end

make_list <- function( left, right ){

  # left = c("beta0","beta1")
  # right = c(0,0)
  var_list = paste0( "'", left, "' = ",
                     paste0("as.numeric(", right, ")") ,   collapse = "," )
  init_list = eval(parse( text = paste0( "list(", var_list, ")"   )  ))

  return( init_list )
}

fit_jags_model <- function( datalist, iters = 5000, thin = 20
                            , vars = 10 ){
#' fits the jags model
#' @param datalist list of model inputs (n, obs1, obs2, sd1, sd2, ... )
#' @param vars number of variables

  beta_     = paste0("beta", 0 : max(vars) )
  beta_list = paste0( "'", beta_, "' = 0",   collapse = "," )

  init.list     = eval(parse( text = paste0( "list(", beta_list, ")"   )  ))
  init.list$tau = 1

  model.fit = rjags::jags.model(file = "model_JAGS.bug"
              , data = datalist
              , inits = init.list
              , n.chains = 1)

model.samples = rjags::coda.samples( model.fit
      , variable.names = c( paste0("beta", 0 : max(vars) ), "sigma")
      , n.iter = iters
      , thin = thin )

ww = data.frame(model.samples[[1]])

return(ww)
}


##### divnet functions #####
read_divnet_output <- function( csv_out ){
  #' Read data file.
  #' @param csv_out location

  # Replicate 0 is actually the estimates for the real data.
  divnet_rs      <- read.table( csv_out, sep = ",", header = TRUE)
  rep0           <- divnet_rs[divnet_rs$replicate == 0, -1]
  rownames(rep0) <- rep0$sample
  rep0$sample    <- NULL

  return( list( divnet_rs = divnet_rs, rep0  = rep0 ))
}


calculate_value <- function( divnet_rs, rep0
                             , reps = 5, f = "DivNet::shannon_true"
                             , variable_name = "shannon"){
  #' Get mean and variance for data and replicates
  #' @param divnet_rs replicates
  #' @param rep0 divnet rep 0 dataset
  #' @param reps number of replicates
  #' @param f function to apply

  rep0_ = apply(rep0, MARGIN = 1, FUN = eval( str2lang( f ) ) )

  # Now calculate the shannon index for the replicates.
  reps_ = sapply(X = 1: reps, FUN = function (i) {
    d            = divnet_rs[divnet_rs$replicate == i, -1]
    rownames(d)  = d$sample
    d$sample     = NULL
    ou = apply(d, MARGIN = 1, eval( str2lang( f ) ) )

    return(ou)
    }
  )

  # What we want is the variance in the diversity estimates for the replicates.
  reps_error <- t(apply(reps_, 1, function (x) {
    c(var(x), sd(x))
  } ))
  colnames(reps_error) <- c("variance", "sd")

  # joining the datasets
  out = tibble::tibble(names = names(rep0_), value = rep0_) %>%
    dplyr::left_join(
      reps_error %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "names") ) %>%
    dplyr::mutate(
      ciupper = value + 1.96 * sd,
      cilower = value - 1.96 * sd )

  #-- naming
  colnames(out)[ colnames(out) == "value"] = variable_name

  return(out)
}


run_divnet_script <- function( script_out , load_module = FALSE ){
  #' running the script in command line
  #' @param script_out location of script
  #' @param load_module if TRUE, loads OpenBLAS
  if (load_module) system("module load OpenBLAS/0.3.1-GCC-7.3.0-2.30")
  system(paste("divnet-rs", script_out))
  return( TRUE )
}


calculate_bray_curtis_var <- function( divnet_rs ){
  #' calculate bray-curtis distance and variance across replicates
  bcc = NULL # storage
  idz = unique( divnet_rs$sample) # assuming ids are all in the same order
  for( rep in unique(divnet_rs$replicate)){ # rep = 0
    bc = divnet_rs %>%
      dplyr::filter( replicate == rep) %>%
      dplyr::select( dplyr::contains("ASV")) %>%
      dplyr::mutate( dplyr::across( .fns = as.numeric )) %>%
      as.matrix() %>%
      DivNet::bray_curtis_true( )

    colnames(bc) = idz # assuming ids are in this order

    #-- turning into a long dataset
    bc2 = bc %>% data.frame() %>%
      dplyr::mutate( Sample1 = idz ) %>%
      tidyr::pivot_longer( cols = -Sample1, names_to = "Sample2" ) %>%
      dplyr::mutate( replicate = rep )

    bcc = rbind.data.frame(bcc, bc2 )
    print( paste("finished replicate:", rep))
  }# end replicate loop
  return(bcc)
}# end function

calculate_bray <- function( divnet_rs,  df , variable = c("pays","sex"), bcc = NULL ){
  #' calculates bray-curtis distance on samples
  #' @param divnet_rs
  #' @Param rep0
  #' @note borrowing heavily from DivNet::simplifyBeta code

# lee
  if ( is.null(bcc)){
    bcc = calculate_bray_curtis_var( divnet_rs )
  }# end bcc not provided

  samp = df@sam_data %>%
    data.frame() %>%
    tibble::rownames_to_column( var = "id") %>%
    dplyr::select( dplyr::all_of( c("id", variable )))

  #-- variables
  vars  = NULL
  p     =  length(variable)
  #-- list of the variables

  varss = samp %>%
    dplyr::mutate(
      dplyr::across( .fns = as.character ))

  #-- identifying variables used (storing all variables in one column)
  vars$covar = apply( data.frame(varss[, -1]), MARGIN = 1, paste0, collapse = ";" )
  vars$ id    = varss$id
  vars = data.frame( id = vars$id, covar = vars$covar)

  #--- removing '0' values of Bray-Curtis distance
  bc = bcc %>%
    dplyr::filter( value > 1e-14) %>%
    dplyr::arrange( Sample1, Sample2 ) %>%
    dplyr::group_by( Sample1, Sample2) %>%
    dplyr::summarize(  beta_est = mean(value[ replicate == 0 ])
                       , beta_var = var( value[ replicate != 0 ]) ) %>%
    #-- adding covariate info for Sample1, Sample2
    dplyr::mutate(     Covar1 = vars[ vars$id %in% Sample1, "covar"]
                     , Covar2 = vars[ vars$id %in% Sample2, "covar"] ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Covar1, Covar2, beta_est, beta_var) %>%
    unique() %>%
    dplyr::distinct(beta_est, .keep_all = TRUE) %>%
    dplyr::mutate(    lower = pmax(0, beta_est - 1.96 * sqrt(beta_var))
                    , upper = pmax(0, beta_est + 1.96 * sqrt(beta_var))
    )

  #--- setting a reference level for contrasts
  bases = unlist( data.frame(unique( varss[,-1]))[1, ])

  #--- specifying a contrast between variables
  bc2 = bc %>%
    tidyr::separate( col = Covar1, sep = ";", into = paste0("x", 1:p)) %>%
    tidyr::separate( col = Covar2, sep = ";", into = paste0("y", 1:p)) %>%
    dplyr::mutate( dplyr::across( .cols = c( paste0("x", 1:p)
                                           , paste0("y", 1:p))
                                  , .names = "_{.col}"
                                  , .fns = function(x)
          ifelse( trimws(as.character(x)) %in% as.character(bases), 1, 0 ) ))

  #-- comparing contrasts
  for ( j in 1:p ){
  bc2[, colnames(varss)[-1][j] ] = bc2[, paste0("_x", j)] - bc2[, paste0("_y", j)]
  }

  #-- showing the baseline levels
  print(bases)

  return( list( bray_curtis = bc2, bcc = bcc ))
}


#### breakaway functions ####
wrapper_betta <- function(bray, variables ){
  #' input from calculate_bray function

  #--breakaway analysis
  bw = breakaway::betta( chats  = bray[ , "beta_est"]
                         , ses  = bray[ , "beta_var"] %>% sqrt()
                         , X    =  cbind( intercept = 1, bray[, variables ])
  ) $ table

  return(bw)
}

