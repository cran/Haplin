## ----setup,include=FALSE------------------------------------------------------
knitr::opts_chunk$set( echo = TRUE )
library( Haplin, quietly = TRUE )

## -----------------------------------------------------------------------------
dir.exmpl <- system.file( "extdata", package = "Haplin" )
exemplary.file1 <- paste0( dir.exmpl, "/HAPLIN.trialdata.txt" )

my.gen.data.haplin <- genDataRead( file.in = exemplary.file1, file.out = "trial_data1",
	dir.out = ".", format = "haplin", n.vars = 0 )

exemplary.file3 <- paste0( dir.exmpl, "/exmpl_data.ped" )
my.gen.data <- genDataRead( exemplary.file3, file.out = "ped_data", dir.out = ".",
	format = "ped" )

## ----eval=FALSE---------------------------------------------------------------
#  ?genDataRead
#  
#  example( genDataRead )

## ----add_cov_data_read--------------------------------------------------------
add.cov.file <- paste0( dir.exmpl, "/add_cov_data2.dat" )
my.gen.data.haplin3 <- genDataRead( file.in = exemplary.file1, file.out = "trial_data3",
	dir.out = ".", format = "haplin", n.vars = 0, cov.file.in = add.cov.file )
my.gen.data.haplin3

add.cov.file2 <- paste0( dir.exmpl, "/add_cov_data.dat" )
my.gen.data2 <- genDataRead( exemplary.file3, file.out = "ped_data2", dir.out = ".",
	format = "ped", cov.file.in = add.cov.file2 )
my.gen.data2

## ----show_summary-------------------------------------------------------------
my.gen.data

## ----showPheno_demo-----------------------------------------------------------
# by default - showing first 5 entries:
showPheno( my.gen.data )
# getting all the info:
head( showPheno( my.gen.data, n = "all" ), n = 20 )
showPheno( my.gen.data, from = 4, to = 15 )
# show information about females only:
head( showPheno( my.gen.data, sex = 2, n = "all" ), n = 20 )

## ----showPheno_demo2----------------------------------------------------------
females.pheno <- showPheno( my.gen.data, sex = 2 )
head( females.pheno )

## ----nindiv_demo--------------------------------------------------------------
nindiv( my.gen.data )
nfam( my.gen.data )

## ----nsnps--------------------------------------------------------------------
nsnps( my.gen.data )

## ----showSNPnames-------------------------------------------------------------
showSNPnames( my.gen.data ) # by default - showing only first 5 SNPs
showSNPnames( my.gen.data, from = 12, to = 31 )

## ----showGen_demo-------------------------------------------------------------
showGen( my.gen.data, markers = c( 10,15,121 ) ) # by default - showing first 5 entries
showGen( my.gen.data, from = 31, to = 231 )

## ----showGen_demo2------------------------------------------------------------
subset.genes <- showGen( my.gen.data, from = 31, to = 231, markers = c( 10,15,121 ) )
subset.genes

## ----eval=FALSE---------------------------------------------------------------
#  my.prepared.gen.data <- genDataPreprocess( data.in = my.gen.data, map.file =
#    "my_gen_data.map", design = "triad", file.out = "my_prepared_gen_data",
#    dir.out = "." )

## ----eval=FALSE---------------------------------------------------------------
#  gen.data.subset <- genDataGetPart( data.in = my.gen.data, markers = c( 3:15,22 ),
#    design = "triad", file.out = "my_gen_data_subset", dir.out = "." )

## ----eval=FALSE---------------------------------------------------------------
#  my.prepared.gen.data <- genDataLoad( filename = "my_prepared_gen_data",
#    dir.in = "." )

