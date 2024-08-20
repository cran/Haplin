# some common variables and file names used in testing

# filenames - ped format
my.ped.data <- "exmpl_data.ped"
my.ped.out.data <- "ped_read"

# filenames - haplin format
my.haplin.data.no.cov <- "HAPLIN.trialdata.txt"
my.haplin.no.cov.out <- "haplin_read_no_cov"
my.haplin.data.cov <- "HAPLIN.trialdata2.txt"
my.haplin.cov.out <- "haplin_read_cov"

# filenames - preprocessed data
my.haplin.cov.preproc.out <- "haplin_prep_cov"
my.haplin.no.cov.preproc.out <- "haplin_prep_no_cov"
my.ped.preproc.out <- "ped_prep"
my.ped.preproc.cc.out <- "ped_cc_prep"
my.ped.preproc.env.out <- "ped_env_prep"
my.ped.preproc.cc.env.out <- "ped_cc_env_prep"

# additional covariate columns
add.cov.file.name <- "add_cov_data.dat"
add.cov.no.head.file.name <- "add_cov_data_no_head.dat"

if.overwrite <- TRUE

# tolerance for comparing estimates:
# (due to differences in random numbers generated on Linux vs Windows?)
tolerance_tests <- 0.1

# for running haplin, haplinStrat and haplinSlide
cc.var <- 4
# env.var <- 7

all.args <- list( 
# varying use.missing and response:
 list( design = "triad", use.missing = FALSE, maternal = TRUE, reference = "reciprocal", response = "free" ),
 list( design = "triad", use.missing = FALSE, maternal = TRUE, reference = "ref.cat", response = "mult" ),
 list( design = "triad", use.missing = TRUE, maternal = TRUE, reference = "reciprocal", response = "free" ),
 list( design = "triad", use.missing = TRUE, maternal = TRUE, reference = "ref.cat", response = "mult" ),
# varying maternal:
# no.5:
 list( design = "triad", use.missing = FALSE, maternal = FALSE, reference = "reciprocal", response = "free" ),
 list( design = "triad", use.missing = FALSE, maternal = FALSE, reference = "ref.cat", response = "mult" ),
 list( design = "triad", use.missing = TRUE, maternal = FALSE, reference = "reciprocal", response = "free" ),
 list( design = "triad", use.missing = TRUE, maternal = FALSE, reference = "ref.cat", response = "mult" ),
# varying poo:
 list( design = "triad", use.missing = FALSE, maternal = TRUE, reference = "reciprocal", response = "free", poo = TRUE ),
# no.10:
 list( design = "triad", use.missing = FALSE, maternal = TRUE, reference = "ref.cat", response = "mult", poo = TRUE ),
 list( design = "triad", use.missing = TRUE, maternal = TRUE, reference = "reciprocal", response = "free", poo = TRUE ),
 list( design = "triad", use.missing = TRUE, maternal = TRUE, reference = "ref.cat", response = "mult", poo = TRUE ),
# varying use.missing and response with design = 'cc.triad'
 list( design = "cc.triad", ccvar = cc.var, use.missing = FALSE, maternal = TRUE, reference = "reciprocal", response = "free" ),
 list( design = "cc.triad", ccvar = cc.var, use.missing = FALSE, maternal = TRUE, reference = "ref.cat", response = "mult" ),
# no.15:
 list( design = "cc.triad", ccvar = cc.var, use.missing = TRUE, maternal = TRUE, reference = "reciprocal", response = "free" ),
 list( design = "cc.triad", ccvar = cc.var, use.missing = TRUE, maternal = TRUE, reference = "ref.cat", response = "mult" ),
# varying maternal with design = 'cc.triad'
 list( design = "cc.triad", ccvar = cc.var, use.missing = FALSE, maternal = FALSE, reference = "reciprocal", response = "free" ),
 list( design = "cc.triad", ccvar = cc.var, use.missing = FALSE, maternal = FALSE, reference = "ref.cat", response = "mult" ),
 list( design = "cc.triad", ccvar = cc.var, use.missing = TRUE, maternal = FALSE, reference = "reciprocal", response = "free" ),
# no.20:
 list( design = "cc.triad", ccvar = cc.var, use.missing = TRUE, maternal = FALSE, reference = "ref.cat", response = "mult" ),
# varying use.missing and response with poo with design = 'cc.triad':
 list( design = "cc.triad", ccvar = cc.var, use.missing = FALSE, maternal = TRUE, reference = "reciprocal", response = "free", poo = TRUE ),
 list( design = "cc.triad", ccvar = cc.var, use.missing = FALSE, maternal = TRUE, reference = "ref.cat", response = "mult", poo = TRUE ),
 list( design = "cc.triad", ccvar = cc.var, use.missing = TRUE, maternal = TRUE, reference = "reciprocal", response = "free", poo = TRUE ),
 list( design = "cc.triad", ccvar = cc.var, use.missing = TRUE, maternal = TRUE, reference = "ref.cat", response = "mult", poo = TRUE ),
# varying use.missing and response with design = 'cc'
# no.25:
 list( design = "cc", ccvar = cc.var, use.missing = FALSE, maternal = FALSE, reference = "reciprocal", response = "free" ),
 list( design = "cc", ccvar = cc.var, use.missing = FALSE, maternal = FALSE, reference = "ref.cat", response = "mult" ),
 list( design = "cc", ccvar = cc.var, use.missing = TRUE, maternal = FALSE, reference = "reciprocal", response = "free" ),
 list( design = "cc", ccvar = cc.var, use.missing = TRUE, maternal = FALSE, reference = "ref.cat", response = "mult" ),
# varying use.missing and response with xchrom and comb.sex = "double"
 list( design = "triad", use.missing = FALSE, maternal = TRUE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "double" ),
# no.30:
 list( design = "triad", use.missing = FALSE, maternal = TRUE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "double" ),
 list( design = "triad", use.missing = TRUE, maternal = TRUE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "double" ),
 list( design = "triad", use.missing = TRUE, maternal = TRUE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "double" ),
# varying maternal with xchrom and comb.sex = "double"
 list( design = "triad", use.missing = FALSE, maternal = FALSE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "double" ),
 list( design = "triad", use.missing = FALSE, maternal = FALSE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "double" ),
# no.35:
 list( design = "triad", use.missing = TRUE, maternal = FALSE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "double" ),
 list( design = "triad", use.missing = TRUE, maternal = FALSE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "double" ),
# varying use.missing and response with xchrom and comb.sex = "females"
 list( design = "triad", use.missing = FALSE, maternal = TRUE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "females" ),
 list( design = "triad", use.missing = FALSE, maternal = TRUE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "females" ),
 list( design = "triad", use.missing = TRUE, maternal = TRUE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "females" ),
# no.40:
 list( design = "triad", use.missing = TRUE, maternal = TRUE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "females" ),
# varying maternal with xchrom and comb.sex = "females"
 list( design = "triad", use.missing = FALSE, maternal = FALSE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "females" ),
 list( design = "triad", use.missing = FALSE, maternal = FALSE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "females" ),
 list( design = "triad", use.missing = TRUE, maternal = FALSE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "females" ),
 list( design = "triad", use.missing = TRUE, maternal = FALSE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "females" ),
# varying use.missing and response with xchrom and comb.sex = "males"
 list( design = "triad", use.missing = FALSE, maternal = TRUE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "males" ),
 list( design = "triad", use.missing = FALSE, maternal = TRUE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "males" ),
 list( design = "triad", use.missing = TRUE, maternal = TRUE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "males" ),
 list( design = "triad", use.missing = TRUE, maternal = TRUE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "males" ),
# varying maternal with xchrom and comb.sex = "males"
 list( design = "triad", use.missing = FALSE, maternal = FALSE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "males" ),
 list( design = "triad", use.missing = FALSE, maternal = FALSE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "males" ),
 list( design = "triad", use.missing = TRUE, maternal = FALSE, reference = "reciprocal", response = "free", xchrom = T, comb.sex = "males" ),
 list( design = "triad", use.missing = TRUE, maternal = FALSE, reference = "ref.cat", response = "mult", xchrom = T, comb.sex = "males" )
)

# grouping the arguments
args.cc <- all.args[ sapply( all.args, function(x){ !is.null( x$ccvar ) } ) ]
args.cc.cc <- args.cc[ sapply( args.cc, function(x){ x$design == "cc" } ) ]
args.cc.triad <- args.cc[ sapply( args.cc, function(x){ x$design != "cc" } ) ]

args.no.cc <- all.args[ sapply( all.args, function(x){ is.null( x$ccvar ) } ) ]
args.xchrom <- args.no.cc[ sapply( args.no.cc, function(x){ !is.null( x$xchrom ) } ) ]
args.simple <- args.no.cc[ sapply( args.no.cc, function(x){ is.null( x$xchrom ) } ) ]

# names of "aux" list:
aux.list.names <- c( "info", "class", "marker.names" )
