# this function is called after loading the library and before sealing the internal package environment
# here, set some variables that will be used by functions within the package
.onLoad <- function( libname, pkgname ){
	assign( ".haplinEnv", new.env(), envir = parent.env( environment() ) )

# 	# the chunk size that would still fit in the memory (used e.g. in readGenData.R):
# 	assign( ".nb.lines.per.chunk", 100, envir = .haplinEnv )
# 	
	# number of columns of the ffmatrix in one chunk (used e.g. in readGenData.R):
	assign( ".nb.cols.per.chunk", 10000, envir = .haplinEnv )
	
	# list of currently available study designs
	assign( ".design.list", c( "triad", "cc.triad", "cc" ), envir = .haplinEnv )
	
	# base name for the objects holding the gen.data chunks
	assign( ".gen.cols.name", "gen.cols", envir = .haplinEnv )
	assign( ".env.cols.name", "env.cols", envir = .haplinEnv )
	
	# names of columns for the covariate data of PED-formatted data
	assign( ".cov.data.colnames", c( "id.fam", "id.c", "id.f", "id.m", "sex", "cc" ), envir = .haplinEnv )
	
	# names of the elements in the haplin-type object list
	assign( ".haplin.data.names", c( "cov.data", "gen.data", "aux" ), envir = .haplinEnv )
	
	# name of the class (S3-type, only for marking the type of data in the object)
	# assign( ".all.haplin.classes", c( get( ".class.data.read.ped", envir = .haplinEnv ), "haplin.ready",	"prep.data", "haplin.data" ),  envir = .haplinEnv )
	
	# mode for creating the ff objects:
	# first, for the raw genetic data (i.e., alleles)
	assign( ".vmode.gen.data", "byte", envir = .haplinEnv )
	# then, for the frequency counts (transformed by f.data, just before main calculations)
	assign( ".vmode.freq.data", "single", envir = .haplinEnv )
	# for the pre-processed data (integers)
	assign( ".vmode.proc.data", "integer", envir = .haplinEnv )
	# and for the output of haplin:
	assign( ".vmode.hap.res.dbl", "double", envir = .haplinEnv )
	assign( ".vmode.hap.res.int", "integer", envir = .haplinEnv )
	
	# possible parallel mechanisms:
	assign( ".para.env", c( "parallel", "Rmpi" ), envir = .haplinEnv )
	assign( ".default.para.env", "parallel", envir = .haplinEnv )
	
	# column names with the covariate data:
	# sex
	assign( "sex.column", "sex.c", envir = .haplinEnv )
	assign( "cc.column", "cc.c", envir = .haplinEnv )
	
	# flag indicating whether R slaves are spawned by Rmpi
	assign( "Rmpi.spawned", FALSE, envir = .haplinEnv )

	# set the ff finalizer to delete the files after an ff object is removed
	options( fffinalizer = "delete" )
	options( fffinonexit = TRUE )
}
