# Haplin 7.1.0

(uploaded 18 April, 2019)

* New helper functions by Julia Romanowska, to ease data handling: `showPheno`, `showGen`, `nsnps`, `nindiv`, and `nfam`.
* New function by Miriam Gjerdevik, `hapRelEff`, with its own vignette, to compute relative efficiencies (sample sizes) for various designs.
* Miscellaneous minor fixes.

# Haplin 7.0.0

(uploaded May 25, 2018)

* Unfortunately, the very popular `GenABEL` package was taken down from CRAN on 24 May 2018 because it’s not being maintained any longer. Haplin relied for a long time on GenABEL to store and handle GWAS data. Thanks mostly to the programming of Julia Romanowska, Haplin Version 7.1.0 now adapts its data storage to the `ff package`. There are thus many updates to the data handling in Haplin, in particular the new functions `genDataRead` and `genDataPreprocess`. For more details, see the documentation.

# Haplin 6.2.1

(uploaded Oct. 18, 2017)

* There is a bug in 6.2.0 that in some instances may cause strata to be given incorrect names in haptables from `haplinSlide` with the strata argument specified. Version 6.2.1 fixes that.

# Haplin 6.2.0

(uploaded Feb 27, 2017)

* The function `hapPowerAsymp` has been added (thanks mostly to Miriam Gjerdevik). `hapPowerAsymp` provides fast power calculations for a wide range of scenarios, without needing simulations.

# Haplin 6.0.1

(uploaded May 26, 2016)

* Fixed a small bug in `hapRun` from 6.0

# Haplin 6.0

(uploaded May 24, 2016)

* Gene-environment interactions implemented in the function gxe

* hapPower does power simulations for a range of designs, haplotype combinations, as well as gxe analyses (thanks to Miriam Gjerdevik)

# Haplin 5.5

(uploaded April 24, 2015)

* P-values for individual effects are now computed with more decimals to facilitate sorting on p-values

* Two new functions snpPower and snpSamplesize, which compute power and sample size for single SNP analyses.

* Two new functions cbindFiles and rbindFiles. cbindFiles takes a number of files and combines them column-wise (side-by-side), i.e. reads each file line by line, pastes corresponding lines, then writes to outfile. rbindFiles takes a number of files and combines them by rows, without reading the full files into memory. The purpose is to facilitate easy merging of very large text files that may be too big to read into R memory.

* Plus a number of minor fixes and cosmetic changes.

# Haplin 5.3

(uploaded May 22, 2013)

* haplin and haplinSlide will now use actual snp names taken from the map file rather than just number them in sequence.

* Plus a number of minor fixes and cosmetic changes.

# Haplin 5.2

(uploaded March 18, 2013)

* Two new functions, lineByLine and convertPed, have been added. The purpose is to facilitate easy conversion/recoding of large text files that may be too big to read into R memory. lineByLine can modify general text files, typically recode values, select lines, columns, etc., whereas convertPed is designed to handle ped files in particular.

# Haplin 5.0.2

(uploaded March 18, 2013)

* Minor fixes and some cosmetic changes.

# Haplin 5.0.1

(uploaded March 8, 2013)

* Fixed a bug (only present in 5.0) that caused haplin to crash if using GWAS data through GenABEL with use.missing = T.

* In haplinSlide, the default of the argument table.output is now set to TRUE, to reduce the risk of running out of memory when running large analyses.

* In haplinSlide, the argument slaveOutfile can now be used also when cpus = 1.

# Haplin 5.0

(uploaded February 7, 2013)

* Haplin now allows parent-of-origin effects through the argument 'poo = T'. The estimated fetal effects will then be split into the effects of haplotypes inherited from the mother and from the father. In addition, a ratio of the two relative risks thus obtained is also printed, with a corresponding p-value testing whether the two effects differ. Parent-of-origin effects can be combined with estimating maternal haplotype effects through 'maternal = T'. Parent-of-origin effects can also be computed on the X-chromosome.

* Previous versions of Haplin did not allow X-chromosome effects to be estimated in the pure case-control design. This is now implemented.

* Plus a number of minor fixes and tidy-ups.

# Haplin 4.1

(uploaded March 15, 2012)

* The function haplinSlide now provides feedback as it runs. A new argument, slaveOutfile, can be used to redirect the feedback to a file.

* A minor bug in the GWAS data concaused the conto crash for some types of family structures. This has now been fixed.

# Haplin 4.0

(uploaded Jan 26, 2012)

* The main new feature is that Haplin now provides full support for GWAS data (through the GenABEL library data structure). Parallel processing is integrated in the sliding window estimation and is activated by simply choosing the number of cpus available (implemented using the snowfall library).

* A new argument comb.sex to be used with X-chromosome analyses. It allows analyses of males and females separately, and different ways of combining males and females in a joint analysis, corresponding to dose-response models or X-inactivation.

# Haplin 3.5

(uploaded May 26, 2010)

* A simple GUI for generating syntax available at haplin.fhi.no. NOTE: implements most but not all current Haplin features.

* Haplin now accepts data from from the X chromosome, using the arguments 'xchrom = T' (and 'sex' to specify column for sex variable).

* Using the 'data.out' argument, haplin can return a processed of the original data file, with all haplotypes estimated/imputed for all families.

* A new function 'haplinSlide' runs haplin on a consecutive series of overlapping windows, useful to scan genes with many SNPs.

* A new function 'haptable' creates a useful summary table of all the important estimation results from a Haplin run.

* A new function 'haplinTDT' that performs various TDT tests for comparison. Written by Øivind Skare.

* A new function 'output' which saves tables and figures from a Haplin run to disk.

# Haplin 3.0.2

(uploaded Jan 11, 2010)

* Changes in R release 2.10 relative to 2.9 caused Haplin to crash when n.vars > 0. This has now been fixed.

* Some minor fixes in documentation etc.

# Haplin 3.0.1

(uploaded May 27, 2009)

* Haplin prepared for upload to the CRAN repository

* Only minor fixes in documentation etc.

# Haplin 3.0

(uploaded April 29, 2009)

This is an extensive update, some of the most important changes are:

* Implemented as a proper R library. S-Plus terminated.

* Haplin now handles the combined case-parent triad + control-parent triad design, and also the standard case-control design.

* Updated output

* Can choose a multiplicative dose-response model if estimation of separate single- and double dose effects is not required

* Killed a number of minor bugs.

# Haplin 2.1.1

(uploaded June 12, 2006)

* Minor cosmetic changes in printout and in graph appearance.

* Got rid of (most of?) the annoying graphics warnings that appeared when plotting result

# Haplin 2.1

(uploaded May 29, 2006)

* Haplin now prints tests of Hardy-Weinberg equilibrium for each marker separately.

* When using a reference category (not the default reciprocal reference) Haplin 2.0 forgot to print the double dose estimate for the reference category, even though it was computed and plotted in the plot. This has now been fixed. (Thanks to S. Hussain for pointing this out).

# Haplin 2.0

(uploaded Jan 13, 2006)

* A change in the format function in the very latest R release (2.2.0, Oct 06, 2005) caused the R of Haplin 1.0 to exit prematurely. This has been fixed in Haplin 2.0.

* This (finally!) has an option to run with triads containing missing data through the EM algorithm, by setting use.missing = T.

* A computation of the standard errors and individual p-values that corrects for missing data and ambiguous haplotypes. The consequence is that the jackknife resampling used in 1.0 is no longer necessary, and only one run of Haplin is needed. Thus, much less patience is required of the user....

* The initial data handling and sorting is considerably improved, with the result that most computations should be faster, and Haplin can now handle more SNPs in a run. The exact number of SNPs possible to run depends a good deal on the data and population structure, but 6-7 SNPs should usually run fine.

# Haplin 1.0

(uploaded June 1, 2005)

* Haplin now uses "reciprocal" as the default reference category. This is slightly different from the previous "population" (which still can be used by setting reference = "population"). See the document on parametrization and reference for more details

* A new "markers" argument that allows the user to select the markers to be used in the analysis. The user no longer has to produce separate data files for each marker selection. See the document on Haplin data format for more details

* Much improved handling of rare haplotypes. In the first (non-numbered) this selection was rather rough. In 1.0 the haplotype frequencies are estimated fairly precisely before removing the rare haplotypes, leading to less loss of data and more precise estimates. The argument "threshold" that determines the frequency limit for retaining haplotypes is now followed much more closely.

* Haplin now prints p-values for all specific effects, and an overall p-value for the locus

* Haplin now reports the amount of  data removed due to missing genotypes and which families that are removed due to Mendelian inconsistencies

* Jackknife resampling procedure (by setting resampling = "jackknife") adjusts the standard error to compensate for the missing information due to unknown haplotypes


