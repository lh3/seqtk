Release 1.5-r133 (1 June 2025)
------------------------------

Notable changes:

 * Improvement: support chromosomes longer than 2Gb (#192), provided @c-zhou

 * New feature: added option `-R` to the seq command to output sequences on
   both strands.

 * New feature: added option `-P` to telo to output scores

(1.5: 1 June 2025, r133)



Release 1.4-r122 (19 May 2023)
------------------------------

Notable changes:

 * Improvement: faster FASTX parsing (#123)

 * New feature: added the `telo` command to output telomere regions.

 * New feature: added the `size` command to count the number of sequences and
   the number of bases. Lighter and thus faster than `comp`.

 * New feature: added the `hpc` command to compress homopolymers in input
   sequences.

 * New feature: added the `split` command to split a large input file into
   multiple smaller files.

 * New feature: added the `gap` command to output non-ACGT regions in the input
   file.

 * New feature: added option `-s` to command `subseq` to support the strand
   field in BED. For the moment, this option does not work with other subseq
   options.

(1.4: 19 May 2023, r122)
