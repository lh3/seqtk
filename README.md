Introduction
------------

Seqtk is a fast and lightweight tool for processing sequences in the FASTA or
FASTQ format. It seamlessly parses both FASTA and FASTQ files which can also be
optionally compressed by gzip.

Tabtk is a tool for processing TAB/SPACE-delimited data (examples at the end of
this page).

Seqtk Examples
--------------

* Convert FASTQ to FASTA:

        seqtk seq -a in.fq.gz > out.fa

* Convert ILLUMINA 1.3+ FASTQ to FASTA and mask bases with quality lower than 20 to lowercases (the 1st command line) or to `N` (the 2nd):

        seqtk seq -aQ64 -q20 in.fq > out.fa
        seqtk seq -aQ64 -q20 -n N in.fq > out.fa

* Fold long FASTA/Q lines and remove FASTA/Q comments:

        seqtk seq -Cl60 in.fa > out.fa

* Convert multi-line FASTQ to 4-line FASTQ:

        seqtk seq -l0 in.fq > out.fq

* Reverse complement FASTA/Q:

        seqtk seq -r in.fq > out.fq

* Extract sequences with names in file `name.lst`, one sequence name per line:

        seqtk subseq in.fq name.lst > out.fq

* Extract sequences in regions contained in file `reg.bed`:

        seqtk subseq in.fa reg.bed > out.fa

* Mask regions in `reg.bed` to lowercases:

        seqtk seq -M reg.bed in.fa > out.fa

* Subsample 10000 read pairs from two large paired FASTQ files (remember to use the same random seed to keep pairing):

        seqtk sample -s100 read1.fq 10000 > sub1.fq
        seqtk sample -s100 read2.fq 10000 > sub2.fq

* Trim low-quality bases from both ends using the Phred algorithm:

        seqtk trimfq in.fq > out.fq

* Trim 5bp from the left end of each read and 10bp from the right end:

        seqtk trimfq -b 5 -e 10 in.fa > out.fa

Tabtk Examples
--------------

* Basic Unix `cut` (duplicated columns ignored):

        tabtk cut -f 5,1-3,6,6- file.txt

* Reorder columns:

        tabtk cut -rf 5,1-3,6 file.txt

* Duplicate columns (duplicated columns not ignored with option `-r`):

        tabtk cut -rf 1,1,1 file.txt

* Use both SPACE and TAB as the delimitor:

        tabtk cut -d isspace -f 1-3 file.txt

