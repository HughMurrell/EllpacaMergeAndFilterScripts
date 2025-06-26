# EllpacaMergeAndFilterScripts

Here we store Julia scripts to merge and filter Ellpaca datasets output by PORPIDpipeline

The PORPIDpipeline outputs `fasta` files in directories for each pool/facility combination.
Each `fasta` file contains sequences from a `sample` from a particular `donor`.

The `MergePools` script merges these data donor/sample files into one fasta files for each donor.

The `FunctionalFilter` script then filters each donor file extracting functional sequences by
comparing them to a reference sequence.

This is done in two stages: 

1) a *matching segment* is extracted from the consensus of the earliest sample 
by searching for the longest open reading frame. If that fails search for starting
and ending matches to the reference.

2) this *matching consensus segment* is used as a reference for extraction 
of matching segments from all sequences in the file. This time the matches are obtained
by searching for starting and ending matches to the ORF of the consensus sequence.

