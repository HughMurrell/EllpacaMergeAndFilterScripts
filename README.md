# EllpacaMergeAndFilterScripts

Here we store Julia scripts to merge and filter Ellpaca datasets output by PORPIDpipeline

The `PORPIDpipeline` outputs `fasta` files in directories for each pool/facility combination.
Each `fasta` file contains sequences from a `sample` from a particular `donor`.

`PORPIDpipeline` has its own directory/file/sequence naming conventions.

These scripts switch to a UCT directory/file/sequence naming conventions and then perform
an NT codon aware functional filter followed by a Merge and Collapse and finally Translate 
to AA and re-align

The scripts should be run in the following order:

```
# get rid of any previous attempt at generating the renamed collection of data
rm -rf ellpaca_collection  
# rename all Ellpaca sample files using  UCT conventions
julia Rename.jl  
# filter non-functionals from each sample file (requires ConsensusC)
julia FilterNonFunctionals.jl  # about 18 minutes on a mac
# merge sample files into donor files and collapse identical sequences
julia MergeAndCollapse.jl  
# translate donor files to AA and re-align using mafft (requires ConsensusC and HXB2)
julia TranslateAlign.jl  # about 35 minutes on a mac
```

Note that some of the scripts make use of the HXB2 reference and a ConsensusC alignment guide
both of which were generated from LANL sequences and are stored in the `pannel` folder.
The `ConsensusMaker.jl` script maybe useful for regenerating these sequences.

Downstream analysis such as the `EllpacaAtlas` project can be performed by copying the
`env_aa_aligned_hxb2` folder to the alignments folder of the Atlas project and then running
the Atlas scripts.
