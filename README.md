# EllpacaMergeAndFilterScripts

The repository https://github.com/HughMurrell/EllpacaMergeAndFilterScripts.git is a collection
of scripts to **rename**, **filter**, **merge and collapse**, **translate and realign** 
the Ellpaca sample fasta files output by PORPIDpipeline

The `PORPIDpipeline` outputs `fasta` files in directories for each pool/facility combination.
Each `fasta` file contains sequences from a `sample` from a particular `patient`.

`PORPIDpipeline` has its own directory/file/sequence naming conventions.

These scripts switch to a UCT directory/file/sequence naming convention and then perform
an NT codon aware functional filter followed by a Merge and Collapse and finally Translate 
to AA and re-align

Before proceeding, the ellpaca output should be collected in one directory called `data_in`.
This can be dome by copying each postproc archive from the PORPIDpipeline runs to the
`data_in` directory and unpacking them there. Below is a screenshot of the current
collection of postproc-archives.

![data_in screenshot](data_in_screenshot.png "postproc-archive collection")

The scripts should be run as follows:

1) get rid of any previous attempt at generating the renamed collection of data
```
rm -rf ellpaca_collection  
```
2) rename all Ellpaca sample files and sequences using  UCT conventions (note that some sample filenames
have been designated **suspect** and this is hardcoded in the script. This designation should be
replaced by an internal renaming mapping once the suspect nature is resolved.
```
julia Rename.jl  
```
3) filter non-functionals from each sample file (requires `ConsensusC_env_nt.fasta` from panels)
```
julia FilterNonFunctionals.jl  # about 18 minutes on a mac
```
4) merge sample files into donor files and collapse identical sequences
```
julia MergeAndCollapse.jl  
```
5) translate donor files to AA and re-align using mafft (requires `ConsensusC_env_aa.fasta` and `hxb2_env_aa.fasta`)
```
julia TranslateAlign.jl  # about 35 minutes on a mac
```
6) finally create a gzipped archive (just 55mb for the ellpaca collection)
```
tar -cvf ellpaca_collection.tar ellpaca_collection
gzip ellpaca_collection.tar
```

Note that some of the scripts make use of the `hxb2` reference and a `ConsensusC` alignment guide
both of which were generated from LANL sequences and are stored in the `pannel` folder.
The `ConsensusMaker.jl` script maybe useful for regenerating these sequences.

Downstream analysis such as the `EllpacaAtlas` project can be performed by copying the
`env_aa_aligned_hxb2` folder to the alignments folder of the Atlas project and then running
the Atlas scripts.
