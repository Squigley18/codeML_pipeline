# codeML_pipeline
CodeML pipeline used in PhD thesis: Genome architecture and  rearrangement: Interplay of  genomic, transcriptomic, and  epigenetic features in  mammals. 

Following steps outlined in the steps outlined in https://academic.oup.com/mbe/article/40/4/msad041/7140562 

# Folders in the repository 

| FILE / FOLDER | PURPOSE |
| :--- | :--- |
| `head.stats.table` | Headers for final stats.table output
| `selection.pipeline.wrapper.sh` | Wrapper script for entire codeML pipeline
| `species.tree` | Species tree used in my analysis - newick format
| `Selection` | Files and scripts for selection.pipeline.wrapper.sh
| `codeml.environment.yml` | Conda environment yml file


# Files and scripts in selection folder 

| FILE / FOLDER | PURPOSE |
| :--- | :--- |
| `codeML_branch_ctrl` | codeML branch control file from https://github.com/abacus-gene/paml-tutorial
| `codeML_M0_ctrl` | codeML null control file from https://github.com/abacus-gene/paml-tutorial
| `FASTAtoPHYL.pl` | FASTA to phyllip format conversion perl script from https://github.com/abacus-gene/paml-tutorial
| `filter.fasta.R` | Script to extract FASTA of given transcript ID from CDS gtf file 
| `one_line_fasta.pl` | FASTA to one line perl script from https://github.com/abacus-gene/paml-tutorial
| `pal2nal.pl` | alignment perl script from https://github.com/abacus-gene/paml-tutorial
| `translatorX.pl` | alignment using translatorX perl script from https://github.com/abacus-gene/paml-tutorial
| `pipeline.sh` | full codeML pipeline called within selection.pipeline.wrapper.sh
| `species.tree` | Species tree used in my analysis - newick format
| `stats.calc.R` | R script to calculate P-value and difference between Lnl values
