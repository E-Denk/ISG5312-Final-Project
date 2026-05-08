This is my final project for ISG5312. 

I am reanalyzing single-cell RNASeq from Single-nucleus transcriptional and chromatin accessibility analyses of maturing mouse Achilles tendon uncover the molecular landscape of tendon stem/progenitor cells.

The main takeaway of the paper is that Cd55 and Cd248 are novel markers of TSPCs (tendon stem/progenitor cells).



Project structure:

-data

-metadata

-README.md

-results

-scripts

  -01_downloadandqc

  -02_cellranger

  -03_ranalysis

  -logs


data is where the fastq files will be downloaded to

scripts is where I will store the scripts to run Cell Ranger and the RStudio analysis

logs is for .err and .out

data and results directories are generated after running the first scripts.

A count_out directory will be made in the results as the output folder for Cell Ranger count.


Citation:
Tsutsumi, H., Chiba, T., Fujii, Y., Matsushima, T., Kimura, T., Kanai, A., Kishida, A., Suzuki, Y., & Asahara, H. (2026). Single-nucleus transcriptional and chromatin accessibility analyses of maturing mouse Achilles tendon uncover the molecular landscape of tendon stem/progenitor cells. eLife, 14. https://doi.org/10.7554/eLife.104768.3

