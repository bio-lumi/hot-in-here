**Project name**: It's Getting Hot in Here: Temperature-Driven Shifts in Microbial Diversity and Community Composition
**Brief description**: This project utilises long-read PacBio 16s rRNA sequencing to investigate experimental microbial community responses to elevated temperatures. This repository contains a bash script to process demultiplexed PacBio circular consensus sequences, as well as R scripts for downstream processing, including network analysis. This project was submitted for the partial fulfilment of the MRes degree in Ecology, Evolution, and Conservation at Imperial College London (2023/2024).
**Languages**: Bash, R version 4.4.1 (Race For Your Life). 
**Dependencies**: For the bioinformatics pipeline, you will need to have already installed QIIME2 (https://qiime2.org/). All dependencies in R are listed in the R scripts. 
**Project structure and usage**: This repository contains three scripts. *loveline-PacBio.sh* contains a Bash script that you can run from the terminal to process demultiplexed PacBio circular consensus sequences using QIIME2. *gluc-net-met.R* contains an R script for downstream sequencing analysis: alpha-diversity, beta-diversity, network construction, and random network generation. *div-net-plots.R* contains an R script for plotting results from *gluc-net-met.R*. 
**Author name and contact**: Loveline Martin, loveline.ubung@gmail.com

