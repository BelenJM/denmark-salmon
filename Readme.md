# The genomic history of the Atlantic salmon in Denmark
Scripts and analysis of the study on the genomic history of the Danish population of Atlantic salmon (article in prep.)

This folder contains the automatic pipeline and postprocessing files used for the project. The file `workflow.py` contains the pipeline commands, that were in turn submitted to the SLURM backend of a computing cluster by the python package `gwf`. You can find a table with the CPU and RAM usage of each command of the pipeline in the file utilization.gwf. This pipelineproject pipeline and postprocessing files includes the following steps and calculations:

* Genomic pipeline to prepare the data: from raw reads to vcf files
* SNP files generated with ANGSD
* SNP files generated with GATK
* Population genetics using R packages adegenet, ...
* Analysis of the mithocondrial genome
* Admixture analysis
* Fst values
* 1D Site Frequency Spectrum
* 2D Site Frequency Spectrum
* Allele-specific error estimates for each individual in the analysis
* Calculation of the D-statistic (ABBA-BABA) test and application of the error correction on it

The file `modules.py` contains the template code for each of the targets called in the workflow. A target is nothing more than a specific sets of commands running together as a job - a target might depend on previous targets being completed, and might have other targets depending on it. For example, the target for the abbababa test between 4 populations is dependent on other targets calculating the allele-specific errors in these populations, and counting the ABBA and BABA combinations of interest in the bam files.

More explanations about each target are contained in the python scripts. By cloning this repository on your folder of a server, you should be able to rerun the whole analysis by simply running `gwf -b slurm -b workflow.py run`. Check the `gwf` github repository and manual to know more about this very useful and intuitive pipeline tool.

Analysis in R and python is included in the R scripts and in the jupyter notebooks in the folder .............. .
