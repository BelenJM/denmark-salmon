###workflow gwf
from gwf import Workflow
import os, sys
import numpy as np
from modules import *

gwf = Workflow()

ancFile = '../steps/admixture/PIKEmappedToLaks.sorted.bam'

######################SAF ANALYSIS

mypath = './filelists'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]

for FILE in filelists:
    runName = FILE#os.path.split( os.path.basename(dataf) )[1]
    gwf.target_from_template( 'dosaf_'+str(FILE) , dosaf(filelist=mypath+'/'+FILE,
                                                         anc=ancFile,
                                                         outFolder='results_all',
                                                         cores=18,
                                                         memory='64g',
                                                         walltime='48:00:00'))


mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
    
for FILE in filelists:
    runName = FILE#os.path.split( os.path.basename(dataf) )[1]
    gwf.target_from_template( 'dosaf_5e6_'+str(FILE) , dosaf(filelist=mypath+'/'+FILE,
                                                         anc=ancFile,
                                                         outFolder='results_5e6',
                                                         regions='pruning/snps_5e6_formatted.tsv',
                                                         cores=18,
                                                         memory='64g',
                                                         walltime='48:00:00'))

######1D sfs    
######REMEMBER: do  `angsd sites index your.file` on the -sites file before running all this
######          so that the SNP list is indexed
######################################SFS1D analysis
###### Here you need -sites and not -rf as option (-rf is only for reading from bam files)


mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]

for i in filelists:
    gwf.target_from_template( 'sfs1d_5e6_'+str(i) , sfs1d(fileSaf=i,
                                                      memory='200g',
                                                      outFolder='results_5e6',
                                                      sites='pruning/snps_5e6.tsv',
                                                      cores=8 ))

##############################2D sfs   
mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
count_i=0
count_j=0
for i in filelists:
    count_j=0
    for j in filelists:
        if count_j>count_i:
            runName = i+'_'+j #os.path.split( os.path.basename(dataf) )[1]
            gwf.target_from_template( 'sfs2d_5e6_'+str(runName) , sfs2d(file1=i,
                                                                file2=j,
                                                                outFolder='results_5e6',
                                                                sites='pruning/snps_5e6.tsv',    
                                                                memory='200g',
                                                                cores=8 ))
        count_j += 1
    count_i += 1


######### Fst


mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
count_i=0
count_j=0
for i in filelists:
    count_j=0
    for j in filelists:
        if count_j>count_i:
            runName = i+'_'+j
            gwf.target_from_template( 'fst_prep_5e6_'+str(runName) , Fst_prepare(file1=i,
                                                                file2=j,
                                                                outFolder='results_5e6',
                                                                #sites='pruning/snps_5e6.tsv',             
                                                                walltime='8:00:00',
                                                                memory='100g',
                                                                cores=8))
            gwf.target_from_template( 'fst_calc_5e6_'+str(runName) , Fst_estimate(file1=i,
                                                                              file2=j,
                                                                              outFolder='results_5e6'))
        count_j += 1
    count_i += 1


######## Calculate thetas for each site. Step 2 at the link
######## http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests#Full_command_list_for_below_examples

mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]

for i in filelists:
    gwf.target_from_template( 'Theta_5e6_'+str(i) , Theta_estimate(file1=i,
                                                                outFolder='results_5e6',
                                                                walltime='2:00:00',
                                                                memory='20g',
                                                                cores=1))


mypath = './filelists'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]

count_i=0
count_j=0
for i in filelists:
    count_j=0
    for j in filelists:
        if count_j>count_i:
            runName = i+'_'+j
            gwf.target_from_template( 'fst_calc_masked_'+str(runName) , Fst_estimate_masked(file1=i,
                                                                file2=j,
                                                                siteFile='baits_350bp_3cols.bed',
                                                                outFolder='results_all'))
        count_j += 1
    count_i += 1



####################### 4 population test analysis
#######################      with D-statistic
#######################           ANGSD

###generate fasta files for the PIKE reference mapped to SALMON and for a modern individual of high quality.
###They are used for error rate estimation in all the genomes.

gwf.target_from_template( 'dofasta_pike' , doFasta(bamFile='../steps/admixture/PIKEmappedToLaks.sorted',
                                                    mode='EBD',
                                                    account='salmon',
                                                    sites='baits_350bp_3cols.bed',
                                                    memory='48g',
                                                    walltime='36:00:00',
                                                    cores=8 ))

gwf.target_from_template( 'dofasta_modern' , doFasta(bamFile='../steps/dedup/VA_16_31_HFH73BBXY.nodup.sorted',
                                                    mode='EBD',
                                                    out='./modern_reference_VA_16_31_HFH73BBXY',
                                                    account='salmon',
                                                    sites='baits_350bp_3cols.bed',
                                                    memory='48g',
                                                    walltime='36:00:00',
                                                    cores=8 ))  


for i in filelists:
    gwf.target_from_template( 'doAncError_modern_'+i , ancError(bamFiles=i,
                                                         anc_genome='../steps/admixture/PIKEmappedToLaks.sorted.fa',
                                                         ref_genome='./modern_reference_VA_16_31_HFH73BBXY.fa',
                                                         regions='pruning/snps_5e6_formatted.tsv',
                                                         out='out_ancError_modernRef'))


###Run the abbababa analysis excluding the genomes that were outliers in the PCA and calculation of genetic quantities.
    
mypath = './filelists_flt_pcaoutlier'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
    
count_i=0
count_j=0
count_k=0
for i in filelists:
    count_j=0
    for j in filelists:
        count_k=0
        for k in filelists:
            if (count_j>count_i)&(count_k>count_j)&(count_k>count_i):
                runName = i+'_'+j+'_'+k
                gwf.target_from_template( 'abbababa_pcaoutliers_'+runName , abbababa(file_1=i,
                                                                         file_2=j,
                                                                         file_3=k,
                                                                         out = 'uncorrected_masked_pcaoutliers/',
                                                                         anc_genome=ancFile,
                                                                         cores=8,
                                                                         regions='pruning/snps_5e6_formatted.tsv',
                                                                         memory='64g',
                                                                         walltime='72:00:00'))
                gwf.target_from_template( 'Dstat_pcaout_modern_'+runName , Dstat_calc(file_1=i,
                                                                         file_2=j,
                                                                         file_3=k,
                                                                         abbababa_folder = 'uncorrected_masked_pcaoutliers/',
                                                                         out = 'modernRefCorrected',
                                                                         err_folder = 'out_ancError_modernRef'))
            count_k += 1
        count_j += 1
    count_i += 1
