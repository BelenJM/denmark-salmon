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
    gwf.target_from_template( 'dosaf_5e5_'+str(FILE) , dosaf(filelist=mypath+'/'+FILE,
                                                         anc=ancFile,
                                                         outFolder='results_5e5',
                                                         regions='pruning/snps_5e5_formatted.tsv',
                                                         cores=18,
                                                         memory='64g',
                                                         walltime='48:00:00'))
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

mypath = './filelists'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]

for i in filelists:
    gwf.target_from_template( 'sfs1d_'+str(i) , sfs1d(fileSaf=i,
                                                      memory='200g',
                                                      outFolder='results_all',
                                                      cores=8 ))


mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]

for i in filelists:
    gwf.target_from_template( 'sfs1d_5e5_'+str(i) , sfs1d(fileSaf=i,
                                                      memory='200g',
                                                      outFolder='results_5e5',
                                                      sites='pruning/snps_5e5.tsv',
                                                      cores=8 ))
    gwf.target_from_template( 'sfs1d_5e6_'+str(i) , sfs1d(fileSaf=i,
                                                      memory='200g',
                                                      outFolder='results_5e6',
                                                      sites='pruning/snps_5e6.tsv',
                                                      cores=8 ))

##############################2D sfs   

mypath = './filelists'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
count_i=0
count_j=0
for i in filelists:
    count_j=0
    for j in filelists:
        if count_j>count_i:
            runName = i+'_'+j #os.path.split( os.path.basename(dataf) )[1]
            gwf.target_from_template( 'sfs2d_'+str(runName) , sfs2d(file1=i,
                                                                file2=j,
                                                                outFolder='results_all',
                                                                memory='200g',
                                                                cores=8 ))
        count_j += 1
    count_i += 1


mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
count_i=0
count_j=0
for i in filelists:
    count_j=0
    for j in filelists:
        if count_j>count_i:
            runName = i+'_'+j #os.path.split( os.path.basename(dataf) )[1]
            gwf.target_from_template( 'sfs2d_5e5_'+str(runName) , sfs2d(file1=i,
                                                                file2=j,
                                                                outFolder='results_5e5',
                                                                sites='pruning/snps_5e5.tsv',    
                                                                memory='200g',
                                                                cores=8 ))
            gwf.target_from_template( 'sfs2d_5e6_'+str(runName) , sfs2d(file1=i,
                                                                file2=j,
                                                                outFolder='results_5e6',
                                                                sites='pruning/snps_5e6.tsv',    
                                                                memory='200g',
                                                                cores=8 ))
        count_j += 1
    count_i += 1


######### Fst

mypath = './filelists'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
count_i=0
count_j=0
for i in filelists:
    count_j=0
    for j in filelists:
        if count_j>count_i:
            runName = i+'_'+j
            gwf.target_from_template( 'fst_prep_'+str(runName) , Fst_prepare(file1=i,
                                                                file2=j,
                                                                outFolder='results_all',
                                                                #sites='pruning/snps_5e6.tsv',             
                                                                walltime='8:00:00',
                                                                memory='100g',
                                                                cores=8))
            gwf.target_from_template( 'fst_calc_'+str(runName) , Fst_estimate(file1=i,
                                                                              file2=j,
                                                                              outFolder='results_all'))
        count_j += 1
    count_i += 1

mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
count_i=0
count_j=0
for i in filelists:
    count_j=0
    for j in filelists:
        if count_j>count_i:
            runName = i+'_'+j
            gwf.target_from_template( 'fst_prep_5e5_'+str(runName) , Fst_prepare(file1=i,
                                                                file2=j,
                                                                outFolder='results_5e5',
                                                                #sites='pruning/snps_5e5.tsv',             
                                                                walltime='8:00:00',
                                                                memory='100g',
                                                                cores=8))
            gwf.target_from_template( 'fst_calc_5e5_'+str(runName) , Fst_estimate(file1=i,
                                                                              file2=j,
                                                                              outFolder='results_5e5'))
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
count_i=0
count_j=0
for i in filelists:
    gwf.target_from_template( 'Theta_5e5_'+str(i) , Theta_estimate(file1=i,
                                                                outFolder='results_5e5',
                                                                walltime='2:00:00',
                                                                memory='20g',
                                                                cores=1))
    gwf.target_from_template( 'Theta_5e6_'+str(i) , Theta_estimate(file1=i,
                                                                outFolder='results_5e6',
                                                                walltime='2:00:00',
                                                                memory='20g',
                                                                cores=1))








##################OLD estimates with bait coordinates
##################QUITE useless but keep them here now
def Fst_estimate_masked(file1, file2, outFolder='results', siteFile='sites.txt', account='salmon', memory='20g', walltime='4:00:00', cores=1):
    inputs = [outFolder+'/out_fst/'+file1+'.'+file2+'.fst.idx']
    outputs = [outFolder+'/out_fst/masked/'+file1+'.'+file2+'.fst.txt']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }
        
  
    spec = '''
    source activate salmon
    mkdir -p {outFolder}/out_fst/masked
    #./angsd/angsd sites index {siteFile}
    ./angsd/misc/realSFS fst stats {outFolder}/out_fst/{file1}.{file2}.fst.idx -sites {siteFile} > {outFolder}/out_fst/masked/{file1}.{file2}.fst.txt
    '''.format(file1=file1, file2=file2, outFolder=outFolder, siteFile=siteFile)

    return inputs, outputs, options, spec

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

######################Probably NEEDS RESTRUCTURING OF OUTPUT FOLDERS
def ancError(bamFiles, anc_genome, ref_genome, out='out_ancError', regions=None, account='salmon', memory='80g', walltime='8:00:00', cores=4):
    inputs = ['./filelists_flt_pcaoutlier/'+bamFiles, ref_genome, anc_genome, f'{ref_genome}.fai', f'{anc_genome}.fai']
    outputs = [out+'/'+bamFiles+'.ancError']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }
        
    if regions is None:
        spec = '''
        source activate salmon
        mkdir -p ./{out}
    
        ./angsd/angsd -doAncError 1 -anc {anc_genome} -ref {ref_genome} -out ./{out}/{bamFiles} -bam ./filelists_flt_pcaoutlier/{bamFiles} -minMapQ 20 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -trim 3 -p {cores}
        '''.format(bamFiles=bamFiles, anc_genome=anc_genome,  ref_genome=ref_genome, out=out, cores=cores)
    else:
        spec = '''
        source activate salmon
        mkdir -p ./{out}
    
        ./angsd/angsd -doAncError 1 -anc {anc_genome} -ref {ref_genome} -out ./{out}/{bamFiles} -bam ./filelists_flt_pcaoutlier/{bamFiles} -minMapQ 20 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -trim 3 -p {cores} -rf {regions}
        '''.format(bamFiles=bamFiles, anc_genome=anc_genome,  ref_genome=ref_genome, out=out, cores=cores, regions=regions)

    return inputs, outputs, options, spec

   

for i in filelists:
    gwf.target_from_template( 'doAncError_'+i , ancError(bamFiles=i,
                                                         anc_genome='../steps/admixture/PIKEmappedToLaks.sorted.fa',
                                                         ref_genome='./GCF_000233375.1_ICSASG_v2_genomic.fa',
                                                         regions='pruning/snps_5e6_formatted.tsv',
                                                         out='out_ancError'))
for i in filelists:
    gwf.target_from_template( 'doAncError_modern_'+i , ancError(bamFiles=i,
                                                         anc_genome='../steps/admixture/PIKEmappedToLaks.sorted.fa',
                                                         ref_genome='./modern_reference_VA_16_31_HFH73BBXY.fa',
                                                         regions='pruning/snps_5e6_formatted.tsv',
                                                         out='out_ancError_modernRef'))



def doFasta(bamFile, mode='EBD', account='salmon',sites=None, out=None, memory='48g', walltime='36:00:00', cores=8):
    inputs = [f'{bamFile}.bam']
    if out==None:
        out=bamFile    
    outputs = [f'{out}.fa',f'{out}.fa.fai']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }



    
    #Three possible modes for running fasta file creation
    if mode=='random':
        mode='1'
    if mode=='most common':
        mode='2 -doCounts 1'
    if mode=='EBD':
        mode='3'
        
    if sites is None:    
        spec = '''
        source activate salmon

        ./angsd/angsd -i {bamFile}.bam -doFasta {mode} -p {cores} -out {out}
        gunzip {out}.fa.gz
        samtools faidx {out}.fa
        '''.format(mode=mode, bamFile=bamFile, cores=cores, out=out)
    else:
        spec = '''
        source activate salmon

        ./angsd/angsd -i {bamFile}.bam -doFasta {mode} -p {cores} -sites {sites} -out {out}
        gunzip {out}.fa.gz
        samtools faidx {out}.fa
        '''.format(mode=mode, bamFile=bamFile, cores=cores, sites=sites, out=out)

    return inputs, outputs, options, spec



def abbababa(file_1, file_2, file_3, anc_genome, out, account='salmon', regions=None, memory='80g', walltime='6:00:00', cores=8):
    inputs = ['./filelists/'+file_1, './filelists/'+file_2, './filelists/'+file_3]
    outputs = ['./out_4pop/'+out+'/'+file_1+'.'+file_2+'.'+file_3+'.abbababa2']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }

    if regions is None:
        spec = '''
        source activate salmon
        mkdir -p ./out_4pop/{out}

        cat ./filelists/{file_1} ./filelists/{file_2} ./filelists/{file_3} > ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.filelist
        echo {anc_genome} >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.filelist
        cat ./filelists/{file_1} | wc -l > ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.size
        cat ./filelists/{file_2} | wc -l >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.size
        cat ./filelists/{file_3} | wc -l >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.size
        echo 1 >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.size
        echo {file_1} > ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.pops
        echo {file_2} >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.pops
        echo {file_3} >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.pops
        echo "SalmoTrutta" >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.pops

        ./angsd/angsd -doAbbababa2 1 -bam ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.filelist -sizeFile ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.size -doCounts 1 -out ./out_4pop/{out}/{file_1}.{file_2}.{file_3} -useLast 1 -minMapQ 20 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -trim 3 -p {cores}
        '''.format(file_1=file_1, file_2=file_2, file_3=file_3, anc_genome=anc_genome, out=out, cores=cores)
    else:
        spec = '''
        source activate salmon
        mkdir -p ./out_4pop/{out}

        cat ./filelists/{file_1} ./filelists/{file_2} ./filelists/{file_3} > ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.filelist
        echo {anc_genome} >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.filelist
        cat ./filelists/{file_1} | wc -l > ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.size
        cat ./filelists/{file_2} | wc -l >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.size
        cat ./filelists/{file_3} | wc -l >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.size
        echo 1 >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.size
        echo {file_1} > ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.pops
        echo {file_2} >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.pops
        echo {file_3} >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.pops
        echo "SalmoTrutta" >> ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.pops

        ./angsd/angsd -doAbbababa2 1 -bam ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.filelist -sizeFile ./out_4pop/{out}/{file_1}.{file_2}.{file_3}.size -doCounts 1 -out ./out_4pop/{out}/{file_1}.{file_2}.{file_3} -useLast 1 -minMapQ 20 -minQ 20 -remove_bads 1 -rf {regions} -only_proper_pairs 1 -trim 3 -p {cores}

        #mkdir -p ./out_4pop/{out}/result/
        #Rscript ./angsd/R/estAvgError.R angsdFile="./out_4pop/{out}/{file_1}.{file_2}.{file_3}" out="./out_4pop/{out}/result/{file_1}.{file_2}.{file_3}" sizeFile="./out_4pop/{out}/{file_1}.{file_2}.{file_3}.size" errFile="./out_4pop/{out}/{file_1}.{file_2}.{file_3}.error" nameFile="./out_4pop/{out}/{file_1}.{file_2}.{file_3}.pops"
        '''.format(file_1=file_1, file_2=file_2, file_3=file_3, anc_genome=anc_genome, out=out, cores=cores, regions=regions)

    return inputs, outputs, options, spec



def Dstat_calc(file_1, file_2, file_3, abbababa_folder, out, err_folder=None, account='salmon', memory='48g', walltime='6:00:00', cores=1):
    inputs = ['./out_4pop/'+abbababa_folder+'/'+file_1+'.'+file_2+'.'+file_3+'.abbababa2']
    if err_folder is not None:
        inputs = ['./out_4pop/'+abbababa_folder+'/'+file_1+'.'+file_2+'.'+file_3+'.abbababa2',
                  err_folder + '/' + file_1 + '.ancError',
                  err_folder + '/' + file_2 + '.ancError',
                  err_folder + '/' + file_3 + '.ancError']
        
    outputs = [f'./results_5e6/Dstat_{out}/'+file_1+'.'+file_2+'.'+file_3+'.Observed.txt',
               f'./results_5e6/Dstat_{out}/'+file_1+'.'+file_2+'.'+file_3+'.TransRem.txt']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }

        
    if err_folder is None:
        spec = '''
        source activate salmon
        mkdir -p ./results_5e6/Dstat_{out}/
        Rscript ./angsd/R/estAvgError.R angsdFile="./out_4pop/{abbababa_folder}/{file_1}.{file_2}.{file_3}" out="./results_5e6/Dstat_{out}/{file_1}.{file_2}.{file_3}" sizeFile="./out_4pop/{abbababa_folder}/{file_1}.{file_2}.{file_3}.size" nameFile="./out_4pop/{abbababa_folder}/{file_1}.{file_2}.{file_3}.pops"
        '''.format(file_1=file_1, file_2=file_2, file_3=file_3, out=out, abbababa_folder=abbababa_folder)
    else:
        spec = '''
        source activate salmon
        mkdir -p ./results_5e6/Dstat_{out}/

        echo {err_folder}/{file_1}.ancError > ./out_4pop/{abbababa_folder}/{file_1}.{file_2}.{file_3}.error
        echo {err_folder}/{file_2}.ancError >> ./out_4pop/{abbababa_folder}/{file_1}.{file_2}.{file_3}.error
        echo {err_folder}/{file_3}.ancError >> ./out_4pop/{abbababa_folder}/{file_1}.{file_2}.{file_3}.error
        echo NA >> ./out_4pop/{abbababa_folder}/{file_1}.{file_2}.{file_3}.error

        Rscript ./angsd/R/estAvgError.R angsdFile="./out_4pop/{abbababa_folder}/{file_1}.{file_2}.{file_3}" out="./results_5e6/Dstat_{out}/{file_1}.{file_2}.{file_3}" sizeFile="./out_4pop/{abbababa_folder}/{file_1}.{file_2}.{file_3}.size" nameFile="./out_4pop/{abbababa_folder}/{file_1}.{file_2}.{file_3}.pops" errFile="./out_4pop/{abbababa_folder}/{file_1}.{file_2}.{file_3}.error"
        '''.format(file_1=file_1, file_2=file_2, file_3=file_3, out=out, err_folder=err_folder, abbababa_folder=abbababa_folder)
        
    return inputs, outputs, options, spec






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
                gwf.target_from_template( 'Dstat_pcaoutliers_'+runName , Dstat_calc(file_1=i,
                                                                         file_2=j,
                                                                         file_3=k,
                                                                         abbababa_folder = 'uncorrected_masked_pcaoutliers/',
                                                                         out = 'salmonRefCorrected',
                                                                         err_folder = 'out_ancError'))
                gwf.target_from_template( 'Dstat_pcaout_modern_'+runName , Dstat_calc(file_1=i,
                                                                         file_2=j,
                                                                         file_3=k,
                                                                         abbababa_folder = 'uncorrected_masked_pcaoutliers/',
                                                                         out = 'modernRefCorrected',
                                                                         err_folder = 'out_ancError_modernRef'))
            count_k += 1
        count_j += 1
    count_i += 1

    
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


#######
#Remove trouts from filelists !!!!!
#Do abbababa when getting an aligned trout
#try abbababa with removed trouts
#SFS 1d with bait mask
#SFS 1d with new selected SNPs when they arrive
#Calculate Fst with new selected SNPs when they arrive
###Also keep the pvalues and puth them in the upper part of the matrix
