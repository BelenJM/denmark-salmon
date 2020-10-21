from gwf import Workflow
import os, sys
import numpy as np

######################SAF ANALYSIS
def dosaf(filelist, anc=None, regions=None, outFolder='results', account='salmon', memory='48g', walltime='48:00:00', cores=8):
    inputFile = filelist
    mypath = os.path.dirname(filelist)
    filelist = os.path.basename(filelist)
    inputs = [inputFile, anc]
    outputs = [outFolder+'/out_saf/'+filelist+'.saf.idx']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }
        
    if regions is None:
        spec = '''
        mkdir -p {outFolder}/out_saf
        source activate salmon
        echo "ciao 1"
        a=(`wc -l {inputFile}`)
        echo "ciao 2"
        minInd=$(($a*8/10))
        ./angsd/angsd -b {inputFile} -anc {anc} -ref {anc} -dosaf 1 -P {cores} -GL 2 -minMapQ 20 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -trim 3 -minInd $minInd -baq 1 -C 50 -out {outFolder}/out_saf/{filelist}
        '''.format(inputFile=inputFile, anc=anc, cores=cores, filelist=filelist, outFolder=outFolder)
    else:
        spec = '''
        mkdir -p {outFolder}/out_saf
        source activate salmon
        echo "ciao 1"
        a=(`wc -l {inputFile}`)
        echo "ciao 2"
        minInd=$(($a*8/10))
        ./angsd/angsd -b {inputFile} -anc {anc} -ref {anc} -rf {regions} -dosaf 1 -P {cores} -GL 2 -minMapQ 20 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -trim 3 -minInd $minInd -baq 1 -C 50 -out {outFolder}/out_saf/{filelist}
        '''.format(inputFile=inputFile, anc=anc, cores=cores, filelist=filelist, outFolder=outFolder, regions=regions)

    return inputs, outputs, options, spec


#Try without filterings ??? Just a copypaste to  be modified
def dosaf_unfilt(filelist, anc=None, regions=None, outFolder='results', account='salmon', memory='48g', walltime='48:00:00', cores=8):
    inputFile = filelist
    mypath = os.path.dirname(filelist)
    filelist = os.path.basename(filelist)
    inputs = [inputFile, anc]
    outputs = [outFolder+'/out_saf/'+filelist+'.saf.idx']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }
        
    if regions is None:
        spec = '''
        mkdir -p {outFolder}/out_saf
        source activate salmon
        echo "ciao 1"
        a=(`wc -l {inputFile}`)
        echo "ciao 2"
        minInd=$(($a*8/10))
        ./angsd/angsd -b {inputFile} -anc {anc} -ref {anc} -dosaf 1 -P {cores} -GL 2 -minMapQ 20 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -trim 3 -minInd $minInd -baq 1 -C 50 -out {outFolder}/out_saf/{filelist}
        '''.format(inputFile=inputFile, anc=anc, cores=cores, filelist=filelist, outFolder=outFolder)
    else:
        spec = '''
        mkdir -p {outFolder}/out_saf
        source activate salmon
        echo "ciao 1"
        a=(`wc -l {inputFile}`)
        echo "ciao 2"
        minInd=$(($a*8/10))
        ./angsd/angsd -b {inputFile} -anc {anc} -ref {anc} -rf {regions} -dosaf 1 -P {cores} -GL 2 -minMapQ 20 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -trim 3 -minInd $minInd -baq 1 -C 50 -out {outFolder}/out_saf/{filelist}
        '''.format(inputFile=inputFile, anc=anc, cores=cores, filelist=filelist, outFolder=outFolder, regions=regions)

    return inputs, outputs, options, spec




###### 1D sfs
######REMEMBER: do  `angsd sites index your.file` on the -sites file before running all this
######          so that the SNP list is indexed
######################################SFS1D analysis
###### Here you need -sites and not -rf as option (-rf is only for reading from bam files)
def sfs1d(fileSaf, sites=None, outFolder='results', account='salmon', memory='48g', walltime='12:00:00', cores=8):
    inputFile1 = outFolder+'/out_saf/'+fileSaf+'.saf.idx'
    inputs = [inputFile1]
    outputs = [outFolder+'/out_1dsfs/'+fileSaf+'.sfs']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }
        
    if sites is None:
        spec = '''
        mkdir -p {outFolder}/out_1dsfs
        source activate salmon
        ./angsd/misc/realSFS {inputFile1} -P {cores} > {outFolder}/out_1dsfs/{fileSaf}.sfs
        '''.format(inputFile1=inputFile1, fileSaf=fileSaf, cores=cores, outFolder=outFolder)
    else:
        spec = '''
        mkdir -p {outFolder}/out_1dsfs
        source activate salmon
        ./angsd/misc/realSFS {inputFile1} -P {cores} -sites {sites} > {outFolder}/out_1dsfs/{fileSaf}.sfs
        '''.format(inputFile1=inputFile1, fileSaf=fileSaf, cores=cores, sites=sites, outFolder=outFolder)
        
    return inputs, outputs, options, spec

##############################2D sfs    
def sfs2d(file1, file2, sites=None, outFolder='results', account='salmon', memory='48g', walltime='12:00:00', cores=8):
    inputFile1 = outFolder+'/out_saf/'+file1+'.saf.idx'
    inputFile2 = outFolder+'/out_saf/'+file2+'.saf.idx'
    inputs = [inputFile1, inputFile2]
    outputs = [outFolder+'/out_2dsfs/'+file1+'.'+file2+'.2dsfs.sfs']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }
        
    if sites is None:
        spec = '''
        mkdir -p {outFolder}/out_2dsfs
        source activate salmon
        ./angsd/misc/realSFS {inputFile1} {inputFile2} -P {cores} -fold 1 > {outFolder}/out_2dsfs/{file1}.{file2}.2dsfs.sfs
        '''.format(inputFile1=inputFile1, inputFile2=inputFile2, file1=file1, file2=file2, cores=cores, outFolder=outFolder)
    else:
        spec = '''
        mkdir -p {outFolder}/out_2dsfs
        source activate salmon
        ./angsd/misc/realSFS {inputFile1} {inputFile2} -P {cores} -sites {sites} -fold 1 > {outFolder}/out_2dsfs/{file1}.{file2}.2dsfs.sfs
        '''.format(inputFile1=inputFile1, inputFile2=inputFile2, file1=file1, file2=file2, cores=cores, sites=sites, outFolder=outFolder)
    return inputs, outputs, options, spec


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

######### Fst
def Fst_prepare(file1, file2, sites=None, outFolder='results', account='salmon', memory='30g', walltime='4:00:00', cores=4):
    inputFile1 = outFolder+'/out_saf/'+file1+'.saf.idx'
    inputFile2 = outFolder+'/out_saf/'+file2+'.saf.idx'
    inputs = [inputFile1, inputFile2, outFolder+'/out_2dsfs/'+file1+'.'+file2+'.2dsfs.sfs']
    outputs = [outFolder+'/out_fst/'+file1+'.'+file2+'.fst.idx']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }
        
    if sites is None:
        spec = '''
        mkdir -p {outFolder}/out_fst
        #source activate salmon
        ./angsd/misc/realSFS fst index {inputFile1} {inputFile2} -sfs {outFolder}/out_2dsfs/{file1}.{file2}.2dsfs.sfs -P {cores} -fold 1 -fstout {outFolder}/out_fst/{file1}.{file2}
        '''.format(inputFile1=inputFile1, inputFile2=inputFile2, file1=file1, file2=file2, cores=cores, outFolder=outFolder)
    else:
        spec = '''
        mkdir -p {outFolder}/out_fst
        #source activate salmon
        ./angsd/misc/realSFS fst index {inputFile1} {inputFile2} -sfs {outFolder}/out_2dsfs/{file1}.{file2}.2dsfs.sfs -sites {sites} -P {cores} -fold 1 -fstout {outFolder}/out_fst/{file1}.{file2}
        '''.format(inputFile1=inputFile1, inputFile2=inputFile2, sites=sites, outFolder=outFolder, file1=file1, file2=file2, cores=cores)
  
    return inputs, outputs, options, spec

def Fst_estimate(file1, file2, sites=None, outFolder='results', account='salmon', memory='20g', walltime='4:00:00', cores=1):
    inputs = [outFolder+'/out_fst/'+file1+'.'+file2+'.fst.idx']
    outputs = [outFolder+'/out_fst/'+file1+'.'+file2+'.fst.txt']
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }
        

    spec = '''
    source activate salmon
    ./angsd/misc/realSFS fst stats {outFolder}/out_fst/{file1}.{file2}.fst.idx > {outFolder}/out_fst/{file1}.{file2}.fst.txt
    '''.format(file1=file1, file2=file2, outFolder=outFolder)
    return inputs, outputs, options, spec

######## Calculate thetas for each site. Step 2 at the link
######## http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests#Full_command_list_for_below_examples
def Theta_estimate(file1, sites=None, outFolder='results', account='salmon', memory='20g', walltime='4:00:00', cores=1):
    safFile = outFolder+'/out_saf/'+file1+'.saf.idx'
    sfsFile = outFolder+'/out_1dsfs/'+file1+'.sfs'
    inputs = [safFile, sfsFile]
    outputs = [outFolder+'/out_theta/'+file1+'.thetas.gz',
               outFolder+'/out_theta/'+file1+'.thetas.idx',
               outFolder+'/out_theta/'+file1+'.logTheta.txt',
               outFolder+'/out_theta/'+file1+'.logTheta.Chr.txt'
               ]
    options = {
        'cores': cores,
        'memory': memory,
        'walltime': walltime,
        'account':account
    }
      
    spec = '''
    source activate salmon
    mkdir -p {outFolder}/out_theta
    ./angsd/misc/realSFS saf2theta {safFile} -sfs {sfsFile} -outname {outFolder}/out_theta/{file1}
    echo "print thetaStat"
    ./angsd/misc/thetaStat print {outFolder}/out_theta/{file1}.thetas.idx 2>/dev/null > {outFolder}/out_theta/{file1}.logTheta.txt
    echo "print thetaStat by chromosome (it works because chromosomes are short and cover an entire window)"
    ./angsd/misc/thetaStat do_stat {outFolder}/out_theta/{file1}.thetas.idx -outnames {outFolder}/out_theta/{file1}.logTheta.Chr
    mv {outFolder}/out_theta/{file1}.logTheta.Chr.pestPG {outFolder}/out_theta/{file1}.logTheta.Chr.txt
    '''.format(file1=file1, outFolder=outFolder, safFile=safFile, sfsFile=sfsFile)
    return inputs, outputs, options, spec


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



