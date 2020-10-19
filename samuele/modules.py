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

