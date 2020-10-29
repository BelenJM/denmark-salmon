##################################################
##### templates for the first part of the
##### analysis (genomic pipeline)
##### author: Belen JM based on genomic pipeline
#####         developed by Alice Manuzzi
##### date of final version: september 2019
#################################################


# unzipping file
def unzip(inputfile, outputfile):
    """A template for unzipping files."""
    inputs = [inputfile]
    outputs = [outputfile]
    options = {
        'cores': 1,
        'memory': '100mb',
    }

    spec = '''
    zcat {} > {}
    '''.format(inputfile, outputfile)

    return inputs, outputs, options, spec

# creating contamination baseline
def create_kraken_lib(folder):
    inputs = []
    outputs = ['data/Kraken_library/bowtie/human.1.bt2',
		'data/Kraken_library/bowtie/human.2.bt2',
		'data/Kraken_library/bowtie/human.3.bt2',
		'data/Kraken_library/bowtie/human.4.bt2',
		'data/Kraken_library/bowtie/human.rev.1.bt2',
		'data/Kraken_library/bowtie/human.rev.2.bt2']
    options = {
        'cores': 27,
        'memory': '200g'
    }

    spec = '''
    # activate the environment with kraken2 installed
    # conda activate salmon-contamination
    # kraken2-build below does not work so I copied files from the Computerome cluster to proceed
    # command followed: 	
    # kraken2-build --standard --db data/Kraken_library/library_contami --use-ftp
    # scp -r bmen@computerome.cbs.dtu.dk:/home/projects/dp_00007/data/Kraken_library/library_contami data/Kraken_library/.
    # kraken2-build --build --threads 27 --db data/Kraken_library/library_contami
    eval "$(conda shell.bash hook)"
    conda activate salmon
    #create index from human genome
    bowtie2-build -f data/Kraken_library/library_contami/library/human/library.fna data/Kraken_library/bowtie/human
'''.format()
    return inputs, outputs, options, spec

def kraken_report(r1, r2, ind2):
    """Template to check for contamination"""
    inputs = ['steps/filtering/{}'.format(r1),
                'steps/filtering/{}'.format(r2),
		'data/Kraken_library/library_contami/hash.k2d',
                'data/Kraken_library/library_contami/library',
                'data/Kraken_library/library_contami/opts.k2d',
                'data/Kraken_library/library_contami/seqid2taxid.map',
                'data/Kraken_library/library_contami/taxo.k2d',
                'data/Kraken_library/library_contami/taxonomy',
                'data/Kraken_library/library_contami/unmapped.txt']
    outputs = ['steps/check_contami/{}_kraken'.format(ind2),
		'steps/check_contami/{}_report'.format(ind2)]
    options = {'cores': 20,
			'memory': '150g',
			'walltime':'00:30:00'
}
    spec = '''
    eval "$(conda shell.bash hook)"
    conda activate salmon-contamination
    kraken2 --threads 20 --db data/Kraken_library/library_contami/ --gzip-compressed --paired steps/filtering/{} steps/filtering/{} --output steps/check_contami/{}_kraken --use-names --report steps/check_contami/{}_report
    '''.format(r1, r2, ind2, ind2)
    return inputs, outputs, options, spec


def visualize_kraken(ind2):
    """Template to visualize kraken report"""
    inputs = ['steps/check_contami/{}_kraken'.format(ind2)]
    outputs = ['steps/check_contami/{}_output_id'.format(ind2),
		'steps/check_contami/{}_output_name'.format(ind2),
		'steps/check_contami/{}_output_merged'.format(ind2),
		'steps/check_contami/{}_output_out.html'.format(ind2)]
    options = {'cores': 1,
			'memory': '2g',
                        'walltime':'00:30:00'}

    spec = """
        eval "$(conda shell.bash hook)"
    	conda activate salmon-contamination

	python scripts/input_for_krona.py steps/check_contami/{}_kraken steps/check_contami/{}_output_id
        cat steps/check_contami/{}_kraken | cut -f2 > steps/check_contami/{}_output_name
        paste steps/check_contami/{}_output_name steps/check_contami/{}_output_id > steps/check_contami/{}_output_merged
        # krona needs first to upload taxonomy:
	ktUpdateTaxonomy.sh
	# and now krona creates output
	ktImportTaxonomy steps/check_contami/{}_output_merged -o steps/check_contami/{}_output_out.html

        #download the report files to desktop to visualize the comparison in R and the html file to see pie charts of krona
        #or to create a plot showing the different percentages
        #find . -name "*_report" | xargs grep -n -H "Homo sapiens" >> homo-sapiens.txt #for every mapped species

    """.format(ind2, ind2, ind2, ind2,ind2,ind2,ind2,ind2,ind2)
    return inputs, outputs, options, spec

def remove_contamination(r1, r2, ind2, organism, path_reads):
    """Template to identify and remove human/bacteria contamination identified in kraken2"""
    inputs = ['{}/{}'.format(path_reads, r1),
					'{}/{}'.format(path_reads, r2),
					'data/Kraken_library/bowtie/']
    outputs = ['steps/check_contami/mapping-{}/{}.mapping_{}.sam'.format(organism, ind2, organism),
					'steps/check_contami/mapping-{}/{}.count.{}'.format(organism, ind2, organism),
					'steps/check_contami/mapping-{}/{}.{}.RemoveList.txt'.format(organism, ind2, organism),
					'steps/filter_humanBact/{}.no{}_1.fastq.gz'.format(ind2, organism),
					'steps/filter_humanBact/{}.no{}_2.fastq.gz'.format(ind2, organism)
					]
    options = {'cores': 10,
			'memory': '30g',
                        'walltime':'00:30:00'}

    spec = """
        eval "$(conda shell.bash hook)"
        conda activate salmon-contamination
	
	# map to the genome
        bowtie2 -q --phred33 --sensitive -I 0 -X 200 --fr -x data/Kraken_library/bowtie/{} -1 {}/{} -2 {}/{} -S steps/check_contami/mapping-{}/{}.mapping_{}.sam

        # Count reads that mapped to the contaminant database and should be removed - bacteria_map_SINGLE/PAIRED.sh
        samtools view -S -c -F 4 steps/check_contami/mapping-{}/{}.mapping_{}.sam > steps/check_contami/mapping-{}/{}.count.{}

        # Make a list of reads that mapped to the contaminant database and should be removed (-F 4 mapped, -f 4 unmapped)
        samtools view -S -F 4 steps/check_contami/mapping-{}/{}.mapping_{}.sam | cut -f 1 | uniq > steps/check_contami/mapping-{}/{}.{}.RemoveList.txt

        ###removal of sequences using bbmap
        filterbyname.sh in={}/{} in2={}/{} out=steps/filter_humanBact/{}.no{}_1.fastq.gz out2=steps/filter_humanBact/{}.no{}_2.fastq.gz names=steps/check_contami/mapping-{}/{}.{}.RemoveList.txt include=f
""".format(organism, path_reads, r1, path_reads, r2, organism, ind2, organism,
			organism, ind2, organism, organism, ind2,organism,
			organism, ind2, organism, organism, ind2, organism,
			path_reads, r1, path_reads, r2, ind2, organism, ind2, organism, organism, ind2, organism)
    return inputs, outputs, options, spec

# add control contamination
# for i in $(cat ../list-names-final.txt); do zcat $i".noHuman_C.fastq.gz" | echo $i $((`wc -l`/4)) >> number_of_reads_cleanedH_C.txt; done

def get_bacteria_genome():
    """Template to collect all bacteria names that have been found in all samples and build a library"""

    inputs = ['steps/check_contami/{}_report', # to do, from all individuals
		'steps/filter_humanBact/{}.noHuman_1.fastq.gz'.format(ind2),
		'steps/filter_humanBact/{}.noHuman_2.fastq.gz'.format(ind2),
		'data/Kraken_library/bacteria/library.fna']
    outputs = ['steps/check_contami/bacteria',
		'steps/check_contami/final-bacteria-list',
		'steps/filter_humanBact/{}.noHuman_1.fastq.gz'.format(ind2),
		'steps/filter_humanBact/{}.noHuman_2.fastq.gz'.format(ind2)]
    options = {'cores': 10,
		'memory': '30g'}

    spec = """
	#identify the list of bacteria found in the samples and build a subset of the library (temporary file)
        cat steps/check_contami/{}_report" | grep "bacteria" | awk '$1 != 0.00' | cut -f6 >> steps/check_contami/bacteria

        cat steps/check_contami/bacteria | sort -u > steps/check_contami/final-bacteria-list

        cat data/Kraken_library/bacteria/library.fna | grep "$i" >> steps/check_contami/bacteria/subset-library-bacteria.txt
        sed 's/>//' steps/check_contami/bacteria/subset-library-bacteria.txt > steps/check_contami/bacteria/subset-library-bacteria-2.txt

        #Extract sequences with names in file name.lst, one sequence name per line:
        seqtk subseq data/Kraken_library/bacteria/library.fna steps/check_contami/bacteria/subset-library-bacteria-2.txt > steps/check_contami/bacteria/subset-bacteria-temporary.fa
        bbmap/36.49/filterbyname.sh in=data/Kraken_library/bacteria/library.fna out=steps/check_contami/bacteria/subset-bacteria-temporary.fa names=steps/check_contami/bacteria/subset-library-bacteria.txt include=t
        #bowtie2-build --wrapper basic-0 -f --large-index -f library.fna bacteria3
        bowtie2-build -f steps/check_contami/bacteria/subset-bacteria-temporary.fa data/Kraken_library/bacteria/bacteria1st

""".format(r1, r2, ind2,
			ind2, ind2,
			ind2, ind2,
			r1, r2, ind2, ind2, ind2,
			ind2)
    return inputs, outputs, options, spec


def check_bacteria():
    """Template to identify and remove bacterial genome identified in kraken2"""
    inputs = ['steps/filter_humanBact/{}.noHuman_1.fastq.gz'.format(ind2),
					'steps/filter_humanBact/{}.noHuman_2.fastq.gz'.format(ind2)]
    outputs = ['steps/filter_humanBact/{}.noHumanBact_1.fastq.gz'.format(ind2),
					'steps/filter_humanBact/{}.noHumanBact_2.fastq.gz'.format(ind2)]
    options = {'cores': 10,
                        'memory': '30g',
                        'walltime':'00:30:00'}

    spec = """

        bowtie2 -q --phred33 --sensitive -I 0 -X 200 --fr -x data/Kraken_library/bacteria/bacteria1st -1 steps/filter_humanBact/{}.noHuman_1.fastq.gz -2 steps/filter_humanBact/{}.noHuman_2.fastq.gz -S steps/check_contami/mapping-bacteria/{}.mapping_bacteria.sam
        samtools view -S -c -F 4 steps/check_contami/mapping-bacteria/{}.mapping_bacteria.sam > steps/check_contami/mapping-bacteria/{}.count.bacteria
        samtools view -S -F 4 steps/check_contami/mapping-bacteria/{}.mapping_bacteria.sam | cut -f 1 | uniq > steps/check_contami/mapping-bacteria/{}.bacteria.RemoveList.txt

        bbmap/36.49/filterbyname.sh in=steps/filter_humanBact/{}.noHuman_1.fastq.gz in2=steps/filter_humanBact/{}.noHuman_2.fastq.gz out=steps/filter_humanBact/{}.noHumanBact_1.fastq.gz out2=steps/filter_humanBact/{}.noHumanBact_2.fastq.gz names=steps/check_contami/mapping-bacteria/{}.bacteria.RemoveList.txt include=f

""".format()
    return inputs, outputs, options, spec

# indexing the genome
def bwa_index(ref_genome):
    """Template for indexing a genome with `bwa index`."""
    inputs = ['data/genome/{}.fa.gz'.format(ref_genome)]
    outputs = ['data/genome/{}.fa.gz.amb'.format(ref_genome),
               'data/genome/{}.fa.gz.ann'.format(ref_genome),
               'data/genome/{}.fa.gz.pac'.format(ref_genome),
               'data/genome/{}.fa.gz.bwt'.format(ref_genome),
               'data/genome/{}.fa.gz.sa'.format(ref_genome),
                'data/genome/{}.dict'.format(ref_genome),
                'data/genome/{}.fa.gz.fai'.format(ref_genome),
                'data/genome/output.fa.gz'.format(),
                'data/genome/output.fa.gz.gzi'.format()]

    options = {
        'cores': 26,
        'memory': '10g',
	'walltime':'01:00:00'
    }

    spec = """
    
    eval "$(conda shell.bash hook)"
    conda activate salmon
    
    cd data/genome/

    bwa index -a bwtsw {ref_genome}.fa.gz

    picard CreateSequenceDictionary R={ref_genome}.fa.gz \
           O={ref_genome}.dict
    zcat {ref_genome}.fa.gz | bgzip -c > output.fa.gz

    samtools faidx output.fa.gz

    mv output.fa.gz.fai {ref_genome}.fa.gz.fai
""".format(ref_genome=ref_genome)

    return inputs, outputs, options, spec

# mapping reads to genome
def bwa_map(ref_genome, r1, r2, bamfile):
    """Template for mapping reads to a reference genome with `bwa` and `samtools`."""
    inputs = ['steps/filter_humanBact/{}'.format(r1),
                'steps/filter_humanBact/{}'.format(r2),
              'data/genome/{}.fa.gz.amb'.format(ref_genome),
              'data/genome/{}.fa.gz.ann'.format(ref_genome),
              'data/genome/{}.fa.gz.pac'.format(ref_genome),
             ]
    outputs = ['steps/mapping/{}.mapped.sam'.format(bamfile),
                'steps/mapping/{}.unsorted.bam'.format(bamfile),
                'steps/mapping/{}.sorted.bam'.format(bamfile),
                'steps/mapping/{}.sorted.bam.bai'.format(bamfile)]
    options = {
        'cores': 24,
        'memory': '50g',
	'walltime': '03:00:00'
    }

    spec = """
    eval "$(conda shell.bash hook)"
    conda activate salmon

    bwa mem -t 24 data/genome/{ref_genome}.fa.gz steps/filter_humanBact/{r1} steps/filter_humanBact/{r2} > steps/mapping/{bamfile}.mapped.sam
    samtools view -bS -q 10 steps/mapping/{bamfile}.mapped.sam > steps/mapping/{bamfile}.unsorted.bam
    samtools sort -o steps/mapping/{bamfile}.sorted.bam steps/mapping/{bamfile}.unsorted.bam
    samtools index steps/mapping/{bamfile}.sorted.bam
    """.format(ref_genome=ref_genome, r1=r1, r2=r2, bamfile=bamfile)

    return inputs, outputs, options, spec


# changing names
def create_RG(bamfile, ind_number, lib_number, ind_name):
    """Template for changing read group to reads"""
    inputs=["steps/mapping/{}.sorted.bam".format(bamfile)]
    outputs=["steps/RG/{}.RG.bam".format(bamfile),
                "steps/RG/{}.RG.sorted.bam".format(bamfile),
                "steps/RG/{}.RG.sorted.bam.bai".format(bamfile)]
    options = {
        'cores': 4,
        'memory': '5g',
        'walltime': '00:15:00'
}

    spec = """
    eval "$(conda shell.bash hook)"
    conda activate salmon

    picard AddOrReplaceReadGroups I=steps/mapping/{bamfile}.sorted.bam O=steps/RG/{bamfile}.RG.bam RGID={ind_number} RGLB={lib_number} RGPL=Illumina RGPU=HiSeq4000 RGSM={ind_name} VALIDATION_STRINGENCY=LENIENT
    samtools view -b steps/RG/{bamfile}.RG.bam > steps/RG/{bamfile}.RG.unsorted.bam
    samtools sort -o steps/RG/{bamfile}.RG.sorted.bam steps/RG/{bamfile}.RG.unsorted.bam
    samtools index steps/RG/{bamfile}.RG.sorted.bam

    """.format(bamfile=bamfile, ind_number= ind_number,lib_number=lib_number, ind_name=ind_name)

    return inputs, outputs,options,spec

# mark duplicates
def mark_dupli(bamfile):
    """Template for marking duplicates"""
    inputs=["steps/RG/{}.RG.sorted.bam".format(bamfile)]
    outputs=["steps/dedup/{}.dedup.unsorted.bam".format(bamfile),
                "steps/dedup/{}.dedup.metrix.txt".format(bamfile),
                "steps/dedup/{}.dedup.sorted.bam".format(bamfile),
                "steps/dedup/{}.dedup.sorted.bam.bai".format(bamfile),
                "steps/dedup/{}.nodup.metrix.txt".format(bamfile),
                "steps/dedup/{}.nodup.unsorted.bam".format(bamfile),
                "steps/dedup/{}.nodup.sorted.bam".format(bamfile),
                "steps/dedup/{}.nodup.sorted.bam.bai".format(bamfile)]
    options = {
        'cores':6,
        'memory': '100g',
        'walltime' : '01:00:00'
}

    spec = """
        eval "$(conda shell.bash hook)"
        conda activate salmon

	picard -Xmx50g MarkDuplicates I=steps/RG/{bamfile}.RG.sorted.bam O=steps/dedup/{bamfile}.dedup.unsorted.bam M=steps/dedup/{bamfile}.dedup.metrix.txt VALIDATION_STRINGENCY=SILENT TAGGING_POLICY=All REMOVE_DUPLICATES=FALSE
	#sort and index the file
        samtools sort -o steps/dedup/{bamfile}.dedup.sorted.bam steps/dedup/{bamfile}.dedup.unsorted.bam
        samtools index steps/dedup/{bamfile}.dedup.sorted.bam
        
	picard -Xmx50g MarkDuplicates I=steps/dedup/{bamfile}.dedup.sorted.bam O=steps/dedup/{bamfile}.nodup.unsorted.bam M=steps/dedup/{bamfile}.nodup.metrix.txt VALIDATION_STRINGENCY=SILENT TAGGING_POLICY=All REMOVE_DUPLICATES=TRUE
        #sort and index the file
        samtools sort -o steps/dedup/{bamfile}.nodup.sorted.bam steps/dedup/{bamfile}.nodup.unsorted.bam
        samtools index steps/dedup/{bamfile}.nodup.sorted.bam
        """.format(bamfile=bamfile)
    return inputs, outputs, options, spec

def indel_realign(bamfile, ref_genome):
    """Target to first create realignment"""
    inputs=["steps/dedup/{}.nodup.sorted.bam".format(bamfile),
                "data/genome/{}.fa".format(ref_genome)]
    outputs=["steps/realign/{}.abra.bam".format(bamfile),
                "steps/realign/{}.abra.log".format(bamfile),
                "steps/realign/{}.abra.filtered.bam".format(bamfile),
                "steps/realign/{}.abra.filtered.bam.bai".format(bamfile)]
    options={'cores': 27,
                'memory':'300g',
                'walltime': '100:00:00'
}
    spec = """

        # in a previous version GATK RealignerTargetCreator+IndelRealigner were used, now ABRA2 is being used
        # java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R data/genome/{ref_genome}.fa -I steps/dedup/{bamfile}.dedup.sorted.bam -o steps/realign/{bamfile}.target
        eval "$(conda shell.bash hook)"
        conda activate salmon
	
        JAVA_TOOL_OPTIONS="-Xmx300G" abra2 --in steps/dedup/{bamfile}.dedup.sorted.bam --out steps/realign/{bamfile}.abra.bam --ref data/genome/{ref_genome}.fa --threads 27 --tmpdir steps/realign/temp > steps/realign/{bamfile}.abra.log
        samtools index steps/realign/{bamfile}.abra.bam

        # create a second file with only properly paired, without duplicates and q>10, to compare
	samtools view -f 2 -F 1024 -q 10 > steps/realign/{bamfile}.abra.filtered.bam
        samtools index steps/realign/{bamfile}.abra.filtered.bam
        """.format(bamfile=bamfile, ref_genome=ref_genome)
    return inputs, outputs, options, spec

def indel_realignGATK(bamfile, ref_genome):
    """Target to first create realignment"""
    inputs=["steps/dedup/{}.nodup.sorted.bam".format(bamfile),
                "data/genome/{}.fa".format(ref_genome)]
    outputs=["steps/realign2/{}.intervals".format(bamfile)]
    options={'cores': 3,
                'memory':'4g',
                'walltime': '15:00:00'
}
    spec = """
        eval "$(conda shell.bash hook)"
        conda activate salmon

        # in a previous version GATK RealignerTargetCreator+IndelRealigner were used, now ABRA2 is being used
        gatk -T RealignerTargetCreator -R data/genome/{ref_genome}.fa -I steps/dedup/{bamfile}.nodup.sorted.bam -o steps/realign2/{bamfile}.intervals

	#JAVA_TOOL_OPTIONS="-Xmx300G" abra2 --in steps/dedup/{bamfile}.dedup.sorted.bam --out steps/realign/{bamfile}.abra.bam --ref data/genome/{ref_genome}.fa --threads 27 --tmpdir steps/realign/temp > steps/realign/{bamfile}.abra.log
        #samtools index steps/realign/{bamfile}.abra.bam
        """.format(bamfile=bamfile, ref_genome=ref_genome)
    return inputs, outputs, options, spec

def indel_realignGATK2(bamfile, ref_genome):
    """Target to first create realignment"""
    inputs=["steps/realign2/{}.intervals".format(bamfile),
                "steps/dedup/{}.nodup.sorted.bam".format(bamfile),
                "data/genome/{}.fa".format(ref_genome) ]
    outputs=["steps/realign2/{}.realigned.bam".format(bamfile),
                "steps/realign2/{}.realigned.sorted.bam".format(bamfile),
                "steps/realign2/{}.realigned.sorted.bai".format(bamfile) ]
    options={'cores': 10,
                'memory':'120g',
                'walltime': '10:00:00'
}
    spec = """
        eval "$(conda shell.bash hook)"
        conda activate salmon

        # in a previous version GATK RealignerTargetCreator+IndelRealigner were used, now ABRA2 is being used
        JAVA_TOOL_OPTIONS="-Xmx120G" gatk -T IndelRealigner -R data/genome/{ref_genome}.fa -I steps/dedup/{bamfile}.nodup.sorted.bam -targetIntervals steps/realign2/{bamfile}.intervals -o steps/realign2/{bamfile}.realigned.bam
 
	samtools sort -o steps/realign2/{bamfile}.realigned.sorted.bam steps/realign2/{bamfile}.realigned.bam
        samtools index steps/realign2/{bamfile}.realigned.sorted.bam
        """.format(bamfile=bamfile, ref_genome=ref_genome)
    return inputs, outputs, options, spec

def snp_callGATK4(bamfile, ref_genome):
    """SNP caller combined with indel realignment"""
    inputs=["steps/dedup/{}.nodup.sorted.bam".format(bamfile),
                "data/genome/{}.fa".format(ref_genome) ]
    outputs=["steps/gatk4/{}.gatk4.g.vcf".format(bamfile)]
    options={'cores': 28,
                'memory':'120g',
                'walltime': '20:00:00'
}
    spec = """
        eval "$(conda shell.bash hook)"
        conda activate gatk

        #gatk4 -R data/genome/{ref_genome}.fa -T HaplotypeCaller -I steps/dedup/{bamfile}.dedup.sorted.bam --dbsnp steps/realign2/{bamfile}.dbSNP.vcf -o steps/realign2/{bamfile}.realigned.gatk4.vcf
	JAVA_TOOL_OPTIONS="-Xmx120G -XX:+UseParallelGC -XX:ParallelGCThreads=28" gatk HaplotypeCaller -R data/genome/{ref_genome}.fa -I steps/dedup/{bamfile}.nodup.sorted.bam -pairHMM FASTEST_AVAILABLE --native-pair-hmm-threads 8 -O steps/gatk4/{bamfile}.gatk4.g.vcf -ERC GVCF -stand-call-conf 10
      """.format(bamfile=bamfile, ref_genome=ref_genome)
    return inputs, outputs, options, spec


def mapDamage(bamfile, ref_genome):
    """Target to run deamination visualization"""
    inputs=["steps/realign/{}.abra.bam".format(bamfile),
                "data/genome/{}.fa".format(ref_genome)]
    outputs=[]  # missing names here !!!
    options={'cores': 2,
                'memory':'100g',
                'walltime': '08:00:00'
}
    spec = """
        eval "$(conda shell.bash hook)"
        conda activate salmon-contamination
        
        #samtools view -b steps/realign/{bamfile}.abra.bam > steps/realign/{bamfile}_corr.abra.bam
        mapDamage -i steps/realign/{bamfile}.abra.bam -r data/genome/{ref_genome}.fa --rescale --merge-reference-sequences -d steps/mapdamage/
       """.format(bamfile=bamfile, ref_genome=ref_genome)
    return inputs, outputs, options, spec










##################################################
##### templates for the second part of the
##### analysis (population genetics)
##### author: Samuele Soraggi
##### date of final version: october 2020
#################################################

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



