#### PIPELINE - FIRST PART FOR GENOMIC LIBRARY
#### Date: September 2019
#### Organism: Atlantic salmon (Salmo salar)
#### Authors: Belen JM based on genomic pipeline developed by Alice Manuzzi
#### Description:
#### This genomic pipeline analyses data from Illumina sequencing platform
#### and: filter adaptors, detects/removes contamination, maps to the genome,
#### identify duplicates, realign indels, measures degradation.
#### Input: raw data (PE reads), reference genome, list of samples
#### Output: vcf file with SNPs
#### Notes: you need the python software gwf to run this pipeline with the command `gwf run`

#Import gwf and target templates
from gwf import Workflow
import os

from templates import *

gwf = Workflow()

# unzipping the genome
gwf.target_from_template('UnzipGenome',
                         unzip(inputfile='data/genome/GCF_000233375.1_ICSASG_v2_genomic.fa.gz',
                               outputfile='data/genome/GCF_000233375.1_ICSASG_v2_genomic.fa'))

# indexing the genome with bwa
gwf.target_from_template('IndexGenome',
                         bwa_index(ref_genome='GCF_000233375.1_ICSASG_v2_genomic'))


# creating kraken library
if not os.path.exists("data/Kraken_library"):
	os.mkdir("data/Kraken_library")
if not os.path.exists("data/Kraken_library/library_contami"):
	os.mkdir("data/Kraken_library/library_contami")

gwf.target_from_template('Build_Kraken_lib',
                         create_kraken_lib(folder='data/Kraken_library/library_contami'))


# loop over samples:
# (1) filtering
# (2) mapping
# (3) mark duplicates
# (4) realignment

stats_inputs = []

with open('data/mysamples') as mysamples:
        ind_number=0
        for line in mysamples:
                ind_number += 1
                print(line)
                fullpath = line.strip()
                partpath = fullpath.split("/")

                print("IND:", partpath)

                folder = partpath[0]
                indiv = partpath[1]
                print("Folder:", folder)
                print("Indiv:", indiv)

                name_target = indiv.split("_")[0:4] # shorten the name of individual
                print("Name:", name_target)
                sep = "_"
                ind2 = str(sep.join(name_target)) # final short name of the individual:
                                                  # name: river+year+sampleID+libID
                lib_name = name_target[3]
                print(lib_name)
                print(ind2)
                #ind_number += 1
                #ind = line.strip()
                #print(ind)
                #name_target = ind.split("_")[0:4] # shorten the name of individual
                #print(name_target)
                #sep = "_"
                #ind2 = str(sep.join(name_target)) # final short name of the individual:
                                                  # name: river+year+sampleID+libID
                #lib_name = name_target[3]
                #print(lib_name)
                #print(ind2)
                #gwf.target("Filtering_{}".format(name_target2), inputs=[], outputs=[]) << """
#cat {}.""".format(name_target2)

                f1 = gwf.target(
                        name="Filtering1_{}".format(ind2),
                        cores=7,
                        memory="700mb",
                        walltime="02:00:00",
                        inputs=[
                                "data/{}/{}_R1.fastq.gz".format(folder, indiv),
                                "data/{}/{}_R2.fastq.gz".format(folder, indiv)
                        ],
                        outputs=[
                                "steps/filtering/{}.trim1.pair1.truncated.gz".format(ind2),
                                "steps/filtering/{}.trim1.pair2.truncated.gz".format(ind2),
                                "steps/filtering/{}.trim1.discarded.gz".format(ind2),
                                "steps/filtering/{}.trim1.singleton.truncated.gz".format(ind2),
                                "steps/filtering/{}.trim1.pair1.truncated_fastqc.html".format(ind2),
                                "steps/filtering/{}.trim1.pair1.truncated_fastqc.zip".format(ind2),
                                "steps/filtering/{}.trim1.pair2.truncated_fastqc.html".format(ind2),
                                "steps/filtering/{}.trim1.pair2.truncated_fastqc.zip".format(ind2),
                                "steps/filtering/{}.trim1.settings".format(ind2)
                        ]
                ) << """
                cd data/{}/

                AdapterRemoval --file1 {}_R1.fastq.gz --file2 {}_R2.fastq.gz --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --basename {}.trim1 --minlength 30 --trimns --trimqualities --minquality 30 --gzip
                fastqc {}.trim1.pair1.truncated.gz
                fastqc {}.trim1.pair2.truncated.gz

                mv {}.t* ../../steps/filtering/.

                """.format(folder, indiv, indiv, ind2,
                                        ind2,
                                        ind2,
                                        ind2)


                # gwf.target_from_template("Filtering_{}".format(name_target2), adapter_removal(data_dir="data/lib01", ind=ind, output_dir="steps/filtering", adapters=("AGAATTTT", "TTTGTAAAA")))

                f2 = gwf.target(
                        name="Filtering2_{}".format(ind2),
                        cores=4,
                        memory="700mb",
                        walltime="02:00:00",
                        inputs=f1.outputs,
                        outputs=["steps/filtering/{}.trim2.pair1.truncated.fastq.gz".format(ind2), # the name is actually "*truncated.gz"
                                                                                                   # but needs to be the final output
                                "steps/filtering/{}.trim2.pair2.truncated.fastq.gz".format(ind2),
                                "steps/filtering/{}.trim2.discarded.gz".format(ind2),
                                "steps/filtering/{}.trim2.singleton.truncated.gz".format(ind2),
                                "steps/filtering/{}.trim2.pair1.truncated_fastqc.html".format(ind2),
                                "steps/filtering/{}.trim2.pair1.truncated_fastqc.zip".format(ind2),
                                "steps/filtering/{}.trim2.pair2.truncated_fastqc.html".format(ind2),
                                "steps/filtering/{}.trim2.pair2.truncated_fastqc.zip".format(ind2),
                                "steps/filtering/{}.trim2.settings".format(ind2)
                        ]
                ) << """
                cd steps/filtering/

                AdapterRemoval --file1 {}.trim1.pair1.truncated.gz --file2 {}.trim1.pair2.truncated.gz --adapter1  GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG --basename {}.trim2 --minlength 30 --trimns --trimqualities --minquality 30 --gzip

                fastqc {}.trim2.pair1.truncated.gz
                fastqc {}.trim2.pair2.truncated.gz

                """.format(ind2,ind2,ind2,
                        ind2,
                        ind2)

                f3 = gwf.target(
                        name="Filtering3_{}".format(ind2),
                        cores=4,
                        memory="700mb",
                        walltime="02:00:00",
                        inputs=f2.outputs,
                        outputs=[
                                "steps/filtering/{}.trim3.pair1.truncated.fastq.gz".format(ind2),
                                "steps/filtering/{}.trim3.pair2.truncated.fastq.gz".format(ind2),
                                "steps/filtering/{}.trim3.pair1.truncated_fastqc.html".format(ind2),
                                "steps/filtering/{}.trim3.pair1.truncated_fastqc.zip".format(ind2),
                                "steps/filtering/{}.trim3.pair2.truncated_fastqc.html".format(ind2),
                                "steps/filtering/{}.trim3.pair2.truncated_fastqc.zip".format(ind2),
                                ]
                ) << """
                cd steps/filtering/

                #cutadapt raises error if the file is not in fastq.gz format
                mv {}.trim2.pair1.truncated.gz {}.trim2.pair1.truncated.fastq.gz
                mv {}.trim2.pair2.truncated.gz {}.trim2.pair2.truncated.fastq.gz

                cutadapt --nextseq-trim=20 --minimum-length 30 -o {}.trim3.pair1.truncated.fastq.gz -p {}.trim3.pair2.truncated.fastq.gz {}.trim2.pair1.truncated.fastq.gz {}.trim2.pair2.truncated.fastq.gz

                fastqc {}.trim3.pair1.truncated.fastq.gz
                fastqc {}.trim3.pair2.truncated.fastq.gz

                """.format(ind2, ind2,
                                ind2, ind2,
                                ind2,ind2,ind2,ind2,
                                ind2,
                                ind2)
		
                f4 = gwf.target(
			name="Filtering4_{}".format(ind2),
			cores=4,
                        memory="700mb",
                        walltime="02:00:00",
                        inputs=f3.outputs,
                        outputs=[
                                "steps/filtering/{}.trim4.pair1.truncated.gz".format(ind2),
                                "steps/filtering/{}.trim4.pair2.truncated.gz".format(ind2),
                                "steps/filtering/{}.trim4.discarded.gz".format(ind2),
                                "steps/filtering/{}.trim4.singleton.truncated.gz".format(ind2),
                                "steps/filtering/{}.trim4.pair1.truncated_fastqc.html".format(ind2),
                                "steps/filtering/{}.trim4.pair1.truncated_fastqc.zip".format(ind2),
                                "steps/filtering/{}.trim4.pair2.truncated_fastqc.html".format(ind2),
                                "steps/filtering/{}.trim4.pair2.truncated_fastqc.zip".format(ind2),
                                "steps/filtering/{}.trim4.settings".format(ind2)
                                ]
                ) << """
                cd steps/filtering/

                AdapterRemoval --file1 {}.trim3.pair1.truncated.fastq.gz --file2 {}.trim3.pair2.truncated.fastq.gz --adapter1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACACGTTACCATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA --basename {}.trim4 --minlength 30 --trimns --trimqualities --minquality 30 --gzip

                fastqc {}.trim4.pair1.truncated.gz
                fastqc {}.trim4.pair2.truncated.gz


                """.format(ind2, ind2, ind2,
                                ind2,
                                ind2)
		
		# filtering for contamination
                if not os.path.exists('steps/check_contami/'):
                    os.mkdir('steps/check_contami/')
                if not os.path.exists("steps/filter_humanBact/"):
                    os.mkdir('steps/filter_humanBact/')

                #gwf.target_from_template(
    		#	"kraken_report_{}".format(ind2),
     		#	kraken_report(r1='{}.trim4.pair1.truncated.gz'.format(ind2),
     		#		r2='{}.trim4.pair2.truncated.gz'.format(ind2),
		#		ind2='{}'.format(ind2)))

                #gwf.target_from_template(
		#	"check_kraken_{}".format(ind2),
		#	visualize_kraken(ind2='{}'.format(ind2)))

                if not os.path.exists('steps/check_contami/mapping-human'):
                    os.mkdir('steps/check_contami/mapping-human')

                human = gwf.target_from_template(
                        "Remove_human_{}".format(ind2), 
                        remove_contamination(r1='{}.trim4.pair1.truncated.gz'.format(ind2),
                        r2='{}.trim4.pair2.truncated.gz'.format(ind2), 
                        ind2='{}'.format(ind2),
                        organism = "human".format(),
                        path_reads = "steps/filtering".format()))
                
                if not os.path.exists('steps/check_contami/mapping-bacteria'):
                    os.mkdir('steps/check_contami/mapping-bacteria')

                bacteria = gwf.target_from_template(
                        "Remove_bacteria_{}".format(ind2),
                        remove_contamination(r1='{}.nohuman_1.fastq.gz'.format(ind2),
                        r2='{}.nohuman_2.fastq.gz'.format(ind2),
                        ind2='{}'.format(ind2),
                        organism = "bacteria".format(),
                        path_reads = "steps/filter_humanBact".format()))

                # add target for control contamination
		# only start when all individuals are finished
                #control_contamination_files.extend(human.outputs)
                
#gwf.target(
#      name="Control_Contamination", 
#      cores = 2,
#      memory="10g",
#      inputs=control_contamination_files, 
#      outputs=['steps/check_contami/bacteria',
#      'steps/check_contami/final-bacteria-list',
#      'steps/check_contami/subset-library-bacteria.txt',
#      'steps/check_contami/subset-library-bacteria-2.txt',
#      'steps/check_contami/subset-bacteria-temporary.fa']
#) << """
#
#eval "$(conda shell.bash hook)"
#conda activate salmon-contamination

# RUN CONTAMINATION CHECK HERE
#identify the list of bacteria found in the samples and build a subset of the library (temporary file)
#cat steps/check_contami/*_report | grep "bacteria" | awk '$1 != 0.00' | cut -f6 >> steps/check_contami/bacteria
#echo "hola"
#cat steps/check_contami/bacteria | sort -u > steps/check_contami/final-bacteria-list
#for i in $(cat steps/check_contami/final-bacteria-list)
#do 
#grep "$i" data/Kraken_library/library_contami/library/bacteria/library.fna >> steps/check_contami/subset-library-bacteria.txt
#done
#sed 's/>//' steps/check_contami/subset-library-bacteria.txt > steps/check_contami/subset-library-bacteria-2.txt
#Extract sequences with names in file name.lst, one sequence name per line:
#seqtk subseq data/Kraken_library/library_contami/library/bacteria/library.fna steps/check_contami/subset-library-bacteria-2.txt > steps/check_contami/subset-bacteria-temporary.fa
#filterbyname.sh in=data/Kraken_library/library_contami/library/bacteria/library.fa out=steps/check_contami/subset-bacteria-temporary.fa names=steps/check_contami/subset-library-bacteria.txt include=t
#bowtie2-build --wrapper basic-0 -f --large-index -f library.fna bacteria3
#bowtie2-build -f steps/check_contami/subset-bacteria-temporary.fa data/Kraken_library/library_contamin/library/bacteria/bacteria1st
#""".format(ind2)

                # mapping reads
                if not os.path.exists("steps/mapping"):
                    os.mkdir("steps/mapping")

                gwf.target_from_template(
                    "MapReads_{}".format(ind2), 
                    bwa_map(ref_genome='GCF_000233375.1_ICSASG_v2_genomic',
                           r1='{}.nobacteria_1.fastq.gz'.format(ind2),
                           r2='{}.nobacteria_2.fastq.gz'.format(ind2),
                           bamfile='{}'.format(ind2)))

                # ReadGroup
                if not os.path.exists("steps/RG"):
                    os.mkdir("steps/RG")

                gwf.target_from_template(
                        "CreateRG_{}".format(ind2),
                        create_RG(bamfile='{}'.format(ind2),
                                ind_number='{}'.format(ind_number),
                                lib_number='{}'.format(lib_name),
                                ind_name='{}'.format(ind2)))

                # mark duplicates
                if not os.path.exists("steps/dedup"):
                        os.mkdir("steps/dedup")
                duplicates = gwf.target_from_template(
                        "MarkDup_{}".format(ind2),
                        mark_dupli(bamfile='{}'.format(ind2)))

                # indel realignment
                #if not os.path.exists("steps/realign"):
                #        os.mkdir("steps/realign")

                #gwf.target_from_template(
                #        "IndelRealign_{}".format(ind2),
                #        indel_realign(bamfile='{}'.format(ind2),
                #                ref_genome='GCF_000233375.1_ICSASG_v2_genomic'))
                #gwf.target_from_template(
                #        "IndelRealign2_{}".format(ind2),
                #        indel_realignGATK(bamfile='{}'.format(ind2),
                #        ref_genome='GCF_000233375.1_ICSASG_v2_genomic'))
                #gwf.target_from_template(
                #        "IndelRealign3_{}".format(ind2),
                #        indel_realignGATK2(bamfile='{}'.format(ind2),
                #        ref_genome='GCF_000233375.1_ICSASG_v2_genomic'))
                gwf.target_from_template(
                        "SNP_GATK4_{}".format(ind2),
                        snp_callGATK4(bamfile='{}'.format(ind2),
                        ref_genome='GCF_000233375.1_ICSASG_v2_genomic'))
                # MapDamage
                #if not os.path.exists("steps/mapdamage"):
                #    os.mkdir("steps/mapdamage")

                #gwf.target_from_template(
                #        "mapDamage_{}".format(ind2),
                #         mapDamage(bamfile='{}'.format(ind2),
                #               ref_genome='GCF_000233375.1_ICSASG_v2_genomic'))
                stats_inputs.extend(duplicates.outputs)

if not os.path.exists("steps/stats"):
    os.mkdir("steps/stats")


gwf.target(
name="stats",
cores=4,
memory="16g",
walltime="08:00:00",
inputs=stats_inputs,
outputs=["steps/stats/output_bioinfo_Dec.txt"]
) << """
# creating the 350bp-intervals around the bait areas (run only once)
#bedtools slop -i data/baits/baits_BED_exact_withoutSexBaits_sorted.bed -g data/genome/chrom_size2 -b 350 > data/baits/baits_350bp.bed
    
   for i in $(cat data/id_names2)
    do echo "$i"

    # no of reads sequenced
    reads_prefilt=$(zcat steps/filtering/$i.trim4.pair1.truncated.gz | wc -l | awk '{print $1 / 4}')

    # no of reads after after-adaptor filter
    reads_after_adaptor=$(zcat steps/filtering/${i}.trim4.pair1.truncated.gz | wc -l | awk '{print $1 / 4}')
                
    # no of reads after-contamination1 filter
    reads_after_contam_human=$(zcat steps/filter_humanBact/${i}.nohuman_1.fastq.gz | wc -l | awk '{print $1 / 4}')
                
    # no of reads after-contamination2 filter
    reads_after_contam_bact=$(zcat steps/filter_humanBact/${i}.nobacteria_1.fastq.gz | wc -l | awk '{print $1 / 4}')

    # no of reads to map
    reads_map=$(samtools flagstat steps/mapping/${i}.sorted.bam | awk 'NR == 1 {print $1}')

    # no of reads mapping back to genome:
    reads_mapped=$(samtools flagstat steps/mapping/${i}.sorted.bam | awk 'NR == 5 {print $1}')

    # no of reads pair1 mapped:
    reads_pair1=$(samtools flagstat steps/mapping/${i}.sorted.bam | awk 'NR == 7 {print $1}')	

    # no of reads pair1 mapped:
    reads_pair2=$(samtools flagstat steps/mapping/${i}.sorted.bam | awk 'NR == 8 {print $1}')
    
    # no of reads mapping back to genome properly paired
    pairs_map=$(samtools flagstat steps/mapping/${i}.sorted.bam | awk 'NR == 9 {print $1}')

    # singletons:
    singletons=$(samtools flagstat steps/mapping/${i}.sorted.bam | awk 'NR == 11 {print $1}')

    # diff chromiosome:
    diff_chrom=$(samtools flagstat steps/mapping/${i}.sorted.bam | awk 'NR == 12 {print $1}')

    # percentage of reads mapping back to genome
    perc_reads_map=$(samtools flagstat steps/mapping/${i}.sorted.bam | awk -F "[(|%]" 'NR == 5 {print $2}')

    # no of reads that are duplicates, and will be removed
    duplic=$(samtools flagstat steps/dedup/{i}.dedup.sorted.bam | awk 'NR == 4 {print $1}')

    # no of total reads proper pairs after removing duplicates
    final_reads=$(samtools flagstat steps/dedup/{i}.nodup.sorted.bam | awk 'NR == 1 {print $1}')
	
    # no of reads mapping back to baits
    bedtools intersect -c -a data/baits/baits_BED_exact_withoutSexBaits_sorted.bed -b steps/dedup/${i}.nodup.sorted.bam > steps/stats_baits/${i}.intersect
    reads_baits=$(cat steps/stats_baits/${i}.intersect | awk '$10 >= 1' | wc -l)

    # mean of reads per bait
    mean_reads_baits=$(cat steps/stats_baits/${i}.intersect | awk '$10>=1' | awk '{total += $10} END { print total/NR }')

    # no of reads mapping back to baits+350bp
    bedtools intersect -c -a data/baits/baits_350bp.bed -b steps/dedup/${i}.nodup.sorted.bam > steps/stats_baits/${i}.intersect350bp
    reads_baits350=$(cat steps/stats_baits/${i}.intersect350bp | awk '$10 >= 1' | wc -l)

    # mean of reads per bait+350bp
    mean_reads_baits350=$(cat steps/stats_baits/${i}.intersect350bp | awk '$10>=1' | awk '{total += $10} END { print total/NR}')

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$i" "$reads_prefilt" "$reads_after_adaptor" "$reads_after_contam_human" "$reads_after_contam_bact" "$reads_map" "$reads_mapped" "$reads_pair1" "$reads_pair2" "$pairs_map" "$singletons" "$diff_chrom" "$perc_reads_map" "$duplic" "$final_reads" "$reads_baits" "$mean_reads_baits" "$reads_baits350" "$mean_reads_baits350" >> steps/stats/output_bioinfo_Dec.txt

done
"""



#### PIPELINE - SECOND PART FOR POPULATION GENETICS
#### Date: October 2020
#### Organism: Atlantic salmon (Salmo salar)
#### Authors: Samuele Soraggi
#### Description:
#### This genomic pipeline analyses data from Illumina sequencing platformThis pipeline runs using primarily the software ANGSD for low- and mid-coverage NGS data. The specific version of ANGSD used for this analysis is 0.931-14-gb9c8ddd together with htslib version 1.10.2-16-g2b820c2. A number of R libraries and python packages have been used to generate plots (python code in form of jupyter notebooks in the folder `notebooks`). You can install all these softwares by using conda and installing an environment starting from the file part_2_environment.yml

#path for the PIKE mapped to the SALMON genome
ancFile = '../steps/admixture/PIKEmappedToLaks.sorted.bam'

#####Run SAF ANALYSIS on each file
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

######Run 1D sfs    
######REMEMBER: do  `angsd sites index your.file` on the -sites file before running all this
######          so that the SNP list is indexed. this has already been done for you and the
#####           file is in the repository folder already formatted
######
###### Here you need -sites and not -rf as option (-rf is only for reading from bam files)
mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]

for i in filelists:
    gwf.target_from_template( 'sfs1d_5e6_'+str(i) , sfs1d(fileSaf=i,
                                                      memory='200g',
                                                      outFolder='results_5e6',
                                                      sites='pruning/snps_5e6.tsv',
                                                      cores=8 ))

   
######Run 2D sfs    
######REMEMBER: do  `angsd sites index your.file` on the -sites file before running all this
######          so that the SNP list is indexed. this has already been done for you and the
#####           file is in the repository folder already formatted
######
###### Here you need -sites and not -rf as option (-rf is only for reading from bam files)
mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
count_i=0
count_j=0
for i in filelists:
    count_j=0
    for j in filelists:
        if count_j>count_i:
            runName = i+'_'+j
            gwf.target_from_template( 'sfs2d_5e6_'+str(runName) , sfs2d(file1=i,
                                                                file2=j,
                                                                outFolder='results_5e6',
                                                                sites='pruning/snps_5e6.tsv',    
                                                                memory='200g',
                                                                cores=8 ))
        count_j += 1
    count_i += 1


#####Run calculation of Fst
mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
count_i=0
count_j=0
for i in filelists:
    count_j=0
    for j in filelists:
        if count_j>count_i:
            runName = i+'_'+j
            #Prepare each file for Fst calculation
            gwf.target_from_template( 'fst_prep_5e6_'+str(runName) , Fst_prepare(file1=i,
                                                                file2=j,
                                                                outFolder='results_5e6',
                                                                walltime='8:00:00',
                                                                memory='100g',
                                                                cores=8))
            #Run calculation
            gwf.target_from_template( 'fst_calc_5e6_'+str(runName) , Fst_estimate(file1=i,
                                                                              file2=j,
                                                                              outFolder='results_5e6'))
        count_j += 1
    count_i += 1


#####Calculate various thetas each individual. More at the link
##### http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests#Full_command_list_for_below_examples

mypath = './filelists_flt'
filelists = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]

for i in filelists:
    gwf.target_from_template( 'Theta_5e6_'+str(i) , Theta_estimate(file1=i,
                                                                outFolder='results_5e6',
                                                                walltime='2:00:00',
                                                                memory='20g',
                                                                cores=1))



####################### 4 population test analysis
#######################      with D-statistic
#######################           ANGSD

#####generate fasta files for the PIKE reference mapped to SALMON and for a modern individual of high quality.
#####They are used for error rate estimation in all the genomes.

#reference PIKE file mapped to the SALMON genome. This was sorted with samtools to generate the fasta file used. The extension must be .fa and not .fasta
gwf.target_from_template( 'dofasta_pike' , doFasta(bamFile='../steps/admixture/PIKEmappedToLaks.sorted',
                                                    mode='EBD',
                                                    account='salmon',
                                                    sites='baits_350bp_3cols.bed',
                                                    memory='48g',
                                                    walltime='36:00:00',
                                                    cores=8 ))
#A second file used as an example of high quality genome without many errors. We choose a modern individual sequenced at high depth. This was sorted with samtools to generate the fasta file used. The extension must be .fa and not .fasta
gwf.target_from_template( 'dofasta_modern' , doFasta(bamFile='../steps/dedup/VA_16_31_HFH73BBXY.nodup.sorted',
                                                    mode='EBD',
                                                    out='./modern_reference_VA_16_31_HFH73BBXY',
                                                    account='salmon',
                                                    sites='baits_350bp_3cols.bed',
                                                    memory='48g',
                                                    walltime='36:00:00',
                                                    cores=8 ))  

#####Determine error rates for each genome. For thee estimation method consult the D-statistic paper from Soraggi et al 2018.
for i in filelists:
    gwf.target_from_template( 'doAncError_modern_'+i , ancError(bamFiles=i,
                                                         anc_genome='../steps/admixture/PIKEmappedToLaks.sorted.fa',
                                                         ref_genome='./modern_reference_VA_16_31_HFH73BBXY.fa',
                                                         regions='pruning/snps_5e6_formatted.tsv',
                                                         out='out_ancError_modernRef'))


###Run the abbababa (D-statistic) analysis excluding few outliers genomes that were excluded looking at the PCA and calculation of genetic quantities. Probably they were too low quality or catalogued uncorrectly as salmons. For more information consult the D-statistic paper from Soraggi et al 2018.
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
                #calculate ABBA and BABA combinations
                gwf.target_from_template( 'abbababa_pcaoutliers_'+runName , abbababa(file_1=i,
                                                                         file_2=j,
                                                                         file_3=k,
                                                                         out = 'uncorrected_masked_pcaoutliers/',
                                                                         anc_genome=ancFile,
                                                                         cores=8,
                                                                         regions='pruning/snps_5e6_formatted.tsv',
                                                                         memory='64g',
                                                                         walltime='72:00:00'))
                #calculate D-statistic and apply error correction
                gwf.target_from_template( 'Dstat_pcaout_modern_'+runName , Dstat_calc(file_1=i,
                                                                         file_2=j,
                                                                         file_3=k,
                                                                         abbababa_folder = 'uncorrected_masked_pcaoutliers/',
                                                                         out = 'modernRefCorrected',
                                                                         err_folder = 'out_ancError_modernRef'))
            count_k += 1
        count_j += 1
    count_i += 1
