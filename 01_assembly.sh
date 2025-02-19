#!/bin/bash
#enable the conda activate in bash
#code from https://stackoverflow.com/questions/60303997/activating-conda-environment-from-bash-script
eval "$(conda shell.bash hook)"
source ~/miniconda3/etc/profile.d/conda.sh

input_dir=$1    #path to directory with fastq files 
output_name=$2  #output directory name
ref_dir=$3      #path to diamond database folder
SPades=$4       #path of SPAdes  

#####
conda activate General_env
for sample_ in $(ls ${input_dir}/*);
do
    #fl_0 -> path_name
    out_dir=$(echo $sample_ | rev | cut -d '/' -f1 | rev | cut -d '.' -f1)
    path_name=$(echo ${output_name}/${out_dir})
    mkdir -p $path_name
 
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #1. Sort reads into bucket
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    zcat $sample_ |  NanoFilt --headcrop 50 --tailcrop 50 | gzip -9 > ${path_name}/trimmed_raw_reads.fq.gz
    diamond blastx --threads 32 -d ${ref_dir}/ref_AA -q ${path_name}/trimmed_raw_reads.fq.gz -o ${path_name}/match.tsv
    cat ${path_name}/match.tsv  | awk '{ if ($3 >= 85 && $4 >= 100) print $1 }' | sort | uniq > ${path_name}/match_read.txt
    seqkit grep -f ${path_name}/match_read.txt $sample_ | gzip -9 > ${path_name}/viral_readT.fq.gz
    read_count=$(seqkit fx2tab ${path_name}/viral_readT.fq.gz | wc -l)

    if [[ $read_count -gt 30000 ]]; then

        ## "downsampling is performed here"
        #1.1 map to reference DNA 
        minimap2 -a ${ref_dir}/ref_DNA.fasta  ${path_name}/viral_readT.fq.gz >  ${path_name}/aln.sam
        samtools sort  ${path_name}/aln.sam -o  ${path_name}/aln.bam
        samtools index  ${path_name}/aln.bam

        #1.2 select the mapped reference with highest frequency of reads
        best_index=$(samtools idxstat  ${path_name}/aln.bam | sort -k3  -nr | head -n +1 | cut -f1)
        samtools view -b ${path_name}/aln.bam $best_index > ${path_name}/sel.bam
        samtools index ${path_name}/sel.bam
        samtools view -F 256 -bo ${path_name}/sel_filter.bam ${path_name}/sel.bam
        samtools index ${path_name}/sel_filter.bam

        #1.3 downsampling by usning rasusa databse (bbnorm could be used in the future version)
        rasusa aln --coverage 1000 ${path_name}/sel_filter.bam | samtools sort -o  ${path_name}/downsample.bam
        samtools fastq  ${path_name}/downsample.bam | gzip -9 >  ${path_name}/viral_read.fq.gz
        seqkit stat ${path_name}/viral_read.fq.gz

    else
        ## if reads are less than 30k, then do nothing.
        mv ${path_name}/viral_readT.fq.gz ${path_name}/viral_read.fq.gz
    
    fi
   
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #2. Denovo assembly with different parameters
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    conda activate SPades_env
    for kmer in 33 55 77 99 111 127;
    do      
        for length in 0 100 300 500;
        do
           for qual in 10 15 20;
            do     

                    #2.1 create dir for SPades results
                    fl=$(echo ${output_name}/${out_dir}/${length}_${kmer}_${qual})
                    mkdir -p $fl
                   
                    #2.2 QC before assembling
                    zcat ${path_name}/viral_read.fq.gz | NanoFilt -q ${qual} -l ${length} | gzip -9 > $fl/read_qc.fq.gz
                   
                    #2.3 perform genome assembly
                    $SPades -s $fl/read_qc.fq.gz -k $kmer -o $fl/assembly --careful --threads 8 --memory 8
                    full_name=$(echo ${Bname}:${length}:${qual}:${kmer})
                    sed 's/NODE/'$full_name'/g' $fl/assembly/scaffolds.fasta | cut -d '_' -f1 | seqkit rename > $fl/assembly/scaffold_rename.fa
                   
              done  
          done
    done

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #3. Curated results
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    conda activate General_env #move back to general enviroment
    #3.1 select the best mapped reference information
    ensemble_dir=$(echo ${output_name}/${out_dir}/ensemble)
    mkdir -p $ensemble_dir

    ##3.1.2 get the reference id with highest frequency and get the sequence in fasta
    ### Amino acid sequence
    ref_name=$(cat ${fl_0}/match.tsv |  awk '{ if ($3 >= 85 && $4 >= 100) print $2 }' | sort | uniq -c | \
        sort -k1 -n | tail -n1 | awk '{ print $2 }')
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  end {printf("\n");}' ${ref_dir}/ref_AA.fasta | \
        grep $ref_name -A1 > ${ensemble_dir}/ref_AA_best.fa
    
    front_ref=$(cat ${ensemble_dir}/ref_AA_best.fa | seqkit subseq -r 1:4 | grep -v ">")    
    back_ref=$(cat ${ensemble_dir}/ref_AA_best.fa | seqkit subseq -r -4:-1 | grep -v ">")
    length_ref=$(cat ${ensemble_dir}/ref_AA_best.fa | seqkit stat | awk '{ print $5 }' | tail -n1 | sed 's\,\\g')

    ### DNA sequence
    best_AA_id=$(grep ">" ${ensemble_dir}/ref_AA_best.fa | cut -d ' ' -f1 | sed 's\>\\g')
    best_DNA_id=$(grep ${best_AA_id} ${ref_dir}/cross_id.tsv | cut -f2)
    seqkit fx2tab ${ref_dir}/ref_DNA.fasta | grep ${best_DNA_id} | seqkit tab2fx > ${ensemble_dir}/ref_DNA_best.fa

    
    #3.2 combine contigs from all parameter sets
    cat ${output_name}/${out_dir}/*/assembly/scaffold_rename.fa | \
        seqkit rename | grep "\S" > ${ensemble_dir}/all_contigs.fa
   
    #3.3 select the contig with full length of AA 
    seqkit translate -f6 -s ${ensemble_dir}/all_contigs.fa | seqkit seq -m $length_ref > ${ensemble_dir}/contigs_full_length_AA.fa
    common_seq=$(seqkit fx2tab ${ensemble_dir}/contigs_full_length_AA.fa |\
                 cut -f2 | sort |uniq -c | grep $front_ref | grep $back_ref | sort -k1 -n | tail -n1 | awk '{ print $2 }')


    seqkit fx2tab ${ensemble_dir}/contigs_full_length_AA.fa | grep $common_seq | seqkit tab2fx > ${ensemble_dir}/QC_contigs.fa
    cat ${ensemble_dir}/QC_contigs.fa | grep ">" | cut -d"_" -f1 | sed 's\>\\g' > ${ensemble_dir}/QC_contigs_id.fa

    ## Plot the AA sequence in each full-length-AA contigs
    Rscript ./02_plot_AA.R ${ensemble_dir}
   
    #3.4 curate UTR region of sequences from (3.3)
    seqkit grep -f ${ensemble_dir}/QC_contigs_id.fa ${ensemble_dir}/all_contigs.fa > ${ensemble_dir}/filter_contigs.fa
    Rscript ./03_flip_strand.R ${ensemble_dir} ${out_dir} ${ref_dir}
   
    #3.5 calculate read depth
    if [[ $read_count -gt 30000 ]]; then

        #depth based on downsampling
        minimap2 -a $ensemble_dir/final_consensus_DNA.fa ${fl_0}/viral_read.fq.gz > $ensemble_dir/map_${out_dir}.sam
        samtools sort $ensemble_dir/map_${out_dir}.sam -o $ensemble_dir/map_${out_dir}.bam
        samtools index $ensemble_dir/map_${out_dir}.bam
        samtools depth $ensemble_dir/map_${out_dir}.bam > $ensemble_dir/map_${out_dir}_depth.txt

        
        #depth based on total read
        minimap2 -a $ensemble_dir/final_consensus_DNA.fa ${fl_0}/viral_readT.fq.gz > $ensemble_dir/Tmap_${out_dir}.sam
        samtools sort $ensemble_dir/Tmap_${out_dir}.sam -o $ensemble_dir/Tmap_${out_dir}.bam
        samtools index $ensemble_dir/Tmap_${out_dir}.bam
        samtools depth $ensemble_dir/Tmap_${out_dir}.bam > $ensemble_dir/Tmap_${out_dir}_depth.txt
    
      
    else

        #depth based on total read
        minimap2 -a $ensemble_dir/final_consensus_DNA.fa ${fl_0}/viral_read.fq.gz > $ensemble_dir/Tmap_${out_dir}.sam
        samtools sort $ensemble_dir/Tmap_${out_dir}.sam -o $ensemble_dir/Tmap_${out_dir}.bam
        samtools index $ensemble_dir/Tmap_${out_dir}.bam
        samtools depth $ensemble_dir/Tmap_${out_dir}.bam > $ensemble_dir/Tmap_${out_dir}_depth.txt
    
    fi
    
    
          
done



