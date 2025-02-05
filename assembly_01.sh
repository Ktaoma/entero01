#!/bin/bash

input_dir=$1    #directory with fastq file 
output_name=$2  #output directory name
ref_dir=$3      #path of RefSeq file 

#####

for sample_ in $(ls ${input_dir}/*);
do
   
    echo $sample_
    out_dir=$(echo $sample_ | rev | cut -d '/' -f1 | rev | cut -d '.' -f1)
    fl_0=$(echo ${output_name}/${out_dir})
    mkdir -p $fl_0
 
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #1. sort reads into bucket
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    zcat $sample_ |  NanoFilt --headcrop 50 --tailcrop 50 | gzip -9 > ${fl_0}/trimmed_raw_reads.fq.gz
    diamond blastx --threads 32 -d ${ref_dir}/ref_AA -q ${fl_0}/trimmed_raw_reads.fq.gz -o ${fl_0}/match.tsv
    cat ${fl_0}/match.tsv  | awk -v sim="$similarity" '{ if ($3 >= 85 && $4 >= 100) print $1 }' | sort | uniq > ${fl_0}/match_read.txt
    seqkit grep -f ${fl_0}/match_read.txt $sample_ | gzip -9 > ${fl_0}/viral_readT.fq.gz
    read_count=$(seqkit fx2tab ${fl_0}/viral_readT.fq.gz | wc -l)

    if [[ $read_count -gt 30000 ]]; then

        echo "downsampling is performed here"
        minimap2 -a ${ref_dir}/ref_DNA.fasta  ${fl_0}/viral_readT.fq.gz >  ${fl_0}/aln.sam
        samtools sort  ${fl_0}/aln.sam -o  ${fl_0}/aln.bam
        samtools index  ${fl_0}/aln.bam

        best_index=$(samtools idxstat  ${fl_0}/aln.bam | sort -k3  -nr | head -n +1 | cut -f1)
        samtools view -b ${fl_0}/aln.bam $best_index > ${fl_0}/sel.bam
        samtools index ${fl_0}/sel.bam
        samtools view -F 256 -bo ${fl_0}/sel_filter.bam ${fl_0}/sel.bam
        samtools index ${fl_0}/sel_filter.bam
        
        rasusa aln --coverage 1000 ${fl_0}/sel_filter.bam | samtools sort -o  ${fl_0}/downsample.bam
        samtools fastq  ${fl_0}/downsample.bam | gzip -9 >  ${fl_0}/viral_read.fq.gz
        seqkit stat ${fl_0}/viral_read.fq.gz

    else

        mv ${fl_0}/viral_readT.fq.gz ${fl_0}/viral_read.fq.gz
    
    fi
   
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #2. denovo with diferent paramter
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    for kmer in 33 55 77 99 111 127;
    do      
        for length in 0 100 300 500;
        do
           for qual in 10 15 20;
            do      
             
                    echo #Prep dir
                    fl=$(echo ${output_name}/${out_dir}/${length}_${kmer}_${qual})
                    mkdir -p $fl
                   
                    echo #Quality control
                    zcat ${fl_0}/viral_read.fq.gz | NanoFilt -q ${qual} -l ${length} | gzip -9 > $fl/read_qc.fq.gz
                   
                    echo #Denovo Assembly
                    ../../setup/SPAdes-4.0.0-Linux/bin/spades.py -s $fl/read_qc.fq.gz -k $kmer -o $fl/assembly --careful --threads 8 --memory 8
                    full_name=$(echo ${Bname}:${length}:${qual}:${kmer})
                    sed 's/NODE/'$full_name'/g' $fl/assembly/scaffolds.fasta | cut -d '_' -f1 | seqkit rename > $fl/assembly/scaffold_rename.fa
                   
              done  
          done
    done

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ###3. first checkpoint: ensemble part
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ensemble_dir=$(echo ${output_name}/${out_dir}/ensemble)
    mkdir -p $ensemble_dir
   
    #3.1 select the best mapped reference information
    ref_name=$(cat ${fl_0}/match.tsv |  awk '{ if ($3 >= 85) print $2 }' | sort | uniq -c | \
        sort -k1 -n | tail -n1 | awk '{ print $2 }')
       
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  end {printf("\n");}' ${ref_dir}/ref_AA.fasta | \
        grep $ref_name -A1 > ${ensemble_dir}/ref_best.fa

    best_AA_id=$(grep ">" ${ensemble_dir}/ref_best.fa | cut -d ' ' -f1 | sed 's\>\\g')
    best_DNA_id=$(grep ${best_AA_id} ${ref_dir}/cross_id.tsv | cut -f2)
    seqkit fx2tab ${ref_dir}/ref_DNA.fasta | grep ${best_DNA_id} | seqkit tab2fx > ${ensemble_dir}/ref_DNA_best.fa

     
    front_ref=$(cat ${ensemble_dir}/ref_best.fa | seqkit subseq -r 1:4 | grep -v ">")    
    back_ref=$(cat ${ensemble_dir}/ref_best.fa | seqkit subseq -r -4:-1 | grep -v ">")
    length_ref=$(cat ${ensemble_dir}/ref_best.fa | seqkit stat | awk '{ print $5 }' | tail -n1 | sed 's\,\\g')

    ##3.2 combine contigs from all parameter sets
    cat ${output_name}/${out_dir}/*/assembly/scaffold_rename.fa | \
        seqkit rename | grep "\S" > ${ensemble_dir}/all_dna.fa
   
    #3.3 select the contig with full length of AA
    seqkit translate -f6 -s ${ensemble_dir}/all_dna.fa | seqkit seq -m $length_ref > ${ensemble_dir}/all_dna_full_length_AA.fa

    #3.4 pass the first checkpoint (if not go to second or third)
    success_=$(awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  end {printf("\n");}' ${ensemble_dir}/all_dna_full_length_AA.fa | grep -v ">" | grep $front_ref | grep $back_ref | wc -l)

    if [[ "$success_" -gt 0 ]]; then
        # pass -> trim UTR
        #3.5 filter AA
        #3.5.1 plot
        Rscript ./filter_AA.R ${ensemble_dir}
       
        #3.5.2 select common one
        common_seq=$(seqkit fx2tab ${ensemble_dir}/all_dna_full_length_AA.fa |\
                     cut -f2 | sort |uniq -c | grep $front_ref | grep $back_ref | sort -k1 -n | tail -n1 | awk '{ print $2 }')
                     
        seqkit fx2tab ${ensemble_dir}/all_dna_full_length_AA.fa | grep $common_seq | seqkit tab2fx > ${ensemble_dir}/optimal_read_AA.fa
        cat ${ensemble_dir}/optimal_read_AA.fa | grep ">" | cut -d"_" -f1 | sed 's\>\\g' > ${ensemble_dir}/optimal_read_id.fa
       

        #3.6 filter UTR        
        #3.6.1: merge+flip strand
        seqkit grep -f ${ensemble_dir}/optimal_read_id.fa ${ensemble_dir}/all_dna.fa > ${ensemble_dir}/optimal_read_DNA.fa
        echo "----------------------------"
        Rscript ./flip_strand.R ${ensemble_dir} ${out_dir} ${ref_dir}
       
        #3.6.3: depth
        #depth down
        minimap2 -a $ensemble_dir/final_consensus_DNA.fa ${fl_0}/viral_read.fq.gz > $ensemble_dir/map_${out_dir}.sam
        samtools sort $ensemble_dir/map_${out_dir}.sam -o $ensemble_dir/map_${out_dir}.bam
        samtools index $ensemble_dir/map_${out_dir}.bam
        samtools depth $ensemble_dir/map_${out_dir}.bam > $ensemble_dir/map_${out_dir}_depth.txt
       
        #depth total
        minimap2 -a $ensemble_dir/final_consensus_DNA.fa ${fl_0}/viral_readT.fq.gz > $ensemble_dir/Tmap_${out_dir}.sam
        samtools sort $ensemble_dir/Tmap_${out_dir}.sam -o $ensemble_dir/Tmap_${out_dir}.bam
        samtools index $ensemble_dir/Tmap_${out_dir}.bam
        samtools depth $ensemble_dir/Tmap_${out_dir}.bam > $ensemble_dir/Tmap_${out_dir}_depth.txt
       
    else

        echo "fail -> fix frame"
        mkdir -p ${ensemble_dir}/checkpoint02

        zcat ${fl_0}/viral_read* | seqkit rmdup | seqkit fq2fa > ${ensemble_dir}/checkpoint02/rep.fa
        cat ${ensemble_dir}/checkpoint02/rep.fa ${ensemble_dir}/all_dna.fa  > ${ensemble_dir}/checkpoint02/merge_final.fa
        minimap2 -a ${ref_dir}/ref_DNA.fasta ${ensemble_dir}/checkpoint02/merge_final.fa > ${ensemble_dir}/checkpoint02/map_ref.sam
        samtools sort ${ensemble_dir}/checkpoint02/map_ref.sam -o ${ensemble_dir}/checkpoint02/map_ref.bam
        samtools index ${ensemble_dir}/checkpoint02/map_ref.bam
        
        id_=$(samtools idxstats ${ensemble_dir}/checkpoint02/map_ref.bam | sort -k3 -nr | head -n +1 | cut -f1)
        samtools consensus --mode simple --min-depth 1 ${ensemble_dir}/checkpoint02/map_ref.bam -o ${ensemble_dir}/checkpoint02/map_ref_consensus.fa
        
        #
        seqkit fx2tab ${ensemble_dir}/checkpoint02/map_ref_consensus.fa | grep ${id_} | seqkit tab2fx > ${ensemble_dir}/checkpoint02/map_ref_consensus2.fa
        seqkit fx2tab ${ref_dir}/ref_DNA.fasta | grep ${id_} | seqkit tab2fx > ${ensemble_dir}/checkpoint02/ref_DNA2.fasta
        #
        
        Rscript ./fix.R ${ensemble_dir}/checkpoint02/ref_DNA2.fasta ${ensemble_dir}/checkpoint02/map_ref_consensus2.fa ${ensemble_dir}/checkpoint02/DNA_fix.fa     
    fi
   
done



