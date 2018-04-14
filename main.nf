#!/usr/bin/env nextflow

// Edit nextflow.configuration!

aux=config.aux_location
data=config.data_location
output=config.output_location

large_core=config.large_core
small_core=config.small_core


// Fetch fqs; alternative suffixes
fq_set = Channel.fromPath(data + "reads/*.fastq.gz")
                .map { n -> [ n.getName(), n ] }

// Other files and parameters
adapters = file("auxillary/TruSeq3-SE.fa")

// ** - Fetch reference genome (fa.gz) and gene annotation file (gtf.gz)
release="WBPS9"
species="ascaris_suum"
prjn="PRJNA62057"
prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${release}/species/${species}/${prjn}"

process fetch_reference {

    publishDir "${output}/reference/", mode: 'copy'
    
    output:
        file("reference.fa.gz") into reference

    """
        echo '${prefix}'
        curl ${prefix}/${species}.${prjn}.${release}.genomic.fa.gz > reference.fa.gz

    """
}

log.info """\
         Parasite smRNA-Seq pipeline:
         transcriptome: ${ prefix }
         fqs          : ${ fq_set }
         """
         .stripIndent()

//TRIM READS
process trimmomatic {
    cpus large_core
    tag { name }

    publishDir "output/", mode: 'copy', pattern: '_trimout.txt'

    input:
        set val(name), file(reads) from fq_set

    output:
        file(name_out) into fq_trim
        file("*_trimout.txt") into trim_log

    script:
    name_out = name.replace('.fastq.gz', '_trim.fq.gz')

    """
        trimmomatic SE -phred33 -threads ${large_core} ${reads} ${name_out} ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 &> ${reads}_trimout.txt
    """
}


//INDEX GENOME - BOWTIE
process build_bowtie_index {

    publishDir "${output}/reference/", mode: 'copy'

    cpus large_core

    input:
        file("reference.fa.gz") from reference

    output:
        file "*.ebwt" into bowtie_indices

    """
        zcat reference.fa.gz > reference.fa
        bowtie-build reference.fa ref_bowtie
    """
}






// Get pre-indexed hisat2 reference for mouse genome - DONE MANUALLY
//mm_HS2_genome="ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38.tar.gz"
//mm_HS2_genometran="ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38_tran.tar.gz"
//manually run in ${data}/reference:
//curl ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38_tran.tar.gz > grcm38_tran.tar.gz
//tar -zxvf grcm38_tran.tar.gz

//load indexes
// hs2_indices = Channel.fromPath(data + '/reference/grcm38_tran/*.ht2') //.println()

// ** - Fetch reference genome (fa.gz) and gene annotation file (gtf.gz) - DONE MANUALLY
//mm_genome="ftp://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
//mm_gtf="ftp://ftp.ensembl.org/pub/release-91/gtf/mus_musculus/Mus_musculus.GRCm38.91.gtf.gz"
//curl ${mm_gtf} > geneset.gtf.gz

// mm_gtf = file("${data}/reference/geneset.gtf.gz")


//** - ALIGNMENT AND STRINGTIE (combined)
// process align_stringtie {

//     publishDir "${output}/expression", mode: 'copy'

//     cpus large_core

//     tag { id }

//     input:
//         set val(id), file(forward), file(reverse) from read_pairs
//         file("geneset.gtf.gz") from mm_gtf
//         //file("genome_tran.5.ht2"), file("genome_tran.3.ht2"), file("genome_tran.4.ht2"), file("genome_tran.6.ht2"), file("genome_tran.1.ht2"), file("genome_tran.8.ht2"), file("genome_tran.2.ht2"), file("genome_tran.7.ht2") from hs2_indices

//     output:
//         file "${id}.hisat2_log.txt" into alignment_logs
//         file("${id}/*") into stringtie_exp

//     //script:
//      //   index_base = hs2_indices[0].toString() - ~/.\d.ht2/

//     """
//         hisat2 -p ${large_core} -x '/home/BIOTECH/zamanian/GitHub/AsELV_RNAseq-nf/data/reference/grcm38_tran/genome_tran' -1 ${forward} -2 ${reverse} -S ${id}.sam --rg-id "${id}" --rg "SM:${id}" --rg "PL:ILLUMINA" 2> ${id}.hisat2_log.txt
//         samtools view -bS ${id}.sam > ${id}.unsorted.bam
//         rm *.sam
//         samtools flagstat ${id}.unsorted.bam
//         samtools sort -@ ${large_core} -o ${id}.bam ${id}.unsorted.bam
//         rm *.unsorted.bam
//         samtools index -b ${id}.bam
//         zcat geneset.gtf.gz > geneset.gtf
//         stringtie ${id}.bam -p ${large_core} -G geneset.gtf -A ${id}/${id}_abund.tab -e -B -o ${id}/${id}_expressed.gtf 
//         rm *.bam
//         rm *.bam.bai
//         rm *.gtf
//     """
// }



//comment out until all else finished
// prepDE = file("${aux}/scripts/prepDE.py")

// process stringtie_table_counts {

//     echo true

//     publishDir "${output}/diffexp", mode: 'copy'

//     cpus small_core

//     output:
//         file ("gene_count_matrix.csv") into gene_count_matrix
//         file ("transcript_count_matrix.csv") into transcript_count_matrix

//     """
//         python ${prepDE} -i ${output}/expression -l 140 -g gene_count_matrix.csv -t transcript_count_matrix.csv

//     """
// }








// // // // // // // // // // // // // // // // // IGNORE BELOW

// ** - Recurse through subdirectories to get all fastqs
// fq_set = Channel.fromPath(data + "fq/*.fastq.gz")
//                 .map { n -> [ n.getName(), n ] }

// SKIP TRIMMING (READS ARE ALREADY TRIMMED)
// process trim {

//     tag { fq_id }

//     publishDir "${data}/fq_trim/", mode: 'move'

//     input:
//         set fq_id, file(forward), file(reverse) from read_pairs

//     output:
//         set file("${fq_id}_1P.fq.gz"), file("${fq_id}_2P.fq.gz") into trim_output
//         //file "${fq_id}.trim_log.txt" into trim_logs

//     """
//     trimmomatic PE -threads ${large_core} $forward $reverse -baseout ${fq_id}.fq.gz ILLUMINACLIP:/home/linuxbrew/.linuxbrew/Cellar/trimmomatic/0.36/share/trimmomatic/adapters/TruSeq3-PE.fa:2:80:10 MINLEN:75 
//     rm ${fq_id}_1U.fq.gz
//     rm ${fq_id}_2U.fq.gz
//     """

// }


// process trimmomatic {

//     cpus small_core

//     tag { name }

//     input:
//         set val(name), file(reads) from fq_set

//     output:
//         file(name_out) into trimmed_reads

//     script:
//     name_out = name.replace('.fastq.gz', '_trim.fq.gz')

//     """
//         trimmomatic SE -phred33 -threads ${small_core} ${reads} ${name_out} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 
//     """
// }







// process stringtie_counts {

//     publishDir "output/expression", mode: 'copy'

//     cpus small_core

//     tag { srid }

//     input:
//         set val(srid), file(bam), file(bai) from hisat2_bams
//         file("geneset.gtf.gz") from geneset_stringtie.first()

//     output:
//         file("${srid}/*") into stringtie_exp

//     """ 
//         zcat geneset.gtf.gz > geneset.gtf
//         stringtie -p ${small_core} -G geneset.gtf -A ${srid}/${srid}_abund.tab -e -B -o ${srid}/${srid}_expressed.gtf ${bam}
//     """
// }



// // ** - ALIGNMENT
// process align {

//     cpus small_core

//     tag { srid }

//     input:
//         set val(srid), file(forward), file(reverse) from read_pairs
//         file hs2_indices from hs2_indices.first()

//     output:
//         set val(srid), file("${srid}.bam"), file("${srid}.bam.bai") into hisat2_bams
//         file "${srid}.hisat2_log.txt" into alignment_logs

//     script:
//         index_base = hs2_indices[0].toString() - ~/.\d.ht2/

//     """
//         hisat2 -p ${small_core} -x $index_base -1 ${forward} -2 ${reverse} -S ${srid}.sam --rg-id "${srid}" --rg "SM:${srid}" --rg "PL:ILLUMINA" 2> ${srid}.hisat2_log.txt
//         samtools view -bS ${srid}.sam > ${srid}.unsorted.bam
//         samtools flagstat ${srid}.unsorted.bam
//         samtools sort -@ ${small_core} -o ${srid}.bam ${srid}.unsorted.bam
//         samtools index -b ${srid}.bam
//         rm *sam
//         rm *unsorted.bam

//     """
// }



// process stringtie_counts {

//     publishDir "output/expression", mode: 'copy'

//     cpus small_core

//     tag { srid }

//     input:
//         set val(srid), file(bam), file(bai) from hisat2_bams
//         file("geneset.gtf.gz") from geneset_stringtie.first()

//     output:
//         file("${srid}/*") into stringtie_exp

//     """ 
//         zcat geneset.gtf.gz > geneset.gtf
//         stringtie -p ${small_core} -G geneset.gtf -A ${srid}/${srid}_abund.tab -e -B -o ${srid}/${srid}_expressed.gtf ${bam}
//     """
// }



// prepDE = file("auxillary/scripts/prepDE.py")

// process stringtie_table_counts {

//     echo true

//     publishDir "output/diffexp", mode: 'copy'

//     cpus small_core

//     tag { sample_id }

//     input:
//         val(sample_file) from stringtie_exp.toSortedList()

//     output:
//         file ("gene_count_matrix.csv") into gene_count_matrix
//         file ("transcript_count_matrix.csv") into transcript_count_matrix

//     """
//         for i in ${sample_file.flatten().join(" ")}; do
//             bn=`basename \${i}`
//             full_path=`dirname \${i}`
//             sample_name=\${full_path##*/}
//             echo "\${sample_name} \${i}"
//             mkdir -p expression/\${sample_name}
//             ln -s \${i} expression/\${sample_name}/\${bn}
//         done;
//         python ${prepDE} -i expression -l 50 -g gene_count_matrix.csv -t transcript_count_matrix.csv

//     """
// }
