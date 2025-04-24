#!/usr/bin/env nextflow

process Decode {
    publishDir "${params.outdir}/${index_step}_demux", mode: 'copy', overwrite: true 
    label 'high'
    //module 'PyTorch/1.12.0-foss-2022a-CUDA-11.7.0'
    tag "${sample.sample}" 
    
    input:
    val (index_step)
    tuple val(sample), path(pe_reads)
    
    output:
    tuple val(sample), file("*.fastq.gz")
    
    script:
    if ("$sample.type" == 42000)
        """
        Demux.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N42}  --out-prefix ${sample.sample} 
        """
    else if ("$sample.type" == 21000)
        """
        Demux.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N21} --out-prefix ${sample.sample} 
        """
    else
        """
        Demux.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} 
        """

}

process Decode_batch {
    publishDir "${params.outdir}/${index_step}_demux", mode: 'copy', overwrite: true 
    label 'high'
    //module 'PyTorch/1.12.0-foss-2022a-CUDA-11.7.0'
    //conda "${moduleDir}/pytorch_cuda_env.yml" to do, still not working
    tag "${sample.sample}" 
    
    input:
    val (index_step)
    tuple val(sample), path(pe_reads)
    
    output:
    tuple val(sample), file("*_demux_{R1,R2}.fastq.gz")
    
    script:
    if ("$sample.type" == 42000)
        """
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N42} --out-prefix ${sample.sample} -C 0 --Seg_id 1/4 &
            pid0=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N42} --out-prefix ${sample.sample} -C 1 --Seg_id 2/4 &
            pid1=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N42} --out-prefix ${sample.sample} -C 0 --Seg_id 3/4 &
            pid2=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N42} --out-prefix ${sample.sample} -C 1 --Seg_id 4/4 &
            pid3=\$!
            wait \$pid0
            wait \$pid1
            wait \$pid2
            wait \$pid3
            cat *_R1.fastq.gz > ${sample.sample}_demux_R1.fastq.gz
            cat *_R2.fastq.gz > ${sample.sample}_demux_R2.fastq.gz
        """
    else if ("$sample.type" == 21000)
        """
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N21} --out-prefix ${sample.sample} -C 0 --Seg_id 1/4 &
            pid0=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N21} --out-prefix ${sample.sample} -C 1 --Seg_id 2/4 &
            pid1=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N21} --out-prefix ${sample.sample} -C 0 --Seg_id 3/4 &
            pid2=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N21} --out-prefix ${sample.sample} -C 1 --Seg_id 4/4 &
            pid3=\$!
            wait \$pid0
            wait \$pid1
            wait \$pid2
            wait \$pid3
            cat *_R1.fastq.gz > ${sample.sample}_demux_R1.fastq.gz
            cat *_R2.fastq.gz > ${sample.sample}_demux_R2.fastq.gz
        """
    else
        """
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 0 --Seg_id 1/4 &
            pid0=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 1 --Seg_id 2/4 &
            pid1=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 0 --Seg_id 3/4 &
            pid2=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 1 --Seg_id 4/4 &
            pid3=\$!
            wait \$pid0
            wait \$pid1
            wait \$pid2
            wait \$pid3
            cat *_R1.fastq.gz > ${sample.sample}_demux_R1.fastq.gz
            cat *_R2.fastq.gz > ${sample.sample}_demux_R2.fastq.gz
        """
}


process Decode_batch_4{
    publishDir "${params.outdir}/${index_step}_demux", mode: 'copy', overwrite: true 
    label 'high'
    //module 'PyTorch/1.12.0-foss-2022a-CUDA-11.7.0'
    //conda "${moduleDir}/pytorch_cuda_env.yml" to do, still not working
    tag "${sample.sample}" 
    
    input:
    val (index_step)
    tuple val(sample), path(pe_reads)
    
    output:
    tuple val(sample), file("*_demux_{R1,R2}.fastq.gz")
    
    script:
    if ("$sample.type" == 69000)
        """
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N69} --out-prefix ${sample.sample} -C 0 --Seg_id 1/8 &
            pid0=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N69} --out-prefix ${sample.sample} -C 1 --Seg_id 2/8 &
            pid1=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N69} --out-prefix ${sample.sample} -C 2 --Seg_id 3/8 &
            pid2=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N69} --out-prefix ${sample.sample} -C 3 --Seg_id 4/8 &
            pid3=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N69} --out-prefix ${sample.sample} -C 0 --Seg_id 5/8 &
            pid4=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N69} --out-prefix ${sample.sample} -C 1 --Seg_id 6/8 &
            pid5=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N69} --out-prefix ${sample.sample} -C 2 --Seg_id 7/8 &
            pid6=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N69} --out-prefix ${sample.sample} -C 3 --Seg_id 8/8 &
            pid7=\$!
            wait \$pid0
            wait \$pid1
            wait \$pid2
            wait \$pid3
            wait \$pid4
            wait \$pid5
            wait \$pid6
            wait \$pid7
            cat *_R1.fastq.gz > ${sample.sample}_demux_R1.fastq.gz
            cat *_R2.fastq.gz > ${sample.sample}_demux_R2.fastq.gz
        """
    else if ("$sample.type" == 34500)
        """
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N34} --out-prefix ${sample.sample} -C 0 --Seg_id 1/8 &
            pid0=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N34} --out-prefix ${sample.sample} -C 1 --Seg_id 2/8 &
            pid1=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N34} --out-prefix ${sample.sample} -C 2 --Seg_id 3/8 &
            pid2=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N34} --out-prefix ${sample.sample} -C 3 --Seg_id 4/8 &
            pid3=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N34} --out-prefix ${sample.sample} -C 0 --Seg_id 5/8 &
            pid4=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N34} --out-prefix ${sample.sample} -C 1 --Seg_id 6/8 &
            pid5=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N34} --out-prefix ${sample.sample} -C 2 --Seg_id 7/8 &
            pid6=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N34} --out-prefix ${sample.sample} -C 3 --Seg_id 8/8 &
            pid7=\$!
            wait \$pid0
            wait \$pid1
            wait \$pid2
            wait \$pid3
            wait \$pid4
            wait \$pid5
            wait \$pid6
            wait \$pid7
            cat *_R1.fastq.gz > ${sample.sample}_demux_R1.fastq.gz
            cat *_R2.fastq.gz > ${sample.sample}_demux_R2.fastq.gz
        """
    else
        """
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 0 --Seg_id 1/8 &
            pid0=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 1 --Seg_id 2/8 &
            pid1=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 2 --Seg_id 3/8 &
            pid2=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 3 --Seg_id 4/8 &
            pid3=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 0 --Seg_id 5/8 &
            pid4=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 1 --Seg_id 6/8 &
            pid5=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 2 --Seg_id 7/8 &
            pid6=\$!
            Demux_batch.py -R1 ${pe_reads[0]} -R2 ${pe_reads[1]} -N ${sample.type} -M ${params.M} -B ${sample.bar} --Ntriage ${params.Ntriage} --Nthresh ${params.N3} --out-prefix ${sample.sample} -C 3 --Seg_id 8/8 &
            pid7=\$!
            wait \$pid0
            wait \$pid1
            wait \$pid2
            wait \$pid3
            wait \$pid4
            wait \$pid5
            wait \$pid6
            wait \$pid7
            cat *_R1.fastq.gz > ${sample.sample}_demux_R1.fastq.gz
            cat *_R2.fastq.gz > ${sample.sample}_demux_R2.fastq.gz
        """
}




