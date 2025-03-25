//Quality Check FastQC and MultiQC

include {FASTQC; MULTIQC} from "../modules/fastQC"

workflow check_QC {
    take:
        step
        reads_ch

    main:
        FASTQC(step, reads_ch)
        
        input_multiqc = FASTQC.out.collect()
        MULTIQC(step, input_multiqc)
}