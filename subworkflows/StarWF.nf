include {Star_Align_QSP; index_bam_QSP} from "../modules/Star"

// workflow Alignment_QSP {
//     take:
//         index_step
//         reads_ch
//         genome
//     main:
//         Star_Align_QSP(index_step, reads_ch, genome)
//         index_bam_QSP(index_step, Star_Align_QSP.out.align_bam)
//     emit:
//         bam_file = Star_Align_QSP.out.align_bam
//         log_bam = Star_Align_QSP.out.log_files
//         index = index_bam_QSP.out.index
//         mock = index_bam_QSP.out.mock
// }

workflow Alignment_QSP {
    take:
        index_step
        reads_ch
        genome
    main:
        Star_Align_QSP(index_step, reads_ch, genome)
        index_bam_QSP(index_step, Star_Align_QSP.out.align_bam)
    emit:
        bam_file = [Star_Align_QSP.out.align_bam, index_bam_QSP.out.index]
        // log_bam = Star_Align_QSP.out.log_files
        // index = index_bam_QSP.out.index
        // mock = index_bam_QSP.out.mock
}
