/*
===========================================
 * Import processes from modules
===========================================
 */


include { SAMTOOLS_FASTQ    } from '../../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_SORT     } from '../../../modules/nf-core/samtools/sort/main'
include { FASTPLONG         } from '../../../modules/nf-core/fastplong/main'
include { SAMTOOLS_IMPORT   } from '../../../modules/nf-core/samtools/import/main'
include { MODKIT_REPAIR     } from '../../../modules/nf-core/modkit/repair/main'
include { RENAME_FASTQ      } from '../../../modules/local/rename_fastq/main'
include { SAMTOOLS_RESET    } from '../../../modules/local/samtools/reset/main'

/*
===========================================
 * Workflows
===========================================
 */


workflow ONT_TRIM_REPAIR {
    take:
    input

    main:

    versions = Channel.empty()

    // Create samtools sort input
    input
        .map { meta, modbam, _ref -> [meta, modbam] }
        .set { ch_sort_in }


    // Create inputs for ref
    input
        .map { meta, _modbam, ref -> [meta, ref] }
        .set { ch_ref_in }


    if (params.reset) {

        SAMTOOLS_RESET(ch_sort_in)
        SAMTOOLS_RESET.out.unaligned_bam
                            . set { ch_reset_bam }

        versions = versions.mix(SAMTOOLS_RESET.out.versions.first())

        SAMTOOLS_SORT(ch_reset_bam, [[],[]], 'bai')

    } else {

        SAMTOOLS_SORT(ch_sort_in, [[],[]], 'bai')

    }

    versions = versions.mix(SAMTOOLS_SORT.out.versions.first())

    // set input to samtools fastq
    SAMTOOLS_SORT.out.bam
        .map { meta, bam -> [meta, bam] }
        .set { fastq_input }

    SAMTOOLS_FASTQ(fastq_input, [])

    versions = versions.mix(SAMTOOLS_FASTQ.out.versions.first())

    RENAME_FASTQ(SAMTOOLS_FASTQ.out.other)

    // Expect: InputFq, adapters_fa, discard_passed, save_failed
    FASTPLONG(RENAME_FASTQ.out.rename_fastq, "${projectDir}/assets/to_remove.fa", false, true)

    // FIXME: fastplong moved to 'topic' -> bellow not working
    //versions = versions.mix(FASTPLONG.out.versions_fastplong.first())
    FASTPLONG.out.log.set { trim_log }

    SAMTOOLS_IMPORT(FASTPLONG.out.reads)
    // Emit also trimmed FASTQ (not repaired yet), as a uBAM
    SAMTOOLS_IMPORT.out.bam
        .join(ch_ref_in)
        .map { meta, trimmed_bam, ref -> [meta, trimmed_bam, ref] }
        .set { not_repaired }


    versions = versions.mix(SAMTOOLS_IMPORT.out.versions.first())

    // Prepare input for modkit repair
    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_IMPORT.out.bam)
        .map { meta, before_trim, after_trim -> [meta, before_trim, after_trim] }
        .set { ch_repair_in }

    MODKIT_REPAIR(ch_repair_in)

    versions = versions.mix(MODKIT_REPAIR.out.versions.first())

    // Prepare input for alignment step


    MODKIT_REPAIR.out.bam
        .join(ch_ref_in)
        .map { meta, trimmed_bam, ref -> [meta, trimmed_bam, ref] }
        .set { dorado_in }

    emit:
    dorado_in
    not_repaired  // FASTQ produced by trimmer, imported as uBAM
    versions
    trim_log
}
