/*
===========================================
 * Import modules
===========================================
 */

include { DORADO_BASECALLER                        } from '../../modules/local/dorado/basecaller/main'
include { DORADO_DEMUX                             } from '../../modules/local/dorado/demux/main'
include { PIGZ_COMPRESS as COMPRESS_DEMUX_SUMMARY  } from '../../modules/nf-core/pigz/compress/main'
include { MINIMAP2_ALIGN as ONT_ALIGN_NOT_REPAIRED } from '../../modules/nf-core/minimap2/align/main'

/*
===========================================
 * Import subworkflows
===========================================
 */

include { FASTQ_UNZIP                      } from './shared_fastqc_unzip/main'
include { ONT_ALIGN                        } from './ont_align/main'
include { ONT_TRIM_REPAIR                  } from './ont_trim_repair/main'
include { BED2BEDGRAPH                     } from './shared_bed2bedgraph/main'
include { INDEX_MODKIT_PILEUP              } from './shared_modkit_pileup/main'
include { ONT_FIBERSEQ                     } from './ont_fiberseq/main'

/*
===========================================
 * ONT Workflows
===========================================
 */


workflow ONT {
    take:
    ch_input

    main:

    ont_versions = Channel.empty()
    map_stat     = Channel.empty()

    // basecall

    ch_input
        .filter { it[1].toString().endsWith('.pod5') || file(it[1]).isDirectory() }
        .map { meta, pod5, _ref -> [meta, pod5] }
        .set { ch_pod5 }

    DORADO_BASECALLER(ch_pod5, params.dorado_model, params.dorado_modification)

    ont_versions = ont_versions.mix(DORADO_BASECALLER.out.versions.first())

    ch_input
        .join ( DORADO_BASECALLER.out.bam )
        .map { meta, _pod5, ref, modbam -> [meta, modbam, ref] }
        .mix ( ch_input.filter { it[1].toString().endsWith('.bam') } )
        .set { ch_input }

    // fastq and gunzip

    FASTQ_UNZIP(ch_input)

    ont_versions = ont_versions.mix(FASTQ_UNZIP.out.versions)
    map_stat = map_stat.mix(FASTQ_UNZIP.out.fastqc_log.collect { it[1] }.ifEmpty([]))

    FASTQ_UNZIP.out.unzip_input.set{ ch_ont }

    if (params.no_trim) {
        if (params.bedgraph) {

            ONT_ALIGN(ch_ont)

            ch_pile_in = ONT_ALIGN.out.ch_pile_in
            ont_versions = ont_versions.mix(ONT_ALIGN.out.versions)
            map_stat = ONT_ALIGN.out.flagstat_out

            INDEX_MODKIT_PILEUP(ch_pile_in)

            ont_versions = ont_versions.mix(INDEX_MODKIT_PILEUP.out.versions)

            BED2BEDGRAPH(INDEX_MODKIT_PILEUP.out.pileup_out)

            ont_versions = ont_versions.mix(BED2BEDGRAPH.out.versions)
        }
        else {

            ONT_ALIGN(ch_ont)

            ch_pile_in = ONT_ALIGN.out.ch_pile_in
            ont_versions = ont_versions.mix(ONT_ALIGN.out.versions)
            map_stat = ONT_ALIGN.out.flagstat_out

            INDEX_MODKIT_PILEUP(ch_pile_in)

            ont_versions = ont_versions.mix(INDEX_MODKIT_PILEUP.out.versions)
        }
    }
    else {

        if (params.dorado_demux_kit) {
            is_demux_sample_sheet = params.dorado_demux_sample_sheet ? params.dorado_demux_sample_sheet : []
            DORADO_DEMUX( ch_ont, [ [id:'kit_info'], params.dorado_demux_kit, is_demux_sample_sheet ] )

            // To correctly pass demux.out to ONT_TRIM_REPAIR
            //      -> 'flatMap' to explode list [meta, [barcode1.bam, barcode2.bam]] -> [meta, barcode1.bam] + [meta, barcode2.bam]
            //      -> Keep old 'meta.id' as 'parent_id' + inherit 'groupName' and other vals to new_meta
            // ENH: Remove 'ref' from out of dorado/demux and use a clean join with 'ch_input' ???
            DORADO_DEMUX.out.bam
                .flatMap { meta, bams, ref -> bams.collect { bam -> [ [id: bam.baseName, parent_id:meta.id, group:meta.group, method:meta.method], bam, ref] } }
                .set { in_trim_demuxed }
            COMPRESS_DEMUX_SUMMARY(DORADO_DEMUX.out.summary)
        }

        if (!params.dorado_demux_kit) {
            in_trim_demuxed = ch_ont
        }

        if (params.bedgraph) {

            ONT_TRIM_REPAIR(in_trim_demuxed)

            ont_versions = ont_versions.mix(ONT_TRIM_REPAIR.out.versions)

            ONT_ALIGN(ONT_TRIM_REPAIR.out.dorado_in)

            // Align not repaired
            ONT_TRIM_REPAIR.out.not_repaired
                .multiMap { meta, modbam, ref ->
                    bam_in: [meta, modbam]
                    ref_in: [meta, ref]
                }
                .set { not_repaired_in }
            ONT_ALIGN_NOT_REPAIRED(not_repaired_in.bam_in, not_repaired_in.ref_in, "bam_format", "bai", [], [])

            ch_pile_in = ONT_ALIGN.out.ch_pile_in
            ont_versions = ont_versions.mix(ONT_ALIGN.out.versions)
            map_stat = ONT_ALIGN.out.flagstat_out

            INDEX_MODKIT_PILEUP(ch_pile_in)

            ont_versions = ont_versions.mix(INDEX_MODKIT_PILEUP.out.versions)

            BED2BEDGRAPH(INDEX_MODKIT_PILEUP.out.pileup_out)

            ont_versions = ont_versions.mix(BED2BEDGRAPH.out.versions)
        }
        else {

            ONT_TRIM_REPAIR(in_trim_demuxed)

            ont_versions = ont_versions.mix(ONT_TRIM_REPAIR.out.versions)

            //Emitted by subwf -> declare empty
            ch_pile_in = Channel.empty()
        }
    }


    // fiberseq

    if (params.fiberseq) {

        ONT_FIBERSEQ(ch_pile_in)

        ont_versions = ont_versions.mix(ONT_FIBERSEQ.out.versions)

    }

    emit:
    ch_pile_in
    ont_versions
    map_stat
}
