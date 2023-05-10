#!/usr/bin/env nextflow

nextflow.enable.dsl=2

VCF_GLOB = params.vcf_glob
CHROMS = (1..22)
chroms_ordered = CHROMS.collect({it -> "chr" + it}).sort()

// https://www.kingrelatedness.com/manual.shtml
// Please do not prune or filter any "good" SNPs that pass QC prior to any KING inference, unless the number of variants is too many to fit the computer memory, e.g., > 100,000,000 as in a WGS study, in which case rare variants can be filtered out. LD pruning is not recommended in KING.

process maf_filter {

    cache 'lenient'
    container 'library://porchard/default/general:20220107'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'

    input:
    tuple val(chrom), path(vcf)

    output:
    path("${chrom}.maf001.bcf")

    """
    bcftools view -f 'PASS,.' --types snps --min-af 0.01:minor -o ${chrom}.maf001.bcf -Ou $vcf
    """

}


process merge_bcfs {

    memory '40 GB'
    time '168h'
    publishDir "${params.results}/bcfs-merged"
    cache 'lenient'
    container 'library://porchard/default/general:20220107'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'

    input:
    path(vcf)

    output:
    path('merged.bcf')

    script:
    tmp = chroms_ordered.collect({x -> x + ".maf001.bcf"}).join(' ')

    """
    bcftools concat -Ou -o merged.bcf $tmp
    """

}


process make_plink {
	
    publishDir "${params.results}/plink"
    container 'library://porchard/default/general:20220107'
    memory '20 GB'
    cache 'lenient'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'

    input:
    path(bcf)

    output:
    tuple path("plink.merged.bed"), path("plink.merged.bim"), path("plink.merged.fam")

    """
    plink --allow-extra-chr 0 --bcf $bcf --make-bed --out plink.merged
    """

}


process get_related {

    publishDir "${params.results}/king/related"
    container 'docker.io/porchard/king:2.2.7'
    memory '100 GB'
    cache 'lenient'
    time '24h'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'
    cpus 20

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path(bed), path(bim), path(fam), path("king.kin0")

    """
    king --cpus 19 -b $bed --degree 5 --related --prefix king
    """

}


process get_unrelated {

    publishDir "${params.results}/king/unrelated"
    container 'library://porchard/default/general:20220107'
    memory '20 GB'
    cache 'lenient'
    time '24h'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'

    input:
    tuple path(bed), path(bim), path(fam), path(king)

    output:
    path('unrelated.txt')

    """
    relatedness-to-unrelated.py $king $fam > unrelated.txt
    """

}


workflow {
    vcf_in = Channel.fromPath(VCF_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}).filter({it -> it[0] != 'chrX'})
    plinks = maf_filter(vcf_in).toSortedList() | merge_bcfs | make_plink | get_related | get_unrelated
}
