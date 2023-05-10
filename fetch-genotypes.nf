#!/usr/bin/env nextflow

nextflow.enable.dsl=2

SAMPLES = params.samples
BCF_GLOB = params.bcf_glob // name: {chrom}.bcf
chroms = (1..22).collect({it -> "chr" + it}) + ['chrX']

process subset_bcf {

    container 'library://porchard/default/general:20220107'
    memory '10 GB'
    time '240h'
    publishDir "${params.results}/bcfs-by-chrom"
    cache 'lenient'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[3-10]'

    input:
    tuple val(chrom), path("in.bcf"), path(s)

    output:
    tuple val(chrom), path("${chrom}.bcf")

    """
    bcftools view --types snps,indels --force-samples --samples-file $s -Ou -i 'TYPE="snp" || ILEN<50' in.bcf | bcftools view --min-ac 1:minor -Ou -o ${chrom}.bcf
    """

}

process update_ids {

    container 'library://porchard/default/general:20220107'
    memory '20 GB'
    time '168h'
    publishDir "${params.results}/vcfs-updated-ids"
    cache 'lenient'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[3-10]'

    input:
    tuple val(chrom), path("in.bcf")

    output:
    tuple val(chrom), path("${chrom}.vcf.gz")

    """
    bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -Oz -o ${chrom}.vcf.gz in.bcf
    """

}

process pass_filter {

    container 'library://porchard/default/general:20220107'
    memory '20 GB'
    time '168h'
    publishDir "${params.results}/vcfs-updated-ids-pass-filter"
    cache 'lenient'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[3-10]'

    input:
    tuple val(chrom), path("in.vcf.gz")

    output:
    tuple val(chrom), path("${chrom}.vcf.gz")

    """
    bcftools view -f 'PASS,.' -Oz -o ${chrom}.vcf.gz in.vcf.gz
    """

}

process index_pass_filter {

    container 'library://porchard/default/general:20220107'
    time '168h'
    publishDir "${params.results}/vcfs-updated-ids-pass-filter"
    cache 'lenient'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[3-10]'

    input:
    tuple val(chrom), path(vcf)

    output:
    path("*.tbi")

    """
    tabix $vcf
    """

}

process make_plink {

    publishDir "${params.results}/plink"
    memory '20 GB'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[3-10]'
    container 'docker.io/porchard/plink2:20221023'

    input:
    tuple val(chrom), path(vcf)

    output:
    tuple path("${chrom}.bed"), path("${chrom}.bim"), path("${chrom}.fam")

    """
    plink2 --vcf $vcf --make-bed --out ${chrom}
    """

}

process make_plink_pass_filter {

    publishDir "${params.results}/plink-pass-filter"
    memory '20 GB'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[3-10]'
    container 'docker.io/porchard/plink2:20221023'

    input:
    tuple val(chrom), path(vcf)

    output:
    tuple path("${chrom}.bed"), path("${chrom}.bim"), path("${chrom}.fam")

    """
    plink2 --vcf $vcf --make-bed --out ${chrom}
    """

}

workflow {
    bcf_in = Channel.fromPath(BCF_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // chrom, bcf
    samples = Channel.fromPath(SAMPLES)

    updated_ids = bcf_in.combine(samples) | subset_bcf | update_ids
    make_plink(updated_ids)
    filtered = pass_filter(updated_ids)
    make_plink_pass_filter(filtered)
    index_pass_filter(filtered)
}
