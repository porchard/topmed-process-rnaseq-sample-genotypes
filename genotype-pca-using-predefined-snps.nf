#!/usr/bin/env nextflow

nextflow.enable.dsl=2

VCF_GLOB = params.vcf_glob
VCF_INDEX_GLOB = params.vcf_index_glob
KEEP_VARIANTS = params.variants
KEEP_SAMPLES = params.samples
ANCESTRIES = params.ancestries

chroms = (1..22).collect({it -> "chr" + it})

process filter_to_variants {

    cache 'lenient'
    tag "${chrom}"
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(chrom), path("in.vcf.gz"), path("in.vcf.gz.tbi")

    output:
    tuple val(chrom), path("${chrom}.vcf.gz"), path("${chrom}.vcf.gz.tbi")

    when:
    chroms.contains(chrom)

    """
    grep -w $chrom $KEEP_VARIANTS | perl -pe 's/\\t/_/g' > keep-variant-ids.1.txt
    # also swap ref and alt as plink doesn't distinguish so we don't know which is correct
    cat keep-variant-ids.1.txt | perl -pe 's/(.*)_(.*)_(.*)_(.*)/\$1_\$2_\$4_\$3/' > keep-variant-ids.2.txt
    cat keep-variant-ids.1.txt keep-variant-ids.2.txt > keep-variant-ids.txt
    rm keep-variant-ids.1.txt keep-variant-ids.2.txt
    bcftools view --samples-file $KEEP_SAMPLES -Ou in.vcf.gz | bcftools filter -i 'ID=@keep-variant-ids.txt' -Oz -o ${chrom}.vcf.gz
    tabix ${chrom}.vcf.gz
    """

}

process merge_vcf {
    
    cache 'lenient'
    publishDir "${params.results}/merged", mode: 'symlink'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'
    container 'library://porchard/default/general:20220107'

    input:
    path(vcfs)
    path(indices)

    output:
    path("merged.sorted.vcf.gz")

    """
    bcftools concat --allow-overlaps ${vcfs.join(' ')} -Ou | bcftools sort -T . -o merged.sorted.vcf.gz -Oz
    """

}

process make_ped_map_from_vcf {

    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'
    container 'library://porchard/default/general:20220107'
    cache 'lenient'
    tag "${chrom}"


    input:
    path(vcf)

    output:
    tuple val('merged'), path("merged.ped"), path("merged.map")

    """
    plink --noweb --keep-allele-order --vcf $vcf --recode --out merged_tmp --make-bed
    cat merged_tmp.fam | awk '\$6=1' > merged_tmp.fam.tmp
    mv merged_tmp.fam.tmp merged_tmp.fam
    plink --noweb --keep-allele-order --bfile merged_tmp --recode --out merged
    """

}

process pca {

    cache 'lenient'
    publishDir "${params.results}/pca", mode: 'symlink'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(plink_prefix_in), path(ped), path(map)

    output:
    path("filtered*")
    tuple path('filtered.pca.evec'), path('filtered.eval'), emit: plot

    """
    cat ${plink_prefix_in}.map | cut -f2 | awk 'length(\$1)>25' > snps-with-long-ids.txt
    plink --noweb --keep-allele-order --file $plink_prefix_in --exclude snps-with-long-ids.txt --recode --out filtered
    smartpca.perl -i filtered.ped -a filtered.map -b filtered.ped -k 15 -m 0 -o filtered.pca -e filtered.eval -p filtered.plot -l filtered.log
    evec2pca-ped.perl 15 filtered.pca.evec filtered.ped filtered.ped.pca
    """

}


process process_pca {

    publishDir "${params.results}/pca", mode: 'symlink'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'
    cache 'lenient'
    container 'docker.io/porchard/general:20220406125608'
    
    input:
    tuple path(pca), path(evalues), path(ancestries) 

    output:
    path("*.png")
    path("*.txt")
    path("PC-scores.txt"), emit: pc_scores

    """
    plot-eigensoft-pca.py $pca $evalues $ancestries
    """

}


workflow {
    ancestries = Channel.fromPath(ANCESTRIES)
    bcf = Channel.fromPath(VCF_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple() // chrom, bcf
    bcf_index = Channel.fromPath(VCF_INDEX_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple() // chrom, bcf_index
    filtered_bcf = filter_to_variants(bcf.combine(bcf_index, by: 0))
    merged = merge_vcf(filtered_bcf.map({it -> it[1]}).toSortedList(), filtered_bcf.map({it -> it[2]}).toSortedList())
    processed = (make_ped_map_from_vcf(merged) | pca).plot.combine(ancestries) | process_pca
}
