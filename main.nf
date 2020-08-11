#!/usr/bin/env nextflow
// Nextflow version of the dragen gvcfgenotyper workflow

Channel
    .fromPath(params.gvcf_list)
    .splitCsv(strip: true)
    .map{ elem -> [elem[0], "${elem[0]}.csi"]}
    .set{ read_regions_from_s3 }
Channel
    .fromPath(params.regions_file)
    .splitCsv(strip: true, limit: params.limit)
    .set{ region }
Channel
    .fromPath(params.reference)
    .into{ first_joint_aggregation_reference; merge_consensus_sites_reference; second_joint_aggregation_reference }

region
    .combine(read_regions_from_s3)
    .groupTuple(size: params.sample_batch_size, remainder: true)
    .into { first_aggregation_extracted_regions; second_aggregation_extracted_regions }


// Parallel for every region and sample batch (region * sample_batch_size)
process first_joint_aggregation {
    cpus 1
    memory "1 GB"

    input:
        each file(reference) from first_joint_aggregation_reference
        tuple val(region), val(first_aggregation_subset), val(index) from first_aggregation_extracted_regions
    output:
        tuple val(region), file("${fixed_region}_first_aggregation.vcf.gz") into sample_consensus_sites
    
    script:   
    fixed_region = region.replaceAll(/[:-]/, "_").replaceAll(/\n/, "")
    """
    set -e

    echo "${first_aggregation_subset.join("\n")}" > ${fixed_region}_gvcf_list

    dragen --sw-mode \
     --enable-gvcf-genotyper=true \
     --enable-map-align=false \
     --gg-extra-params="--drop-genotypes" \
     --gg-enable-indexing=true \
     --output-directory . \
     --output-file-prefix ${fixed_region}_first_aggregation \
     --variant-list ${fixed_region}_gvcf_list \
     --ht-reference ${reference} \
     --gg-max-alternate-alleles $params.max_alleles \
     --gg-regions ${region} \
     --num-threads $params.threads
    """
}

sample_consensus_sites
    .groupTuple()
    .set{ grouped_sample_consensus_sites }

// Parallel for every region    
process merge_consensus_sites {
    cpus 1
    memory "1 GB"

    input:
        tuple val(region), file("*_consensus_sites_to_merge.vcf.gz") from grouped_sample_consensus_sites
        each file(reference) from merge_consensus_sites_reference

    output:
        tuple val(region), file("sites_with_consensus_samples_${fixed_region}.vcf.gz"), file("sites_with_consensus_samples_${fixed_region}.vcf.gz.csi") into consensus_sites_for_second_aggregation
    
    script:
    fixed_region = region.replaceAll(/[:-]/, "_").replaceAll(/\n/, "")
    
    """
    set -e
    ls -lh

    ls -1 *_consensus_sites_to_merge.vcf.gz > ${fixed_region}_consensus_site_file_list

    bcftools concat -a \
     -f ${fixed_region}_consensus_site_file_list \
     -Ov \
    | bcftools norm -m+any \
     -f ${reference} \
     -d none \
     -Oz \
     -o sites_with_consensus_samples_${fixed_region}.vcf.gz

    bcftools index sites_with_consensus_samples_${fixed_region}.vcf.gz
    """   
}

second_aggregation_items = second_aggregation_extracted_regions
    .combine(consensus_sites_for_second_aggregation, by: 0)

// Parallel for every region and sample batch (region * sample_batch_size)
process second_joint_aggregation {
    cpus 4
    memory "1 GB"

    input:
        tuple val(region), val(second_aggregation_subset), val(second_aggregation_subset_index), file(merged_consensus_sites), file(merged_consensus_sites_index) from second_aggregation_items
        each file(reference) from second_joint_aggregation_reference

    output:
        tuple val(region), file("${fixed_region}_second_aggregation.vcf.gz"), file("${fixed_region}_second_aggregation.vcf.gz.tbi") into merge_samples

    script:
    fixed_region = region.replaceAll(/[:-]/, "_").replaceAll(/\n/, "")
    """
    set -e

    echo "${second_aggregation_subset.join("\n")}" > ${fixed_region}_gvcf_list

    dragen --sw-mode \
     --enable-gvcf-genotyper=true \
     --enable-map-align=false \
     --gg-enable-indexing=true \
     --variant-list ${fixed_region}_gvcf_list \
     --gg-extra-params="--allele-list ${merged_consensus_sites}" \
     --gg-regions ${region} \
     --ht-reference ${reference} \
     --gg-max-alternate-alleles $params.max_alleles \
     --num-threads $params.threads \
     --output-directory . \
     --output-file-prefix ${fixed_region}_second_aggregation \
     --gg-enable-concat=false
    """
}

// Parallel for every region  
process merge_chunks {
    cpus 1
    memory "1 GB"

    publishDir "./results", mode: "copy"

    input:
        tuple val(region), file("*_second_aggregation_subset.vcf.gz"), file("*_second_aggregation_subset.vcf.gz.tbi") from merge_samples.groupTuple()
    output:
        tuple file("all_chunks_merged_${fixed_region}.vcf.gz"), file("all_chunks_merged_${fixed_region}.vcf.gz.tbi")
    
    script:
    fixed_region = region.replaceAll(/[:-]/, "_").replaceAll(/\n/, "") 
    """
    set -e

    ls -1 *_second_aggregation_subset.vcf.gz > "${fixed_region}_second_joint_aggregation_file_list"

    bcftools merge -m all \
     -l ${fixed_region}_second_joint_aggregation_file_list \
     -Oz \
     -o all_chunks_merged_${fixed_region}.vcf.gz

    bcftools index -t all_chunks_merged_${fixed_region}.vcf.gz
    """
}
