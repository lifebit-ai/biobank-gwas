#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/biobank-gwas
========================================================================================
 lifebit-ai/biobank-gwas GWAS pipeline built for Genomics England using SAIGE
 #### Homepage / Documentation
 https://github.com/lifebit-ai/biobank-gwas
----------------------------------------------------------------------------------------
*/

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/

ch_pheno = params.pheno_data ? Channel.value(file(params.pheno_data)) : Channel.empty()
(phenoCh_gwas_filtering, phenoCh) = ch_pheno.into(2)

Channel
  .fromFilePairs("${params.grm_plink_input}",size:3, flat : true)
  .ifEmpty { exit 1, "PLINK files not found: ${params.grm_plink_input}.\nPlease specify a valid --grm_plink_input value. eg. testdata/*.{bed,bim,fam}" }
  .set { plinkCh }
Channel
  .fromPath(params.plink_keep_pheno)
  .set {plink_keep_pheno_ch}
Channel
  .fromPath(params.vcfs_list)
  .ifEmpty { exit 1, "Cannot find CSV VCFs file : ${params.vcfs_list}" }
  .splitCsv(skip:1)
  .map { chr, vcf, index -> [file(vcf).simpleName, chr, file(vcf), file(index)] }
  .set { vcfsCh }
Channel
  .fromPath(params.gwas_cat)
  .ifEmpty { exit 1, "Cannot find GWAS catalogue CSV  file : ${params.gwas_cat}" }
  .set { ch_gwas_cat }



  /*--------------------------------------------------
  Pre-GWAS filtering - download, filter and convert VCFs
  ---------------------------------------------------*/
if (params.trait_type == 'binary'){
  process gwas_filtering_bin {
    tag "$name"
    publishDir "${params.outdir}/gwas_filtering", mode: 'copy'

    input:
    set val(name), val(chr), file(vcf), file(index) from vcfsCh
    each file(phe_file) from phenoCh_gwas_filtering
    each file(plink_keep_file) from plink_keep_pheno_ch

    output:
    set val(name), val(chr), file("${name}.filtered_final.vcf.gz"), file("${name}.filtered_final.vcf.gz.csi") into filteredVcfsCh
    
    script:
    // TODO: (High priority) Only extract needed individuals from VCF files with `bcftools -S samples.txt` - get from samples file?
    // TODO: (Not required) `bcftools -T sites_to_extract.txt`
    // Optional parameters
    extra_plink_filter_missingness_options = params.plink_keep_pheno != "s3://lifebit-featured-datasets/pipelines/biobank-gwas/testdata/nodata" ? "--keep ${plink_keep_file}" : ""
    """
    # Download, filter and convert (bcf or vcf.gz) -> vcf.gz
    bcftools view -q ${params.q_filter} $vcf -Oz -o ${name}_filtered.vcf.gz
    bcftools index ${name}_filtered.vcf.gz
  
    # Create PLINK binary from vcf.gz
    plink2 \
      --make-bed \
      --set-missing-var-ids @:#,\\\$r,\\\$a \
      --vcf ${name}_filtered.vcf.gz \
      --out ${name}_filtered \
      --vcf-half-call m \
      --double-id \
      --set-hh-missing \
      --new-id-max-allele-len 60 missing

    #Filter missingness
    plink \
      --bfile ${name}_filtered \
      --pheno $phe_file \
      --pheno-name PHE \
      --allow-no-sex \
      --test-missing midp \
      --out ${name} \
      --1 \
      --keep-allele-order \
      ${extra_plink_filter_missingness_options}

    awk '\$5 < ${params.thres_m} {print}' ${name}.missing > ${name}.missing_FAIL

    #Filter HWE
    plink \
      --bfile ${name}_filtered \
      --pheno $phe_file \
      --pheno-name PHE \
      --allow-no-sex \
      --hwe ${params.thres_HWE} midp \
      --out ${name}.misHWEfiltered \
      --make-just-bim \
      --exclude ${name}.missing_FAIL \
      --1 \
      --keep-allele-order \
      ${extra_plink_filter_missingness_options}

    bcftools view ${name}_filtered.vcf.gz | awk -F '\\t' 'NR==FNR{c[\$1\$4\$6\$5]++;next}; c[\$1\$2\$4\$5] > 0' ${name}.misHWEfiltered.bim - | bgzip > ${name}.filtered_temp.vcf.gz
    bcftools view -h ${name}_filtered.vcf.gz -Oz -o ${name}_filtered.header.vcf.gz
    cat ${name}_filtered.header.vcf.gz ${name}.filtered_temp.vcf.gz > ${name}.filtered_final.vcf.gz
    bcftools index ${name}.filtered_final.vcf.gz
    """
  }
}


