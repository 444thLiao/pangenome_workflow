from MA_pipelines import config, base_luigi_task,luigi,valid_path,run_cmd

from pre import MarkDuplicate
#########7

class HaplotypeCaller(base_luigi_task):
    def requires(self):
        return MarkDuplicate(infodict=self.infodict, dry_run=self.dry_run)
    
    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.bam',
                                                           '.raw_variants.vcf'))
        
    def run(self):
        valid_path(self.output().path, check_ofile=1)
        
        ref=self.infodict.get('REF','')
        input=self.input().path
        output=self.output().path
        extra_str=extra_str
        gatk4=config.gatk_pro
        
        cmdline = f"{gatk4} HaplotypeCaller --java-options '-Xmx4g' --native-pair-hmm-threads 10 --reference {ref} --input {input} --sample-ploidy 1 -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {output}"
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())

#########9
class SelectVariants(base_luigi_task):
    def requires(self):
        return HaplotypeCaller(infodict=self.infodict,
                               dry_run=self.dry_run)

    def output(self):
        if self.object_type == "snp":
            ofile_name = '.raw_snps.vcf'
        elif self.object_type == "indel":
            ofile_name = '.raw_indels.vcf'
        else:
            raise Exception

        return luigi.LocalTarget(self.input().path.replace('.raw_variants.vcf',
                                                           ofile_name))
    def run(self):
        valid_path(self.output().path, check_ofile=1)
        if self.object_type == "snp":
            selecttype = "SNP"
        elif self.object_type == "indel":
            selecttype = "INDEL"
        else:
            raise Exception
        gatk4=config.gatk_pro
        REF=self.infodict.get('REF','')
        input_f=self.input().path
        output_f=self.output().path
        selecttype=selecttype
        cmdline = f"{gatk4} SelectVariants -R {REF} -V {input_f} -select-type {selecttype} -O {output_f}"
        
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())


#########10
class VariantFiltration(base_luigi_task):
    def requires(self):
        return SelectVariants(infodict=self.infodict,
                              dry_run=self.dry_run,
                              object_type=self.object_type)
    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.raw_',
                                                           '.filter_'))
    def output(self):
        odir = self.infodict.get("odir", '')
        project_name = self.infodict.get("project_name", '')
        sample_name = self.infodict.get("SampleID", '')
        return luigi.LocalTarget(config.output_fmt.format(
            path=odir,
            PN=project_name,
            SN=sample_name) + '.sam')
    def run(self):
        valid_path(self.output().path, check_ofile=1)
        SNP_QUAL=100.0
        SNP_MQ=40.0
        SNP_SOR=3.0
        INDEL_QUAL=100.0
        INDEL_MQ=40.0
        INDEL_SOR=10.0
        if self.object_type == "snp":
            filterExpression = f"QUAL < {SNP_QUAL} || MQ < {SNP_MQ} || SOR > {SNP_SOR}"
            # filterExpression = "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        elif self.object_type == "indel":
            filterExpression = f"QUAL < {INDEL_QUAL} || MQ < {INDEL_MQ} || SOR > {INDEL_SOR}"
            # filterExpression = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
        else:
            raise Exception
        
        gatk4=config.gatk_pro
        REF=self.infodict.get('REF','')
        input_f=self.input().path
        output_f=self.output().path
        object_type=self.object_type
        
        cmdline = f"""{gatk4} VariantFiltration -R {REF} -V {input_f} --filter-expression "{filterExpression}" --filter-name \"luolab_{object_type}_filter\" -O {output_f}"""
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())


#########13
class CombineVariants(base_luigi_task):
    
    def requires(self):
        required_task = {ot: VariantFiltration(infodict=self.infodict,
                                               dry_run=self.dry_run,
                                               object_type=ot)
                         for ot in ["snp", "indel"]}
        return required_task
    
    def output(self):
        return luigi.LocalTarget(self.input()["snp"].path.replace('.filter_snps.vcf',
                                                                    '.merged.vcf'))
    def run(self):
        valid_path(self.output().path, check_ofile=1)
        
        input_indel=self.input()["indel"].path
        input_snp=self.input()["snp"].path
        output_f=self.output().path
        
        cmdline = f"""{config.gatk_pro} MergeVcfs --java-options "-Xmx4g" -R {self.infodict.get('REF','')} --INPUT {input_indel} --INPUT {input_snp} --OUTPUT {output_f}"""
        
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())


if __name__ == '__main__':
    luigi.run()

    #
