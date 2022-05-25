from MA_pipelines import config, base_luigi_task,luigi,valid_path,run_cmd
# from MA_pipelines import HaploCaller

from pre import MarkDuplicate
from HaploCaller import *


class BaseRecalibrator(base_luigi_task):
    def requires(self):
        return [MarkDuplicate(infodict=self.infodict,
                              dry_run=self.dry_run),
                CombineVariants(infodict=self.infodict,
                              dry_run=self.dry_run)
                ]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam',
                                                              '.recal.table'))
    def run(self):
        output_f=self.output().path
        input_f=self.input()[0].path
        mergedvcf=self.input()[1].path
        REF = self.infodict.get('REF','')
        
        cmdline = f"{config.gatk_pro} BaseRecalibrator --java-options {config.java_option}  --input {input_f} --output {output_f} --reference {REF} --known-sites {mergedvcf} "
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)

class ApplyBQSR(base_luigi_task):
    def requires(self):
        return [BaseRecalibrator(infodict=self.infodict,
                                 dry_run=self.dry_run),
                MarkDuplicate(infodict=self.infodict,
                              dry_run=self.dry_run),
                ]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace(".recal.table",
                                                              ".recal.bam"))

    def run(self):
        output_f=self.output().path
        REF = self.infodict.get('REF','')
        recal_base=self.input()[0].path
        bam_input=self.input()[1].path
        output_f=self.output().path
        cmdline = f"{config.gatk_pro} ApplyBQSR -R {REF} -I {bam_input} --bqsr-recal-file {recal_base} -O {output_f}"
        
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        
        cmdline = '{samtools_pth} index {ofile}'.format(samtools_pth=config.samtools_pro,
                                                        ofile=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)

## 2nd call SNP
class HC_2nd(HaplotypeCaller):
    def requires(self):
        return ApplyBQSR(infodict=self.infodict, dry_run=self.dry_run)
    def output(self):        
        ofile = self.input().path.replace('.recal.bam',f'.raw_variants.rnd2.vcf')
        return luigi.LocalTarget(ofile)
    
class SV_2nd(SelectVariants):
    def requires(self):
        return HC_2nd(infodict=self.infodict,
                      dry_run=self.dry_run)

class VF_2nd(VariantFiltration):
    def requires(self):
        return SV_2nd(infodict=self.infodict,
                              dry_run=self.dry_run,
                              object_type=self.object_type)
class CV_2nd(CombineVariants):
    def requires(self):
        required_task = {ot: VF_2nd(infodict=self.infodict,
                                               dry_run=self.dry_run,
                                               object_type=ot)
                         for ot in ["snp", "indel"]}
        return required_task
    
class BR_2nd(BaseRecalibrator):
    def requires(self):
        return [ApplyBQSR(infodict=self.infodict,
                              dry_run=self.dry_run),
                CV_2nd(infodict=self.infodict,
                              dry_run=self.dry_run)]
    def output(self):
        ofile = self.input()[0].path.replace('.recal.bam','.recal.rnd2.table')
        return luigi.LocalTarget(ofile)

class AB_2nd(ApplyBQSR):
    def requires(self):
        return [BR_2nd(infodict=self.infodict,
                                 dry_run=self.dry_run),
                ApplyBQSR(infodict=self.infodict,
                              dry_run=self.dry_run)]
    def output(self):
        ofile = self.input()[0].path.replace(".table",".bam")
        # suffix is .recal.rnd2.bam
        return luigi.LocalTarget(ofile)
    
## 3nd call SNP
class HC_3nd(HaplotypeCaller):
    def requires(self):
        return AB_2nd(infodict=self.infodict, dry_run=self.dry_run)
    def output(self):
        ofile = self.input().path.replace('.recal.rnd2.bam',f'.raw_variants.rnd3.vcf')
        return luigi.LocalTarget(ofile)
    
class SV_3nd(SelectVariants):
    def requires(self):
        return HC_3nd(infodict=self.infodict,
                      dry_run=self.dry_run)

class VF_3nd(VariantFiltration):
    def requires(self):
        return SV_3nd(infodict=self.infodict,
                              dry_run=self.dry_run,
                              object_type=self.object_type)

class CV_3nd(CombineVariants):
    def requires(self):
        required_task = {ot: VF_3nd(infodict=self.infodict,
                                    dry_run=self.dry_run,
                                    object_type=ot)
                         for ot in ["snp", "indel"]}
        return required_task


class BR_3nd(BaseRecalibrator):
    def requires(self):
        return [AB_2nd(infodict=self.infodict,
                              dry_run=self.dry_run),
                CV_3nd(infodict=self.infodict,
                              dry_run=self.dry_run)]
    def output(self):
        ofile = self.input()[0].path.replace('.recal.rnd2.bam','.recal.rnd3.table')
        # suffix is .recal.rnd3.table
        return luigi.LocalTarget(ofile)

class AB_3nd(ApplyBQSR):
    def requires(self):
        return [BR_3nd(infodict=self.infodict,
                                 dry_run=self.dry_run),
                AB_2nd(infodict=self.infodict,
                       dry_run=self.dry_run)]
    def output(self):
        ofile = self.input()[0].path.replace(".recal.rnd3.table",".recal.rnd3.bam")
        # suffix is .recal.rnd3.bam
        return luigi.LocalTarget(ofile)
    
## 4nd call SNP
class HC_4nd(HaplotypeCaller):
    def requires(self):
        return AB_3nd(infodict=self.infodict, dry_run=self.dry_run)
    def output(self):
        ofile = self.input().path.replace('.recal.rnd3.bam',f'.raw_variants.rnd4.vcf')
        return luigi.LocalTarget(ofile)
    
class SV_4nd(SelectVariants):
    def requires(self):
        return HC_4nd(infodict=self.infodict,
                      dry_run=self.dry_run)

class VF_4nd(VariantFiltration):
    def requires(self):
        return SV_4nd(infodict=self.infodict,
                              dry_run=self.dry_run,
                              object_type=self.object_type)

class CV_4nd(CombineVariants):
    def requires(self):
        required_task = {ot: VF_4nd(infodict=self.infodict,
                                    dry_run=self.dry_run,
                                    object_type=ot)
                         for ot in ["snp", "indel"]}
        return required_task

