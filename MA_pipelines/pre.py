import os

from MA_pipelines import config,run_cmd, valid_path,base_luigi_task,luigi
from pipelines.pangenome_pipelines import fastp


class GenerateSam_pair(base_luigi_task):
    def requires(self):
        input1 = self.infodict.get("R1", "")
        input2 = self.infodict.get("R2", "")
        sample_name = self.infodict.get("SampleID", '')
        return fastp(PE1=input1,
                     PE2=input2,
                     sample_name=sample_name)
    def output(self):
        odir = self.infodict.get("odir", '')
        sample_name = self.infodict.get("SampleID", '')
        return luigi.LocalTarget(config.output_fmt.format(
            path=odir,
            SN=sample_name) + '.sam')
        
    def run(self):
        valid_path(self.output().path, check_ofile=1)
        sample_name = self.infodict.get("SampleID", '')
        input_file1 = self.input()[0].path
        if len(self.input()) == 1:
            input_file2 = ''
        else:
            input_file2 = self.input()[1].path
        cmdline = "bwa mem -M -t {bwa_thread} -o {ofile} -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1} {i2} ".format(
            bwa_thread=config.bwa_thread,
            SN=sample_name,
            REF=config.REF_file_path,
            i1=input_file1,
            i2=input_file2,
            ofile=self.output().path)
        if self.dry_run:
            run_cmd("touch %s" % self.output().path,
                    dry_run=False)
        run_cmd(cmdline,
                dry_run=self.dry_run,
                log_file=self.get_log_path())


class Convertbam(base_luigi_task):
    
    def requires(self):
        return GenerateSam_pair(infodict=self.infodict,
                                dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.sam', '.bam'))

    def run(self):
        if config.samtools_version > 1:
            cmdline = "samtools view -G 100 -bu %s -o %s" % (self.input().path, self.output().path)
        else:
            cmdline = "samtools view -F 0x100 -bSu %s -o %s" % (self.input().path, self.output().path)
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)

        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())


class sorted_bam(base_luigi_task):
    def requires(self):
        return Convertbam(infodict=self.infodict,
                          dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.bam', '_sorted.bam'))

    def run(self):
        if config.samtools_version > 1:
            cmdline = "{samtools_pth} sort -m {sort_sam_ram} -@ {sort_sam_thread} -o {output_path} {input_path} ".format(
                samtools_pth=config.samtools_pro,
                sort_sam_ram=config.sort_sam_ram,
                sort_sam_thread=config.sort_sam_thread,
                output_path=self.output().path,
                input_path=self.input().path)
        else:
            cmdline = "{samtools_pth} sort -m {sort_sam_ram} -f -@ {sort_sam_thread} {input_path} {output_path} ".format(
                samtools_pth=config.samtools_pro,
                sort_sam_ram=config.sort_sam_ram,
                sort_sam_thread=config.sort_sam_thread,
                output_path=self.output().path,
                input_path=self.input().path)

        run_cmd(cmdline,
                dry_run=self.dry_run, log_file=self.get_log_path())
        cmdline = '{samtools_pth} index {ofile}'.format(samtools_pth=config.samtools_pro,
                                                        ofile=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)

#########2
class MarkDuplicate(base_luigi_task):
    def requires(self):
        return sorted_bam(infodict=self.infodict,
                          dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('_sorted.bam', '.dedup.bam'))

    def run(self):
        valid_path(self.output().path,
                   check_ofile=1)
        if config.PCR_ON:
            cmdline = "touch %s" % self.output().path
        else:
            input=self.input().path
            output=self.output().path
            odir=os.path.dirname(self.output().path)
            SN = self.infodict.get("SampleID", '')
            cmdline = f"gatk MarkDuplicates --java-options '{config.java_option}' --INPUT {input} --OUTPUT {output} --METRICS_FILE {odir}/{SN}_dedup_metrics.txt --CREATE_INDEX true --REMOVE_DUPLICATES true -AS true"
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


# class cal_coverage_info(base_luigi_task):
#     """
#     generate 4 files actually even thought there are only two files detected at self.output()
#     """
#     def requires(self):
#         return [sorted_bam(infodict=self.infodict, dry_run=self.dry_run)]

#     def output(self):
#         return [luigi.LocalTarget(self.input()[1].path.replace("_sorted.bam",
#                                                                ".sorted_cov.info"))
#                 ]

#     def run(self):
#         if self.dry_run:
#             for _o in self.output():
#                 run_cmd("touch %s" % _o.path, dry_run=False)

#         for _o, _i in zip(self.output(), self.input()):
#             if not self.dry_run:
#                 bam2info(bam_path=_i.path,
#                          output_cov=_o.path,
#                          bed_file=config.bed_file_path,
#                          REF_file=config.REF_file_path
#                          )
#                 summarize_covinfo(_o.path,
#                                   output_f=_o.path.replace('cov.info',
#                                                            'cov_summary.info'))
#                 # summaize the cov info with fixed format
#                 # todo: change the fixed format?? does it needed??
#             else:
#                 run_cmd("run bam2info for %s" % _i.path, dry_run=self.dry_run, log_file=self.get_log_path())


# ############################################################
# class quality_assessment(base_luigi_task):
#     """
#     calculating the coverage of given bam files.
#     Embedded into WES pipelines. Add it into main entry because it is independent.
#     """
#     tab_file = luigi.Parameter()
#     odir = luigi.Parameter()

#     def requires(self):
#         if "Normal" in self.infodict:
#             # if given a pair dict contains both normal and tumor
#             return [cal_coverage_info(infodict=self.infodict["Normal"],
#                                      dry_run=self.dry_run),
#                     cal_coverage_info(infodict=self.infodict["Tumor"],
#                                       dry_run=self.dry_run)]
#         else:
#             # just a single sample dict.
#             return cal_coverage_info(infodict=self.infodict,
#                                  dry_run=self.dry_run)

#     def output(self):
#         return luigi.LocalTarget(os.path.join(str(self.odir),
#                                               'quality_accessment_raw.csv'))

#     def run(self):
#         # summarize python script
#         # it will iterate all samples contains at `self.tab_file`
#         py_file = os.path.join(config.project_root_path,
#                                "api",
#                                "quality_accessment.py")

#         run_cmd("python3 {pyfile} -i {input} -o {output}".format(
#             pyfile=py_file,
#             input=self.tab_file,
#             output=self.odir),
#             dry_run=self.dry_run,
#             log_file=self.get_log_path())
#         if self.dry_run:
#             run_cmd("touch %s" % self.output().path, dry_run=False)
