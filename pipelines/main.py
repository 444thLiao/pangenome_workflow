import sys
from os.path import dirname

sys.path.insert(0, dirname(dirname(__file__)))
from pipelines import *


class workflow(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)
    log_path = luigi.Parameter(default=None)

    def output(self):
        return luigi.LocalTarget(os.path.join(str(self.odir),
                                              "pipelines_summary"))

    def requires(self):
        inputdata = fileparser(self.tab)
        # it will validate the input df
        require_tasks = {}
        pairreads = inputdata.get_PE_info()
        singlereads = inputdata.get_SE_info()
        other_info = inputdata.get_full_PE_info()
        ############################################################
        valid_path(self.odir, check_odir=True)
        unify_kwargs = dict(odir=self.odir,
                            dry_run=self.dry_run,
                            log_path=self.log_path,
                            PE_data=pairreads, )

        require_tasks["fastqc_before"] = multiqc(status='before',
                                                 **unify_kwargs
                                                 )
        require_tasks["fastqc_after"] = multiqc(status='after',
                                                **unify_kwargs
                                                )
        require_tasks["fastqc_quast"] = multiqc(status='quast',
                                                other_info=json.dumps(other_info),
                                                **unify_kwargs
                                                )
        require_tasks["abricate"] = abricate(SE_data=singlereads,
                                             **unify_kwargs)
        require_tasks["fasttree"] = fasttree(SE_data=singlereads,
                                             **unify_kwargs)
        require_tasks["pandoo"] = pandoo(SE_data=singlereads,
                                         **unify_kwargs)
        require_tasks["ISEscan_summary"] = ISEscan_summary(SE_data=singlereads,
                                                           **unify_kwargs)

        require_tasks["detect_plasmid"] = detect_plasmid(**unify_kwargs)
        require_tasks["detect_prophage"] = phigaro_summary(SE_data=singlereads,
                                                           **unify_kwargs)
        return require_tasks

    def run(self):
        # post pipelines
        post_analysis(self)


if __name__ == '__main__':
    luigi.run()

    # python -m luigi --module pipelines.luigi_pipelines workflow --tab /home/liaoth/project/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/project/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # local cmd

    # python -m luigi --module pipelines.luigi_pipelines workflow --tab /home/liaoth/tools/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/tools/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # server cmd
