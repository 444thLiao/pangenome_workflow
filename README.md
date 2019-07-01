# Pangenome analysis workflow

For collecting all required tasks into a single pipelines, this project had been started.

## installation

This is a project composed with a lot of other tools or software. Some software also contain its own database. 

There is no way could easily handle all these software at one-click. 
Here we provide a `environment.yml` for easy create a environment with **ananconda**.

Here is a list of necessary software

* trimmomatic
* fastqc
* multiqc
* shovill
* prokka
* fasttree
* ISEscan
* abricate
* Spades
* mlst
* kraken2
* seqtk

## config
After installing all these stuff, you must fulfill a `config` file which located at `pipelines/soft_db_path.py`.

## testing
```bash
python3 pipelines/main.py testdata -o output_dir
```

## QuickStart

After all above things, you could finally start using this pipelines for your own data.

Following the header and separator of `toolkit/data_input.template`, fulfill a new `data_input.tab`.

With this tab, you could run:
```bash
python3 pipelines/main.py run -- workflow --tab data_input.tab --odir output_dir --workers 2 --log-path output_dir/cmd_log.txt
```


Besides the params `--tab`, `--odir`, `--analysis-type`, `--log-path`, other params are luigi implemented. 

Here describe a little bit about these params. For more detailed, you should check the documentation of luigi at [luigi doc](https://luigi.readthedocs.io/en/stable/)

* `--tab`: given a path(could be relative/absolute) of `input_data.tab`
* `--odir`: jus the path of output directory. a little be need to say is that, different pipelines like `otu, deblur, dada2`, it will separately located the final output of different pipelines. So **don't worry using same output dir will confuse the result**.
* `--analysis-type`: for now, three options including *otu, deblur, dada2* could be selected, if you want to perform all at once. You could pass `all` param to it. Because there are a lot of overlapped tasks among three different pipelines, it would save a lot of time than running these separately with different `odir`. 
* `--log-path`: it just record the cmd history.*(optional)*
* `--workers`: it could control how many tasks could be parallel.



## about the `data_input.tab`

If you look at the `toolkit/data_input.template`, there are only three header. 

`Tab` is taken as separator of this `data_input.tab` for better handle some weird filename.

Inside this `data_input.tab`, you could append more columns besides the necessary `five columns(sample ID	R1	R2  ref gff)`. This pipelines only check you have these three instead of only have these five.

Besides that, if you don't know which `ref` or which `gff` you need to choose, you could just left blank and it will bypass this step.

It just for quality accessment for now.

## Feedback

Please file questions, bugs or ideas 
to the [Issue Tracker](https://github.com/444thLiao/pangenome_workflow)

## Authors
Tianhua Liao
email: l0404th@gmail.com

 