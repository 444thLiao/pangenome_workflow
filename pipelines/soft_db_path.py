############################################################
# exe path
############################################################
fastqc_path = "/home/liaoth/tools/FastQC/fastqc"

multiqc_path = "/usr/local/bin/multiqc"

trimmomatic_dir = "/home/liaoth/tools/Trimmomatic-0.36"
trimmomatic_setting = "ILLUMINACLIP:%s/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:5" % trimmomatic_dir
trimmomatic_jar_path = "%s/trimmomatic-0.36.jar" % trimmomatic_dir

shovill_path = "/home/liaoth/tools/anaconda3/envs/pangenome_pipelines/bin/shovill"
prokka_path = "/home/liaoth/tools/prokka/bin//prokka"
roary_path = "/usr/bin/roary"
roary_plot_path = "/tools/Roary/contrib/roary_plots/roary_plots.py"
quast_path = "/usr/local/bin/quast.py"
pandoo_path = "/home/liaoth/tools/anaconda3/envs/pangenome_pipelines/bin/pandoo"

fasttree_path = "/usr/bin/fasttreeMP"

ISEscan_path = "/home/liaoth/tools/ISEScan/isescan.py"
abricate_path = "/home/liaoth/tools/anaconda3/envs/pandoo/bin/abricate"
abricate_py_path = "/home/liaoth/project/genome_pipelines/toolkit/batch_abricate.py"
plasmid_detect_path = "/home/liaoth/project/genome_pipelines/toolkit/process_plasmid.py"

phigaro_path = "/usr/local/bin/phigaro"
gubbins_path = "/tools/anaconda3/envs/pangenome_pipelines/bin/run_gubbins.py"
mlst_path = '/home/liaoth/tools/mlst/bin/mlst'
kraken2_path='/usr/bin/kraken2'
seqtk_path ='/home/liaoth/tools/anaconda3/envs/pangenome_pipelines/bin/seqtk'

############################################################
# database path
############################################################
abricate_db = "/tools/abricate/db/"

ariba_db = "/home/db_public/ariba_db"
ariba_db_list = ["CARD", "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_full", "virulencefinder"]
ariba_str = str({k: "{db}/{name}".format(db=ariba_db, name=k.lower()) for k in ariba_db_list})

mlst_db = '/home/liaoth/tools/mlst/db/pubmlst/'
phigaro_config = "/home/liaoth/.phigaro/config.yml"
kraken2_db = "/home/liaoth/data2/kraken2_db/minikraken2_v2_8GB_201904_UPDATE"