############################################################
# exe path
############################################################
fastqc_path = "/tools/anaconda3/envs/pangenome_pipelines/bin/fastqc"

multiqc_path = "/tools/anaconda3/envs/pangenome_pipelines/bin/multiqc"

trimmomatic_dir = "/home/liaoth/tools/Trimmomatic-0.36"
trimmomatic_setting = "ILLUMINACLIP:%s/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:5" % trimmomatic_dir
trimmomatic_jar_path = "%s/trimmomatic-0.36.jar" % trimmomatic_dir

shovill_path = "/tools/anaconda3/envs/pangenome_pipelines/bin/shovill"
prokka_path = "/tools/prokka/bin//prokka"
roary_path = "/tools/anaconda3/envs/pangenome_pipelines/bin/roary"
quast_path = "/tools/anaconda3/envs/pangenome_pipelines/bin/quast.py"
pandoo_path = "/tools/anaconda3/envs/pangenome_pipelines/bin/pandoo"
fasttree_path = "/usr/bin/fasttreeMP"

ISEscan_path = "/tools/ISEScan/isescan.py"
abricate_path = "/tools/abricate/bin/abricate"
abricate_py_path = "/home/liaoth/tools/genome_pipelines/toolkit/batch_abricate.py"
plasmid_detect_path = "/home/liaoth/tools/genome_pipelines/toolkit/process_plasmid.py"

phigaro_path = "/tools/anaconda3/envs/pangenome_pipelines/bin/phigaro"
gubbins_path = "/tools/anaconda3/envs/pangenome_pipelines/bin/run_gubbins.py"

mlst_path = '/tools/mlst/bin/mlst'
kraken2_path='/tools/kraken2/kraken2'
seqtk_path ='/usr/bin/seqtk'

############################################################
# database path
############################################################
abricate_db = "/tools/abricate/db/"

ariba_db = "/home/db_public/ariba_db"
ariba_db_list = ["CARD", "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_full", "virulencefinder"]
ariba_str = str({k: "{db}/{name}".format(db=ariba_db, name=k.lower()) for k in ariba_db_list})

mlst_db = '/tools/mlst/db/pubmlst/'
phigaro_config = "/home/db_public/phigaro_db/config.yml"
kraken2_db = "/home/db_public/kraken2_db/"
