#################################################################################
#### File for storage all static path of exe or command line(unformatted)
####
#### If you want to use the pipelines at your computer/workstation/server, Please first check the exe pth first. Or run `#todo: `blabla test` to check whether any exe could find
############################################################
# exe path
############################################################
fastqc_path = "/home/liaoth/tools/FastQC//fastqc"

multiqc_path = "/usr/local/bin/multiqc"

trimmomatic_dir = "/home/liaoth/tools/Trimmomatic-0.36"
trimmomatic_setting = "ILLUMINACLIP:%s/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:5" % trimmomatic_dir
trimmomatic_jar_path = "%s/trimmomatic-0.36.jar" % trimmomatic_dir

shovill_path = "/home/liaoth/tools/anaconda3/envs/pangenome_pipelines/bin/shovill"
prokka_path = "/home/liaoth/tools/prokka/bin//prokka"
roary_path = "/usr/bin/roary"
quast_path = "/usr/local/bin/quast.py"
pandoo_path = "/usr/local/bin/pandoo"
fasttree_path = "/usr/bin/fasttreeMP"

ISEscan_path = "/home/liaoth/tools/ISEScan/util/batch4hmp.py"
abricate_path = "/home/liaoth/tools/anaconda3/envs/pandoo/bin/abricate"
abricate_py_path = "/home/liaoth/project/genome_pipelines/toolkit/batch_abricate.py"
############################################################
# database path
############################################################
abricate_db = "/tools/abricate/db/"

ariba_db = "/home/db_public/ariba_db"
ariba_db_list = ["CARD", "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_full", "virulencefinder"]
ariba_str = str({k: "{db}/{name}".format(db=ariba_db, name=k.lower()) for k in ariba_db_list})
############################################################
# unformatted command line
############################################################
fastqc_cmd = "{exe_path} {in_files} -t 2 -o {odir} --quiet"
multiqc_cmd = "{exe_path} {indir} --outdir {odir} --filename {fn} --force -q"
trimmomatic_cmd = """java -jar {exe_path} PE -threads {threads} {R1} {R2} -trimlog {log} {clean_r1} {unpaired_r1} {clean_r2} {unpaired_r2} {params}"""

shovill_cmd = """{exe_path} --outdir {odir} --ram {ram} --R1 {R1} --R2 {R2} --depth {depth} --cpus {thread} --minlen 500 --force"""
# force otherwise it will exit because of pre-created the directory.
prokka_cmd = "{exe_path} {infile} --outdir {odir} --prefix {sn} --force --quiet"
roary_cmd = "{exe_path} -r -v -e --mafft -p {thread} -f {odir} {gff_pattern}"

quast_cmd = """{exe_path} {contig} --circos --gene-finding -1 {R1} -2 {R2} --threads {threads} --no-check -o {odir} {extra_str}"""

pandoo_cmd = """{exe_path} run -i {input} -o {odir} -t -r -b "{ariba_str}" -c {thread}"""
# -t: --infer_tree_on; -r: --abaricate_on; -b: --ariba_dbs;

fasttree_cmd = "{exe_path} -nt -gtr {in_aln} > {o_newick}"

isescan_cmd = "{exe_path} {in_list} {odir}"
abricate_cmd = "python {py_path} -i {indir} -r {roary_dir} -o {odir} -db {db} -mc {mincov} --abricate_path {exe_path} --threads {thread}"