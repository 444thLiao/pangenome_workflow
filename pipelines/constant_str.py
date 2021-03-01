#################################################################################
#### File for storage all static path of exe or command line(unformatted)
####
#### If you want to use the pipelines at your computer/workstation/server, Please first check the exe pth first.
# Or run `#todo: `blabla test` to check whether any exe could find

## mlst scheme
specific_species = "Acinetobacter baumannii"

############################################################
# unformatted command line
############################################################
fastqc_cmd = "{exe_path} {in_files} -t 2 -o {odir} --quiet"
multiqc_cmd = "{exe_path} {indir} --outdir {odir} --filename {fn} --force -q {extra_str}"
trimmomatic_cmd = """java -jar {exe_path} PE -threads {threads} {R1} {R2} -trimlog {log} {clean_r1} {unpaired_r1} {clean_r2} {unpaired_r2} {params}"""

shovill_cmd = """{exe_path} --outdir {odir} --ram {ram} --R1 {R1} --R2 {R2} --depth {depth} --cpus {thread} --minlen 500 --force --tmpdir `realpath {odir}` """
# force otherwise it will exit because of pre-created the directory.
prokka_cmd = """{exe_path} {infile} --outdir {odir} --prefix {sn} --locustag {sn} --cpus {thread} --force --quiet"""
roary_cmd = "rm -r {odir}* ;{exe_path} -r -v -e -g 100000 --mafft -p {thread} -f {odir} {gff_pattern} "

quast_cmd = """{exe_path} {contig} -1 {R1} -2 {R2} --threads {threads} --no-check -o {odir} {extra_str}"""

pandoo_cmd = """{exe_path} run -i {input} -o {odir} -t -r -b "{ariba_str}" -c {thread}"""
# -t: --infer_tree_on; -r: --ariba_on; -b: --ariba_dbs;

fasttree_cmd = "{exe_path} -nt -gtr {in_aln} > {o_newick}"

isescan_cmd = """python3 {exe_path} {infile} {proteome_dir} {hmm_dir} -odir {odir} -sn {sn}"""
abricate_cmd = """python3 {py_path} -i {indir} -o {odir} -db {db} -mc {mincov} --abricate_path {exe_path} --threads {thread} {extra_str}"""
plasmid_detect_cmd = """python3 "{py_path}" -i {indir} -o {ofile} """

phigaro_cmd = "{exe_path} -f {infile} -e txt -t {thread} -o {ofile} -p -c '{phigaro_config}' "   # For run this software, we need to remove the y/N

gubbins_cmd = "{exe_path} {infile} --iterations 10 -c {thread} -p {oprefix}"
#################################################################################
#### some output dir
#################################################################################
summary_dir = "summary_output"
#################################################################################
####  parameter
#################################################################################
from multiprocessing import cpu_count
import psutil

total_thread = 20
available_ram = int(psutil.virtual_memory().available / (1024 ** 3))


# p_trimmomatic = int(total_thread-1)
# p_abricate = int(total_thread-1)
mincov_abricate = 80
# p_quast = int(total_thread-1)
# p_shovill = int(total_thread-1)
ram_shovill = 40
# p_roary = int(total_thread-1)
# p_pandoo = int(total_thread-1)
# p_phigaro = int(total_thread - 1)
# p_kraken2 = int(total_thread - 1)

############################################################
# from ncbi genome refseq header
assembly_summary_header = ['assembly_accession',
                           'bioproject',
                           'biosample',
                           'wgs_master',
                           'refseq_category',
                           'taxid',
                           'species_taxid',
                           'organism_name',
                           'infraspecific_name',
                           'isolate',
                           'version_status',
                           'assembly_level',
                           'release_type',
                           'genome_rep',
                           'seq_rel_date',
                           'asm_name',
                           'submitter',
                           'gbrs_paired_asm',
                           'paired_asm_comp',
                           'ftp_path',
                           'excluded_from_refseq',
                           'relation_to_type_material']
############################################################
# from mlst

# From https://github.com/tseemann/mlst/blob/master/db/species_scheme_map.tab
FORCE_MLST_SCHEME = {"Acinetobacter baumannii": ("abaumannii_2",  # i.e., Pasteur
                                                 "abaumannii"),
                     "Achromobacter": "achromobacter",
                     "Aeromonas": "aeromonas",
                     "Aspergillus afumigatus": "afumigatus",
                     "Anaplasma aphagocytophilum": "aphagocytophilum",
                     "Arcobacter": "arcobacter",
                     "Borrelia burgdorferi": "bburgdorferi",
                     "Burkholderia cepacia": "bcc",
                     "Bacillus cereus": "bcereus",
                     "Brachyspira hyodysenteriae": "bhyodysenteriae",
                     "Bifidobacterium bifidobacterium": "bifidobacterium",
                     "Brachyspira intermedia": "bintermedia",
                     "Bacillus licheniformis": "blicheniformis",
                     "Bordetella pertussis": "bordetella",
                     "Brachyspira pilosicoli": "bpilosicoli",
                     "Burkholderia pseudomallei": "bpseudomallei",
                     "Brachyspira": "brachyspira",
                     "Candida albicans": "calbicans",
                     "Campylobacter jejuni": "campylobacter",
                     "Campylobacter coli": "campylobacter",
                     "Clostridium botulinum": "cbotulinum",
                     "Campylobacter concisus": "cconcisus",
                     "Peptoclostridium difficile": "cdifficile",
                     "Clostridium difficile": ("cdifficile",
                                               "cdifficile_2"),
                     "Corynebacterium diphtheriae": "cdiphtheriae",
                     "Campylobacter fetus": "cfetus",
                     "Candida glabrata": "cglabrata",
                     "Campylobacter helveticus": "chelveticus",
                     "Chlamydia": "chlamydiales",
                     "Campylobacter hyointestinalis": "chyointestinalis",
                     "Campylobacter insulaenigrae": "cinsulaenigrae",
                     "Candida krusei": "ckrusei",
                     "Campylobacter lanienae": "clanienae",
                     "Campylobacter lari": "clari",
                     "Cryptococcus neoformans": "cneoformans",
                     "Cronobacter": "cronobacter",
                     "Clostridium septicum": "csepticum",
                     "Clonorchis sinensis": "csinensis",
                     "Campylobacter sputorum": "csputorum",
                     "Candida tropicalis": "ctropicalis",
                     "Campylobacter upsaliensis": "cupsaliensis",
                     "Enterobacter cloacae": "ecloacae",
                     "Escherichia": ["ecoli", "ecoli_2"],
                     "Shigella": "ecoli",
                     "Enterococcus faecalis": "efaecalis",
                     "Enterococcus faecium": "efaecium",
                     "Flavobacterium psychrophilum": "fpsychrophilum",
                     "Haemophilus": "haemophilus",
                     "Helicobacter cinaedi": "hcinaedi",
                     "Haemophilus parasuis": "hparasuis",
                     "Helicobacter pylori": "hpylori",
                     "Haematopinus suis": "hsuis",
                     "Klebsiella oxytoca": "koxytoca",
                     "Klebsiella pneumoniae": "kpneumoniae",
                     "Lactobacillus casei": "lcasei",
                     # "Legionella": "legionella", #STOP. Omp locus problem
                     "Leptospira": ("leptospira", "leptospira_2", "leptospira_3"),
                     "Listeria monocytogenes": "lmonocytogenes",
                     "Lactobacillus salivarius": "lsalivarius",
                     "Mycobacterium abscessus": "mabscessus",
                     "Mycoplasma agalactiae": "magalactiae",
                     "Moraxells catarrhalis": "mcatarrhalis",
                     "Mannheimia haemolytica": "mhaemolytica",
                     "Mycoplasma hyorhinis": "mhyorhinis",
                     "Mycobacterium massiliense": "mmassiliense",
                     "Melissococcus plutonius": "mplutonius",
                     "Neisseria": "neisseria",
                     "Propionibacterium acnes": "pacnes",
                     "Pseudomonas aeruginosa": "paeruginosa",
                     "Pantoea agglomerans": "pagglomerans",
                     "Pseudomonas fluorescens": "pfluorescens",
                     "Propionibacterium freudenreichii": "pfreudenreichii",
                     "Porphyromonas gingivalis": "pgingivalis",
                     "Pasteurella multocida": "pmultocida_multihost",
                     "Pediococcus pentosaceus": "ppentosaceus",
                     "Plesiomonas shigelloides": "pshigelloides",
                     "Streptococcus agalactiae": "sagalactiae",
                     "Staphylococcus aureus": "saureus",
                     "Streptococcus canis": "scanis",
                     "Streptococcus dysgalactiae": "sdysgalactiae",
                     "Salmonella enterica": "senterica",
                     "Staphylococcus epidermidis": "sepidermidis",
                     "Streptococcus gallolyticus": "sgallolyticus",
                     "Sinorhizobium": "sinorhizobium",
                     "Stenotrophomonas maltophilia": "smaltophilia",
                     "Streptococcus oralis": "soralis",
                     "Streptococcus pneumoniae": "spneumoniae",
                     "Staphylococcus pseudintermedius": "spseudintermedius",
                     "Streptococcus pyogenes": "spyogenes",
                     "Streptococcus suis": "ssuis",
                     "Streptococcus thermophilus": ["sthermophilus", "sthermophilus_2"],
                     "Streptomyces": "streptomyces",
                     "Streptococcus uberis": "suberis",
                     "Streptococcus equi": "szooepidemicus",
                     "Taylorella": "taylorella",
                     "Vibrio cholerae": "vcholerae",
                     "Vibrio": "vibrio",
                     "Vibrio parahaemolyticus": "vparahaemolyticus",
                     "Vibrio tapetis": "vtapetis",
                     "Vibrio vulnificus": "vvulnificus",
                     "Wolbachia": "wolbachia",
                     "Xylella fastidiosa": "xfastidiosa",
                     "Yersinia": "yersinia",
                     "Yersinia pseudotuberculosis": "ypseudotuberculosis",
                     "Yersinia ruckeri": "yruckeri"}
