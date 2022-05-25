# setting file , we normally just verify the path of executor.

##############################################################
##
## Usually you need to change part.
##############################################################

##Noramlly need to change part
###    Change Change
from os.path import dirname, join

###### server programe path
bcftools_path = ''
vt_pro = ''
vep_pro = ''
bgzip_pro = ''
tabix_pro = ''
gemini_pro = '' # (optional, could stay default)
annovar_pro = '/home/liaoth/tools/annovar'
gatk_pro = '/home-user/thliao/software/gatk-4.2.6.1/gatk'
samtools_pro = '/home-user/thliao/anaconda3/bin/samtools'
samtools_version = 1.12  
fastp = ''
trimmomatic_jar = '/home/liaoth/tools/Trimmomatic-0.36/trimmomatic-0.36.jar'
gatkv36_path = "/home/liaoth/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"

java_option = '-Xmx4g'
##############################################################
##
## Usually you don't need to change part.
##############################################################

## file structure of output result, Normaly don't need to change.
output_fmt = '{path}/{SN}/{SN}'

bip = ''
max_memory = 10240

## software params
trimmomatic_thread = 5
sort_sam_ram = "4G"
sort_sam_thread = 2

gatk_thread = 5
annovar_thread = 1  # be carefull this... each one will take a lot of memory
gemini_thread = 5
vep_thread = 5
java_option = "-Xmx4g"
bwa_thread = 5

## DB files, Normaly don't need to change.
####For cal_fun
cal_sample_name = ''
###### For format file. Normaly don't need to change.
columns_need_to_extract = ['A',
                           'C',
                           'G',
                           'T',
                           'A_Rate',
                           'C_Rate',
                           'G_Rate',
                           'T_Rate']

column_retained = ['Chr',
                   'Start',
                   'End',
                   'Ref',
                   'Alt',
                   'Func.refGene',
                   'Gene.refGene',
                   'GeneDetail.refGene',
                   'ExonicFunc.refGene',
                   'AAChange.refGene',
                   'A',
                   'C',
                   'G',
                   'T',
                   'A_Rate',
                   'C_Rate',
                   'G_Rate',
                   'T_Rate',
                   'snp138',
                   'CLINSIG',
                   'CLNDBN',
                   'CLNACC',
                   'CLNDSDB',
                   'CLNDSDBID',
                   'GT', 'AD', 'AF', 'Filter']

blastn_header = ['qseid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
                 'evalue', 'bitscore']

bed_file_fmt = 'chr\tstart\tend\tgene_name\tstrand\n'
