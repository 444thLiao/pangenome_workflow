


python3 /home-user/thliao/software/pangenome_workflow/toolkit/extra_toolkit_luolab/extract_16s.py ./pipelines_o/prokka_o ./16S.fasta
# python3 /home-user/thliao/script/EzbioCloudMatcher/run.py "1155136557@link.cuhk.edu.hk" "12345678" ./16S.fasta
python3 /home-user/thliao/script/EzbioCloudMatcher/run.py  ./16S.fasta
#python3 $HOME/script/EzbioCloudMatcher/run.py "1155136557@link.cuhk.edu.hk" "12345678" ./16S.fasta
python3 /home-user/thliao/software/pangenome_workflow/toolkit/extra_toolkit_luolab/modify_prokka_o.py ./pipelines_o/prokka_o ./new_proteins ./matchResults.csv

#python3 $HOME/software/pangenome_workflow/toolkit/extra_toolkit_luolab/modify_prokka_o.py ./pipelines_o/prokka_o ./new_proteins ./matchResults.csv