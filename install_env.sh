# todo: finish it and test it
conda create -n pangenome_workflow --file environment.yml
export toolspath=/media/liaoth/tools

cd $toolspath
git clone https://github.com/tseemann/mlst.git
export mlst_path=$toolspath/mlst



