CoNIFER
=======

CoNIFER (Copy Number Inference From Exome Reads) can discover CNVs from exome read-depth data. See Krumm et al., 2012 (Genome Research) for more information.

===============
SETUP Instructions for MAC: (work in progress)

export PYENV="conifer_env"

mkdir ~/Envs
pip install virtualenv
pip install virtualenvwrapper
#add to .bash_profile: export WORKON_HOME=~/Envs
#add to .bash_profile: source /usr/local/bin/virtualenvwrapper.sh
mkvirtualenv $PYENV
pip install numpy
pip install numexpr
pip install Cython
pip install tables
cd ~/Envs/$PYENV/
git clone https://github.com/matplotlib/matplotlib.git && cd matplotlib
python setup.py build
python setup.py install
pip install scipy
pip install pysam

#Example run

#Collect some count datoa
./pull_data.sh

#Run first step
python conifer.py analyze --probes targeted_data/BROCA_v6_0465501.bed --rpkm_dir targeted_data/ --output targeted_data/test.hdf5 --svd 0 --write_svals targeted_data/singular_values.txt --plot_scree targeted_data/screeplot.png --write_sd targeted_data/sd_values.txt
