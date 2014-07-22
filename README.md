CoNIFER
=======

CoNIFER (Copy Number Inference From Exome Reads) can discover CNVs from exome read-depth data. See Krumm et al., 2012 (Genome Research) for more information.

SETUP:
virtualenv conifer-env
source conifer-env/bin/activate
git clone https://github.com/matplotlib/matplotlib.git && cd matplotlib
python setup.py build
python setup.py install
pip install -r requirements.txt
pip install -e .