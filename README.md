CoNIFER
=======

CoNIFER (Copy Number Inference From Exome Reads) can discover CNVs from exome read-depth data. See Krumm et al., 2012 (Genome Research) for more information.

===============
SETUP:
#TODO
# Try to get install via Noah's wheel stuff.
# Add his HDF5 code back

#ZLIB
cd ~/src
wget http://zlib.net/zlib-1.2.8.tar.gz && tar -xf zlib-1.2.8.tar.gz
cd zlib-1.2.8 && ./configure --prefix=$HOME/local && make test && make install

#SZIP
cd ~/src
wget http://www.hdfgroup.org/ftp/lib-external/szip/2.1/src/szip-2.1.tar.gz && tar -xf szip-2.1.tar.gz
cd szip-2.1 && ./configure --prefix=$HOME/local && make test && make install

# HDF5
cd ~/src
wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.13.tar.gz && tar -xf hdf5-1.8.13.tar.gz && cd hdf5-1.8.13
export LDFLAGS=-L/usr/local/opt/readline/lib
export CPPFLAGS=-I/usr/local/opt/readline/include
./configure --prefix=$HOME/local/hdf5 --enable-fortran --enable-cxx --with-szlib=$HOME/local/lib
make && make check && make install && make check-install

virtualenv conifer-env
source conifer-env/bin/activate
git clone https://github.com/matplotlib/matplotlib.git && cd matplotlib
python setup.py build
python setup.py install

#Not working, so using individual installs
#pip install -r requirements.txt
#pip install -e .

pip install numpy
pip install numexpr
pip install Cython
export HDF5_DIR=$HOME/local/hdf5
pip install tables
pip install pysam
pip install scipy

#Example run
export LD_LIBRARY_PATH=/home/local/AMC/hermands/local/hdf5/lib

#Collect some count data
dev/pull_data.sh

#Run first step
bpython conifer.py analyze --probes targeted_data/BROCA_v6_0465501.bed --rpkm_dir targeted_data/ --output targeted_data/test.hdf5 --svd 0 --write_svals targeted_data/singular_values.txt --plot_scree targeted_data/screeplot.png --write_sd targeted_data/sd_values.txt
