#!/bin/bash

p=$PWD
## sudo apt-get -y install python3-pyqt4 python-qt4-dev python3-sip-dev libqt4-dev

# yum install -y git python34 python34-devel python34-numpy

pyvenv --without-pip orange3env
source orange3env/bin/activate

Install pip3.4
curl https://bootstrap.pypa.io/get-pip.py | python3.4

echo "/usr/lib/python3/dist-packages/" > "orange3env/lib/python3.4/site-packages/0.pth"
# Update numpy
pip3.4 install --upgrade numpy

# Install scipy
pip3.4 install scipy

git clone https://github.com/biolab/orange3
cd orange3
pip3.4 install -r requirements.txt
#pip3.4 install -r requirements-gui.txt
#pip3.4 install -r requirements-sql.txt
python3.4 setup.py develop
cd ..
git clone https://github.com/biolab/orange-bio
cd orange-bio
python3.4 setup.py develop
