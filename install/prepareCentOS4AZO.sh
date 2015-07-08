echo "Installing dependent packages..."
yum install -y epel-release
yum install -y tcsh
yum install -y openssh-clients
yum groupinstall -y "Development Tools"
yum install -y gcc
yum install -y gcc-c++
yum install -y gcc-gfortran
yum install -y nspr
yum install -y mpich  
#yum install -y libcr-dev
yum install -y swig	
yum install -y atlas-devel
yum install -y python-devel
yum install -y numpy
yum install -y sip
yum install -y python-networkx 
yum install -y PyQt4
yum install -y lapack-devel
#yum install -y python-qwt5-qt4
#yum install -y git-core gitosis
yum install -y git
yum install -y mercurial
yum install -y nspr-devel
yum install -y tbb-devel
yum install -y libjpeg-turbo-devel
yum install -y libpng-devel
yum install -y libtiff-devel
yum install -y jasper-devel
yum install -y libdc1394-devel
yum install -y libv4l-devel
yum install -y PyQt4-webkit-devel
# hg clone ssh://hg@bitbucket.org/biolab/orange
# hg clone https://bitbucket.org/biolab/orange

# ==========  Needed by cinfony ===========
#Openbabel
yum install -y openbabel
yum install -y python-openbabel
#RDKit
yum install -y cmake
yum install -y bison
yum install -y flex
yum install -y sqlite
yum install -y sqlite-devel
#yum install -y libboost-all-dev
yum install -y boost-devel
#CDK
yum install -y python-jpype
# ==========  Needed by cinfony ===========

#fminer
yum install -y openbabel-devel
yum install -y gsl-devel

sudo ln -s /usr/include/libv4l1-videodev.h   /usr/include/linux/videodev.h 

echo "Updating db..."
sudo updatedb
echo ""
echo "============================================================="
echo "             Finished preparation of CentOS"
echo "-------------------------------------------------------------"
echo "Configure Git if you plan to check in changes:"
echo '      git config --global user.name "YouGitUserName"'
echo "      git config --global user.email YourEmail@ServerX"
echo ""

