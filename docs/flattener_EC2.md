# flattener.py on EC2

Starting up EC2
----------------
Startup instance, using Amazon Linux 2 AMI, lattice-R-jupyter security group and lattice_ec2.pem


Anaconda installation
----------------
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
bash Anaconda3-2021.11-Linux-x86_64.sh
source .bashrc
```

Jupyter Notebook installation
----------------
Create and activate scanpy notebook environment (https://dataschool.com/data-modeling-101/running-jupyter-notebook-on-an-ec2-server/):
```
conda create -n scanpy python=3.7
conda activate scanpy
pip install scanpy
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph leidenalg notebook
```

Configure Jupyter notebook settings, first create config file and run IPython:
```
jupyter notebook --generate-config
ipython
```
Enter IPython the following commandline. Enter your password when prompted. COPY HASH AND SAVE IT FOR LATER
```
from IPython.lib import passwd
passwd()
exit
```

Open jupyter config file:
```
cd .jupyter
vim jupyter_notebook_config.py
```
Paste the following, with your hash replaced (including 'sha1:'):
```
conf = get_config()
conf.NotebookApp.open_browser = False
conf.NotebookApp.ip = '0.0.0.0'
conf.NotebookApp.password = u'YOUR PASSWORD HASH'
conf.NotebookApp.port = 8888
```

Run notebook:
```
jupyter notebook
```
Go to website: http://[your Public DNS]:8888/, and use password that you has originally put in.


Python and R installation requirements
----------------
Create and activate lattice\_submit environment as documented on https://github.com/Lattice-Data/lattice-tools. Additional python library to install is:
```
pip install rpy2
pip install boto3
```
Install R and Seurat:
```
sudo amazon-linux-extras install R4
sudo yum -y install libcurl-devel
sudo yum install openssl-devel
sudo R -e "install.packages('curl', repos='http://cran.us.r-project.org')"
sudo R -e "install.packages('openssl', repos='http://cran.us.r-project.org')"
sudo R -e "install.packages('httr', repos='http://cran.us.r-project.org')"
sudo R -e "install.packages('Seurat', repos='http://cran.us.r-project.org', Ncpus=4)"
sudo amazon-linux-extras install epel -y
sudo yum install libxml2-devel
```

Install correct version of hdf5 fortran wrappers and SeuratDisk
```
wget https://cbs.centos.org/kojifiles/packages/hdf5/1.8.13/7.el7/x86_64/hdf5-1.8.13-7.el7.x86_64.rpm
wget ftp://ftp.pbone.net/mirror/ftp.centos.org/7.9.2009/cloud/x86_64/openstack-queens/Packages/h/hdf5-devel-1.8.13-7.el7.x86_64.rpm
sudo yum localinstall hdf5-1.8.13-7.el7.x86_64.rpm
sudo yum localinstall hdf5-devel-1.8.13-7.el7.x86_64.rpm
sudo R -e "install.packages('remotes', repos='http://cran.us.r-project.org')"
sudo R -e "remotes::install_github('mojaveazure/seurat-disk')"
```

R Studio installation
----------------
```
wget https://download2.rstudio.org/server/centos7/x86_64/rstudio-server-rhel-2021.09.1-372-x86_64.rpm
sudo yum install rstudio-server-rhel-2021.09.1-372-x86_64.rpm
sudo useradd <username>
sudo passwd <username>
```
Add 'www-port=80' to this file: /etc/rstudio/rserver.conf:
```
sudo vi /etc/rstudio/rserver.conf
```
Restart R studio:
```
sudo rstudio-server restart
```
Go to website: http://[your Public DNS]:80

AWS credentials
----------------
From your local computer, copy files to ec2:
```
scp -i lattice_ec2.pem flattener.py lattice.py ~/.aws/credentials ec2-user@<Public DNS>:~/
```

On ec2, create .aws directory, and move credentials there:
```
mkdir ~/.aws
mv credentials ~/.aws/
```







