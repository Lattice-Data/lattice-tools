# lattice-tools
External scripts used to interact with the Lattice Database

Environment configuration
---------------- 
Create a new environment and activate it
```
$ conda create -n lattice_submit python=3.7
```
```
$ conda activate lattice_submit
```
*Note: the examples call the environment `lattice_submit` but you can name it anything as long as it is clearly distinguishable from the enviroment you use to launch the encoded app*

Install the following packages
```
$ pip install python-magic requests openpyxl Pillow gspread gspread_formatting oauth2client scanpy
$ pip install google-cloud-storage google-auth-httplib2
$ conda install anndata -c conda-forge
$ conda install -c conda-forge pint
$ conda install pandas
```
Define variables in your environment based on the various servers you might submit to based on an alias for each server
(`ALIAS_KEY`, `ALIAS_SECRET`, `ALIAS_SERVER`). For example, when submitting to a local instance of the app, you might call this `local`.  
So you'd define the following three variables. 
```
$ conda env config vars set LOCAL_KEY=<key>
```
```
$ conda env config vars set LOCAL_SECRET=<secret>
```
```
$ conda env config vars set LOCAL_SERVER=http://localhost:6543
```
After defining those, you'll need to reactivate your environment
```
$ conda activate lattice_submit
```
You can then confirm that they are defined
```
$ conda env config vars list
```
