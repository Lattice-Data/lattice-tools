Schema 5.1.0 QA Testing
----------------
QA validation tests for moving towards schema 5.1.0 migration. Almost all tests will be through pytest; `validate_notebook_workflow.ipynb` checks that validation through the curator notebook workflow does not break.

This iteration utilizes a Makefile to help with environment setup, updating the `cellxgene-schema` cli tool, and running all the pytests.

This Makefile can create both venv and conda envs depending on your preference.

General process
---------------- 
Setup test env, run pytest through the command line, run the notebook.


Installation
---------------- 
To get proper testing environment: 
- Local, up-to-date [CZI single-cell-curation repo](https://github.com/chanzuckerberg/single-cell-curation)
- `pytest` is part of test env
- `cellxgene-schema` pip package is installed/updated to at least 5.0.1 from this repo.

Create a seperate test environment since schema will be > 5.0.0. This can be done by cloning the lattice env or creating a new one.

Running `make` will display a helpful, abbreviated outline of the following steps:

### 1. Create testing virtual venv:
```
make venv
```
OR
#### Create a new conda env named `cxg51testing`
```
make conda
```
OR
#### Clone a conda env
```
make conda-clone
```
By default, this clones a conda env named lattice. Use the following to clone another local conda env:
```
make conda-clone CONDA_ENV=name_of_conda_env
```
### 2. Activate testing env
Activate your env with `source venv/bin/activate` for venv or `conda activate name_of_conda_env`
### 3. Install rest of env dependencies
```
make venv-install
```
OR
```
make conda-install
```
This step will also build `cellxgene-schema` off the `main` branch of the `single-cell-curation` repo

### 4. To update to the latest version of `cellxgene-schema`
```
make cxg-schema
```
This uninstalls and reinstalls `cellxgene-schema` to get the latest version from the `main` branch of the `single-cell-curation` repo. There will likely be a prompt confirming with y/n that you want to uninstall `cellxgene-schema`. You can also use any valid git reference to install off of a specific branch, commit, tag, etc:
```
make cxg-schema GIT_REF=your_git_ref_of_choice
```
The tests assume the standard location of cloned repos that Lattice uses:
```
~/GitClones/CZI/single-cell-curation/
```

Running tests
---------------- 
Make sure the test env is activated.
If immediately following the above directions to create a new conda env,
the env might need to be deactivated and reactivated to make sure pytest uses
the correct testing env.
```
make tests
```
This will uninstall and reinstall the latest version of `cellxgene-schema` and then run all tests with the `-vvv` flag.
Navigate to this directory and run pytest from command line with desired flags.
You can also use the normal cli pytest commands:
```
pytest -vv
```
Cleanup
---------------- 
```
make clean
```
This will `rm -rf` the `venv` directory and any `__pycache__` files.

Run `validate_notebook_workflow.ipynb` to make sure CLI interface of validator still works.