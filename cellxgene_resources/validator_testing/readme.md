Schema QA Testing
----------------
QA validation tests any given migration. Almost all tests will be through pytest

There is a  Makefile to help with environment setup, updating the `cellxgene-schema` cli tool, and running all the pytests.

This Makefile can create both venv and conda envs depending on your preference.

For each schema update, create a new directory in `validator_testing/` and make use of the common `fixtures/` directory

General process
---------------- 
Setup test env, run pytest through the command line, run the notebook.


Installation
---------------- 
To get proper testing environment: 
- Local, up-to-date [CZI single-cell-curation repo](https://github.com/chanzuckerberg/single-cell-curation) - this is useful for making changes to build your own version of `cellxgene-schema`. Otherwise, the validators in the local python env will be used.
- `pytest` is part of test env
- `pytest-xdist` plugin for parallel test running
- `cellxgene-schema` pip package is installed/updated to at least the latest published version from this repo.

Create a seperate test environment from the normal curation env. This can be done by cloning the lattice env or creating a new one.

Running `make` will display a helpful, abbreviated outline of the following steps:

### 1. Create testing virtual venv:
```
make venv
```
OR
#### Create a new conda env named `cxg53testing`
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


Running tests
---------------- 
Make sure the test env is activated.
If immediately following the above directions to create a new conda env,
the env might need to be deactivated and reactivated to make sure pytest uses
the correct testing env.

IMPORTANT:

For fixture imports to work correctly, pytest MUST be run from the `validator_testing/` directory. 
```
pytest 5.3.0/*multi*
```
will work and collect all 5.3 test files with multi in the name.

Changing into the `5.3.0` directory and running `pytest *multi*` WILL NOT WORK
```
make tests
```
This will uninstall and reinstall the latest version of `cellxgene-schema` and then run ALL TESTS (~4500 from 5.0 to 5.3) in parallel with the `-vvv` flag.

Navigate to this directory and run pytest from command line with desired flags.
You can also use the normal cli pytest commands:
```
pytest 5.3.0/*multi* -vv
```
By providing the -n arument to `pytest`, tests will be run in parallel with the `pytest-xdist` plugin:
```
pytest 5.3.0/*atac* -vvv -n auto
```
By providing the -k argument to pytest, you can specify which test function or test class to run
```
pytest 5.3.0/*atac* -vvv -n auto -k TestFragmentCol1Chr
```
Cleanup
---------------- 
```
make clean
```
This will `rm -rf` the `venv` directory and any `__pycache__` files.