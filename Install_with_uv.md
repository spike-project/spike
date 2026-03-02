# installing Spike using `uv`

We'll use `uv` to install and update `Spike`. `uv`  is a powerfull package and project manager for Python, which replaces tools like `pip`,  `poetry`, `pyenv`, `virtualenv` etc...    
Check [docs.astral.sh/uv](https://docs.astral.sh/uv/)

We will use the terminal to do it... *(using the MacOs or Linux syntax)*

**Remark** if you already have a `conda` installation (from `anaconda`), you should first deactivate it:

	conda deactivate


### install `uv` if not installed yet.
The simpler is to follow instructions here : [docs.astral.sh/uv/getting-started/installation/](https://docs.astral.sh/uv/getting-started/installation/)    
For instance using the `pip` or the *standalone* methods, but any should do.

### create a specific environment
create a folder, and move into it

	mkdir Spike; cd Spike

### create a virtual environment
This will isolate from the whole computer every thing we install in the environment.
We'll use an uv anonymous environment.

	uv venv
and activate it :

	source .venv/bin/activate
	
### insure that you are using python 3.12 of higher
	uv python install '>= 3.12'

### then install Spike and dependencies
First dependencies

	uv pip install jupyter lab nbclassic scipy tables pandas matplotlib ipympl threadpoolctl
then `Spike` itself,

- either from pip repository

		uv pip install spike-py

- or if you prefer, downloading the source from github :

		uv pip install "git+https://github.com/spike-project/spike"
	
if you have the DATA tests files, you can launch the test suite:

	python -m spike.Tests  -D  xx-your-dir-xx/DATA_test

## using the program
go the the installation directory, 
deactivate conda if needed

	conda deactivate
activate the uv installation

	source .venv/bin/activate
that's it




mardi, 24. février 2026 04:11 
