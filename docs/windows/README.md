# Development on Windows

To develop on Windows you will have to install a number of dependencies. I have not quite gotten everything to work, but the following should at least get you started.

## 1. Install conda for Windows

[conda][] is a cross-platform package and environment management software. It lets you install packages (that have been build for conda) as well as dependenies. It also lets you manage these packages in different environments which can be loaded and unloaded. This lets you maintain separate versions of packages for projects where they may conflict (e.g., Python 3.8 for one project and Python 3.9 for another project).

The easiest way to install is to download the [miniconda][] version of **conda** (the **anaconda**/full version just has more packages built-in by default). You can download the version from the [miniconda][] page (you can download the latest 64-bit Python version since you can install different versions of Python in conda).

## 2. Launch conda for Windows

Once installed you should be able to find an **Anaconda Prompt (Miniconda3)** application which will open up a terminal command-line prompt with conda loaded up for you to use.

![conda-prompt.png][]

The conda command is `conda`. You can run `conda info --envs` to see a list of environments available for you to use. At the beginning there will just be one called **base** but we will install a new one with the necessary dependencies to work with the Genomics Data Index software.

![conda-terminal.png][]

## 3. Install dependencies in conda

To install a new environment in **conda** we can use the `conda env create` command.

```bash
conda env create --name gdi --file environment.yml
```

Here the `--name gdi` is the name of the environment and `--file environment.yml` lists the file containing a list of dependencies to install. This can be found at [environment.yml][].

Conda should download and install all the necessary packages (listed in `environment.yml`). This may take a while to complete.

## 4. Activate environment

Once all the dependencies are installed, you can activate the environment with `conda activate [environment name]`.

```bash
conda activate gdi
```

This will change the prefix of the string shown when you run commands in the terminal to show you are in the **gdi** environment:

![conda-prefix.png][]

## 5. Install `genomics-data-index`

Finally, let's install the Python `genomics-data-index` package.

```bash
pip install --no-deps genomics-data-index
```

The `--no-deps` option tells **pip** to ignore any dependencies required by `genomics-data-index` (most of these are installed by conda, but some I haven't gotten to work in Windows so I want to avoid trying to install them for now).

## 6. Test out `genomics-data-index`

Once you've installed `genomics-data-index` you can test out like:

```bash
gdi --version
```


[conda]: https://docs.conda.io/projects/conda/en/latest/index.html
[miniconda]: https://docs.conda.io/en/latest/miniconda.html
[conda-prompt.png]: images/conda-prompt.png
[conda-terminal.png]: images/conda-terminal.png
[environment.yml]: environment.yml
