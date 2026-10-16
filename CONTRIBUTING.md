# Contributing to ASALIPY

Thank you for contributing to ASALIPY! This is largely a volunteer effort, and small fixes to documentation, examples and tests are very welcome.

ASALIPY is a Python package for chemical reactor models. It is part of a small ecosystem:

- [ASALI](https://github.com/srebughini/ASALI): the GUI and the thermodynamic and transport property library in several languages. Contributions to those belong there; see its `CONTRIBUTING.md`.
- [JASALI](https://github.com/srebughini/JASALI): the JavaScript version of the property library.

## Set up

ASALIPY needs [Anaconda](https://www.anaconda.com/), because stable versions of Cantera and Assimulo are not available on PyPI (see `README.md`). To work on the source code, use the environment file from the repository root, as described in `README.md`:

```bash
conda env create -f environment.yml
conda activate asali
```

`environment.yml` asks for Python `>=3.7,<=3.12`, Cantera, Assimulo, termcolor, matplotlib, pyyaml, sympy, pandas and openpyxl.

## Run the tests

Run the commands from the repository root. The tests read their input files from the `tests/` folder.

```bash
python run_tests.py                                   # all tests
python run_tests.py --class BatchReactor              # one class
python run_tests.py --class BatchReactor --local      # see below
```

## Check the documentation

If you change `README.md` or `EXAMPLES.md`, check that the links to the examples point to existing files. From the repository root:

```bash
grep -o -h 'blob/main/examples/[A-Za-z0-9_./-]*' EXAMPLES.md README.md | sort -u | while read l; do
  f=${l#blob/main/}; [ -e "$f" ] && echo "OK      $f" || echo "BROKEN  $f"
done
```

The example scripts in `examples/` use paths such as `os.path.join('files', ...)`, which are relative to the `examples/` folder.

## Find something to work on

Browse the [open issues](https://github.com/srebughini/ASALIPY/issues), and look for the labels `good first issue` and `help wanted`. Before you start, check the comments of the issue and the [open pull requests](https://github.com/srebughini/ASALIPY/pulls), and leave a comment saying that you would like to work on it.

## Reporting bugs and asking questions

Open a [GitHub issue](https://github.com/srebughini/ASALIPY/issues). Please include:

* The exact steps that reproduce the problem, and the example or input file you used.
* What you see, and what you expected to see.
* Your operating system, your Python and Cantera versions, and how you installed ASALIPY (conda package or `environment.yml`).

## Submitting changes

* Fork the repository, make your changes in a branch of your fork, and open a [pull request](https://github.com/srebughini/ASALIPY/pulls) against the `main` branch.
* Keep each pull request focused on one issue, and include the issue number in the description.
* Describe the problem and the solution, and give the exact commands you ran to check your change and their results.
* Use the same commit message and pull request conventions as ASALI: present tense and imperative mood ("Add feature"), and reference issues liberally.

The project has not defined a code style, a formatter or a review process. Please follow the style of the file you are editing.
