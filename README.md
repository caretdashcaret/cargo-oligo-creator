# CARGO Olgio Creator

CARGO Oligo Creator is a script to generate oligos for gRNAs for the [CARGO cloning method](http://science.sciencemag.org/content/early/2018/01/24/science.aao3136/tab-figures-data).

## Setup

Set up the virtualenv:

```shell
$ python3 -m venv .cargo-oligo-creator
```

Activate

```shell
$ source .cargo-oligo-creator/bin/activate
```

Install requirements

```shell
$ pip3 install -r requirements.txt
```

### Run

To run, pass in your guides (gRNA) as the the input `-i`.

Example:

```shell
$ python3 run.py -i GGGCGAGGAGCTGTTCACCG GCTGCACGCCGTAGGTCAGGG GGTGAACCGCATCGAGCTGA GGTGTTCTGCTGGTAGTGGT
```

Results will be printed to standard out.

If you don't want the `|` separator and instead want the raw sequence, you can pass in `-r`.

Example:

```shell
$ python3 run.py -r -i GGGCGAGGAGCTGTTCACCG GCTGCACGCCGTAGGTCAGGG GGTGAACCGCATCGAGCTGA GGTGTTCTGCTGGTAGTGGT
```

### Tests

To run tests:

```shell
pytest
```
