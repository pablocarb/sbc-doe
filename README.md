# SYNBIOCHEM DOE - Design of Experiments tools

## Requirements

1. Create `doe` environment and activate it.

```
conda create --name doe --file doe.conda
conda activate doe
```
2. Install `R` packages `numbers`, `crossdes`, `combinat`, `R.utils`, `planor`,  `DoE.base`.

```
R -e install.packages(c('numbers','crossdes','combinat',\
'R.utils','planor','DoE.base'),\
dependencies=TRUE, repos='http://rforge.net/')
```

## Workflow:

1. Define design specifications in Excel sheet (see templates).

2. Use JMP for the statistical design (optional).

  2.1. Generate JMP scripts:

  ```bash
  python doe2jmp.py inputSpecifications librarySize
  ```

  2.2. Compute DoE in JMP.

3. Generate the design:

```bash
python doepy.py inputSpecifications name -o -b -r -V -j jmpDesign -v "Title" -bro
```

Current internal designs:

* Full factorial;
* Regular fractional factorial;
* Orthogonal arrays.

## Full options


```
usage: doepy.py [-h] [-p] [-s S] [-i] [-r] [-o] [-x [X]] [-O O] [-b] [-g] [-V]
                [-c] [-v V] [-j J] [-w] [-G G] [-k] [-nolab] [-blankPromoter]
                [-bro]
                f id

SBC-DeO. Pablo Carbonell, SYNBIOCHEM, 2016

positional arguments:
  f               Input file with specifications (excel or txt format)
  id              Design id

optional arguments:
  -h, --help      show this help message and exit
  -p              Full positional permutation (default: random latin square)
  -s S            Excel sheet number (default 1)
  -i              Ignore segment calculations based on promoters
  -r              No regular fractional factorial design
  -o              No orthogonal array design
  -x [X]          Random seed (default 100) [or pick random number] for oa
                  design
  -O O            Output path
  -b              Do not generate sbol file
  -g              Generate pigeon cad image
  -V              Generate viscad diagram
  -c              Generate construct fasta files
  -v V            Project description
  -j J            DoE from JMP
  -w              DoE from json (web version)
  -G G            Regenerate pigeon from file and exit
  -k              Keep pigeon files
  -nolab          Do not use labels in pigeon figures
  -blankPromoter  Add blank promoter even if not explicitly given
  -bro            Add file with full list of bridging oligos
```

```
usage: doe2jmp.py [-h] [-O O] [-s S] [-o O] [-r] [-x X] [-n N] [-l L]
                  inputFile libSize

Prepare JMP Files. Pablo Carbonell, SYNBIOCHEM, 2018

positional arguments:
  inputFile   Input file with specifications (excel format)
  libSize     Library size

optional arguments:
  -h, --help  show this help message and exit
  -O O        Output folder (default: same as input)
  -s S        Excel sheet number (default 1)
  -o O        Output file (default: same as input file with jmp extension)
  -r          Overwrite if exists
  -x X        Seed (default random)
  -n N        Number of starts (default 100)
  -l L        Log file
```
