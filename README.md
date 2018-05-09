# FlexiDot: Highly customizable, ambiguity-aware dotplots for visual sequence analyses

![alt text](https://github.com/molbio-dresden/flexidot/blob/master/images/Selfdotplots_banner4.png "FlexiDot self dotplots")

FlexiDot is a cross-platform dotplot suite generating high quality self, pairwise and all-against-all visualizations. To improve dotplot suitability for comparison of consensus and error-prone sequences, FlexiDot harbors routines for strict and relaxed handling of mismatches and ambiguous residues. The custom shading modules facilitate dotplot interpretation and motif identification by adding information on sequence annotations and sequence similarities to the images. Combined with collage-like outputs, FlexiDot supports simultaneous visual screening of a large sequence sets, allowing dotplot use for routine screening.


## Citation

Kathrin M. Seibt, Thomas Schmidt, and Tony Heitkam (2018) "FlexiDot: Highly customizable, ambiguity-aware dotplots for visual sequence analyses". Bioinformatics, in press, doi 10.1093/bioinformatics/bty395


## Documentation

* [in depth documentation](https://github.com/molbio-dresden/flexidot/blob/master/documentation/SupplementaryData.pdf) 
* [parameter cheat sheet](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage.pdf)
* [artificial test sequences used for the examples](https://github.com/molbio-dresden/flexidot/tree/master/test-data)


## Implementation

FlexiDot is implemented in [Python 2.7](https://www.python.org/), using 

* [numpy](https://pypi.python.org/pypi/numpy)
* [matplotlib](https://pypi.python.org/pypi/matplotlib)
* [biopython](https://pypi.python.org/pypi/biopython)
* [regex](https://pypi.python.org/pypi/regex)
* [colormap](https://pypi.python.org/pypi/colormap)
* [easydev](https://pypi.python.org/pypi/easydev) (required for colormap)
* [colour](https://pypi.python.org/pypi/colour)

Upon **first starting FlexiDot**, the program calls all needed modules. If absent, it installs them automatically using Python’s install manager pip. If this fails, please try again with **administrator** privileges. 

Please note, that the dependency **Biopython** requires a C compiler. In case of errors during Biopython installation, installing Microsoft Visual C++ Compiler (Windows), GCC (Linux) or Apple’s XCode suite (Mac OS) may help.


## Use FlexiDot

Download the [FlexiDot script](https://github.com/molbio-dresden/flexidot/blob/master/code/flexidot.py).

To run FlexiDot, [**Python 2.7**](https://www.python.org/download/releases/2.7/) must be installed on the machine. 
FlexiDot is started via **command line** in the console. For a brief introduction to the command line interface, check out this nice [tutorial](https://tutorial.djangogirls.org/en/intro_to_command_line/). 

In brief, the console can be started the following way:

* **Windows** 
     * start console: WINDOWS key + type `CMD` + ENTER (Shift + ENTER starts console as administrator)
     * prepare directory
          * select directory and add python script "flexidot.py" and sequence files   
          * copy userpath from address bar (e.g.: C:\Users\Documents\Test)
     * navigate to directory in console: type `cd userpath` + ENTER (paste userpath using right click)
     * start Flexidot with the command below (with your specific fasta file name)
* **Linux/MacOS**
     * start console: Applications → Utilities [Linux] or Accessories [MacOS] → Terminal
     * prepare directory (see above, e.g. /Users/Documents/Test)
     * navigate to directory in console: type `cd userpath` + ENTER (paste userpath using right click)
     * start Flexidot with the command below (with your specific fasta file name)
     
The general FlexiDot command depends on whether one or multiple fasta files are used as input via:

```
# use individual fasta file (can contain multiple sequences)
python flexidot.py -i input.fas [optional arguments]

# use multiple fasta files
python flexidot.py -i input1.fas,input2.fas [optional arguments]

# use all fasta files in current directory
python flexidot.py -a [optional arguments]
```

Optional arguments are explained below and in detail in the [**usage**](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage.pdf). Importantly, `-k` defines the word size (e.g. `-k 10`) and `-t` specifies the sequence type (`-t y` for DNA [default]; `-t n` for proteins). The plotting mode is chosen via `-p` and described below.



## Plotting modes

FlexiDot allows sequence investigation in three run modes via the option `-p/--plotting_mode`: 

`-p 0`    self sequence comparison 
`-p 1`    pairwise sequence comparison
`-p 2`    all-to-all sequence comparison


### Self dotplots

with `-p/--plotting_mode 0`

In **self** dotplot mode, each sequence is compared with itself. The resulting dotplots can be combined to form a **collage** [default] or written to separate files.

![alt text](https://github.com/molbio-dresden/flexidot/blob/master/images/Selfdotplots_banner.png "FlexiDot self dotplots")

```
python flexidot.py -i test-seqs.fas -p 0 -D y -f 1 -k 10 -w y -r y -x n -m 6 -P 15 -g example.gff3 -G gff_color.config
```


### Pairwise comparisons

with `-p/--plotting_mode 1`

For **pairwise** dotplots, the collage output is recommended for larger numbers of sequences. The collage output of the 15 pairwise dotplots for the test sequences is shown below. By default, dotplot images are in square format (panel A). This maximizes the visibility of matches, if the compared sequences differ drastically in length. To enable scaling according to the respective sequence lengths, the FlexiDot scaling feature is callable via option `-L/--length_scaling` (panel B). If scaling is enabled, a red line indicates the end of the shorter sequence in the collage output. 

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/pairwise_low_res.png" width="600">

```
Panel A$ python flexidot.py -i test-seqs.fas -p 1 -D y -f 0 -k 10 -w y -r y -m 5 -c y -L n 
Panel B$ python flexidot.py -i test-seqs.fas -p 1 -D y -f 0 -k 10 -w y -r y -m 5 -c y -L y
```


### All-against-all comparisons

with `-p/--plotting_mode 2`

In **all-against-all** mode, FlexiDot compares each pair from a set of input sequences. To enable the identification of long shared subsequences at a glance, FlexiDot offers similarity shading (switched on/off via option `-x/--lcs_shading`) based on the LCS length in all-against-all comparisons (see below). 

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/all_against_all.png" width="500">

```
python flexidot.py -i test-seqs.fas -p 2 -D y -f 0 -t y -k 10 -w y -r y -x y -y 0
```


## Major features

### Mismatch and ambiguity handling

In diverged or distantly related sequences matches may be interrupted by mismatches or residues might be represented as ambiguities to refer to frequent variants or mutations. Similarly, relaxed matching is helpful when analyzing error-prone sequences like SMRT reads. The achieved relaxation of the matching conditions thus increases sensitivity, while decreasing specificity. 

Firstly, FlexiDot handles **ambiguous residues**, often found in consensus sequences. This allows the comparison of species-specific representations of multigene or repeat families as well as common variants or sequence subfamilies. The ambiguity handling is controlled via`-w/--wobble_conversion Y/N`.

Secondly, a defined number of **mismatches** within the window can be allowed with `-S/--substitution_count [number of allowed mismatches (substitutions)]`. This is even less stringent than the ambiguity handling. Please note, that only substitution mutations are allowed but not indels. 

Lastly, both mismatch and ambiguity handling can be combined for the analysis.

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/Fig-Suppl-MismatchesWobbles.png" width="600">

```
Panel tl$ python flexidot.py -i Seq4.fas,Seq1.fas -p 1 -D n -f 0 -c n -k 10 -w n -r y -x n
Panel tm$ python flexidot.py -i Seq4.fas,Seq1.fas -p 1 -D n -f 0 -c n -k 10 -w n -r y -x n -S 1
Panel tr$ python flexidot.py -i Seq4.fas,Seq1.fas -p 1 -D n -f 0 -c n -k 10 -w n -r y -x n -S 2
Panel bl$ python flexidot.py -i Seq4.fas,Seq1.fas -p 1 -D n -f 0 -c n -k 10 -w y -r y -x n
Panel bm$ python flexidot.py -i Seq4.fas,Seq1.fas -p 1 -D n -f 0 -c n -k 10 -w y -r y -x n -S 1
Panel br$ python flexidot.py -i Seq4.fas,Seq1.fas -p 1 -D n -f 0 -c n -k 10 -w y -r y -x n -S 2
```


### Annotation-based shading

In FlexiDot self dotplots, annotated sequence regions can be highlighted by **shading** to allow clear assignment of dotplot regions to specific sequence contexts (see Seq2 in self dotplots). The underlying **annotation** information must be provided in general feature format (**gff3**), either as individual file or file list via the `-g/--input_gff_files` option. To customize GFF-based shading, a user-defined configuration file can be provided via the `-G/--gff_color_config option`. Example files are provided in the test-data directory. Please note, that a legend is generated in a separate file.

If you wish to find out more on the gff3 file format used here, Ensembl provides a [good overview](https://www.ensembl.org/info/website/upload/gff3.html).

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/Selfdotplot_shaded.png" width="500">

```
python flexidot.py -i Seq2.fas -p 0 -D y -f 0 -k 10 -w y -r y -x n -m 12 -P 5 -g example.gff3 -G gff_color.config
```


### Similarity shading

In all-against-all mode, FlexiDot compares each pair from a set of input sequences. To enable the identification of long shared subsequences at a glance, FlexiDot offers similarity shading (switched on/off via option `-x/--lcs_shading`) based on the **LCS length** (longest common subsequence, or longest match if mismatches are considered) in all-against-all comparisons. Longer matches are represented by darker background shading. A separate shading **legend** output file is created written according to mathematical interval notation, where interval boundaries are represented by a pair of numbers. Consequently, the symbols “(” or “)” represent exclusion, whereas “[” or “]” represent inclusion of the respective number.

FlexiDot similarity shading is highly customizable with the following parameters, explained in depth in the documentation:
* Reference for shading (option `-y/--lcs_shading_ref`)
* Number of shading intervals (option `-X/--lcs_shading_num`)
* Shading based on sequence orientation (option `-z/--lcs_shading_ori`)

Shading examples based on sequence orientation (forward, panel A; reverse, panel B; both, panel C) are shown:

![alt text](https://github.com/molbio-dresden/flexidot/blob/master/images/all_against_all_shaded_orientation2.png "FlexiDot shaded dotplots")

```
Panel A$ python flexidot.py -i test-seqs.fas -p 2 -D y -f 0 -t y -k 10 -w n -r y -x y -y 0 -z 0
Panel B$ python flexidot.py -i test-seqs.fas -p 2 -D y -f 0 -t y -k 10 -w n -r y -x y -y 0 -z 1
Panel C$ python flexidot.py -i test-seqs.fas -p 2 -D y -f 0 -t y -k 10 -w n -r y -x y -y 0 -z 2
```

### Custom matrix shading

When comparing related sequences, multiple sequence alignments are frequently applied. The resulting pairwise **sequence similarities** can be integrated in the FlexiDot images by providing a **matrix file** via `-u/--input_user_matrix_file <matrix.txt>`. This allows a shading of the upper right triangle according to the matrix (here orange). With `-U/--user_matrix_print y` the matrix values can be printed into the respective fields. Besides, also **text** information can be provided in the matrix, but then shading is suppressed.

In the example, LCS and matrix shading are combined to visualize the relationships between different members of a repeat family. 

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/Beetle_matrix_shading.png" width="750">

```
python flexidot.py -i Beetle.fas -p 2 -x y -k 10 -S 1 -r n -u custom_matrix.txt -U y
```

