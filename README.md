# FlexiDot: Highly customizable ambiguity-aware dotplots for visual sequence analyses

![alt text](https://github.com/molbio-dresden/flexidot/blob/master/images/Selfdotplots_banner4.png "FlexiDot self dotplots")

FlexiDot is a cross-platform dotplot suite generating high quality self, pairwise and all-against-all visualizations. To improve dotplot suitability for comparison of consensus and error-prone sequences, FlexiDot harbors routines for strict and relaxed handling of mismatches and ambiguous residues. The custom shading modules facilitate dotplot interpretation and motif identification by adding information on sequence annotations and sequence similarities to the images. Combined with collage-like outputs, FlexiDot supports simultaneous visual screening of a large sequence sets, allowing dotplot use for routine screening.


## Citation

Kathrin M. Seibt, Thomas Schmidt, and Tony Heitkam "FlexiDot: Highly customizable ambiguity-aware dotplots for visual sequence analyses". *in prep.*


## Documentation

* in depth documentation (will follow) 
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

Upon first starting FlexiDot, the program calls all needed modules. If absent, it installs them automatically using the Python’s install manager pip. If this fails, please try again with administrator privileges.

## General FlexiDot command

```
python flexidot.py -i input.fas [optional arguments]
```

## Plotting modes

FlexiDot allows sequence investigation in three run modes via the option `-p/--plotting_modes`: 

`-p 0`    self sequence comparison 
`-p 1`    pairwise sequence comparison
`-p 2`    all-to-all sequence comparison


### Self dotplots

with `-p/--plotting_mode 0`

In self dotplot mode, each sequence is compared with itself. The resulting dotplots can be combined to form a collage [default] or written to separate files.

![alt text](https://github.com/molbio-dresden/flexidot/blob/master/images/Selfdotplots_banner.png "FlexiDot self dotplots")

```
python flexidot.py -i test-seqs.fas -p 0 -D y -f 1 -k 10 -w y -r y -x n -m 6 -P 15 -g example.gff3 -G gff_color.config
```


### Pairwise comparisons

with `-p/--plotting_mode 1`

Similar to self dotplots, pairwise dotplots can be combined to a collage or saved to separate output files. The collage output of the 15 pairwise dotplots for the test sequences is shown. By default, dotplot images are in square format (panel A). This maximizes the visibility of matches, if the compared sequences differ drastically in length. To enable scaling according to the respective sequence lengths, the FlexiDot scaling feature is callable via option `-L/--length_scaling` (panel B). If scaling is enabled, a red line indicates the end of the shorter sequence in the collage output. 

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/pairwise_low_res.png" width="600">

```
Panel A$ python flexidot.py -i test-seqs.fas -p 1 -D y -f 0 -k 10 -w y -r y -m 5 -c y -L n 
Panel B$ python flexidot.py -i test-seqs.fas -p 1 -D y -f 0 -k 10 -w y -r y -m 5 -c y -L y
```


### All-against-all comparisons

with `-p/--plotting_mode 2`

In all-against-all mode, FlexiDot compares each pair from a set of input sequences. To enable the identification of long shared subsequences at a glance, FlexiDot offers similarity shading (switched on/off via option `-x/--lcs_shading`) based on the LCS length in all-against-all comparisons. Longer matches are represented by darker background shading. 

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/all_against_all.png" width="500">

```
python flexidot.py -i test-seqs.fas -p 2 -D y -f 0 -t y -k 10 -w y -r y -x y -y 0
```


## Major features

### Mismatch and ambiguity handling

#### Mismatches 

with `-S/--substitution_count [number of allowed mismatches (substitutions)]`

FlexiDot allows toleration of a specified mismatch number per window, thus increasing specificity, yet decreasing specificity.


#### Ambiguities

with `-w/--wobble_conversion Y/N`

FlexiDot handles base ambiguities, often found in consensus sequences. This allows the comparison of species-specific representations of multigene or repeat families as well as common variants or sequence subfamilies. 

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/ambiguities.png" width="500">

```
Panel A$ python flexidot.py -i Seq4.fas -p 0 -D y -f 0 -t y -k 10 -w n -r y -x n
Panel B$ python flexidot.py -i Seq4.fas -p 0 -D y -f 0 -t y -k 10 -w y -r y -x n
```


### Annotation-based shading

In FlexiDot self dotplots, annotated sequence regions can be highlighted by shading to allow clear assignment of dotplot regions to specific sequence contexts (see Seq2 in self dotplots). The underlying annotation information must be provided in general feature format (gff3), either as individual file or file list via the `-g/--input_gff_files` option. To customize GFF-based shading, a user-defined configuration file can be provided via the `-G/--gff_color_config option`. Example files are provided in the test-data directory. Please note, that a legend is generated in a separate file.

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/Selfdotplot_shaded.png" width="500">

```
python flexidot.py -i Seq2.fas -p 0 -D y -f 0 -k 10 -w y -r y -x n -m 12 -P 5 -g example.gff3 -G gff_color.config
```


### Similarity shading

In all-against-all mode, FlexiDot compares each pair from a set of input sequences. To enable the identification of long shared subsequences at a glance, FlexiDot offers similarity shading (switched on/off via option `-x/--lcs_shading`) based on the LCS length in all-against-all comparisons. Longer matches are represented by darker background shading. A separate shading legend output file is created written according to mathematical interval notation, where interval boundaries are represented by a pair of numbers. Consequently, the symbols “(” or “)” represent exclusion, whereas “[” or “]” represent inclusion of the respective number.

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

with `-u/--input_user_matrix_file <matrix.txt>` and `-U/--user_matrix_print y`

coming soon...

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/all_against_all_shaded_orientation_custom_matrix.png" width="680">

```
python flexidot.py -i test-seqs.fas -p 2 -D y -f 2 -t y -k 10 -w y -r y -x y -y 0 -u custom_matrix.txt -U y
```

