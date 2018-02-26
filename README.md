# FlexiDot: Highly customizable ambiguity-aware dotplots for visual sequence investigation

![alt text](https://github.com/molbio-dresden/flexidot/blob/master/images/Selfdotplots_banner.png "FlexiDot self dotplots")



FlexiDot is a cross-platform dotplot suite generating high quality self, pairwise and all-against-all visualizations by exact matching. To improve dotplot suitability for comparison of consensus and error-prone sequences, FlexiDot harbors routines for strict and relaxed handling of ambiguous residues. Our two custom shading modules facilitate dotplot interpretation and motif identification by adding information on sequence annotations and sequence similarities to the images. Combined with collage-like outputs, FlexiDot supports simultaneous visual screening of a large sequence sets, allowing dotplot use for routine screening.

## Implementation

FlexiDot is implemented in Python 2.7, using 

* [numpy](https://pypi.python.org/pypi/numpy)
* [matplotlib](https://pypi.python.org/pypi/matplotlib)
* [biopython](https://pypi.python.org/pypi/biopython)
* [colour](https://pypi.python.org/pypi/colour)

Upon first starting FlexiDot, the program tries to call all needed modules. If absent, it installs them automatically via the Python’s install manager pip. If this fails, please try again with administrator privileges.

## General FlexiDot command

`python flexidot.py -i input.fas [optional arguments]`


