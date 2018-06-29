# FlexiDot version changes

![alt text](https://github.com/molbio-dresden/flexidot/blob/master/images/Selfdotplots_banner4.png "FlexiDot self dotplots")

## Version 1.04
*29.06.2018* 

* [new parameter cheat sheet v1.04](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage_v1.04.pdf) 
* [new FlexiDot script v1.04](https://github.com/molbio-dresden/flexidot/blob/master/code/flexidot_v1.04.py)


**[New feature] Graphic formatting options for all-against-all dotplots:**  
On request we added two parameters: With `-M/--mirror` it is now possible to mirror the middle diagonal.

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/all_against_all_Flexi1.04_Para_Mirror.png" width="700">

Basic commands:
```
python flexidot.py -i test-seqs.fas -p 2 -M n
python flexidot.py -i test-seqs.fas -p 2 -M y
```
Command plus aesthetics as shown here (as described in version update 1.03):
```
python flexidot.py -i test-seqs.fas -p 2 -M n -g example2.gff3 -G gff_color.config -x y -k 10 -F 0.06 -A 1.5
python flexidot.py -i test-seqs.fas -p 2 -M y -g example2.gff3 -G gff_color.config -x y -k 10 -F 0.06 -A 1.5
```

The `-R/--representation` parameter allows partial dotplotting, either printing the complete `-R 0`, the top `-R 1`, or the bottom dotplot `-R 2`. 

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/all_against_all_Flexi1.04_Para_Representation.png" width="900">

Basic commands:
```
python flexidot.py -i test-seqs.fas -p 2 -R 0
python flexidot.py -i test-seqs.fas -p 2 -R 1
python flexidot.py -i test-seqs.fas -p 2 -R 2
```
Command plus aesthetics as shown here (as described in version update 1.03):
```
python flexidot.py -i test-seqs.fas -p 2 -R 0 -g example2.gff3 -G gff_color.config -x y -k 10 -F 0.06 -A 1.5
python flexidot.py -i test-seqs.fas -p 2 -R 1 -g example2.gff3 -G gff_color.config -x y -k 10 -F 0.06 -A 1.5
python flexidot.py -i test-seqs.fas -p 2 -R 2 -g example2.gff3 -G gff_color.config -x y -k 10 -F 0.06 -A 1.5
```

**[Bugfix]:**  

We also fixed a distortion issue in `-p/--plotting_mode 0` (self dotplots).


## Version 1.03
*17.06.2018* 

* [new parameter cheat sheet v1.03](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage_v1.03.pdf) 
* [new FlexiDot script v1.03](https://github.com/molbio-dresden/flexidot/blob/master/code/flexidot_v1.03.py)


**[New feature] Annotation-based shading also available for all-against-all dotplots:**   
Previously only available for self dotplots, we added annotation-based shading to all-against-all dotplots, allowing for many new visualizations. As before, annotation information is provided as general feature file (GFF3). These features are added to the middle diagonal (see our example below).

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/all_against_all_annotation_based_shading_cool.png" width="700">

Basic command:
```
python flexidot.py -i test-seqs.fas -g example2.gff3 -G gff_color.config -p 2
```

Command plus aesthetics as shown here (+ LCS shading, wordsize 10, change of subplot spacing and line width):
```
python flexidot.py -i test-seqs.fas -g example2.gff3 -G gff_color.config -p 2 -x y -k 10 -F 0.06 -A 1.5
```

The test files used here are provided:
* [test-seqs.fas](https://github.com/molbio-dresden/flexidot/blob/master/test-data/test-seqs.fas)
* [example2.gff3](https://github.com/molbio-dresden/flexidot/blob/master/test-data/example2.gff3)
* [gff_color.config](https://github.com/molbio-dresden/flexidot/blob/master/test-data/gff_color.config)



## Version 1.02 
*09.05.2018*

* [new parameter cheat sheet v1.02](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage_v1.02.pdf) 
* [new FlexiDot script v1.02](https://github.com/molbio-dresden/flexidot/blob/master/code/flexidot_v1.02.py)

Changed handling of `-T` parameter: The character count of the sequence titles has been limited to `20` by default. This limit can be changed with `-T`. If an `E` (end) is added to the limit, the last characters are chosen instead of the first. 

```
-T 20  (the first 20 characters)     
-T 20E (the last 20 characters)
```


## Version 1.01 
*21.04.2018*


* [parameter cheat sheet v1.01](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage_v1.01.pdf)
* [FlexiDot script v1.01](https://github.com/molbio-dresden/flexidot/blob/master/code/flexidot_v1.01.py)

minor bugfixing



## Version 1.00 
*21.03.2018*

first FlexiDot release
