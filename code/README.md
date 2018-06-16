# FlexiDot version changes

![alt text](https://github.com/molbio-dresden/flexidot/blob/master/images/Selfdotplots_banner4.png "FlexiDot self dotplots")

## Version 1.03 (stub, uploads will follow)
*xx.06.2018* 

* [new parameter cheat sheet v1.03](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage_v1.03.pdf) 

**New feature**: Annotation-based shading is also available for all-to-all dotplots, allowing for many new visualizations. Annotation information will be added to the middle diagonal as in our example below.

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/all_against_all_annotation_based_shading_cool.png" width="700">

```
command will follow
```


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
