# Boundary Neural Fields Globalization

This is the code for Boundary Neural Fields globalization method. The technical report of the method can be found at http://arxiv.org/pdf/1511.02674v1.pdf 

If you use this software please cite our CVPR 2016 paper:

@InProceedings{gberta_2016_CVPR, 
author = {Gedas Bertasius and Jianbo Shi and Lorenzo Torresani},
title = {Semantic Segmentation with Boundary Neural Fields},
booktitle = {The IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
month = {June},
year = {2016}
}

@InProceedings{gberta_2015_ICCV,
author = {Gedas Bertasius and Jianbo Shi and Lorenzo Torresani},
title = {High-for-Low and Low-for-High: Efficient Boundary Detection from Deep Object Features and its Applications to High-Level Vision},
booktitle = {The IEEE International Conference on Computer Vision (ICCV)},
month = {December},
year = {2015}
}


## Installation

1. VL Feat

	Compile the VL Feat library in the folder 'libs/'

2. Normalized Cuts

	Compile Normalized Cuts code in the directory 'Ncuts_9/''

## Usage

Check out files 'BNF_binary_class_demo.m' and 'BNF_multi_class_demo.m' for binary and multi class segmentations respectively.


## Notes

This version of the code is slightly different than the one presented in the technical report.
