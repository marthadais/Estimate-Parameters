Copyright 2017 Martha Dais Ferreira, Rodrigo Fernandes de Mello

This file is part of ImageFalseNearest.

ImageFalseNearest is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ImageFalseNearest is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ImageFalseNearest.  If not, see <http://www.gnu.org/licenses/>.

See the file "COPYING" for the text of the license.

Contact: 
	Martha D. Ferreira: daismf@icmc.usp.br
	Rodrigo Fernandes de Mello: mello@icmc.usp.br

----------------------------------------------
# Estimate CNN Parameters

This source codes are used to estimate the most adequate convolution mask sizes and number of units for a convolution layer of CNN. Based on the False Nearest Neighbors (FNN), a well-known tool from the area of Dynamical Systems, this method helps estimating CNN architectures that are both less complex and more effective.

## Getting Started

These instructions of how to use package ImageFalseNearest come with a copy of this package, which contains our FNN approach that supports the estimation of CNN architectures. All source codes and instructions permit to run the experiments on your local machine.

## Requirements

Before starting, it is necessary to have the R Statistical Software installed, some R packages, and then package ImageFalseNearest on your local machine. 

 1. R language:

```
sudo apt-get install r-base
```

 2. R packages (need be in the R environment -- 'sudo R'):

```
install.packages(c('png', 'jpeg','pixmap'))
```
	
 3. Package ImageFalseNearest needs to be installed. This can be done following the README provided with this project (folder 'imagefalsenearest'):

```
cd imagefalsenearest
sudo R CMD build ImageFalseNearest
sudo R CMD INSTALL ImageFalseNearest_0.1-01.tar.gz
```

 4. If you want run the CNN training, it is necessary to install the Caffe deep learning framework, according the instructions available at <http://caffe.berkeleyvision.org/installation.html>.

## Estimating Parameters for CMU Face Images dataset

This step is used to estimate the convolutional mask size and the number of units for the first convolutional layer of the CNN architecture. The source code used to estimate the convolutional mask size requires the name of the dataset and a file with the path to the input images. Source codes maskEstimation.r and unitsEstimation.r are prepared to read PGM files, however it is possible to modify both codes to deal with other file extension. 

If you want make the estimation over PNG or JPG extension files, you need to modify files maskEstimation.r and unitsEstimation.r. In maskEstimation.r, you need comment lines 137,186,190 and uncomment lines 135,184,188 for PNG extension or lines 136,185,189 for JPG extension. In unitsEstimation.r, you need comment lines 95,143,147 and uncomment lines 93,141,145 for PNG extension or lines 94,142,146 for JPG extension.

Lines:

```
.
.
.
# img = readPNG(as.character(Files[[1]][1])) #png
# img = as.matrix(readJPEG(as.character(Files[[1]][1]))) #jpg
img = read.pnm(file = as.character(Files[[1]][1]),cellres=1)@grey #pgm
.
.
.
if(channels==3){
# img = readPNG(as.character(files[i]))[,,color] #png
# img = as.matrix(readJPEG(as.character(files[i]))) #jpg
img = read.pnm(file = as.character(files[i]),cellres=1)@grey #pgm
} else {
# img = readPNG(as.character(files[i])) #png
# img = as.matrix(readJPEG(as.character(files[i]))) #jpg
img = read.pnm(file = as.character(files[i]),cellres=1)@grey #pgm
}
```

The example to estimate the parameters is provided in R language, thus the following commands need to be running in the R environment. In this context, we suggest the download of the CMU Face Images dataset (available at <https://archive.ics.uci.edu/ml/datasets/CMU+Face+Images>), in which the image files are in PGM extension.

### Convolutional Mask Sizes Estimation

While in the R environment (terminal command 'R'), we load the mask evaluation file (maskEstimation.r) that will read the images and apply the ImageFalseNearest method to each individual image and compute the stop criterion. This evaluation is applied on each color channel of every image under analysis. The mask size estimation requires the name of the dataset, and a file containing the path for images.

Steps for evaluating the convolutional mask sizes:

```
source("maskEstimation.r")

#Running evaluation of mask size
maskRes = mask_evaluation("CMU","cmu-list.txt")
```

This command will take a bit longer to finish, and a folder will be create to save the results. During the mask evaluation, output messages will appear, indicating image names that are under analysis, and the error produced after it. Some part of output is present below:

```
-------Compute FNN for mask estimation.-------
Analysing  440  files
Maximum Mask Selected:  15 x 15 

Analysing color  GRAY 
file[ 1 ]= bpm_right_neutral_sunglasses_4 
Calculating mask size...
Frobenius Norm:  2.915476 
file[ 2 ]= at33_right_angry_open_4 
Calculating mask size...
Frobenius Norm:  0.3535534 
file[ 3 ]= cheyer_up_happy_sunglasses_4 
Calculating mask size...
Frobenius Norm:  0.4859127 
file[ 4 ]= phoebe_right_neutral_open_4 
Calculating mask size...
	.
	.
	.
file[ 104 ]= sz24_up_angry_open_4 
	Calculating mask size...
	Frobenius Norm:  0.03357871 
file[ 105 ]= saavik_up_neutral_open_4 
	Calculating mask size...
	Frobenius Norm:  0.01904011 
	
Best Mask:
row col
[1,]   4   4
[2,]   3   5
[3,]   5   3
[4,]   6   2
[5,]   2   8	
```

At the end of the process, folder 'CMU' will have all results obtained by our FNN method. This folder will have another folder named as 'fnnMaskEstimate', which contains the average FNN matrix, histograms, this information for each channel color, and the five best masks resulted. Variable 'maskRes' will contain the dataset name ('maskRes$dataset'), the folder in which the results are saved ('maskRes$folder_initial'), the histogram matrix ('maskRes$bestMasks$histogram'), the sorted mask size list according to the histogram ('maskRes$bestMasks$maskList'), and the five best mask sizes with the counting number ('maskRes$bestMasks$maskKList').

### Number of Convolutional Units Estimation

Inside the R environment, we load the units evaluation file (unitsEstimation.r) that will read the images and apply function ImageFalseNearest on each individual image, while assessing the stop criterion. This evaluation is also applied on each color channel of every image under analysis. Units estimation requires the name of the dataset, a file containing the path of images, and the result obtained by code maskEstimation.r.

Steps for assessing the number of convolutional units:

```
source("unitsEstimation.r")
#Reading masks found
bMasks = as.matrix(read.table("./CMU/fnnMaskEstimate/bestMasks.dat"))
nkEst = NK_evaluation("CMU","cmu-list.txt", bMasks[,1], bMasks[,2])
```
This commands will take a bit longer to finish, and the results will be saved in the same folder created by maskEstimation.r. After reading file 'bestMasks,dat', variable 'bMasks' will contain the five best mask sizes:

```
     V1 V2
[1,]  4  4
[2,]  3  5
[3,]  5  3
[4,]  6  2
[5,]  2  8
```

During the units evaluation, output messages will appear, indicating the mask size and image names that are under analysis, as well as the error produced after it. Some part of output is presented below:

```
-------Compute FNN for number of units estimation.-------
mkdir: cannot create directory “CMU”: Permission denied
Analysing  440  files

Analysing color  GRAY  for mask  4 x 4 
	file[ 1 ]= kk49_left_angry_open_4 
		Calculating number of kernels...
		Euclidian Norm:  1.021914 
	file[ 2 ]= glickman_right_happy_open_4 
		Calculating number of kernels...
		Euclidian Norm:  0.005726963 
	file[ 3 ]= at33_up_neutral_sunglasses_4 
		Calculating number of kernels...
		Euclidian Norm:  0.005345669 
	file[ 4 ]= an2i_up_sad_sunglasses_4 
		Calculating number of kernels...
.
.
.
Analysing color  GRAY  for mask  3 x 5 
	file[ 1 ]= kk49_left_angry_open_4 
		Calculating number of kernels...
		Euclidian Norm:  1.208448 
	file[ 2 ]= glickman_right_happy_open_4 
		Calculating number of kernels...
		Euclidian Norm:  0.02535648 
	file[ 3 ]= at33_up_neutral_sunglasses_4 
		Calculating number of kernels...
		Euclidean Norm:  0.01223306 
	file[ 4 ]= an2i_up_sad_sunglasses_4 
		Calculating number of kernels...
.
.
.
Best Number of Units:
	Mask: 4 x 4 	Number of Units:  10 
	Mask: 3 x 5 	Number of Units:  10 
	Mask: 5 x 3 	Number of Units:  10 
	Mask: 6 x 2 	Number of Units:  10 
	Mask: 2 x 8 	Number of Units:  10 

```

At the end of the process, folder 'CMU' will have all results obtained by our FNN method. This folder will have another subfolder named as 'fnnNKernelEstimate', which contains the results for every five best mask estimated, including the FNN vector for each channel color and its average, and the best number of units. Variable 'nkEst' will contain the dataset name ('nkEst$dataset'), the folder in which results are saved ('nkEst$folder_initial'), and the best number of units for each mask size estimated by the maskEstimation.r.

## Training CNN to Validate the Designed Architecture

Before starting, the following instructions require Caffe deep learning framework [1] installed, and a brief knowledge about it. Instructions for installing Caffe are available at <http://caffe.berkeleyvision.org/installation.html>, and the basic tutorials at <http://caffe.berkeleyvision.org/gathered/examples/mnist.html>. In addition, this example requires the CMU Face Images dataset, which is available at <https://archive.ics.uci.edu/ml/datasets/CMU+Face+Images>.

Files for training the CNN architecture are available in folder 'caffe-experiments', that contains the example files. 

### Executing CNN training for CMU Face Images dataset

Our example topology is composed of one convolutional layer with ReLU activation function, followed by a sub-sampling max-pooling layer, a fully-connected layer, and a soft-max classifier. This fully-connected layer has $20$ units, which is the number of labels to classify the person in a given image. However, this dataset can be also classified according to other labels as pose (left, right, up, straight), expression (happy, sad, angry, neutral), and eyes (wearing sunglasses or not). It is worth to mention that we are working with images, while the Caffe tutorial work with LMDB databases.

The files used in this example are example_solver.prototxt and example_topology.prototxt. The solver file uses the SGD to train the CNN with a learning rate equals to 0.001, momentum 0.9, and regularization term 0.004 (see <http://caffe.berkeleyvision.org/tutorial/solver.html> for more information). The topology file has the CNN architecture information, in which the mask convolutional size is 4x4 with 10 units, and the max-pooling has size 3x3 with stride equals to 2x2. Having all files prepared, we can perform the CNN training:

```
./caffe-master/build/tools/caffe train --solver=./example_solver.prototxt
```

The trained models will be saved in folder 'models'. They permit to execute just a test in CNN, as follows:

```
./caffe-master/build/tools/caffe test -model example_topology.prototxt -weights models/cmu_example_train_iter_10000.caffemodel 
```

### References

[1] Jia, Yangqing, et al. "Caffe: Convolutional architecture for fast feature embedding." _Proceedings of the 22nd ACM international conference on Multimedia_. ACM, 2014.

### Citing this Estimate Parameters code 

Please cite this article in your publications if this code helps your research:

    @article{Ferreira2018,
    title = "Designing architectures of convolutional neural networks to solve practical problems",
    journal = "Expert Systems with Applications",
    volume = "94",
    number = "Supplement C",
    pages = "205 - 217",
    year = "2018",
    issn = "0957-4174",
    doi = "https://doi.org/10.1016/j.eswa.2017.10.052"
    author = "Martha Dais Ferreira and Débora Cristina Corrêa and Luis Gustavo Nonato and Rodrigo Fernandes de Mello"}

Or:

[2] Martha Dais Ferreira, Débora Cristina Corrêa, Luis Gustavo Nonato, Rodrigo Fernandes de Mello, "Designing architectures of convolutional neural networks to solve practical problems", _In Expert Systems with Applications_, Volume 94, 2018, Pages 205-217, ISSN 0957-4174, https://doi.org/10.1016/j.eswa.2017.10.052.

