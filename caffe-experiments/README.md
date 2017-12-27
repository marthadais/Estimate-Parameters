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
# Training CNN with Estimated Parameters

This source code is used to train CNN architectures designed after our FNN method, using the Caffe deep learning framework [1]. 

## Getting Started

These instructions are a brief example of how to train the CNN architecture using the estimated parameters (see <ImageFalseNearest>). All source codes and instructions permit you to run the experiments on your local machine.

## Requirements

Before starting, it is necessary to have Caffe deep learning framework installed [1], and a brief knowledge about it:

1. Install the Caffe deep learning framework, according to the instructions available at <http://caffe.berkeleyvision.org/installation.html>.

2. Follow the tutorials available at <http://caffe.berkeleyvision.org/gathered/examples/mnist.html> to have a brief knowledge of its execution.

3. Before proceeding, it is required to download the CMU Face Images dataset (available at <https://archive.ics.uci.edu/ml/datasets/CMU+Face+Images>).

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

Please cite this article in your publications if it helps your research:

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
