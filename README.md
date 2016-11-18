# VBCCA: Variational Bayesian Canonical Correlation Analysis

Developed by Yusuke Fujiwara (email: yureisoul [at] gmail [dot] com), 2013/03/11.

This package provides a Matlab (object-oriented) implementation of Bayesian CCA.
Parameter estimation and prediction of Bayesian CCA are demonstrated using visual image reconstruction data from Miyawaki et al. 2008.

## Directories

- `vbBCCA`
    - Bayesian CCA source code.
- `script`
    - Examples; parameter estimation, visual image reconstruction and identification.
- `data`
    - Visual images and corresponding fMRI data for area V1.

## Files

### vbBCCA

- `BCCAtrainMain.m`
    - Estimate parameters for Bayesian CCA.
- `BCCApredOneWay.m`
    - Predict data1 from data2 or data2 from data1.
- `BCCApredBoth.m`
    - Predict data1 from data2 and data2 from data1.
- `vbBCCA.m`
    - Superclass of BCCAtrain and BCCApred. Interface of data input.
- `BCCAtrain.m`
    - Object for parameter estimation.
- `BCCApred.m`
    - Object for prediction.

### script

- `bcca_trainRandom_testFigure.m`
    - Visual image reconstruction.
- `bcca_Random_identification.m`
    - Visual image identification by comparing brain activity patterns.
- `setfigure.m`
    - Figure setting in bcca_trainRandom_testFigure.

### data/

- `V1_raw_random.mat`
    - Visual image and fMRI data of "random image session."
- `V1_mean_figure.mat`
    - Visual image and fMRI data of "figure image session."

## Citation


Fujiwara Y, Miyawaki Y and Kamitani Y. (2013). Modular encoding and decoding models derived from Bayesian Canonical Correlation Analysis. *Neural Computation* **25**, 979-1005. <http://www.mitpressjournals.org/doi/abs/10.1162/NECO_a_00423>

## References

- Miyawaki Y, Uchida H, et al., and Kamitani Y. (2008). Visual image reconstruction from human brain activity using a combination of multi scale local image decoders. Neuron, 60, 915-929.
- Fujiwara, Y., Miyawaki, Y., and Kamitani, Y. (2009). Estimating image bases for visual image reconstruction from human brain activity. Advances in neural information processing systems, 22, 576-584. (<http://books.nips.cc/papers/files/nips22/NIPS2009_0804.pdf>)

## License

Copyright (c) 2013, Yusuke Fujiwara
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
  in the documentation and/or other materials provided with the distribution.
- Neither the name of the Advanced Telecommunications Research Institute International nor the names of its contributors may be
  used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
