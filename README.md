# Brain2mesh: a one-liner for 3D brain mesh generation

* Version: 0.8
* Project Maintainer: Qianqian Fang <q.fang at neu.edu>
* Contributors: See AUTHORS.txt for details
* Address:
  * Department of Bioengineering
  * Northeastern University
  * 360 Huntington Ave, Boston, MA 02115
* License: GPL v3 or later, see LICENSE.txt
* URL: http://mcx.space/brain2mesh

## Introduction 

The Brain2Mesh toolbox provides a streamlined matlab function to convert a segmented brain 
volumes and surfaces into a high-quality multi-layered tetrahedral brain/full head mesh. 

The details of this toolbox is described in the paper listed in the [Reference](#reference) section.

This tool does not handle the segmentation of MRI scans, but examples of how commonly 
encountered segmented datasets can be used to create meshes are available under the `examples` folder.

The Brain2Mesh toolbox is also extensively dependent on:
1. Iso2Mesh toolbox (http://iso2mesh.sf.net), not included, download at https://github.com/fangq/iso2mesh
2. MATLAB Image-Processing toolbox (such as `imfill`, `imdilate`)
3. `intriangulation.m` (by Adam Aitkenhead and Johannes Korsawe)

## Overview of the functions

The function `brain2mesh` handles the conversion of segmented volumes into high-quality 3D meshes. 
It takes an 4D array as input, with different assumptions as to the number of layers. Typically, the layers
are assumed to contain: white matter (WM), grey matter (GM), cerebrospinal fluid (CSF), bone and scalp.
It also is able to handle inputs missing a bone segmentation, and missing bone+scalp segmentation. 

Patient-specific segmentations can be done using common neuroimaging tools such as FSL, SPM, 
FreeSurfer and BrainSuite. There also exists series of available atlas databases that offer segmented volumes.

In a near future release, scripts will be made available to accomodate the combination of segmented surface meshes
such as the ones produced in FreeSurfer and BrainSuite as part of the input data.

Another function `brain1020` provides an automated interface to compute head landmarks (10-20/10-5 systems 
or user-customizable divisions). Users can either interactively select 5 initial landmarks (nasion, inion, 
left and right ear-lobes, i.e. LPA/RPA, and vertex, i.e. CZ), the function automatically computes all brain
landmarks on the scalp surface using user-specified density.

Your acknowledgement of Brain2Mesh in your publications or presentations 
would be greatly appreciated by the authors of this toolbox. The citation 
information can be found in the Introduction section.

## Reference 

If you use Brain2Mesh or Brain Mesh Library in your publications, the authors of this toolbox 
greatly appreciate if you can cite the below paper

* Anh Phong Tran†, Shijie Yan†, Qianqian Fang*, (2020) "Improving model-based fNIRS analysis using mesh-based anatomical and light-transport models," Neurophotonics, 7(1), 015008, URL: http://dx.doi.org/10.1117/1.NPh.7.1.015008

## Acknowledgement 

This project is funded by the National Institutes of Health (NIH) / National Institute of General 
Medical Sciences (NIGMS) under the grant number R01-GM114365, and NIH/NINID/NIBIB under the grant
number R01-EB026998.

###  Copyright disclaimer for `intriangulation.m`

Copyright (c) 2016, Johannes Korsawe
Copyright (c) 2013, Adam H. Aitkenhead
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution
* Neither the name of Volkswagen AG nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
* Neither the name of The Christie NHS Foundation Trust nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
