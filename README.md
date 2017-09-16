# Brain2mesh: a one-liner for human brain 3D mesh generation

* Authors: Anh Phong Tran* <tran.anh@husky.neu.edu>, Qianqian Fang** <q.fang at neu.edu>
  * *Department of Chemical Engineering and **Department of Bioengineering
  * Northeastern University
  * 360 Huntington Ave, Boston, MA 02115
* Version: 0.5.0
* License: GPL v2 or later 
* URL: http://www.mcx.space/brain2mesh

## Introduction 

The Brain2Mesh toolbox provides a streamlined matlab function to convert a segmented brain 
volumes and surfaces into a high-quality multi-layered tetrahedral brain/full head mesh. 

This tool does not handle the segmentation of MRI scans, but examples of how commonly 
encountered segmented datasets can be used to create meshes are available under the "examples" folder.

The details of this toolbox can be found in the following paper:
Anh Phong Tran and Qianqian Fang, "Fast and high-quality tetrahedral mesh generation \
from neuroanatomical scans,". In: arXiv pre-print (August 2017). arXiv:1708.08954v1 [physics.med-ph]

The Brain2Mesh toolbox is also extensively dependent on:
1. Iso2Mesh toolbox (http://iso2mesh.sf.net)
   Not included - use the latest version at: https://github.com/fangq/iso2mesh
2. MATLAB Image-Processing toolbox (such as imfill)
3. intriangulation.m (by Adam Aitkenhead and Johannes Korsawe)
4. max_filter.m

## Overview of the functions

The function "brain2mesh.m" handles the conversion of segmented volumes into high-quality 3D meshes. 
It takes an 4D array as input, with different assumptions as to the number of layers. Typically, the layers
are assumed to contain: white matter (WM), grey matter (GM), cerebrospinal fluid (CSF), bone and scalp.
It also is able to handle inputs missing a bone segmentation, and missing bone+scalp segmentation. 

Patient-specific segmentations can be done using common neuroimaging tools such as FSL, SPM, 
FreeSurfer and BrainSuite. There also exists series of available atlas databases that offer segmented volumes.

In a near future release, scripts will be made available to accomodate the combination of segmented surface meshes
such as the ones produced in FreeSurfer and BrainSuite as part of the input data.

Your acknowledgement of brain2mesh in your publications or presentations 
would be greatly appreciated by the authors of this toolbox. The citation 
information can be found in the Introduction section.
