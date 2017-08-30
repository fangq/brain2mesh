# Brain2Mesh - a one-liner for human brain 3D mesh generation

The Brain2Mesh toolbox provides a streamlined matlab function to convert a segmented MRI scan into a high-quality multi-layered tetrahedral brain/full head mesh. 

The Brain2Mesh toolbox is extensively dependent on the Iso2Mesh toolbox (http://iso2mesh.sf.net) with additional dependency on MATLAB Image-Processing toolbox (for imfill) and 3rd party open-source toolboxes: intriangulation.m (by Adam Aitkenhead and Johannes Korsawe) and max_filter.m.