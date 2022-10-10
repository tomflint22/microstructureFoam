# CFWeldFoam


Using the PF approach described in this paper:

https://www.nature.com/articles/s41524-021-00524-6




order to do things:

cp -r initial 0
blockMesh
setFields
TesselateFoam
microstructureFoam
