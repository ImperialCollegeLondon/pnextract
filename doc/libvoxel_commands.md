# libvoxel commands

The following are commands will be processed by libvoxel when reading images.  They shall be given in after  ElementDataFile keyword in the .mhd files to avoid incompatibility with Fiji-is-ImageJ and Paraview.  Furthermore, any additional application specific keywords shall be provided after these since libvoxel stops processing the commands once it encounters any unrecognised keyword (i.e. other than comments and keywords below).  

The list of supported commands by libvxel can be obtained from voxelImageProcess app by running:    

		voxelImageProcess ?

Specific commands' usage can be shown by running:

		voxelImageProcess ? command


Note these are in additional to the mhd header data provided in the top of .mhd files:

		ObjectType =  Image
		NDims =       3
		ElementType = MET_UCHAR
		ElementByteOrderMSB = False

		DimSize =    	1000	1000	1000
		ElementSize = 	1.6 	1.6 	1.6
		Offset =      	0   	0   	0

		ElementDataFile = Image.raw.gz

		//HeaderSize = 0 
		//CompressedData = True
		//OutputFormat = .raw.gz // alias for  DefaultImageFormat, not part of MHD specification

Additional MHD rules:

- The order of the first 6 keywords should not be changed for compatibility with third-party software (ImageJ and Paraview).    
- Use `#` or  `//` for commenting and deactivating commads/lines.    
- Each keyword and its data should be given in a single line.   This is because libvoxel currently does not rely on the more user friendly InputFile.h syntax, used by pnextract and pnflow ... codes.  
- Libvoxel detects only basic data/keywords provided in Image.am files and Image.tif files, without the need for a .mhd file. However, the commands bellow can only be provided in a separate .mhd (text) file, in which the name of Image.tif or Image.am file is given as an argument to ElementDataFile keyword (see above).


------------------------------------

## Basic commands


**info**:    
    Show image information/stats
 
		info


**pore**:    
    Pore/selectPore  thresholds image and sets the given range to zero (pore voxel values)  
 
		pore	  0  150

**threshold**:    
    Same as  selectPore   
 
		threshold  0  150

**threshold101**:    
    Same as  selectPore,  101 means inside the range will be set to zero and outside to 1
 
		threshold101  0  150


**direction**:    
    Re-interprets the image direction so that the flow direction becomes in the given  direction

		direction z


**crop**:    
    Crops the image, from a lower I, J, K and before an upper bound I J K values
 
		crop   100 100 0   900 900 1000

**resampleMean**:    
    Resamples the image while setting the new values to the mean of original voxel values in case of refining (ratio<1)   

		resampleMean   2

**resampleMax**:    
    Resample the image while setting the new values to the max of original voxels in case of refining (ratio<1)   
 
		resampleMax   2

**resampleMode**:    
    Resample the image while setting the new values to the mode of original voxels in case of refining (ratio<1)  
 
		resampleMode  2

**resliceZ**:    
    Same as resample but only in z direction (for correcting CT images, limited)  
 
		resliceZ

**replaceRange**:    
    Sets range [min max] to a given value   
 
		replaceRange  200 255      1


**keepLargest0**:    
    Keeps the largest connected region with value of 0, sets the rest to maxT-1 (=254 for 8bit images). 
    Use this to compute effective porosity.

**write**:    
    Write the image into file
 
		write   image_copy.raw
		write   image_copy.am
		write   image_copy.mhd

**write8bit**:    
    Write the image as unsigned char    
 
		write8bit image_8bit.raw 
		write8bit image_8bit.am 
		write8bit image_8bit.mhd 

**modeFilter**:    
    Apply  few iterations (first argument) of mode filter, reject mode if frequency less than  second argument   
 
		modeFilter      1       2    // 1 iterations, 2 for the threshold mode frequency to change a voxel

**medianFilter**:    
    Apply a median filter  
 
		medianFilter   1

**medianX**:    
    Apply few iterations of  median filter  only in x direction ( can help reduce file size of grey images when compressed) 
 
		medianX       1

**FaceMedian06**:    
    Median filter over 6 closest adjacent voxels
 
		FaceMedian06   1 

**PointMedian032**:    
    Median filter over 32 closest adjacent voxels
    
		PointMedian032   1  //<- one iterations

**mapFrom**:    
    Assigns voxel values from another image based on their physical location (X0+ijk*dX).  Image offsets should be assigned correctly.
 
		mapFrom "image_name.mhd"   ;
		mapFrom "image_name.mhd"  0 65536  0  0  ; //<- \[vv_min vv_max\] increase-factor shift  


**operation**:    
    Apply a mathematical unary (!, &, ~) or binary (+, -, \*, b, e).  Here b stands for begin or min and e stands for maximum or end for voxel values.  The second argument shall be empty in case of unary operations, or a number or an image file name in case of binary operators.     
 
		operation   + 10               //<- add 10 to voxel values
		operation   + image2.mhd -50   //<- -50 is a shift value to avoid overfow/capping
		operation   b 1                //<- set voxel values below <10 to 10
		operation   e 250              //<- set voxel values above 250> to 250

**operat**:    
    operat is an alias for `operation` keyword


------------------------------------

## Keywords used for tests/synthetic image generation


**read**:    
    Read an image, overwriting any previously loaded data
 
		read  input_image.tif 
		read  input_image.am 
		read  input_image.mhd 

**reset**:    
    Reset N ( image size) dx (voxel size) or X0 (offset/origin)
 
		reset	  N   100  100  100
		reset	  dx  3e-6 3e-6 3e-6
		reset	  X0  0.  0.  0.
		reset	  Nd0 100  100  100   3e-6 3e-6 3e-6   0. 0. 0.
		reset	  dx 3e-6

**Paint**:    
    Paint a 3D shape into the image    
 
		Paint     s   30 30 30  10    //<- shape-type(sphere)  Centre:Xc(30,30,30) and radius R=10 

**PaintAdd**:    
    Paints over a 3D shape on the image, increasing the previous voxel values (brightness), otherwise same as Paint    
 
		PaintAdd  s  30 30 30  10    //<- shape-type(sphere)  Centre:Xc(30,30,30) and radius R=10


------------------------------------

## Obsolete keywords

The following keywords fate in further libvoxel versions shall be reviewed 

**cropD**:    
    Same as crop, D is a reminder for \[\) meaning beginning of range is included but end is not    

		cropD  100 100 0   900 900 1000

**rescale**:    
    Adjust voxel values brightness /contrast (redundant, use `operat` keyword instead)
 
		rescale  0  10000

**Offset**:    
    Change the image origin 
 
		Offset 0 0 0

**delense032**:    
    Primitive lens artefact remover, obsolete(?)
 
		delense032 

**circleOut**:    
    Set voxels outside a circles (cylinder) of radius R centrrd at Xc Yc    

		circleOut ;   //<- centre set to centre of image and radius to half of image (width+height)/2
		circleOut  500 500   450 ;   //<- centre_X=500 Y=500, radius=450 voxels

**growLabel**:    
    Propagate a voxel value iteratively to its adjacent ones, obsolete
 
		growLabel             0         5 

**fillHoles**:    
    Filter out small features (obsolete feature)    
 
		fillHoles  2


**readAtZ**:    
    Read a new image and assign it at the beginning of the given slice number  (obsolete, see `mapFrom` instead)
   
		readAtZ   1150

**maskWriteFraction**:    
    writes fraction of voxels of two images having the same value (TODO remove?)
 
		maskWriteFraction  "maskname.mhd"  "outName.txt"   0    2  1000000
