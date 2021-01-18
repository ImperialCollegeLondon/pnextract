
import os, sys;  "msRoot" in os.environ or sys.exit("Retry after `source .../src/bashrc`");  from msrc import *; DbgMsg('*** main, ignore ***')


nErrs=0


with open("spack_in.mhd", 'w') as f1:
	f1.write("""DimSize = 20 20 20 
				Offset =      0  0  0
				replaceRange 0 255 0   // everything set to 0
				reset   dx 1
				//---------  shape    centre    radius inside-value
				shapeToVoxel sphere  0   0   0     12     1
				shapeToVoxel sphere  0   20  0     12     1
				shapeToVoxel sphere  0   20  20    12     1
				shapeToVoxel sphere  0   0   20    12     1
				shapeToVoxel sphere  20  0   0     12     1
				shapeToVoxel sphere  20  20  0     12     1
				shapeToVoxel sphere  20  20  20    12     1
				shapeToVoxel sphere  20  0   20    12     1
				reset   dx 5e-6
				""");

runBashScript(".", "voxelImageConvert spack_in.mhd  spack.mhd");
runBashScript(".", "pnextract spack.mhd");
if fileFloatDiffersFrom("pnextract.log","nElems:",5,0.01): nErrs+=1

exit(nErrs)
