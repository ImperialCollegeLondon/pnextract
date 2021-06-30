
import os; ''' ========== set up paths  =========== '''
if not ("msRoot" in os.environ): 
  print("try again after running:\nsource .../src/bashrc"); exit(-1);
from msrc  import *



nErrs=0


with open("spack_in.mhd", 'w') as f1:
	f1.write("""DimSize = 20 20 20 
				Offset =      0  0  0
				replaceRange 0 255 0   //< everything set to 0
				reset   dx 1
				//--- shape     centre   radius inside-value
				Paint sphere  0   0   0     12     1
				Paint sphere  0   20  0     12     1
				Paint sphere  0   20  20    12     1
				Paint sphere  0   0   20    12     1
				Paint sphere  20  0   0     12     1
				Paint sphere  20  20  0     12     1
				Paint sphere  20  20  20    12     1
				Paint sphere  20  0   20    12     1
				reset   dx 5e-6
				""");

runSh('.', "voxelImageProcess spack_in.mhd  spack.mhd");
runSh('.', "pnextract spack.mhd");
if fileFloatDiffersFrom("pnextract.log","nElems:",5.,0.01): nErrs+=1

exit(nErrs)
