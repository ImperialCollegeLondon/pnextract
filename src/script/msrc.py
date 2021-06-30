#!/usr/bin/env python3

# basic python utilities used in testing etc.
# set the PATHs needed to run porefoam... codes from python3
# developed by Ali Q. Raeini (2019), 



# python initializer for testing and launching porefoam and network modelling apps

import math
from math import exp, expm1, floor
import os,sys,copy,time, subprocess
from sys import stdin,stdout,stderr
import re
import io
import inspect



######################  SET SCRIPT DIRECTORIES    ######################

try:    disp(msEnv['msSrc'])
except: msEnv = os.environ.copy(); pass

msSrc = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # script directory
if "msSrc" in msEnv: 
	if msSrc!=msEnv['msSrc'] :disp("Info: not resetting msSrc: "+msEnv['msSrc'])
	msSrc= msEnv['msSrc'];
else :
	disp(("Setting msSrc:",msSrc))
	sep=os.pathsep;
	msRoot = os.path.abspath(msSrc+"/..")
	msEnv["PATH"] = msEnv["PATH"] +sep+ msSrc+"/script"
	msEnv["PATH"] = msEnv["PATH"] +sep+ msSrc+"/porefoam1f/script"
	msEnv["PATH"] = msEnv["PATH"] +sep+ msSrc+"/porefoam2f/script"
	msEnv["PATH"] = msEnv["PATH"] +sep+ msSrc+"/gnm/script/PNM"
	try:msEnv["LD_LIBRARY_PATH"] = msEnv["LD_LIBRARY_PATH"] +sep+msRoot+"/lib"
	except:msEnv["LD_LIBRARY_PATH"] = msRoot+"/lib"
	msEnv["PATH"] = msEnv["PATH"] +sep+ msRoot+"/bin"
	sys.path.append(msSrc+"/pylib")
	msEnv.update({'msSrc':msSrc})
	os.environ.update({'msSrc':msSrc})
	os.environ.update({"PYTHONPATH":msSrc+"/pylib"+sep+os.environ.get("PYTHONPATH","")});

######################  BASIC TEST UTILITIES  ##########################

# Linux commands !
import functools
disp  = functools.partial(print, flush=True)
mkdr = functools.partial(os.makedirs, exist_ok=True)
def pwd(): return os.getcwd()+'/'
def cd(workdir):  
	pwd_main = os.getcwd()
	if len(workdir) and workdir!='.' and pwd_main!=workdir:
		try: os.chdir(workdir)
		except OSError:  DbgMsg('Error cannot get into '+workdir+'/, from '+pwd_main) ;  exit(-1)
	return pwd_main


def DbgMsg(message='', nTrace=4, isError=0, endl='\n',nSkipTrace=1):
	stak = inspect.stack()
	ErTyp="Error in " if isError else ""
	if nTrace==0:  print(f'{ErTyp} {inspect.getframeinfo(stak[nSkipTrace][0]).function},  {message}'); return  -1
	if min(nTrace,len(stak))<nSkipTrace+2:
		er = inspect.getframeinfo(stak[nSkipTrace][0]) # 0 represents this line, 1 represents line at caller
		print(f'File "{er.filename}", line {er.lineno}, in {er.function}  {message}')
	else:
		for ii in range(0,min(nTrace,len(stak)-nSkipTrace)):
			er = inspect.getframeinfo(stak[ii+nSkipTrace][0]) # 0 represents this line, 1 represents line at caller
			print(f'File "{er.filename}", line {er.lineno}, in {er.function}')
		if nTrace<len(stak)-nSkipTrace: 
			er = inspect.getframeinfo(stak[-1][0]) # in main module
			print(f'File "{er.filename}", line {er.lineno}, in {er.function}')
		disp(message,endl)
	return isError

def alert(message='', abort=0, err=None, nTrace=1, nSkipTrace=2): 
	if err: disp(message);  raise Exception(message) from err 
	else: 
		DbgMsg( message, nTrace=nTrace, isError=True, nSkipTrace=nSkipTrace); 
		if abort:  exit(abort);

def ensure(okey, message='', abort=0, nTrace=1, nSkipTrace=2): 
	if not okey:  alert(message, abort, None, nTrace, nSkipTrace)

def runSh(resDir, script, logfile=None, envs={}):
	'''if script starts with space or len(script[0])<=3 then logfile=sys.stdout'''
	if not os.path.exists(resDir):   disp('error directory'+resDir+' not present'); exit(-1)
	oNam=str(script).replace('[','').split(' ')[0]; disp('>> '+str(script)) # 
	if not logfile and len(oNam)>3 and oNam[0].isalnum():  
		logfile = open(resDir+'/'+oNam+'.log','wb')
	myenv = msEnv.copy();	myenv.update(envs) 
	disp(f'Running {script} in {resDir}, envs: {envs or ""} > {oNam or ""}');
	return subprocess.run(script, stdout=logfile, stderr=logfile, shell=True, cwd=resDir, env=myenv)


def grepFloatInStr(lines='',keyword='Kx=', fnamHint=''):
	try:
		for ele in re.split(':|=| |,|\n|\t|;',lines.split(keyword,1)[1],20) :
			if len(ele):
				try: return float(ele) ;   
				except ValueError: 
					DbgMsg('no valid data for '+keyword+' in '+pwd()+fnamHint+', got '+ele);  pass ; 
					return float('NaN')
	except : DbgMsg(keyword+' not in file (?) '+pwd()+fnamHint+" ",6); pass
	return float('NaN')

def grepFloatInFile(inFIle='summary.txt',keyword='Kx='):
	try:     lines = open(inFIle, 'r').read()
	except: DbgMsg('cannot open '+inFIle) ; return 0.
	return grepFloatInStr(lines,keyword, inFIle)

def fileFloatDiffersFrom(inFIle, keyword, val, frac=0.01, delta=1e-32):
	fVal=grepFloatInFile(inFIle,keyword)
	if abs(fVal-val)<frac*abs(val)+delta : 
		print(inFIle+" -> "+keyword+": "+str(fVal) +" ~= "+str(val))
		return 0
	else                                  : 
		print(inFIle+" -> "+keyword+": "+str(fVal) +" != "+str(val))
		return 1
