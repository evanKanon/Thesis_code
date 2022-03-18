#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import os
import argparse
import subprocess
#------------------------------
#- Options
parser=argparse.ArgumentParser(description='DiseaseNet: a disease association network')
parser.add_argument('-i', '--input', type=str, help='file path for input directory')
args=parser.parse_args()
#------------------------------

#- define inputfile path
inputdir=args.input

#- gunzip files which end with gz
for file in os.listdir(inputdir):
	s=os.path.abspath(inputdir)
	if file.endswith('.gz'):
		filename=s+'/'+file
		newfile='.'.join(file.split('.')[:-1])
		newname=s+'/'+newfile
		unzip=subprocess.call('gunzip -c '+filename+' > '+newname,shell=True)
		if unzip==0:
			pass
		else:
			print('Problem, unzipping failed: '+filename)

