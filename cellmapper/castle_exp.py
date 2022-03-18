#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# castle comparison
#------------------------------
import os
import subprocess
import pandas as pd
from slacker import Slacker
#from cellmapfunctions import *
#------------------------------

#- set Slack paramaters
user='@evan'
slack=Slacker("xoxp-56103715248-56114624145-154905509446-72c754c99fabcd5e9e1772648f3ee7c4")

#- set inputdir path
inputdir='/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/Castle_data/'

#- set outputfolder
outputfolder='/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/output/castle_comparison/'

#- save paths of directories included
folders=[x[0] for x in os.walk(inputdir)]

#- run comparison between datafiles
for folder in folders:
	if folder==inputdir:
		pass
	else:
		#- save full file paths of files for analysis from each folder in list so that half matrix can be run using them
		myfiles=[]
		for filename in os.listdir(folder):
			s=os.path.abspath(folder)
			if filename.endswith('stdt.csv'):
				myfiles.append(s+'/'+filename)
		#- run the comparison
		for i in range(len(myfiles)):
			x=max(range(len(myfiles)))
			while x!=i:
				#- using the files one way round
				mapfile=myfiles[x]
				queryfile=myfiles[i]
				#print(mapfile,queryfile)
				cmd='python /mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/scripts/new_version_draft/cellmapper.py -q '+queryfile+' -m '+mapfile+' -o '+outputfolder+' -log -slack -cor'
				comparison=subprocess.call(cmd,shell=True)
				if comparison==0:
					slack.chat.post_message(user,'Comparison between '+mapfile+' and '+queryfile+' succeeded, rejoice')
				else:
					slack.chat.post_message(user,'Comparison between '+mapfile+' and '+queryfile+' failed, check it!')
				#- using the files the other way round
				queryfile=myfiles[x]
				mapfile=myfiles[i]
				#print(mapfile,queryfile)
				cmd='python /mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/scripts/new_version_draft/cellmapper.py -q '+queryfile+' -m '+mapfile+' -o '+outputfolder+' -log -slack -cor'
				comparison=subprocess.call(cmd,shell=True)
				if comparison==0:
					slack.chat.post_message(user,'Comparison between '+mapfile+' and '+queryfile+' succeeded, rejoice')
				else:
					slack.chat.post_message(user,'Comparison between '+mapfile+' and '+queryfile+' failed, check it!')
				x+=(-1)
