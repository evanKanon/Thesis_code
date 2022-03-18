#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import os
import subprocess
#------------------------------

outputdir='/mnt/d/PhD/Results/KFold/'
inputdir='/mnt/d/PhD/input/KFold/'

a=[1,2,3,4,5]

for mya in a:
	eons=subprocess.call('python /home/evangelos/eons/eons.py -hp -s '+outputdir+' -f '+inputdir+'GSE19429_train_'+str(mya)+'.expression -feat 25 -eval -prep -log -class '+inputdir+'GSE19429_train_'+str(mya)+'_disease.classes', shell=True)
	eons=subprocess.call('python /home/evangelos/eons/eons.py -hp -s '+outputdir+' -f '+inputdir+'GSE19429_train_'+str(mya)+'.expression -feat 25 -eval -prep -log -class '+inputdir+'GSE19429_train_'+str(mya)+'_specimen.classes', shell=True)





