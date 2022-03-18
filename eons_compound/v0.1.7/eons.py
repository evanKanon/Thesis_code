#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

#EONS master script 
#Coordinate the running of all other scripts

#------------------------------
#v0.0.0 initial version of the new code, adapted from EONS_main version 0.2
#v0.0.1 added log functionality, improved serv function, added error messages and Slack notification
#v0.0.2 added config file, changed EONS thresh defaults, redefined outputdirs, made single EONS_funct and deleted unused functions (version descr now includes changes in any/all scripts)
#v0.0.3 added extreme network selection functionality and the the relevant option for the distance multiplier selection, removed config hardpath and renamed EONS_main to EONS
#v0.0.4 changed MCL feature selection with PCA based on networktools code
#v0.0.5 changed feature selection to only principal components, deleted unused code
#v0.0.6 merged MCL feature selection and PCA, created multiple feature subsets through different inflation values, made it possible to notify many users simulataneously, added option to run on normalised data as well, automated static deletion of IDENTIFIER column of data, simplified EONS_main code
#v0.0.7 created prep script, changed file names from capital to lowercase, changed networks so that all networks are similarity and not distance networks
#v0.0.8 gave two options for feature selection, either MCL or PCA-based feature selection, added classfile option and config path as well as working and useful eval script
#v0.1.0 BASELINE EONS added failed network functionality to eons_02 making all distances runable, corrected multiplier functionality in eons_03, made eons_eval save proper output, used similarity and not distance networks in eons_funct i.e. 0 is least similar
#v0.1.1 changed make_expression function to accomodate file names other than .soft, added NA dealing option to eons_prep, improved config option for eons_main, made code python 3 compatible
#v0.1.2 established perpendicular point plots
#v0.1.3 improved folder names, added .layout visualisation, updated environment file, deleted unecessary options of eons_eval and turned all p-p networks to similarity networks
#v0.1.4 made running the MCL component faster and allowed user to select inflation values
#v0.1.5 changed naming of folders and certain default values
#v0.1.6 corrected PCA random selection in eons_01, corrected network names in eons_02 and deleted unused inputfile option in eons_03 with the appropriate corrections in eons_funct and _main
#v0.1.7 added correct evaluation measure

#------------------------------
__version__ = '0.1.7'
#------------------------------
import argparse, subprocess, datetime, time, os, sys, json
from slacker import Slacker
#------------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#------------------------------

#- Options
#- General functionality options
parser=argparse.ArgumentParser(description='Exhaustive Observation of Network Space Analysis')
parser.add_argument('-s', '--save', type=str, help='file path for results and figures to be saved', default='0')
parser.add_argument('-f', '--file', type=str, help='file path for input file')
parser.add_argument('-class', '--classfile', type=str, help='file path for class file, i.e. file with desired outcome', default='0')
parser.add_argument('-thr', '--thresholded', default=False, help='run a thresholded version of EONS', action='store_true')
parser.add_argument('-log', '--log_file', default=False, action='store_true', help='create log file containing the argunents used, the date run, the version and the input file name')
parser.add_argument('-slack', '--slack', default=False, action='store_true', help='send messages through Slack platform')
parser.add_argument('-mu', '--mult_users', default=False, action='store_true', help='send Slack messages to multiple users')
parser.add_argument('-eval', '--evaluation', default=False, help='run evaluation of networks based on how well they separate subgroups', action='store_true')
parser.add_argument('-prep', '--preparation', default=False, help='run default script to reformat input data into appropriate input', action='store_true')
#- developed this option so that eons can be run with hard paths provided at input disregarding i/o config filepaths
parser.add_argument('-hp', '--use_hard', default=False, action='store_true', help='use custom config file instead of hardpath for file paths, NOTE all other paths in config will be used')
parser.add_argument('-conf', '--configfile', default=os.path.join(scriptpath,'config.json'), help='custom config file')
#- EONS prep options
parser.add_argument('-na', '--na_comp', choices=[0,1], default=0, type=float, help='choose whether nas will be dropped or imputed to the mean')
#- EONS I options
parser.add_argument('-mcl', '--mcl_selection', default=False, help='use MCL as a way to subset groups', action='store_true')
parser.add_argument('-infl', '--infl_values', nargs='*', type=float, default=1.6, help='set list of inflation values to use for MCL feature selection method')
parser.add_argument('-feat','--features', default=19131261, type=int, help='set number of features for PCA')
parser.add_argument('-chs', '--choose_threshold', default=False, action='store_true', help='use mcxquery to find threshold for the creation of the initial gene-gene network')
parser.add_argument('-norm', '--normalise', action='store_true', default=False, help='normalise input and use that as additional data input')
#- EONS II options 
parser.add_argument('-full', '--full_analysis', default=False, help='run analysis with all corelation types possible', action='store_true')
parser.add_argument('-q', '--quick', default=False, help='quick run of analysis with restricted predefined parameters', action='store_true')
#- EONS III options
parser.add_argument('-abs', '--absolute', default=1, type=float, help='use absolute value for network retention i.e. if cluster contains more than 10 networks keep it')
parser.add_argument('-rel', '--relative', default=0, type=float, help='use relative value for network retention i.e. if cluster contains more than 10 per cent of all networks keep it')
parser.add_argument('-mult', '--multiplier', nargs='*', type=float, help='multiplier value for the distance threshold of the extreme network selection algorithm')
args=parser.parse_args()

#- open config file
with open(args.configfile) as json_data_file:
	config=json.load(json_data_file)

#- SET SLACK PARAMETERS
slack_arg=args.slack
slack=Slacker(config['Slack_env']['address'])
#- make allowance for more than one user to be notified
multiple_users=args.mult_users
#- set the main user who will be notified anyway if slack option is on
main_user=config['Slack_env']['user1']
#- set list of all users to be notified
users=[]
if multiple_users and slack_arg:
	for myuser in config['Slack_env'].keys():
		if myuser.startswith('user'):
			users.append(config['Slack_env'][myuser])
else:
	users.append(main_user)

#- SET MAIN PATHS

#- define inputfile path
if args.use_hard==False:
	inputfile=str(config['File_paths']['input'])+args.file
else:
	try:
		inputfile=args.file
	except:
		print('Error: inputfile needed !')
		#! future improvement
		#raise error : inputfile required
#- define outputdir path
if args.use_hard==False:
	if args.save=='0':
		outputdir=str(config['File_paths']['output'])
	else:
		outputdir=str(config['File_paths']['output'])+args.save
else:
	outputdir=args.save
#- define classfile path
if args.use_hard==False:
	if args.classfile=='0':
		classfile=str(config['File_paths']['classfile'])
	else:
		classfile=str(config['File_paths']['classfile'])+args.classfile
else:
	classfile=args.classfile

#- define filepath for scripts (this is done so that main can be run from anywhere)
filepath=scriptpath+'/'


#- Processing arguments for interprocess functionality
absolute=args.absolute
relative=args.relative

thr=1
if args.thresholded==True:
	thr=0

#- select number of distance metrics to be used
my_opt=None
if args.full_analysis==True:
	my_opt=1
	name_add='_extensive'
elif args.quick==True:
	my_opt=2
	name_add='_quick'
else:
	my_opt=3
	name_add='_default'

#- use mcxquery for feature network threshold selection
choose=0
if args.choose_threshold==True:
	choose=1

#- use MCL for feature selection
mcl=0
name_add_m='_pca'
if args.mcl_selection==True:
	mcl=1
	name_add_m='_mcl'

#- Set number of features for PCA
fet=args.features

#- Activate evaluation
evaluat=args.evaluation

#- Perform normalisation on data
norma=0
if args.normalise==True:
	norma=1

#- set whether missing data will be imputed or dropped (imputation is only for mean at present)
nas=0
if args.na_comp==1:
	nas=1

#- Create overarching output folder
daterun='_'+str(time.strftime("%d%m%Y"))
#analysislabel = os.path.split(inputfile)[1].replace('.expression','')
analysislabel = os.path.split(inputfile)[1].split('.')[0]
if evaluat==True:
	e_end='_'+classfile.split('/')[-1].split('.')[0]
else:
	e_end=''
outputdir = os.path.join(outputdir,analysislabel+daterun+'_EONSv'+str(__version__))+name_add+name_add_m+e_end+'/'

#- Optionally create log file for run
logfile=args.log_file

if logfile==True:
	if not os.path.isdir(outputdir):
		try:
			os.mkdir(outputdir)
		except:
			updir = os.path.split(outputdir)[0]
			os.mkdir(updir)
			os.system("/bin/chmod 755 %s"%(updir))
			os.mkdir(this_dir)
	if outputdir.endswith('/'):
		log=open(outputdir+'logfile.txt','w')
	else:
		log=open(outputdir+'/logfile.txt','w')
	log.write("Date run: "+time.strftime("%d%m%Y")+'\n')
	log.write("On datafile: "+inputfile+'\n')
	log.write('Using EONS version: '+__version__+'\n')
	log.write('With arguments:\n')
	arglist=args
	#print(sys.argv)
	#print(args)
	log.write(str(arglist)+'\n')
	log.write('\n')
	log.write('Full command used was:\n')
	log.write("%s\n"%(' '.join([str(x) for x in sys.argv])))
	log.close()
else:
	pass

#- Set extreme network multiplier/s
m=args.multiplier
if not m:
	m=[1]

multi=" ".join(str(x) for x in m)

#- Set MCL inflation values for EONS_01
inf_v=args.infl_values

#- Inform user that process commences
if slack_arg==True:
	for user in users:
		slack.chat.post_message(user, 'Commencing EONS analysis with '+str(inputfile))
else:
	print('Commencing EONS analysis with '+str(inputfile))

#- Optionally, run preparation script
if args.preparation:
	prep1=subprocess.call('python '+filepath+'eons_prep.py -i'+inputfile+' -o '+outputdir+' -na '+str(nas), shell=True)
	if prep1 == 0:
		if slack_arg==True:
			for user in users:
				slack.chat.post_message(user, '...Preparation of '+inputfile+' completed...')
		else:
			print('...Preparation of '+inputfile+' completed...')
	else:
		if slack_arg==True:
			for user in users:
				slack.chat.post_message(user, '...Preparation of '+inputfile+' failed...')
		else:
			print('...Preparation of '+inputfile+' failed...')
	inputfile=os.path.splitext(inputfile)[0]+'_prepped.expression'
else:
	pass

#- Initiate Phase I
eons1=subprocess.call('python '+filepath+'eons_01.py -s '+outputdir+' -f '+inputfile+' -feat '+str(fet)+' -chs '+str(choose)+' -norm '+str(norma)+' -mcl '+str(mcl)+' -infl '+str(inf_v), shell=True)
if eons1 == 0:
	if slack_arg==True:
		for user in users:
			slack.chat.post_message(user, 'EONS_01 ran, continuing with EONS_02...')
	else:
		print('EONS_01 ran, continuing with EONS_02...')
	#- Initiate Phase II
	eons2=subprocess.call('python '+filepath+'eons_02.py -s '+outputdir+' -f '+inputfile+' -a '+str(my_opt)+' -thr '+str(thr), shell=True)
	if eons2 == 0:
		if slack_arg==True:
			for user in users:
				slack.chat.post_message(user, '...EONS_02 ran, continuing with EONS_03...')
		else:
			print('...EONS_02 ran, continuing with EONS_03...')
		#- Initiate Phase III
		eons3=subprocess.call('python '+filepath+'eons_03.py -s '+outputdir+' -abs '+str(absolute)+' -rel '+str(relative)+' -thr '+str(thr)+' -m '+multi, shell=True)
		if eons3 == 0:
			if slack_arg==True:
				for user in users:
					slack.chat.post_message(user,'......EONS_03 ran, continuing with evaluation if selected, otherwise, ANALYSIS COMPLETE')
			else:
				print('......EONS_03 ran, continuing with evaluation if selected, otherwise, ANALYSIS COMPLETE')
			#- Initiate Evaluation of Networks
			if evaluat==True:
				evaluation=subprocess.call('python '+filepath+'eons_eval.py -s '+outputdir+' -class '+classfile+' -thr '+str(thr), shell=True)
				if evaluation==0:
					if slack_arg==True:
						for user in users:
							slack.chat.post_message(user,'EONS evaluation ran, ANALYSIS COMPLETE')
					else:
						print('EONS evaluation ran, ANALYSIS COMPLETE')
				else:
					if slack_arg==True:
						for user in users:
							slack.chat.post_message(user,'EONS evaluation failed, check problem manually')
					else:
						print('EONS evaluation failed, check problem manually')
			else:
				pass
		else:
			if slack_arg==True:
				for user in users:
					slack.chat.post_message(user, 'EONS_03 failed, check problem manually')
			else:
				print('EONS_03 failed, check problem manually')
	else:
		if slack_arg==True:
			for user in users:
				slack.chat.post_message(user, 'EONS_02 failed, check problem manually')
		else:
			print('EONS_02 failed, check problem manually')
else:
	if slack_arg==True:
		for user in users:
			slack.chat.post_message(user, 'EONS_01 failed, check problem manually')
	else:
		print('EONS_01 failed, check problem manually')
