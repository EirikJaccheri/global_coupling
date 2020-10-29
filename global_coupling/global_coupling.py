import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import copy
import time

#work
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
#home
#sys.path.append("/home/eirik/CERN/beta-Beat/Python_Classes4MAD/")
try:
     from metaclass import *
except:
     from metaclass25 import *

     
madx_path = "no madx_path"
response_path = "no response_path"
C_min_path = "no C_min_path"
FineTuneCoupling_path = "no FineTuneCoupling_path"

reset_response_dict = {"unchanged_key" : "unchanged_value"}
reset_beta_error_dict = {"unchanged_key" : "unchanged_value"}



def set_global_paths(path_dict):
	path_key = "madx_path"
	if path_key in path_dict.keys():
		global madx_path
		error_string = "no file at: " + path_dict[path_key]
		#assert os.path.isfile(path_dict[path_key]), error_string
		madx_path = path_dict[path_key]
	
	path_key = "response_path"
	if path_key in path_dict.keys():
		global response_path
		error_string = "no file at: " + path_dict[path_key]
		assert os.path.isfile(path_dict[path_key]), error_string
		response_path = path_dict[path_key]
	
	path_key = "C_min_path"
	if path_key in path_dict.keys():
		global C_min_path
		error_string = "no file at: " + path_dict[path_key]
		assert os.path.isfile(path_dict[path_key]), error_string
		C_min_path = path_dict[path_key]
	
	path_key = "FineTuneCoupling_path"
	if path_key in path_dict.keys():	
		global FineTuneCoupling_path
		error_string = "no file at: " + path_dict[path_key]
		assert os.path.isfile(path_dict[path_key]), error_string
		FineTuneCoupling_path = path_dict[path_key]
	
	

def set_template(path,change_dict):
	with open(path,'r') as mask:
		template = mask.read()
	
	for key in change_dict:
		template = template.replace(key,change_dict[key])
	
		
	new_path = path.replace(".madx","_temp.madx")
	with open(new_path, "w") as job_file:
		job_file.write(template)
	
	error_str = "% in template " + new_path
	assert "%" not in template, error_str


def run_madx(path,change_dict):
	error_string = "no file at: " + path
	assert os.path.isfile(path), error_string
	path_temp = path.replace(".madx","_temp.madx")
	set_template(path,change_dict)
	os.system(madx_path + path_temp)


def get_twiss(file_path,twiss_path,change_dict):
	run_madx(file_path,change_dict)
	return twiss("output_files/" + twiss_path)

def get_f(file_path,twiss_path,change_dict):
	tw40cm = get_twiss(file_path,twiss_path,change_dict)
	tw40cm.Cmatrix()
	f_R, f_I = np.array(tw40cm.F1001R) , np.array(tw40cm.F1001I)
	f = f_R + 1j * f_I
	return f	
	
def change_value(change_dict,key,value):
	assert key in change_dict, "key not in dictionary"	
	change_dict[key] = value


def set_changes(change_dict,reset_dict):
	for key in reset_dict:
		change_value(change_dict,key,reset_dict[key])

def set_reset_dict(reset_dict,reset_dict_new):
	reset_dict.clear()
	for key in reset_dict_new:
		reset_dict[key] = reset_dict_new[key]
					

def get_responsematrix(change_dict):
	#set all errors to zero
	change_dict_local = copy.deepcopy(change_dict)
	assert "unchanged_key" not in reset_response_dict, "you havent set the reset_response_dict"
	set_changes(change_dict_local,reset_response_dict)
	
	delta_knob = 0.0001
	
	change_value(change_dict_local,"%knob_Re_value",str(delta_knob))
	change_value(change_dict_local,"%knob_Im_value","0.")
	run_madx(response_path,change_dict_local)
	tw40cm = twiss("output_files/twiss.original")
	tw40cm.Cmatrix()
	f_R_1 = np.array(tw40cm.F1001R)
	f_I_1 = np.array(tw40cm.F1001I)
	f_1 = np.concatenate((f_R_1,f_I_1))
	
	change_value(change_dict_local,"%knob_Re_value","0.")
	change_value(change_dict_local,"%knob_Im_value",str(delta_knob))
	run_madx(response_path,change_dict_local)
	tw40cm = twiss("output_files/twiss.original")
	tw40cm.Cmatrix()
	f_R_2 = np.array(tw40cm.F1001R)
	f_I_2 = np.array(tw40cm.F1001I)
	f_2 = np.concatenate((f_R_2,f_I_2))

	R = np.array([f_1/delta_knob,f_2/delta_knob]).T
	R_inverse = np.linalg.pinv(R)
	return R_inverse



def get_response_knobs(R_inverse,change_dict):
	change_dict_local = copy.deepcopy(change_dict)
	run_madx(response_path,change_dict_local)
	tw40cm = twiss("output_files/twiss.original")
	tw40cm.Cmatrix()
	f_R = np.array(tw40cm.F1001R)
	f_I = np.array(tw40cm.F1001I)
	f = np.concatenate((f_R,f_I))
	f_c = f_R + 1j * f_I
	knobs = -1 * np.dot(R_inverse,f)
	knob_Re, knob_Im = knobs[0], knobs[1]
	return knob_Re, knob_Im

def get_madx_knobs(change_dict):
	run_madx(FineTuneCoupling_path,change_dict)
	knob_file = open("output_files/knob.txt")
	lines = knob_file.readlines()
	knob_Re = float(lines[0].split("=")[1][0:-2:])
	knob_Im = float(lines[1].split("=")[1][0:-2:])
	return knob_Re, knob_Im


def get_C_min(change_dict):
	run_madx(C_min_path,change_dict)
	tw40cm = twiss("output_files/twiss.C_min")
	Q1, Q2 = tw40cm.Q1, tw40cm.Q2
	C_min = abs(Q1 - Q2 - 2)
	return C_min

def get_beta_error(change_dict):
	change_dict_local = copy.deepcopy(change_dict)
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm.Cmatrix()
	betax1, betay1 =  np.array(tw40cm.BETX), np.array(tw40cm.BETY)

	assert "unchanged_key" not in reset_beta_error_dict, "you havent set the reset_beta_error_dict"
	set_changes(change_dict_local,reset_beta_error_dict)
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm.Cmatrix()
	betax0 , betay0 = np.array(tw40cm.BETX), np.array(tw40cm.BETY)
	
	N = len(betax0)
	betax_rmsError = np.sum(((betax1 - betax0) / betax0)**2) / N
	betay_rmsError = np.sum(((betay1 - betay0) / betay0)**2) / N
	betax_error = (betax1 - betax0) / betax0
	betay_error = (betay1 - betay0) / betay0
	beta_max_error = max(np.max(betax_error),np.max(betay_error))
	beta_beat = (betax0 - betax1)/betax0
	return beta_max_error, beta_beat, betax0, betax1, betay0,betay1
	
	
def get_mean_strength(change_dict):
	change_dict_local = copy.deepcopy(change_dict)
	change_value(change_dict_local,"%twiss_pattern",".")
	run_madx(response_path,change_dict_local)
	tw40cm = twiss("output_files/twiss.original")	
	K1SL = np.array(tw40cm.K1SL)
	K1SL_nonzero = K1SL[K1SL != 0]
	return len(K1SL_nonzero), np.mean(np.abs(K1SL_nonzero))

def set_knobs(change_dict,knob_Re,knob_Im):
	change_value(change_dict,"%knob_Re_value",str(knob_Re))
	change_value(change_dict,"%knob_Im_value",str(knob_Im))		






