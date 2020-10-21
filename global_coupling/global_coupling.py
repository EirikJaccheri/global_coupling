import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import copy
import time

sys.path.append("/home/eirik/CERN/beta-Beat/Python_Classes4MAD/")
try:
     from metaclass import *
except:
     from metaclass25 import *

     
madx_path = "no madx_path"
response_path = "no response_path"
C_min_path = "no C_min_path"
FineTuneCoupling_path = "no FineTuneCoupling_path"

def set_global_paths(path_dict):
	global madx_path
	path_key = "madx_path"
	error_string = "no file at: " + path_dict[path_key]
	#assert os.path.isfile(path_dict[path_key]), error_string
	madx_path = path_dict[path_key]
	
	global response_path
	path_key = "response_path"
	error_string = "no file at: " + path_dict[path_key]
	assert os.path.isfile(path_dict[path_key]), error_string
	response_path = path_dict[path_key]
	
	global C_min_path
	path_key = "C_min_path"
	error_string = "no file at: " + path_dict[path_key]
	assert os.path.isfile(path_dict[path_key]), error_string
	C_min_path = path_dict[path_key]
	
	global FineTuneCoupling_path
	path_key = "FineTuneCoupling_path"
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

def change_value(change_dict,key,value):
	assert key in change_dict, "key not in dictionary"	
	change_dict[key] = value
	
def get_responsematrix(change_dict):
	#set all errors to zero
	change_dict_local = copy.deepcopy(change_dict)
	change_value(change_dict_local,"%error_strength","0.")
	change_value(change_dict_local,"%measured_C_min_comment","!")
	
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
	change_value(change_dict_local,"%knob_Re_value","0.")
	change_value(change_dict_local,"%knob_Im_value","0.")
	
	run_madx(response_path,change_dict_local)
	tw40cm = twiss("output_files/twiss.original")
	tw40cm.Cmatrix()
	f_R = np.array(tw40cm.F1001R)
	f_I = np.array(tw40cm.F1001I)
	f = np.concatenate((f_R,f_I))
	f_c = f_R + 1j * f_I
	knobs = -1 * np.dot(R_inverse,f)
	knob_Re, knob_Im = knobs[0], knobs[1]
	return tw40cm, knob_Re, knob_Im

def get_madx_knobs(change_dict):
	run_madx(FineTuneCoupling_path,change_dict)
	knob_file = open("knob.txt")
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
	
def get_mean_strength(change_dict):
	change_dict_local = copy.deepcopy(change_dict)
	change_value(change_dict_local,"%twiss_pattern",".")
	run_madx(response_path,change_dict_local)
	tw40cm = twiss("output_files/twiss.original")	
	K1SL = np.array(tw40cm.K1SL)
	K1SL_nonzero = K1SL[K1SL != 0]
	return len(K1SL_nonzero), np.mean(np.abs(K1SL_nonzero))

	
def plot_real_error_correction(change_dict,n_steps,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	measured_C_min_Re = float(change_dict_local["%measured_C_min_Re_value"])
	measured_C_min_Im = float(change_dict_local["%measured_C_min_Im_value"])	
	
	
	R_inverse = get_responsematrix(change_dict_local)
	run_madx(response_path,change_dict_local)
	tw40cm0 = twiss("output_files/twiss.original")
	tw40cm0.Cmatrix()
	f_Re = np.array(tw40cm0.F1001R)
	f_Im = np.array(tw40cm0.F1001I)
	f0 = f_Re + 1j * f_Im
	C_min_0 = get_C_min(change_dict_local)
	
	knob_Re_madx, knob_Im_madx = get_madx_knobs(change_dict_local)
	change_value(change_dict_local,"%knob_Re_value",str(knob_Re_madx))
	change_value(change_dict_local,"%knob_Im_value",str(knob_Im_madx))
	run_madx(response_path,change_dict_local)
	tw40cm1= twiss("output_files/twiss.original")
	tw40cm1.Cmatrix()
	S = np.array(tw40cm1.S)
	f_Re = np.array(tw40cm1.F1001R)
	f_Im = np.array(tw40cm1.F1001I)
	f1 = f_Re + 1j * f_Im
	C_min_1 = get_C_min(change_dict_local)
	
	
	change_value(change_dict_local,"%knob_Re_value","0.")
	change_value(change_dict_local,"%knob_Im_value","0.")
	change_value(change_dict_local,"%measured_C_min_Re_value",str(measured_C_min_Re))
	change_value(change_dict_local,"%measured_C_min_Im_value",str(measured_C_min_Im))
	n0, mean_strength0 = get_mean_strength(change_dict_local)
	
	
	change_value(change_dict_local,"%knob_Re_value",str(knob_Re_madx))
	change_value(change_dict_local,"%knob_Im_value",str(knob_Im_madx))
	change_value(change_dict_local,"%measured_C_min_Re_value","0.")
	change_value(change_dict_local,"%measured_C_min_Im_value","0.")
	n1, mean_strength1 = get_mean_strength(change_dict_local)
	
	
	fig = plt.figure()
	fig.suptitle("$C_-^0$ = " + str(C_min_0) + "	$C_-^1$ = " + str(C_min_1) + "\n" + "$<|K1SL1>_{before}$ = " + "{:.2e}".format(mean_strength0) + "$<|K1SL1>_{correction}$ = " + "{:.2e}".format(mean_strength1))
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(S,abs(f0),label = "before correction")
	ax1.plot(S,abs(f1),label = "corrected with sq knobs")
	ax1.legend()
	ax1.set_xlabel("S")
	ax1.set_ylabel("|F1001|")
	plt.savefig("plots/" + savepath)
	plt.show()

