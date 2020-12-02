import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import copy
import time

#work
#sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
#home
sys.path.append("/home/eirik/CERN/beta-Beat/Python_Classes4MAD/")
try:
     from metaclass import *
except:
     from metaclass25 import *

     
madx_path = "no madx_path"
response_path = "no response_path"
C_min_path = "no C_min_path"
FineTuneCoupling_path = "no FineTuneCoupling_path"

reset_response_dict = {"unchanged_key" : "unchanged_value"}
reset_analyticresponse_dict = {"unchanged_key" : "unchanged_value"}
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
	tw40cm = twiss("output_files/" + twiss_path)
	tw40cm.Cmatrix()
	return tw40cm

def get_f(file_path,twiss_path,change_dict):
	tw40cm = get_twiss(file_path,twiss_path,change_dict)
	tw40cm.Cmatrix()
	f_R, f_I = np.array(tw40cm.F1001R) , np.array(tw40cm.F1001I)
	f = f_R + 1j * f_I
	return f
	
def get_S_f(file_path,twiss_path,change_dict):
	tw40cm = get_twiss(file_path,twiss_path,change_dict)
	S = np.array(tw40cm.S)
	f = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	return S , f
	
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

	R = np.array(
	[f_1/delta_knob,f_2/delta_knob]).T
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
	
def get_beta_error_rms(change_dict):
	beta_max_error, beta_beat, betax0, betax1, betay0,betay1 = get_beta_error(change_dict)
	beta_beatx = (betax0 - betax1)/betax0
	beta_beaty = (betay0 - betay1)/betay0
	betax_rms = np.sqrt(np.mean(np.square(beta_beatx)))
	betay_rms = np.sqrt(np.mean(np.square(beta_beaty)))
	return betax_rms , betay_rms , beta_beatx, beta_beaty
	
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
	
	
#analytical responsematrix code	



def deltaPhi(mu_BPM,mu_error,Q):
	condlist = np.array([mu_BPM >= mu_error, mu_BPM < mu_error])
	choicelist = np.array([mu_BPM - mu_error, mu_BPM - mu_error + Q])
	return np.select(condlist,choicelist)

def B_components(mux_BPM,muy_BPM,Qx,Qy,betax_error,betay_error,mux_error,muy_error):
	n_BPM = len(mux_BPM)
	n_error = len(mux_error)

	mux_error_grid, mux_BPM_grid = np.meshgrid(mux_error,mux_BPM)	
	muy_error_grid, muy_BPM_grid,  = np.meshgrid(muy_error,muy_BPM)
	
	C =  1 / (4*(1 - np.exp(2*np.pi*1j*(Qx - Qy))))
	B_1 = np.tensordot(np.ones(n_BPM),np.sqrt(betax_error*betay_error),axes = 0)
	B_2 = np.exp(2*np.pi*1j*(deltaPhi(mux_BPM_grid,mux_error_grid,Qx) - deltaPhi(muy_BPM_grid,muy_error_grid,Qy)))
	return C , B_1, B_2
def B_matrix(mux_BPM,muy_BPM,Qx,Qy,betax_error,betay_error,mux_error,muy_error):
	C , B_1, B_2 = B_components(mux_BPM,muy_BPM,Qx,Qy,betax_error,betay_error,mux_error,muy_error)
	B = C * B_1 * B_2
	B.real *= -1 #multiplying by -1 since there are different sign connventions in franchi and madx
	return B


def f_1001(mux_BPM,muy_BPM,Qx,Qy,betax_error,betay_error,mux_error,muy_error,KS):
	B = B_matrix(mux_BPM,muy_BPM,Qx,Qy,betax_error,betay_error,mux_error,muy_error)
	return np.dot(B,KS)

	
def get_knob_matrix(change_dict):
	change_dict_local = copy.deepcopy(change_dict)
	change_value(change_dict_local,"%twiss_pattern",".")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	name_l = tw40cm.NAME
	

	N_KQS = 12
	knob_Re_name = change_dict_local["%knob_Re_type"]
	knob_Im_name = change_dict_local["%knob_Im_type"]
	
	with open("/home/eirik/CERN/global_coupling_correction/analytical_test/lhc/lhc_as-built.seq",'r') as read_obj:
		data_lhc = read_obj.readlines()
	
	
	KQS_matrix = []
	KQS_index_l = []
	KQS_name_l = []
	with open("/home/eirik/CERN/global_coupling_correction/global_coupling/knob_matrix.txt",'r') as read_obj:
		for line1 in read_obj:
			line1 = line1.replace('(','')
			line1 = line1.replace(')','')
			line1_l = line1.split(' ')
			KQS = line1_l[0]
			Knob_Re_index = line1_l.index(knob_Re_name) - 2
			Knob_Im_index = line1_l.index(knob_Im_name) - 2
			for line2 in data_lhc:
				if line2.__contains__(KQS.lower()):
					line2_l = line2.split(' ')
					QS = line2_l[2]
					QS = QS.replace(',','')
					KQS_index_l.append(name_l.index(QS))
					KQS_matrix.append([float(line1_l[Knob_Re_index]), float(line1_l[Knob_Im_index])])
					KQS_name_l.append(QS)
	return np.array(KQS_matrix), np.array(KQS_index_l),KQS_name_l,tw40cm

	
def get_responsematrix_components(change_dict,B_scaling = 1):
	change_dict_local = copy.deepcopy(change_dict)
	
	assert "unchanged_key" not in reset_analyticresponse_dict, "you havent set the reset_analyticresponse_dict"
	set_changes(change_dict_local,reset_analyticresponse_dict)
	

	KQS_matrix, KQS_index_l, KQS_name_l, tw40cm = get_knob_matrix(change_dict)[0:4]
	
	
	change_value(change_dict_local,"%twiss_pattern","BPM")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm.Cmatrix()
	f1001_0 = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	betx_BPM = np.array(tw40cm.BETX)
	bety_BPM = np.array(tw40cm.BETY)
	mux_BPM = np.array(tw40cm.MUX)
	muy_BPM = np.array(tw40cm.MUY)
	
	
	change_value(change_dict_local,"%twiss_pattern",".")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	name_l = tw40cm.NAME
	betx = np.array(tw40cm.BETX)
	bety = np.array(tw40cm.BETY)
	mux = np.array(tw40cm.MUX)
	muy = np.array(tw40cm.MUY)
	Qx = tw40cm.Q1
	Qy = tw40cm.Q2
	

	betx_error = np.take(betx,KQS_index_l)
 	bety_error = np.take(bety,KQS_index_l)
	mux_error = np.take(mux,KQS_index_l)
	muy_error = np.take(muy,KQS_index_l)
	
	l = 0.32
	C , B_1 , B_2 = B_components(mux_BPM,muy_BPM,Qx,Qy,betx_error,bety_error,mux_error,muy_error)
	B = B_scaling * B_matrix(mux_BPM,muy_BPM,Qx,Qy,betx_error,bety_error,mux_error,muy_error)
	Z = np.dot(B,KQS_matrix) * l
	Z = np.dot(B,KQS_matrix) * l
	return np.real(Z), np.imag(Z), B , C , B_1 , B_2 , KQS_matrix , KQS_index_l , KQS_name_l , tw40cm


def get_analytic_responsematrix(change_dict,B_scaling = 1):
	Z_Re, Z_Im = get_responsematrix_components(change_dict,B_scaling)[0:2]
	
	R = np.block([[Z_Re],[Z_Im]])
	R_inverse = np.linalg.pinv(R)
	return R_inverse
	


