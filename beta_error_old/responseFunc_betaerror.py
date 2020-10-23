import sys
import os
import numpy as np
import cmath
import matplotlib.pyplot as plt
import time
import copy
from mpl_toolkits import mplot3d


def qmean(num):
    return sqrt(sum(n*n for n in num)/len(num))

#sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
sys.path.append("/home/eirik/CERN/beta-Beat/Python_Classes4MAD/")
try:
     from metaclass import *
except:
     from metaclass25 import *

N_error = 12
madx_path = "/home/eirik/madx"

def delQ_min_new(Q1,Q2,f_R,f_I,S,MUX,MUY):
	N = len(f_I)
	delta = np.abs(Q1%1 - Q2%1)
	f = f_R + 1j * f_I

	deltaS_list = [np.abs(S[i +1] - S[i]) for i in range(N-1)]
	deltaS_list.append(S[-1] - S[-2])
	deltaS = np.array(deltaS_list)
	return 4 * delta / S[-1] * abs(np.sum(deltaS * f * np.exp(-1j * 2*np.pi*(MUY - MUX))))


def delQ_min_old(Q1,Q2,f_R,f_I,S,MUX,MUY):
	N = len(f_I)
	delta = np.abs(Q1%1 - Q2%1)
	f = f_R + 1j * f_I
	return 4 * delta/N * np.sum(np.abs(f))


def set_template(path,knob_Re, knob_Im, error):
	replace_dict = {}
	replace_dict["%knob_Re"] = knob_Re
	replace_dict["%knob_Im"] = knob_Im
	
	replace_dict["%error_component"] = error[1]
	replace_dict["%error_strength"] = error[2]
	replace_dict["%pattern_1"] = error[3]
	replace_dict["%pattern_2"] = error[4]
	replace_dict["%quad_component"] = error[5]
	replace_dict["%quad_pattern_1"] = error[6]
	replace_dict["%quad_strength"] = error[7]	
	
	for i in range(N_error):
		replace_dict["%Var" + str(i)] = error[0][i]

	with open(path,'r') as mask:
		template = mask.read()
	
	squeeze = error[8]
	if squeeze:
		template = template.replace("%CMRS.b1","CMRS.b1_sq")
		template = template.replace("%CMIS.b1","CMIS.b1_sq")
	else:
		template = template.replace("%CMRS.b1","CMRS.b1")
		template = template.replace("%CMIS.b1","CMIS.b1")
	template = template.replace("%colknob1",error[9])
	template = template.replace("%colknob5",error[10])
	for key, value in replace_dict.items():
		template = template.replace(str(key), str(value))
	
	new_path = path.replace(".madx","_temp.madx")
	with open(new_path, "w") as job_file:
		job_file.write(template)

def get_beta_error(error):
	error0 = copy.deepcopy(error)
	set_template("response_betaerror.madx",0.,0.,error0)
	os.system(madx_path + " response_betaerror_temp.madx")
	tw40cm = twiss("twiss.original")
	tw40cm.Cmatrix()
	betax1, betay1 =  np.array(tw40cm.BETX), np.array(tw40cm.BETY)	
	
	
	error1 = copy.deepcopy(error)
	error1[2] = "0."
	error1[7] = "0."
	error1[9] = "0."
	error1[10] = "0."	

	set_template("response_betaerror.madx",0.,0.,error1)
	os.system(madx_path + " response_betaerror_temp.madx")
	tw40cm = twiss("twiss.original")
	tw40cm.Cmatrix()
	betax0 , betay0 = np.array(tw40cm.BETX), np.array(tw40cm.BETY)
	
	N = len(betax0)
	betax_rmsError = np.sum(((betax1 - betax0) / betax0)**2) / N
	betay_rmsError = np.sum(((betay1 - betay0) / betay0)**2) / N
	betax_error = (betax1 - betax0) / betax0
	betay_error = (betay1 - betay0) / betay0
	beta_max_error = max(np.max(betax_error),np.max(betay_error))
	return betax_rmsError,betay_rmsError,beta_max_error

def get_responsematrix(squeeze):
	error_l = np.zeros(N_error)
	error_type = "quadrupole"
	error_size = "0."
	pattern_1 = "."
	pattern_2 = "."
	quad_component = "quadrupole"
	quad_pattern = "."
	quad_strength = "0."
	colknob1 = "0."
	colknob5 = "0."
	error = [error_l,error_type,error_size,pattern_1,pattern_2,quad_component,quad_pattern,quad_strength,squeeze,colknob1,colknob5]
			

	knob_Re_1 = 0.0001
	knob_Im_1 = 0.
	set_template("response_betaerror.madx",knob_Re_1, knob_Im_1,error)
	os.system(madx_path + " response_betaerror_temp.madx")
	tw40cm = twiss("twiss.original")
	tw40cm.Cmatrix()
	f_R_1 = np.array(tw40cm.F1001R)
	f_I_1 = np.array(tw40cm.F1001I)
	f_1 = np.concatenate((f_R_1,f_I_1))

	knob_Re_2 = 0.
	knob_Im_2 = 0.0001
	set_template("response_betaerror.madx",knob_Re_2, knob_Im_2, error)
	os.system(madx_path + " response_betaerror_temp.madx")
	tw40cm = twiss("twiss.original")
	tw40cm.Cmatrix()
	f_R_2 = np.array(tw40cm.F1001R)
	f_I_2 = np.array(tw40cm.F1001I)
	f_2 = np.concatenate((f_R_2,f_I_2))

	R = np.array([f_1/knob_Re_1,f_2/knob_Im_2]).T
	R_inverse = np.linalg.pinv(R)
	return R_inverse

def get_madx_knobs(error):
	set_template("FineTuneCoupling_betaerror.madx",0.,0.,error)
	os.system(madx_path + " FineTuneCoupling_betaerror_temp.madx")
	knob_file = open("knob.txt")
	lines = knob_file.readlines()
	knob_Re = float(lines[0][10:-2])
	knob_Im = float(lines[1][10:-2])
	return knob_Re, knob_Im

def get_response_knobs(R_inverse,error):
	set_template("response_betaerror.madx",0.,0., error)
	os.system(madx_path + " response_betaerror_temp.madx")
	tw40cm = twiss("twiss.original")
	tw40cm.Cmatrix()
	f_R = np.array(tw40cm.F1001R)
	f_I = np.array(tw40cm.F1001I)
	f = np.concatenate((f_R,f_I))
	f_c = f_R + 1j * f_I
	knobs = -1 * np.dot(R_inverse,f)
	knob_Re, knob_Im = knobs[0], knobs[1]
	return tw40cm, knob_Re, knob_Im


def get_C_min(knob_Re,knob_Im,error):
	set_template("exact_C_min_betaerror.madx",knob_Re, knob_Im, error)
	os.system(madx_path + " exact_C_min_betaerror_temp.madx")
	tw40cm = twiss("twiss.C_min")
	tw40cm.Cmatrix()
	Q1, Q2 = tw40cm.Q1, tw40cm.Q2
	C_min = abs(Q1 - Q2 - 2)
	return C_min


def get_twiss_rdt(knob_Re,knob_Im,error):
	set_template("response_betaerror.madx",knob_Re,knob_Im,error)
	os.system(madx_path + " response_betaerror_temp.madx")
	tw40cm = twiss("twiss.original")
	tw40cm.Cmatrix()
	f_R, f_I = np.array(tw40cm.F1001R) , np.array(tw40cm.F1001I)
	return tw40cm, f_R + 1J * f_I
	


def plot_response(error,savepath):
	knob_Re_madx, knob_Im_madx = get_madx_knobs(error)	
	f_madx = get_twiss_rdt(knob_Re_madx,knob_Im_madx,error)[1]	
	C_min_0 = get_C_min(0.,0.,error)

	R_inverse = get_responsematrix(error[8])
	tw40cm0, knob_Re_res, knob_Im_res = get_response_knobs(R_inverse,error)
	f_res = get_twiss_rdt(knob_Re_res,knob_Im_res,error)[1]
	C_min_res = get_C_min(knob_Re_res,knob_Im_res,error)
	
	S = np.array(tw40cm0.S)	
	f_0 = np.array(tw40cm0.F1001R) + 1j * np.array(tw40cm0.F1001I)
	
	
	betax_rmsError,betay_rmsError,beta_max_error = get_beta_error(error)
	fig = plt.figure()
	fig.suptitle("$C_-^0$ =" + str(C_min_0) + " $C_-^{res}$ " + str(C_min_res) + "\n" +
 r"$\frac{\delta\beta_{beta_max}}{beta}$ =" + str(round(beta_max_error,3)) + r"	$\frac{\delta\beta}{beta}_{RMS}$ =" + str(round(betax_rmsError,3)) + "	colknob1 =" + error[9] +"	colknob5 =" + error[10])
	
	ax = fig.add_subplot(1,2,1)
	ax.plot(S,abs(f_0),label = "before response")
	ax.plot(S,abs(f_madx),label = "madx response")
	ax.plot(S,abs(f_res),label = "response matrix")
	ax.set_xlabel("S")
	ax.set_ylabel("F1001")
	ax.legend()

	#test to compare betafunction at IR1 and IR5
	ax2 = fig.add_subplot(1,2,2)	
	ax2.plot(S,tw40cm0.BETX)
	ax2.plot(S,tw40cm0.BETY)

	plt.savefig('plots/' + savepath)
	plt.show()

def plot_colknob_comparison(error,savepath,n_colknob,colknob_step,colknob_type):
	R_inverse = get_responsematrix(error[8])
	error[9] , error[10] = "0." , "0."
	colknob_list = colknob_step*np.arange(n_colknob)
	C_min_madx_list = np.zeros(n_colknob)
	C_min_res_list = np.zeros(n_colknob)
	C_min_0_list = np.zeros(n_colknob)
	beta_max_error_list = np.zeros(n_colknob)
	betax_rms_error_list = np.zeros(n_colknob)
	for i in range(n_colknob):
		if colknob_type == "1":
			error[9] = str(colknob_list[i])
		elif colknob_type == "5":
			error[10] = str(colknob_list[i])
		else:
			print("error")
			break
		
		knob_Re_madx, knob_Im_madx = get_madx_knobs(error)
		tw40cm0, knob_Re_res, knob_Im_res = get_response_knobs(R_inverse,error)
		betax_rms_error,betay_rms_error,beta_max_error = get_beta_error(error)

		f_res = get_twiss_rdt(knob_Re_res,knob_Im_res,error)[1]
		f_madx = get_twiss_rdt(knob_Re_madx,knob_Im_madx,error)[1]	

		C_min_0 = get_C_min(0.,0.,error)
		C_min_res = get_C_min(knob_Re_res,knob_Im_res,error)
		C_min_madx = get_C_min(knob_Re_madx, knob_Im_madx,error)

		C_min_madx_list[i] = C_min_madx
		C_min_res_list[i] = C_min_res
		C_min_0_list[i] = C_min_0
		beta_max_error_list[i] = beta_max_error
		betax_rms_error_list[i] = betax_rms_error
	
	fig = plt.figure()
	ax1 = fig.add_subplot(1,2,1)
	ax1.plot(colknob_list, C_min_0_list,label = "C_min_0")
	ax1.plot(colknob_list, C_min_madx_list,label = "C_min_madx")
	ax1.plot(colknob_list, C_min_res_list,label = "C_min_res")
	ax1.set_xlabel("colknob" + colknob_type)
	ax1.set_ylabel("$C_-$")
	ax1.legend()

	ax2 = fig.add_subplot(1,2,2)
	ax2.plot(colknob_list,beta_max_error_list,label = "max_beta_error")
	ax2.plot(colknob_list,betax_rms_error_list,label = "rms_beta_error")
	ax2.set_xlabel("colknbob"+ colknob_type)
	ax2.set_ylabel("beta_error")
	ax2.legend()

	plt.savefig("plots/" + savepath)
	#plt.show()
		
		
		

error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0."
pattern_1 = "."
pattern_2 = "."
quad_component = "quadrupole"
quad_pattern_1 = "."
#quad_strength = "0.000005 * GAUSS()"
quad_strength = "0."
squeeze = True
colknob1 = "0."
colknob5 = "0."
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,colknob1,colknob5]
#plot_response(error,"noError_colknob1.pdf")
#plot_colknob_comparison(error,"colknob_comparison_5.pdf",10,0.3,"5")
#plot_colknob_comparison(error,"colknob_comparison_1.pdf",10,0.3,"1")


error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0.0001"
pattern_1 = "3[1-4]R"
pattern_2 = "3[1-4]R"
quad_component = "quadrupole"
quad_pattern_1 = "."
#quad_strength = "0.000005 * GAUSS()"
quad_strength = "0."
squeeze = True
colknob1 = "0"
colknob5 = "0"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,colknob1,colknob5]
#plot_colknob_comparison(error,"colknob_comparison_error_5.pdf",10,0.3,"5")
#plot_colknob_comparison(error,"colknob_comparison_error_1.pdf",10,0.3,"1")

error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0."
pattern_1 = "."
pattern_2 = "."
quad_component = "quadrupole"
quad_pattern_1 = "."
#quad_strength = "0.000005 * GAUSS()"
quad_strength = "0."
squeeze = True
colknob1 = "0"
colknob5 = "0"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,colknob1,colknob5]
plot_colknob_comparison(error,"colknob_5_zeroError_fewSteps.pdf",3,0.1,"5")
plot_colknob_comparison(error,"colknob_1_zeroError_fewSteps.pdf",3,0.5,"1")
#plot_response(error,"colknob1_blowup_test")

error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0."
pattern_1 = "."
pattern_2 = "."
quad_component = "quadrupole"
quad_pattern_1 = "."
#quad_strength = "0.000005 * GAUSS()"
quad_strength = "0."
squeeze = True
colknob1 = "3."
colknob5 = "0"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,colknob1,colknob5]
#plot_response(error,"local_coupling_noSkewQuad4.pdf")


error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0.0000005"
pattern_1 = "."
pattern_2 = "."
quad_component = "quadrupole"
quad_pattern_1 = "."
quad_strength = "0.000005 * GAUSS()"
#quad_strength = "0."
squeeze = True
colknob1 = "0."
colknob5 = "0"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,colknob1,colknob5]
#plot_response(error,"local_coupling_colknob_0.pdf")

error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0.0000005"
pattern_1 = "."
pattern_2 = "."
quad_component = "quadrupole"
quad_pattern_1 = "."
quad_strength = "0.000005 * GAUSS()"
#quad_strength = "0."
squeeze = True
colknob1 = "0"
colknob5 = "5"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,colknob1,colknob5]
#plot_response(error,"local_coupling_colknob_05.pdf")
#plot_response(error,"test_random_error")

error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0.0001"
pattern_1 = "3[1-4]R"
pattern_2 = "3[1-4]R"
quad_component = "quadrupole"
quad_pattern_1 = "."
#quad_strength = "0.000005 * GAUSS()"
quad_strength = "0."
squeeze = True
colknob1 = "3"
colknob5 = "0"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,colknob1,colknob5]
#plot_response(error,"local_coupling_colknob5_test.pdf")



