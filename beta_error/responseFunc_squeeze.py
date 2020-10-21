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
	
	template = template.replace("%colknob1","0.")
	template = template.replace("%colknob5","0.")
	template = template.replace("%twiss_pattern",error[9])

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
	return betax_rmsError,betay_rmsError,beta_max_error, (betax0 - betax1)/betax0

def get_responsematrix(squeeze):
	error_l = np.zeros(N_error)
	error_type = "quadrupole"
	error_size = "0."
	pattern_1 = "."
	pattern_2 = "."
	quad_component = "quadrupole"
	quad_pattern = "."
	quad_strength = "0."
	twiss_pattern = "BPM"
	error = [error_l,error_type,error_size,pattern_1,pattern_2,quad_component,quad_pattern,quad_strength,squeeze,twiss_pattern]

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
	knob_Re = float(lines[0].split("=")[1][0:-2:])
	knob_Im = float(lines[1].split("=")[1][0:-2:])
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

def plot_squeeze_comparison(error,savepath):	
	C_min_0 = get_C_min(0.,0.,error)

	error[8] = False
	knob_Re_madx, knob_Im_madx = get_madx_knobs(error)	
	f_madx = get_twiss_rdt(knob_Re_madx,knob_Im_madx,error)[1]
	C_min_madx = get_C_min(knob_Re_madx, knob_Im_madx,error)
	
	R_inverse = get_responsematrix(error[8])
	tw40cm0, knob_Re_res, knob_Im_res = get_response_knobs(R_inverse,error)
	f_res = get_twiss_rdt(knob_Re_res,knob_Im_res,error)[1]
	C_min_res = get_C_min(knob_Re_res,knob_Im_res,error)
	
	error[9] = "."
	tw40cm1 = get_twiss_rdt(knob_Re_res,knob_Im_res,error)[0]
	K1SL = np.array(tw40cm1.K1SL)
	K1SL_nonzero = K1SL[K1SL != 0]
	mean_strength = np.mean(abs(K1SL_nonzero))
	error[9] = "BPM"
	

	error[8] = True
	knob_Re_madx_sq, knob_Im_madx_sq = get_madx_knobs(error)	
	f_madx_sq = get_twiss_rdt(knob_Re_madx,knob_Im_madx,error)[1]
	C_min_madx_sq = get_C_min(knob_Re_madx_sq, knob_Im_madx_sq,error)
	
	R_inverse_sq = get_responsematrix(error[8])
	tw40cm0_sq, knob_Re_res_sq, knob_Im_res_sq = get_response_knobs(R_inverse_sq,error)
	f_res_sq = get_twiss_rdt(knob_Re_res_sq,knob_Im_res_sq,error)[1]
	C_min_res_sq = get_C_min(knob_Re_res_sq,knob_Im_res_sq,error)
	
	error[9] = "."
	tw40cm1_sq = get_twiss_rdt(knob_Re_res_sq,knob_Im_res_sq,error)[0]
	K1SL_sq = np.array(tw40cm1_sq.K1SL)
	K1SL_nonzero_sq = K1SL_sq[K1SL_sq != 0]
	mean_strength_sq = np.mean(abs(K1SL_nonzero_sq))
	error[9] = "BPM"
	
	print(C_min_res_sq)
	print(C_min_res)
	time.sleep(20)
	
	S = np.array(tw40cm0.S)	
	f_0 = np.array(tw40cm0.F1001R) + 1j * np.array(tw40cm0.F1001I)	
	betax_rms_error,betay_rms_error,beta_max_error, beta_TEST = get_beta_error(error)

	fig = plt.figure(figsize=plt.figaspect(0.5))

	ax1 = fig.add_subplot(1,2,1)
	ax1.plot(S,abs(f_0),label = "before response")
	ax1.plot(S,abs(f_madx),label = "madx matching correction")
	ax1.plot(S,abs(f_res),label = "response matrix")
	ax1.set_title("regular knobs \n $C_-^{res}$ =" + "{:.2e}".format(C_min_res) + "	$C_-^{match}$ =" + "{:.2e}".format(C_min_madx) +"\n"+ r"<|K1SL|> =" + "{:.2e}".format(mean_strength))
	ax1.set_xlabel("S")
	ax1.set_ylabel("F1001")	
	ax1.legend()	

	ax2 = fig.add_subplot(1,2,2)
	ax2.plot(S,abs(f_0),label = "before response")
	ax2.plot(S,abs(f_madx_sq),label = "madx matching correction")
	ax2.plot(S,abs(f_res_sq),label = "response matrix")
	ax2.set_title("squeeze knobs \n $C_-^{res sq}$ =" + "{:.2e}".format(C_min_res_sq) + "	$C_-^{match sq}$ =" + "{:.2e}".format(C_min_madx_sq) + "\n" + r"$<|K1SL|>_{sq}$ =" + "{:.2e}".format(mean_strength_sq))
	ax2.set_xlabel("S")
	ax2.set_ylabel("F1001")
	ax2.legend()
	fig.tight_layout(rect=[0, 0.03, 1, 0.95])	
	plt.savefig('plots/' + savepath)
	plt.show()
	#plt.plot(S,beta_TEST)
	#plt.show()


error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0.00004*gauss()"
pattern_1 = "R7"
pattern_2 = "R7"
quad_component = "quadrupole"
quad_pattern_1 = "R3"
quad_strength = "0.00008"
squeeze = True
twiss_pattern = "BPM"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
#plot_squeeze_comparison(error,"collisionComparison_randomLocalSQ_localrQ.pdf")
quad_strength = "0."
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
#plot_squeeze_comparison(error,"collisionComparison_randomLocalSQ_localrQ_noBetabeat.pdf")


error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0.0000003"
pattern_1 = "."
pattern_2 = "."
quad_component = "quadrupole"
quad_pattern_1 = "R3"
quad_strength = "0.00008"
squeeze = True
twiss_pattern = "BPM"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
plot_squeeze_comparison(error,"collisionComparison_UniformSQ_localQ.pdf")
quad_strength = "0."
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
plot_squeeze_comparison(error,"collisionComparison_UniformSQ_localQ_noBetabeat.pdf")

error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0.000012"
pattern_1 = "L7"
pattern_2 = "R8"
quad_component = "quadrupole"
quad_pattern_1 = "R3"
quad_strength = "0.00008"
squeeze = False
twiss_pattern = "BPM"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
#plot_squeeze_comparison(error,"collisionComparison_twoSQ_localQ.pdf")
quad_strength = "0."
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
#plot_squeeze_comparison(error,"collisionComparison_twoSQ_localQ_noBetabeat.pdf")
	
"""
error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0.00004*gauss()"
pattern_1 = "R7"
pattern_2 = "R7"
quad_component = "quadrupole"
quad_pattern_1 = "R3"
quad_strength = "0.00008"
squeeze = True
twiss_pattern = "BPM"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
plot_squeeze_comparison(error,"injectionComparison_randomLocalSQ_localrQ.pdf")
quad_strength = "0."
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
plot_squeeze_comparison(error,"injectionComparison_randomLocalSQ_localrQ_noBetabeat.pdf")


error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0.000002"
pattern_1 = "."
pattern_2 = "."
quad_component = "quadrupole"
quad_pattern_1 = "R3"
quad_strength = "0.00008"
squeeze = True
twiss_pattern = "BPM"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
plot_squeeze_comparison(error,"injectionComparison_UniformSQ_localQ.pdf")
quad_strength = "0."
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
plot_squeeze_comparison(error,"injectionComparison_UniformSQ_localQ_noBetabeat.pdf")

error_l = np.zeros(N_error)
error_component = "quadrupole"
error_strength = "0.00002"
pattern_1 = "L7"
pattern_2 = "R8"
quad_component = "quadrupole"
quad_pattern_1 = "R3"
quad_strength = "0.00008"
squeeze = False
twiss_pattern = "BPM"
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
plot_squeeze_comparison(error,"injectionComparison_twoSQ_localQ.pdf")
quad_strength = "0."
error = [error_l,error_component,error_strength,pattern_1,pattern_2,quad_component,quad_pattern_1,quad_strength,squeeze,twiss_pattern]
plot_squeeze_comparison(error,"injectionComparison_twoSQ_localQ_noBetabeat.pdf")
"""

"""
R_inverse = get_responsematrix(error[8])
tw40cm0, knob_Re_res, knob_Im_res = get_response_knobs(R_inverse,error)
tw40cm1, f_res = get_twiss_rdt(knob_Re_res,knob_Im_res,error)
S = np.array(tw40cm1.S)
K1SL1 = np.array(tw40cm1.K1SL)
error[8] = True
R_inverse = get_responsematrix(error[8])
tw40cm0, knob_Re_res, knob_Im_res = get_response_knobs(R_inverse,error)
tw40cm1, f_res = get_twiss_rdt(knob_Re_res,knob_Im_res,error)
S = np.array(tw40cm1.S)
K1SL2 = np.array(tw40cm1.K1SL)
print(K1SL1[K1SL1 != 0])
print(K1SL2[K1SL2 != 0])
print(np.mean(abs(K1SL1[K1SL1 != 0])))
print(np.mean(abs(K1SL2[K1SL2 != 0])))
#a = np.array([0.,313,0.,0.,0.,321,33,0.000,4,0.3])
#print(a[a != 0])
"""
