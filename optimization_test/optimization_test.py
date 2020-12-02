#-*- coding: utf-8 -*-

import sys
from scipy.optimize import minimize
sys.path.append("../global_coupling/")
from global_coupling import *
      
modulename = "global_coupling"
if modulename in sys.modules:
    print("You have imported the {} module".format(modulename))

N_error = 12

#home
madx_path = "/home/eirik/madx "
folder_path = "/home/eirik/CERN/global_coupling_correction/optimization_test/"
lhc_path = "/home/eirik/CERN/lhc2018/2018"

#work
#madx_path = "/home/ehoydals/madx "
#folder_path = "/home/ehoydals/global_coupling_correction/analytical_test/"
#lhc_path =  "/afs/cern.ch/eng/lhc/optics/runII/2018"

response_path = folder_path + "response_optimization_test.madx"
C_min_path = folder_path + "exact_C_min_optimization_test.madx"
FineTuneCoupling_path = folder_path + "FineTuneCoupling_optimization_test.madx"

path_dict = {}
path_dict["madx_path"] = madx_path
path_dict["response_path"] = response_path
path_dict["C_min_path"] = C_min_path
path_dict["FineTuneCoupling_path"] = FineTuneCoupling_path
set_global_paths(path_dict)



reset_response_dict_new = {
	"%error_strength" : "0.",
	"%quad_strength" : "0.",
	"%colknob1" : "0.",
	"%colknob5" : "0."
}

reset_analyticresponse_dict_new = {
	"%error_strength" : "0.",
	"%quad_strength" : "0.",
	"%colknob1" : "0.",
	"%colknob5" : "0."
}

reset_beta_error_dict_new = {
	"%error_strength" : "0.", 
	"%quad_strength" : "0.",
	"%colknob1" : "0.",
	"%colknob5" : "0."
}


set_reset_dict(reset_response_dict,reset_response_dict_new)
set_reset_dict(reset_analyticresponse_dict,reset_analyticresponse_dict_new)
set_reset_dict(reset_beta_error_dict,reset_beta_error_dict_new)



def test_absolute_response(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	S , f_0 = get_S_f(response_path,"twiss.original",change_dict_local)
	
	R_Re , R_Im = get_responsematrix_components(change_dict_local)[0:2]
	R = ((np.real(f_0) / abs(f_0)) * R_Re.T).T + ((np.imag(f_0) / abs(f_0)) * R_Im.T).T
	R_inverse_abs = np.linalg.pinv(R)
	knobs = - np.dot(R_inverse_abs,abs(f_0))
	knob_Re_abs, knob_Im_abs = knobs[0] , knobs[1]
	set_knobs(change_dict_local,knob_Re_abs,knob_Im_abs)
	f_res_abs = get_f(response_path,"twiss.original",change_dict_local)
	
	set_knobs(change_dict_local,0.,0.)
	R = np.block([[R_Re],[R_Im]])
	R_inverse = np.linalg.pinv(R)
	knob_Re_res , knob_Im_res = get_response_knobs(R_inverse,change_dict_local)
	set_knobs(change_dict_local,knob_Re_res,knob_Im_res)
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	
	fig = plt.figure()
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_res_abs),label = "f_res_abs")
	ax1.plot(S,abs(f_res),label = "f_res")
	ax1.legend()
	plt.show()
	
def test(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	S, f_0 = get_S_f(response_path,"twiss.original",change_dict_local)
	C_min_0 = get_C_min(change_dict_local)
	
	R_inverse = get_analytic_responsematrix(change_dict_local)
	knob_Re,knob_Im = get_response_knobs(R_inverse,change_dict_local)
	R = np.linalg.pinv(R_inverse)
	knobs = np.array([knob_Re,knob_Im])
	f_predicted_L = np.dot(R,-knobs)
	f_Re = f_predicted_L[0:len(f_predicted_L)/2:]
	f_Im = f_predicted_L[len(f_predicted_L)/2:len(f_predicted_L):]
	f_predicted = f_Re + 1j * f_Im
	
	set_knobs(change_dict_local,knob_Re,knob_Im)
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res = get_C_min(change_dict_local)
	
	
	f_res_predicted = f_0 - f_predicted
	f_res_predicted_l = np.concatenate((np.real(f_res_predicted),np.imag(f_res_predicted)))
	knobs = -np.dot(R_inverse,f_res_predicted_l)
	knob_Re_double , knob_Im_double = knobs[0] , knobs[1]
	set_knobs(change_dict_local,knob_Re + knob_Re_double,knob_Im + knob_Im_double)
	f_res_double = get_f(response_path,"twiss.original",change_dict_local)
	
	print(knobs)
	fig = plt.figure()
	fig.suptitle(r"$C_-^0$ = " + "{:.3e}".format(C_min_0) + r"$C_-^{res}$ = " + "{:.3e}".format(C_min_res))
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_predicted),label = "f_predicted")
	ax1.plot(S,abs(f_res_predicted),label = "f_correction_predicted")
	ax1.plot(S,abs(f_res_double),label = "f_res_double")
	ax1.plot(S,abs(f_res),label = "f_res")
	ax1.set_ylim(ymin = 0)
	ax1.legend()
	plt.savefig("plots/" + savepath)
	plt.show()
	
def plot_betafunc_comparison(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	C_min_0 = get_C_min(change_dict_local)
	
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	S = np.array(tw40cm.S)
	betx_perturbed , bety_perturbed = np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	
	change_value(change_dict_local,"%error_strength","0.")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	betx , bety =  np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	
	fig = plt.figure()
	fig.suptitle(str(C_min_0))
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(S,(betx_perturbed - betx) / betx,label = "betx")
	ax1.plot(S,(bety_perturbed - bety) / bety,label = "bety")
	#ax1.plot(S,betx_perturbed,label = "betx_perturbed")
	#ax1.plot(S,bety_perturbed,label = "bety_perturbed")
	ax1.legend()
	plt.savefig("plots/" + savepath)
	plt.show()
	
	
	

	
def optimize_knobs(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	S = np.array(tw40cm.S)
	f_0 = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	Q1 , Q2 = tw40cm.Q1 , tw40cm.Q2
	MUX , MUY = np.array(tw40cm.MUX) , np.array(tw40cm.MUY)
	
	R_Re , R_Im = get_responsematrix_components(change_dict_local)[0:2]
	R_z = R_Re + 1j * R_Im
	R = np.block([[R_Re],[R_Im]])
	R_inverse = np.linalg.pinv(R)
	
	N = len(f_0)
	delta = np.abs(Q1%1 - Q2%1)
	deltaS_list = [np.abs(S[i +1] - S[i]) for i in range(N-1)]
	deltaS_list.append(S[-1] - S[-2])
	deltaS = np.array(deltaS_list)
	#check complex conjugate in dot
	def C_min_local(knobs):
		return  4 * delta / S[-1] * abs(np.sum(deltaS * (f_0 + np.dot(R_z,knobs)) * np.exp(-1j * 2*np.pi*(MUY - MUX))))
	
	
	x0 = np.array([0,0])
	res = minimize(C_min_local, x0, method='nelder-mead',options={'xatol': 1e-8, 'disp': True})
	knob_Re_opt , knob_Im_opt = res.x[0] , res.x[1]
	set_knobs(change_dict_local,knob_Re_opt,knob_Im_opt)
	f_opt = get_f(response_path,"twiss.original",change_dict_local)
	C_min_opt = get_C_min(change_dict_local)
	
	set_knobs(change_dict_local,0.,0.)
	knob_Re_res , knob_Im_res = get_response_knobs(R_inverse,change_dict_local)
	set_knobs(change_dict_local,knob_Re_res,knob_Im_res)
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res = get_C_min(change_dict_local)
	
	
	fig = plt.figure()
	fig.suptitle("C_min_opt = " + "{:.3e}".format(C_min_opt) + "C_min_res = " + "{:.3e}".format(C_min_res))
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_res),label = "f_res")
	ax1.plot(S,abs(f_opt),label = "f_opt")
	ax1.legend()
	plt.savefig("plots/" + savepath)
	plt.show()
	
def delQ_min_new(Q1,Q2,f,S,R,MUX,MUY):
	N = len(f)
	#should this be abs
	delta = Q1 - Q2 - 2

	deltaS_list = [np.abs(S[i +1] - S[i]) for i in range(N-1)]
	deltaS_list.append(S[-1] - S[-2])
	deltaS = np.array(deltaS_list)
	return abs(4 * delta / (2*np.pi*R) * np.sum(deltaS * f * np.exp(-1j * 2*np.pi*(MUY - MUX - S * delta / R))))


def delQ_min_old(Q1,Q2,f,S,MUX,MUY):
	N = len(f)
	delta = np.abs(Q1%1 - Q2%1)
	return 4 * delta/N * np.sum(np.abs(f))	

def projection_test(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)	
	
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	C_min_0 = get_C_min(change_dict_local)
	S = np.array(tw40cm.S)
	f_0 = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	Q1 , Q2 = tw40cm.Q1 , tw40cm.Q2
	delta = Q1 - Q2 - 2
	R = tw40cm.LENGTH / 2 / np.pi
	MUX , MUY = np.array(tw40cm.MUX) , np.array(tw40cm.MUY)
	N = len(S)
	deltaS_list = [np.abs(S[i +1] - S[i]) for i in range(N-1)]
	deltaS_list.append(S[-1] - S[-2])
	deltaS = np.array(deltaS_list)
	
	
	x = deltaS * np.exp(-2*np.pi*1j*(MUY - MUX - S * delta / R))
	x_hat = x / np.sqrt(np.dot(x,x))
	f_projection = np.dot(x_hat,f_0) * x_hat 
	f_ort = f_0 - f_projection
	
	print(C_min_0)
	print(delQ_min_new(Q1,Q2,f_0,S,R,MUX,MUY))
	print(delQ_min_new(Q1,Q2,f_projection,S,R,MUX,MUY))
	print(delQ_min_new(Q1,Q2,f_ort,S,R,MUX,MUY))
	
	
	fig, (ax1,ax2) = plt.subplots(2)
	ax1.plot(S,np.real(f_0),label = "Re{f_0}")
	ax1.plot(S,np.real(f_projection),label = "Re{f_projection}")
	ax1.legend()
	
	ax2.plot(S,np.imag(f_0),label = "Im{f_0}")
	ax2.plot(S,np.imag(f_projection),label = "Im{f_projection}")
	ax2.legend()
	plt.show()
	
	
	plt.plot(S,abs(f_0),label = "f_0")
	plt.plot(S,abs(f_ort),label = "f_ort")
	plt.legend()
	plt.savefig("plots/" + savepath)
	plt.show()
	
	
def double_correction(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)	
	
	S , f_0 = get_S_f(response_path,"twiss.original",change_dict_local)
	C_min_0 = get_C_min(change_dict_local)
	
	R_Re , R_Im = get_responsematrix_components(change_dict_local)[0:2]
	R_z = R_Re + 1j * R_Im
	R = np.block([[R_Re],[R_Im]])
	R_inverse = np.linalg.pinv(R)
	
	knob_Re_res , knob_Im_res = get_response_knobs(R_inverse,change_dict_local)
	set_knobs(change_dict_local,knob_Re_res,knob_Im_res)
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	
	del_knob_Re_res , del_knob_Im_res = get_response_knobs(R_inverse,change_dict_local)
	knob_Re_double_res , knob_Im_double_res = knob_Re_res + del_knob_Re_res, knob_Im_res + del_knob_Im_res
	set_knobs(change_dict_local,knob_Re_double_res,knob_Im_double_res)
	f_double_res = get_f(response_path,"twiss.original",change_dict_local)

	R = np.linalg.pinv(R_inverse)
	knobs_res = np.array([knob_Re_res,knob_Im_res])
	f_res_predicted = f_0 + np.dot(R_z,knobs_res)
	
	del_knobs_res = np.array([del_knob_Re_res,del_knob_Im_res])
	f_double_res_predicted = f_res + np.dot(R_z,del_knobs_res)
	
	fig, ax1 = plt.subplots(1)
	fig.suptitle("$C^-_0 =$" + "{:.3e}".format(C_min_0))		
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_res),label = "f_res")
	ax1.plot(S,abs(f_res_predicted),label = "f_res_predicted")
	ax1.plot(S,abs(f_double_res),label = "f_double_res")
	ax1.plot(S,abs(f_double_res_predicted),label = "f_double_res_predicted")
	ax1.legend()
	plt.savefig("plots/" + savepath)
	plt.show()




change_dict = {}
change_dict["%lhc_path"] = lhc_path
#change_dict["%opticsfile"] = "opticsfile.19"
change_dict["%opticsfile"] =  "opticsfile.1"
change_dict["%knob_Re_value"] = "0."
change_dict["%knob_Im_value"] = "0."
change_dict["%knob_Re_type"] = "CMRS.b1_sq"
change_dict["%knob_Im_type"] = "CMIS.b1_sq"
change_dict["%error_component"] = "quadrupole"
change_dict["%error_strength"] = "0.00003*gauss()"
change_dict["%pattern_1"] = "."
change_dict["%pattern_2"] = "."
change_dict["%quad_component"] = "quadrupole"
change_dict["%quad_pattern_1"] = "R5"
change_dict["%quad_pattern_2"] = "R5"
#change_dict["%quad_strength"] = "0.00018"
change_dict["%quad_strength"] = "0."
change_dict["%twiss_pattern"] = "BPM"
change_dict["%colknob1"] = "0."
change_dict["%colknob5"] = "0."
change_dict["%Qx"] = "62.29"
change_dict["%Qy"] = "60.34"

#test(change_dict,"test.pdf")
#optimize_knobs(change_dict,"test_optimization.pdf")
change_dict_local = copy.deepcopy(change_dict)
#projection_test(change_dict_local,"projection_test.pdf")

change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%quad_strength","0.")
#plot_betafunc_comparison(change_dict_local,"betafunc_comparison.pdf")

change_dict_local = copy.deepcopy(change_dict)
#optimize_knobs(change_dict_local,"optimize_knobs.pdf")

change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%quad_strength","0.00018")
change_value(change_dict_local,"%error_strength","0.00003*gauss()")
change_value(change_dict_local,"%pattern_1","R4")
change_value(change_dict_local,"%pattern_2","L3")
double_correction(change_dict_local,"double_correction_R4L3.pdf")


change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%quad_strength","0.00018")
change_value(change_dict_local,"%error_strength","0.00003*gauss()")
change_value(change_dict_local,"%pattern_1",".")
change_value(change_dict_local,"%pattern_2",".")
double_correction(change_dict_local,"double_correction_all.pdf")




