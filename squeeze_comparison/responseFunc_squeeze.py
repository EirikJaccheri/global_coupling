import sys
sys.path.append("../global_coupling/")
from global_coupling import *
      
modulename = "global_coupling"
if modulename in sys.modules:
    print("You have imported the {} module".format(modulename))

N_error = 12

madx_path = "/home/eirik/madx "

folder_path = "/home/eirik/CERN/global_coupling_correction/squeeze_comparison/"
response_path = folder_path + "response_squeeze.madx"
C_min_path = folder_path + "exact_C_min_squeeze.madx"
FineTuneCoupling_path = folder_path + "FineTuneCoupling_squeeze.madx"

path_dict = {}
path_dict["madx_path"] = madx_path
path_dict["response_path"] = response_path
path_dict["C_min_path"] = C_min_path
path_dict["FineTuneCoupling_path"] = FineTuneCoupling_path
set_global_paths(path_dict)

#lhc_path =  "/afs/cern.ch/eng/lhc/optics/runII/2018" #work
lhc_path = "/home/eirik/CERN/lhc2018/2018" #home

reset_response_dict_new = {"%error_component" : "quadrupole",
	"%error_strength" : "0.",
	"%pattern_1" : ".",
	"%pattern_2" : ".",
	"%quad_component" : "quadrupole",
	"%quad_pattern_1" : ".",
	"%quad_strength" : "0.",
	"%twiss_pattern" : "BPM",
}



reset_beta_error_dict_new = {"%error_strength" : "0.", "%quad_strength" : "0."}


set_reset_dict(reset_response_dict,reset_response_dict_new)
set_reset_dict(reset_beta_error_dict,reset_beta_error_dict_new)

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



def plot_squeeze_comparison(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	change_value(change_dict_local,"%knob_Re_value","0.")
	change_value(change_dict_local,"%knob_Im_value","0.")
		
	C_min_0 = get_C_min(change_dict_local)
	
	
		
	#squeeze knobs
	change_value(change_dict_local,"%knob_Re_type","CMRS.b1_sq")
	change_value(change_dict_local,"%knob_Im_type","CMIS.b1_sq")
	knob_Re_madx_sq, knob_Im_madx_sq = get_madx_knobs(change_dict_local)
	set_knobs(change_dict_local,knob_Re_madx_sq, knob_Im_madx_sq) 
	f_madx_sq =  get_f(response_path,"twiss.original",change_dict_local)
	C_min_madx_sq = get_C_min(change_dict_local)
	set_knobs(change_dict_local,0.,0.)

	
	R_inverse_sq = get_responsematrix(change_dict_local)
	tw40cm0_sq = get_twiss(response_path,"twiss.original",change_dict_local)
	knob_Re_res_sq, knob_Im_res_sq = get_response_knobs(R_inverse_sq,change_dict_local)
	set_knobs(change_dict_local,knob_Re_res_sq, knob_Im_res_sq) 
	f_res_sq = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res_sq = get_C_min(change_dict_local)
	#double correction
	dKnob_Re_res_sq , dKnob_Im_res_sq = get_response_knobs(R_inverse_sq,change_dict_local)
	set_knobs(change_dict_local,knob_Re_res_sq + dKnob_Re_res_sq , knob_Im_res_sq + dKnob_Im_res_sq )
	f_double_res_sq =  get_f(response_path,"twiss.original",change_dict_local)
	C_min_double_res_sq = get_C_min(change_dict_local)
	set_knobs(change_dict_local,0.,0.)
	
	change_value(change_dict_local,"%twiss_pattern",".")
	set_knobs(change_dict_local,knob_Re_res_sq, knob_Im_res_sq) 
	tw40cm1_sq = get_twiss(response_path,"twiss.original",change_dict_local)
	K1SL_sq = np.array(tw40cm1_sq.K1SL)
	K1SL_nonzero_sq = K1SL_sq[K1SL_sq != 0]
	mean_strength_sq = np.mean(abs(K1SL_nonzero_sq))
	set_knobs(change_dict_local,0.,0.)	
	change_value(change_dict_local,"%twiss_pattern","BPM")
	

	#regular knobs
	change_value(change_dict_local,"%knob_Re_type","CMRS.b1")
	change_value(change_dict_local,"%knob_Im_type","CMIS.b1")
	knob_Re_madx, knob_Im_madx = get_madx_knobs(change_dict_local)
	set_knobs(change_dict_local,knob_Re_madx, knob_Im_madx) 
	f_madx =  get_f(response_path,"twiss.original",change_dict_local)
	C_min_madx = get_C_min(change_dict_local)
	set_knobs(change_dict_local,0.,0.)
			
	
	R_inverse = get_responsematrix(change_dict_local)
	tw40cm0 = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm0.Cmatrix()
	knob_Re_res, knob_Im_res = get_response_knobs(R_inverse,change_dict_local)
	set_knobs(change_dict_local,knob_Re_res, knob_Im_res) 
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res = get_C_min(change_dict_local)
	#double correction
	dKnob_Re_res , dKnob_Im_res = get_response_knobs(R_inverse,change_dict_local)
	set_knobs(change_dict_local,knob_Re_res + dKnob_Re_res , knob_Im_res + dKnob_Im_res )
	f_double_res =  get_f(response_path,"twiss.original",change_dict_local)
	C_min_double_res = get_C_min(change_dict_local)
	set_knobs(change_dict_local,0.,0.)
			
	
	change_value(change_dict_local,"%twiss_pattern",".")
	set_knobs(change_dict_local,knob_Re_res, knob_Im_res)
	tw40cm1 = get_twiss(response_path,"twiss.original",change_dict_local)
	K1SL = np.array(tw40cm1.K1SL)
	K1SL_nonzero = K1SL[K1SL != 0]
	mean_strength = np.mean(abs(K1SL_nonzero))
	set_knobs(change_dict_local,0.,0.)
	change_value(change_dict_local,"%twiss_pattern","BPM")


	
	S = np.array(tw40cm0.S)	
	f_0 = np.array(tw40cm0.F1001R) + 1j * np.array(tw40cm0.F1001I)	
	beta_max_error, beta_TEST = get_beta_error(change_dict)[0:2:]
	
	quad_pattern = change_dict["%quad_pattern_1"]
	if change_dict["%quad_pattern_2"] != change_dict["%quad_pattern_1"]:
		quad_pattern += " and " + change_dict["%quad_pattern_2"]  

	fig = plt.figure(figsize=plt.figaspect(0.5))
	fig.suptitle(r"$C_0$ = " + "{:.3e}".format(C_min_0)  + r"	$\beta_{max}$ =" + str(round(beta_max_error,3))+"	quad_section = " + quad_pattern)

	ax1 = fig.add_subplot(1,2,1)
	ax1.plot(S,abs(f_0),label = "before response")
	ax1.plot(S,abs(f_madx),label = "madx matching correction")
	ax1.plot(S,abs(f_res),label = "response matrix")
	ax1.plot(S,abs(f_double_res),label = "response matrix double correction")
	ax1.set_title("regular knobs \n $C_-^{res}$ =" + "{:.2e}".format(C_min_res) +  " $C_-^{Double res}$ =" + "{:.2e}".format(C_min_double_res) +"\n"" $C_-^{match}$ =" + "{:.2e}".format(C_min_madx) + r"<|K1SL|> =" + "{:.2e}".format(mean_strength))
	ax1.set_xlabel("S")
	ax1.set_ylabel("F1001")	
	ax1.legend()	

	ax2 = fig.add_subplot(1,2,2)
	ax2.plot(S,abs(f_0),label = "before response")
	ax2.plot(S,abs(f_madx_sq),label = "madx matching correction")
	ax2.plot(S,abs(f_res_sq),label = "response matrix")
	ax2.plot(S,abs(f_double_res_sq),label = "response matrix double correction")
	ax2.set_title("squeeze knobs \n $C_-^{res sq}$ =" + "{:.2e}".format(C_min_res_sq) +  " $C_-^{Double res sq}$ =" + "{:.2e}".format(C_min_double_res_sq) + "\n" +" $C_-^{match sq}$ =" + "{:.2e}".format(C_min_madx_sq) + r"$<|K1SL|>_{sq}$ =" + "{:.2e}".format(mean_strength_sq))
	ax2.set_xlabel("S")
	ax2.set_ylabel("F1001")
	ax2.legend()
	fig.tight_layout(rect=[0, 0.03, 1, 0.95])	
	plt.savefig('plots/' + savepath)
	#plt.show()
	#plt.plot(S,beta_TEST)
	#plt.show()

def plot_correction_test(changse_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)	
	C_min_0 = get_C_min(change_dict_local)
	tw40cm0 = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm0.Cmatrix()
	S = np.array(tw40cm0.S)
	f_0 = np.array(tw40cm0.F1001R) + 1j * np.array(tw40cm0.F1001I)	
	
	knob_Re_madx, knob_Im_madx = get_madx_knobs(change_dict_local)
	set_knobs(change_dict_local,knob_Re_madx, knob_Im_madx)
	f_madx = get_f(response_path,"twiss.original",change_dict_local)
	C_min_madx = get_C_min(change_dict_local)
	set_knobs(change_dict_local,0., 0.)
	
	R_inverse = get_responsematrix(change_dict_local)
	knob_Re_res, knob_Im_res = get_response_knobs(R_inverse,change_dict_local)
	set_knobs(change_dict_local,knob_Re_res, knob_Im_res)
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res = get_C_min(change_dict_local)
	set_knobs(change_dict_local,0., 0.)

	beta_max_error, beta_TEST = get_beta_error(change_dict)[0:2]
	quad_pattern = change_dict["%quad_pattern_1"]
	if change_dict["%quad_pattern_2"] != change_dict["%quad_pattern_1"]:
		quad_pattern += " and " + change_dict["%quad_pattern_2"]  
	
	fig = plt.figure(figsize=plt.figaspect(1))
	fig.suptitle("$C_0$ = " + "{:.3e}".format(C_min_0) + r" $\beta_{max}$ = " + str(round(beta_max_error,3)) +"	quad_section = " + quad_pattern )
	
	ax1 = fig.add_subplot(1,1,1)
	ax1.set_title("$C_{res}$ = " + "{:.3e}".format(C_min_res) + "	$C_{matching}$ = " + "{:.3e}".format(C_min_madx))
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_madx),label = "f_matching")
	ax1.plot(S,abs(f_res),label = "f_response")
	ax1.set_xlabel("S")
	ax1.set_ylabel("abs(F1001)")
	
	plt.savefig("plots/" + savepath)
	plt.show()
	plt.plot(S,beta_TEST)
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
change_dict["%quad_strength"] = "0.00008"
change_dict["%twiss_pattern"] = "BPM"
change_dict["%colknob1"] = "0."
change_dict["%colknob5"] = "0."



change_dict["%quad_strength"] = "0.00016"
change_dict["%quad_pattern_1"] = "R5"
change_dict["%quad_pattern_2"] = "R5"
change_dict["%error_strength"] = "0.00003*gauss()"
plot_squeeze_comparison(change_dict,"injectionComparison_randomSQ_localQ_gauss3_R5.pdf")
change_dict["%error_strength"] = "0.00002*gauss()"
plot_squeeze_comparison(change_dict,"injectionComparison_randomSQ_localQ_gauss2_R5.pdf")

change_dict["%quad_strength"] = "0.0001"
change_dict["%quad_pattern_1"] = "R8"
change_dict["%quad_pattern_2"] = "R8"
change_dict["%error_strength"] = "0.00003*gauss()"
plot_squeeze_comparison(change_dict,"injectionComparison_randomSQ_localQ_gauss3_R8.pdf")
change_dict["%error_strength"] = "0.00002*gauss()"
plot_squeeze_comparison(change_dict,"injectionComparison_randomSQ_localQ_gauss2_R8.pdf")

change_dict["%quad_strength"] = "0.00005"
change_dict["%quad_pattern_1"] = "R3"
change_dict["%quad_pattern_2"] = "R8"
change_dict["%error_strength"] = "0.00003*gauss()"
plot_squeeze_comparison(change_dict,"injectionComparison_randomSQ_localQ_gauss3_R3R8.pdf")
change_dict["%error_strength"] = "0.00002*gauss()"
plot_squeeze_comparison(change_dict,"injectionComparison_randomSQ_localQ_gauss2_R3R8.pdf")

change_dict["%quad_strength"] = "0."
change_dict["%quad_pattern_1"] = "R3"
change_dict["%quad_pattern_2"] = "R8"
change_dict["%error_strength"] = "0.00003*gauss()"
plot_squeeze_comparison(change_dict,"injectionComparison_randomSQ_zeroQ_gauss3_R3R8.pdf")
change_dict["%error_strength"] = "0.00002*gauss()"
plot_squeeze_comparison(change_dict,"injectionComparison_randomSQ_zeroQ_gauss2_R3R8.pdf")


