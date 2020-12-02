import sys
sys.path.append("../global_coupling/")
from global_coupling import *
      
modulename = "global_coupling"
if modulename in sys.modules:
    print("You have imported the {} module".format(modulename))

N_error = 12

#home
madx_path = "/home/eirik/madx "
folder_path = "/home/eirik/CERN/global_coupling_correction/analytical_test/"
lhc_path = "/home/eirik/CERN/lhc2018/2018"

#work
#madx_path = "/home/ehoydals/madx "
#folder_path = "/home/ehoydals/global_coupling_correction/analytical_test/"
#lhc_path =  "/afs/cern.ch/eng/lhc/optics/runII/2018"

response_path = folder_path + "response_beta_beat_impact.madx"
C_min_path = folder_path + "exact_C_min_analytic_test.madx"
#FineTuneCoupling_path = folder_path + "FineTuneCoupling_squeeze.madx"

path_dict = {}
path_dict["madx_path"] = madx_path
path_dict["response_path"] = response_path
path_dict["C_min_path"] = C_min_path
#path_dict["FineTuneCoupling_path"] = FineTuneCoupling_path
set_global_paths(path_dict)



reset_response_dict_new = {"%error_component" : "quadrupole",
	"%error_strength" : "0.",
	"%pattern_1" : ".",
	"%pattern_2" : ".",
	"%quad_component" : "quadrupole",
	"%quad_pattern_1" : ".",
	"%quad_strength" : "0.",
	"%twiss_pattern" : "BPM",
}


reset_analyticresponse_dict_new = {"%knob_Re_value" : "0.", 
"%knob_Im_value" : "0.",
"%error_strength":"0."
}

reset_beta_error_dict_new = {"%error_strength" : "0.", "%quad_strength" : "0."}


set_reset_dict(reset_response_dict,reset_response_dict_new)
set_reset_dict(reset_analyticresponse_dict,reset_analyticresponse_dict_new)
set_reset_dict(reset_beta_error_dict,reset_beta_error_dict_new)




def plot_f1001_comparison_colknob(change_dict,KS_names,KS,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	change_value(change_dict_local,"%colknob1","0.")
	change_value(change_dict_local,"%colknob5","0.")
	change_value(change_dict_local,"%error_strength","0.")
	change_value(change_dict_local,"%quad_strength","0.")

	change_value(change_dict_local,"%twiss_pattern",".")
	
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	name_l = tw40cm.NAME
	betx = np.array(tw40cm.BETX)
	bety = np.array(tw40cm.BETY)
	mux = np.array(tw40cm.MUX)
	muy = np.array(tw40cm.MUY)
	Qx = tw40cm.Q1
	Qy = tw40cm.Q2
	
	n_error = len(KS_names)
	betx_error = np.zeros(n_error)
 	bety_error = np.zeros(n_error)
	mux_error = np.zeros(n_error)
	muy_error = np.zeros(n_error)
	
	KS_index = np.zeros(n_error) 
	for i in range(n_error):
		error_index = name_l.index(KS_names[i])	
		KS_index[i] = error_index
		
		betx_error[i] = betx[error_index] 
	 	bety_error[i] = bety[error_index]
		mux_error[i] = mux[error_index]
		muy_error[i] = muy[error_index]
	
	change_value(change_dict_local,"%twiss_pattern","BPM")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	betx_BPM = np.array(tw40cm.BETX)
	bety_BPM = np.array(tw40cm.BETY)
	mux_BPM = np.array(tw40cm.MUX)
	muy_BPM = np.array(tw40cm.MUY)
	Qx = tw40cm.Q1
	Qy = tw40cm.Q2
	
	l = 0.223
	S = np.array(tw40cm.S)
	f1001_anal = f_1001(mux_BPM,muy_BPM,Qx,Qy,betx_error,bety_error,mux_error,muy_error,KS * l)
	
	change_value(change_dict_local,"%colknob5",str(KS[1] / 1e-4))
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm.Cmatrix()
	f1001_madx = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	
	fig = plt.figure()
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(S,abs(f1001_anal),label = "analytical")
	ax1.plot(S,abs(f1001_madx),label = "madx")
	ax1.set_ylim(bottom=0)
	ax1.set_ylim(top = 2*max(np.max(abs(f1001_madx)),np.max(abs(f1001_anal))))
	ax1.legend()
	plt.savefig("plots/" + savepath)
	plt.show()

def plot_f1001_knob_comparison(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	knob_Re = float(change_dict_local["%knob_Re_value"])
	knob_Im = float(change_dict_local["%knob_Im_value"])
	knob_arr = np.array([knob_Re,knob_Im])
	
	
	KQS_matrix, KQS_index_l, KQS_name_l = get_knob_matrix(change_dict)[0:3]
	
	
	change_value(change_dict_local,"%twiss_pattern","BPM")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm.Cmatrix()
	S = np.array(tw40cm.S)
	f1001_madx = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	C_min_madx = get_C_min(change_dict_local)
	
	change_value(change_dict_local,"%twiss_pattern","BPM")
	change_value(change_dict_local,"%knob_Re_value","0.")
	change_value(change_dict_local,"%knob_Im_value","0.")
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
	KS = np.dot(KQS_matrix,knob_arr) * l
	f1001_anal = f_1001(mux_BPM,muy_BPM,Qx,Qy,betx_error,bety_error,mux_error,muy_error,KS)
	
	
	fig = plt.figure()
	fig.suptitle("C_min_madx = " + str(C_min_madx))
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(S,abs(f1001_anal),label = "analytical")
	ax1.plot(S,abs(f1001_madx),label = "madx")
	#ax1.plot(S,np.real(f1001_anal),label = "analytical real")
	#ax1.plot(S,np.real(f1001_madx),label = "madx real")
	#ax1.plot(S,np.imag(f1001_anal),label = "analytical imag")
	#ax1.plot(S,np.imag(f1001_madx),label = "madx imag")
	ax1.set_ylim(bottom=0)
	ax1.set_ylim(top = 2*max(np.max(abs(f1001_madx)),np.max(abs(f1001_anal))))
	ax1.set_xlabel("S")
	ax1.set_ylabel("abs(f1001)")
	ax1.legend()
	plt.savefig("plots/"+savepath)
	plt.show()
	

def plot_response_matrix_comparison(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	set_knobs(change_dict_local,"0.", "0.") 
	f_0 = get_f(response_path,"twiss.original",change_dict_local)
	C_min_0 = get_C_min(change_dict_local)
	betax_rms , betay_rms , beta_beatx, beta_beaty = get_beta_error_rms(change_dict_local)
	
	change_value(change_dict_local,"%quad_strength","0.")
	R_inverse_analytic =  get_analytic_responsematrix(change_dict_local)
	change_value(change_dict_local,"%quad_strength",change_dict["%quad_strength"])
	knob_Re_analytic_res, knob_Im_analytic_res = get_response_knobs(R_inverse_analytic,change_dict_local)
	set_knobs(change_dict_local,knob_Re_analytic_res, knob_Im_analytic_res) 
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm.Cmatrix()
	S = np.array(tw40cm.S)
	f_analytic_res = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	C_min_analytic_res = get_C_min(change_dict_local)
	
	set_knobs(change_dict_local,"0.", "0.") 
	R_inverse = get_responsematrix(change_dict_local)
	tw40cm0 = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm0.Cmatrix()
	knob_Re_res, knob_Im_res = get_response_knobs(R_inverse,change_dict_local)
	set_knobs(change_dict_local,knob_Re_res, knob_Im_res) 
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res = get_C_min(change_dict_local)

	fig = plt.figure()
	fig.suptitle(r"$\beta_x^{RMS}$ = " + str(round(betax_rms,3)) + r"	$\beta_y^{RMS}$ = " + str(round(betay_rms,3)) + '\n' + r"$C_0$ = " + "{:.3e}".format(C_min_0) + r"	$C_{analytic-res}$ = " + "{:.3e}".format(C_min_analytic_res) + r"	$C_{res}$ = " + "{:.3e}".format(C_min_res))
	
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_analytic_res),label = "f_analytic_res")
	ax1.plot(S,abs(f_res),label = "f_res")
	ax1.set_xlabel("S")
	ax1.set_ylabel("|f1001|")
	ax1.legend()
	
	fig.tight_layout(rect=[0, 0.03, 1, 0.9])
	plt.savefig('plots/' + savepath)
	plt.show()
	
	plt.plot(S,np.real(f_analytic_res))
	plt.plot(S,np.real(f_res))
	plt.plot(S,np.imag(f_analytic_res))
	plt.plot(S,np.imag(f_res))
	plt.show()
	

def plot_response_betabeating_impact(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	S = np.array(tw40cm.S)
	f_0 = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	C_min_0 = get_C_min(change_dict_local)
	beta_max, beta_beat, betax0, betax1, betay0,betay1 = get_beta_error(change_dict_local)
	beta_beatx = (betax0 - betax1)/betax0
	beta_beaty = (betay0 - betay1)/betay0
	betax_rms = np.sqrt(np.mean(np.square(beta_beatx)))
	betay_rms = np.sqrt(np.mean(np.square(beta_beaty)))
	
	
	R_inverse_betabeating = get_analytic_responsematrix(change_dict_local)
	knob_Re_betabeating, knob_Im_betabeating = get_response_knobs(R_inverse_betabeating,change_dict_local)
	set_knobs(change_dict_local,knob_Re_betabeating,knob_Im_betabeating)
	f_res_betabeating = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res_betabeating = get_C_min(change_dict_local)
	
	set_knobs(change_dict_local,0.,0.)
	change_value(change_dict_local,"%quad_strength","0.")
	R_inverse = get_analytic_responsematrix(change_dict_local)
	change_value(change_dict_local,"%quad_strength",change_dict["%quad_strength"])
	knob_Re, knob_Im = get_response_knobs(R_inverse,change_dict_local)
	set_knobs(change_dict_local,knob_Re,knob_Im)
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res = get_C_min(change_dict_local)
	
	fig = plt.figure(figsize = plt.figaspect(0.3))
	fig.suptitle(r"$\beta_{xRMS}$ = " + str(round(betax_rms,3))+ r"	$\beta_{yRMS}$ = " + str(round(betay_rms,3)) + r"	$C_-^0$ = " + "{:.3e}".format(C_min_0) + '\n' + r"$C_-^{res-betabeat}$ = " + "{:.3e}".format(C_min_res_betabeating)+ r"	$C_-^{reS}$ = " + "{:.3e}".format(C_min_res))
	
	ax1 = fig.add_subplot(1,2,1)
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_res),label = "f_analytic_res")
	ax1.plot(S,abs(f_res_betabeating),label = "f_analytic_ res_betabeating")
	ax1.set_xlabel("S")
	ax1.set_ylabel("abs(f1001)")
	ax1.legend()
	
	ax2 = fig.add_subplot(1,2,2)
	ax2.plot(S,(np.sqrt(betax1*betay1) - np.sqrt(betax0*betay0)) / np.sqrt(betax0*betay0),label = "betabeating")
	ax2.set_xlabel("S")
	ax2.set_ylabel(r"$\frac{\Delta \sqrt{\beta_x\beta_y}}{\sqrt{\beta_{x0}\beta_{y0}}}$")
	ax2.legend()
	
	KQS_matrix, KQS_index_l, KQS_name_l, tw40cm = get_knob_matrix(change_dict_local)
	S_all = np.array(tw40cm.S)
	positions = np.take(S_all,KQS_index_l)
	knobs_res_betabeating = np.array([knob_Re_betabeating, knob_Im_betabeating])
	knobs_res = np.array([knob_Re,knob_Im])
	strengths_betabeating = np.dot(KQS_matrix,knobs_res_betabeating)
	strengths = np.dot(KQS_matrix,knobs_res)
	
	
	ax3 = ax2.twinx()
	w = 50
	max_value = max(np.max(abs(strengths_betabeating)),np.max(abs(strengths)))
	ax3.bar(positions - w / 2 * np.ones(len(positions)),strengths_betabeating,width = w,color = "red",label = "strengths betabeating")
	ax3.bar(positions + w / 2 * np.ones(len(positions)),strengths,width = w,color = "orange",label = "strengths")
	ax3.set_ylim(bottom=-1.5*max_value)
	ax3.set_ylim(top = 1.5*max_value)
	ax3.set_ylabel("strength of skew corrector")
	ax3.legend()
	
	fig.tight_layout(rect=[0, 0.03, 1, 0.9])
	plt.savefig("plots/" + savepath)
	#plt.show()
	
	plt.plot(S,beta_beatx)
	#plt.show()
	
	
def tune_scaling_correction(change_dict,Qx_matrix,Qy_matrix,Qx_correction,Qy_correction,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	#getting response knobs for matrix achieved at different tunes then the correction
	change_value(change_dict_local,"%Qx",str(Qx_matrix))
	change_value(change_dict_local,"%Qy",str(Qy_matrix))
	alpha = (1 - np.exp(2*np.pi*1j*(Qx_matrix - Qy_matrix))) / (1 - np.exp(2*np.pi*1j*(Qx_correction - Qy_correction)))
	R_inverse_scaled = get_analytic_responsematrix(change_dict_local,B_scaling = alpha)
	
	change_value(change_dict_local,"%Qx",str(Qx_correction))
	change_value(change_dict_local,"%Qy",str(Qy_correction))
	knob_Re_scaled, knob_Im_scaled = get_response_knobs(R_inverse_scaled,change_dict_local)
	
	#getting responseknobs for matrix achieved at same tunes
	R_inverse = get_analytic_responsematrix(change_dict_local)
	knob_Re, knob_Im = get_response_knobs(R_inverse,change_dict_local)
	
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	S = np.array(tw40cm.S)
	f_0 = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	C_min_0 = get_C_min(change_dict_local)
	
	set_knobs(change_dict_local,knob_Re,knob_Im)
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res = get_C_min(change_dict_local)
	
	set_knobs(change_dict_local,knob_Re_scaled,knob_Im_scaled)
	f_res_scaled = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res_scaled = get_C_min(change_dict_local)
	
	fig = plt.figure()
	fig.suptitle(r"$C_-^0$ = " + "{:.3e}".format(C_min_0) + '\n' + r"$C_-^{res}$ = " + "{:.3e}".format(C_min_res) +  r"	$C_-^{res-scaled}$ = " + "{:.3e}".format(C_min_res_scaled))
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_res),label = "f_res (Qx:" + str(Qx_correction) + " Qy:" + str(Qy_correction)+")")
	ax1.plot(S,abs(f_res_scaled),label = "f_res_scaled (Qx:" + str(Qx_matrix) + " Qy:" + str(Qy_matrix)+")")
	ax1.set_xlabel("S")
	ax1.set_ylabel("|f1001|")
	ax1.legend()
	
	plt.savefig("plots/" + savepath)
	plt.show()

def test(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	S, f_0 = get_S_f(response_path,"twiss.original",change_dict_local)
	C_min_0 = get_C_min(change_dict_local)
	
	
	change_value(change_dict_local,"%quad_strength","0.")
	R_inverse =  get_analytic_responsematrix(change_dict_local)
	change_value(change_dict_local,"%quad_strength",change_dict["%quad_strength"])
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
change_dict["%quad_strength"] = "0.00018"
change_dict["%twiss_pattern"] = "BPM"
change_dict["%colknob1"] = "0."
change_dict["%colknob5"] = "0."
change_dict["%Qx"] = "62.31"
change_dict["%Qy"] = "60.32"



change_dict_local = copy.deepcopy(change_dict)
#change_value(change_dict_local,"%quad_strength","0.")
#test_absolute_response(change_dict,"test.pdf")
#test(change_dict_local,"test2_normal_tunes.pdf")
change_value(change_dict_local,"%quad_pattern_1","R5")
change_value(change_dict_local,"%quad_pattern_2","R5")
change_value(change_dict_local,"%Qx","62.30")
change_value(change_dict_local,"%Qy","60.33")
#test(change_dict_local,"test2_tuneshift003.pdf")


change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%error_strength","0.")
change_value(change_dict_local,"%quad_strength","0.")
change_value(change_dict_local,"%colknob1","0.0001")
KS_names = np.array(["MQSX.3L5","MQSX.3R5"])
KS = np.array([-0.001,0.001])
#plot_f1001_comparison_colknob(change_dict_local,KS_names,KS,"colknob5_tunesplit001.pdf")
change_value(change_dict_local,"%Qx","62.29")
change_value(change_dict_local,"%Qy","60.34")
#plot_f1001_comparison_colknob(change_dict_local,KS_names,KS,"colknob5_tunesplit005.pdf")


change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%knob_Re_type","CMRS.b1_sq")
change_value(change_dict_local,"%knob_Im_type","CMIS.b1_sq")
change_value(change_dict_local,"%error_strength","0.")
change_value(change_dict_local,"%quad_strength","0.")
change_value(change_dict_local,"%knob_Re_value","0.001")
change_value(change_dict_local,"%knob_Im_value","0.001")
#plot_f1001_knob_comparison(change_dict_local,"knob_comparison0.001_tunesplit001.pdf")
change_value(change_dict_local,"%knob_Re_value","0.005")
change_value(change_dict_local,"%knob_Im_value","0.005")
#plot_f1001_knob_comparison(change_dict_local,"knob_comparison0.005_tunesplit001.pdf")
change_value(change_dict_local,"%Qx","62.29")
change_value(change_dict_local,"%Qy","60.34")
change_value(change_dict_local,"%knob_Re_value","0.001")
change_value(change_dict_local,"%knob_Im_value","0.001")
#plot_f1001_knob_comparison(change_dict_local,"knob_comparison0.001_tunesplit005.pdf")
change_value(change_dict_local,"%knob_Re_value","0.005")
change_value(change_dict_local,"%knob_Im_value","0.005")
#plot_f1001_knob_comparison(change_dict_local,"knob_comparison0.005_tunesplit005.pdf")


#have to go through and check these functions
change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%knob_Re_type","CMRS.b1_sq")
change_value(change_dict_local,"%knob_Im_type","CMIS.b1_sq")
change_value(change_dict_local,"%error_strength","0.00002*gauss()")
change_value(change_dict_local,"%quad_strength","0.00018")
change_value(change_dict_local,"%knob_Re_value","0.")
change_value(change_dict_local,"%knob_Im_value","0.")
change_value(change_dict_local,"%quad_pattern_1","R3")
change_value(change_dict_local,"%quad_pattern_2","R3")
#plot_response_matrix_comparison(change_dict_local,"response_matrix_comparison_R3_localQ.pdf")
change_value(change_dict_local,"%quad_pattern_1","R5")
change_value(change_dict_local,"%quad_pattern_2","R5")
change_value(change_dict_local,"%error_strength","0.00003*gauss()")
#plot_response_matrix_comparison(change_dict_local,"response_matrix_comparison_R5_localQ.pdf")
change_value(change_dict_local,"%error_strength","0.00002*gauss()")
change_value(change_dict_local,"%quad_strength","0.")
change_value(change_dict_local,"%quad_pattern_1","R3")
change_value(change_dict_local,"%quad_pattern_2","R3")
#plot_response_matrix_comparison(change_dict_local,"response_matrix_comparison_R3_zeroQ.pdf")
change_value(change_dict_local,"%error_strength","0.00003*gauss()")
change_value(change_dict_local,"%quad_pattern_1","R5")
change_value(change_dict_local,"%quad_pattern_2","R5")
#plot_response_matrix_comparison(change_dict_local,"response_matrix_comparison_R5_zeroQ.pdf")

change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%knob_Re_type","CMRS.b1_sq")
change_value(change_dict_local,"%knob_Im_type","CMIS.b1_sq")
change_value(change_dict_local,"%error_strength","0.00003*gauss()")
change_value(change_dict_local,"%quad_strength","0.00020")
change_value(change_dict_local,"%knob_Re_value","0.")
change_value(change_dict_local,"%knob_Im_value","0.")
change_value(change_dict_local,"%quad_pattern_1","R2")
change_value(change_dict_local,"%quad_pattern_2","R2")
#plot_response_betabeating_impact(change_dict_local,"betabeat_impact_R2.pdf")
change_value(change_dict_local,"%quad_pattern_1","R8")
change_value(change_dict_local,"%quad_pattern_2","R8")
#plot_response_betabeating_impact(change_dict_local,"betabeat_impact_R8.pdf")
change_value(change_dict_local,"%knob_Re_type","CMRS.b1")
change_value(change_dict_local,"%knob_Im_type","CMIS.b1")
change_value(change_dict_local,"%quad_pattern_1","R2")
change_value(change_dict_local,"%quad_pattern_2","R2")
#plot_response_betabeating_impact(change_dict_local,"betabeat_impact_R2_regular_knobs.pdf")
#change_value(change_dict_local,"%quad_strength","0.0001")
change_value(change_dict_local,"%quad_pattern_1","R8")
change_value(change_dict_local,"%quad_pattern_2","R8")
#plot_response_betabeating_impact(change_dict_local,"betabeat_impact_R8_regular_knobs.pdf")



change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%knob_Re_type","CMRS.b1_sq")
change_value(change_dict_local,"%knob_Im_type","CMIS.b1_sq")
change_value(change_dict_local,"%error_strength","0.00003*gauss()")
change_value(change_dict_local,"%quad_strength","0.00018")
change_value(change_dict_local,"%knob_Re_value","0.")
change_value(change_dict_local,"%knob_Im_value","0.")

Qx_matrix , Qy_matrix = 62.30 , 60.33
Qx_correction , Qy_correction = 62.31 , 60.32
#tune_scaling_correction(change_dict_local,Qx_matrix,Qy_matrix,Qx_correction,Qy_correction,"Rescaling_correction_matrix03_correction01.pdf")

Qx_matrix , Qy_matrix = 62.29 , 60.34
Qx_correction , Qy_correction = 62.31 , 60.32
#tune_scaling_correction(change_dict_local,Qx_matrix,Qy_matrix,Qx_correction,Qy_correction,"Rescaling_correction_matrix05_correction01.pdf")

change_value(change_dict_local,"%quad_strength","0.")
change_value(change_dict_local,"%error_strength","0.000018*gauss()")
Qx_matrix , Qy_matrix = 62.3125 , 60.3175
Qx_correction , Qy_correction = 62.31 , 60.32
#tune_scaling_correction(change_dict_local,Qx_matrix,Qy_matrix,Qx_correction,Qy_correction,"Rescaling_correction_matrix005_correction01_zeroQ.pdf")

