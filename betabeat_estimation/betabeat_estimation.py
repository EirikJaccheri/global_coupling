#-*- coding: utf-8 -*-

import sys
import time

sys.path.append("../global_coupling/")
from global_coupling import *
      
modulename = "global_coupling"
if modulename in sys.modules:
    print("You have imported the {} module".format(modulename))

N_error = 12

#home
madx_path = "/home/eirik/madx "
folder_path = "/home/eirik/CERN/global_coupling_correction/betabeat_estimation/"
lhc_path = "/home/eirik/CERN/lhc2018/2018"

#work
#madx_path = "/home/ehoydals/madx "
#folder_path = "/home/ehoydals/global_coupling_correction/analytical_test/"
#lhc_path =  "/afs/cern.ch/eng/lhc/optics/runII/2018"

response_path = folder_path + "response_betabeat_estimation.madx"
C_min_path = folder_path + "exact_C_min_betabeat_estimation.madx"
FineTuneCoupling_path = folder_path + "FineTuneCoupling_betabeat_estimation.madx"

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


def beta_beating_estimate_oneerror(change_dict,savepath,error_name,error_strength):
	change_dict_local = copy.deepcopy(change_dict)
	
	#quad_strength, skew_strength = zero , zero
	change_value(change_dict_local,"%twiss_pattern",".")
	change_value(change_dict_local,"%quad_strength","0.")
	change_value(change_dict_local,"%error_strength","0.")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	name_l = np.array(tw40cm.NAME)
	func_BPM = np.frompyfunc(lambda x: 'BPM' in x, 1,1)
	func_name = np.frompyfunc(lambda x: error_name in x, 1,1)
	error_index = func_name(name_l).astype(bool)
	BPM_index = func_BPM(name_l).astype(bool)
	
	Qx , Qy = tw40cm.Q1 , tw40cm.Q2
	mux0_all , muy0_all = np.array(tw40cm.MUX) , np.array(tw40cm.MUY)
	betx0_all , bety0_all = np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	
	mux0_error , muy0_error = mux0_all[error_index] , muy0_all[error_index]
	mux0_BPM , muy0_BPM = mux0_all[BPM_index] , muy0_all[BPM_index]
	betx0_error , bety0_error = betx0_all[error_index] , bety0_all[error_index]
	
	
	#quad_strength, skew_strength = nozero , zero
	change_value(change_dict_local,"%twiss_pattern",".")
	change_value(change_dict_local,"%quad_strength",change_dict["%quad_strength"])
	change_value(change_dict_local,"%error_strength","0.")
	
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	betx1_all , bety1_all = np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	S_all = np.array(tw40cm.S)
	
	betx1_error , bety1_error = betx1_all[error_index] , bety1_all[error_index]
	S = S_all[BPM_index]
	f_0_Re = np.array(tw40cm.F1001R)[BPM_index]
	f_0_Im = np.array(tw40cm.F1001I)[BPM_index]
	f_0 = f_0_Re + 1j * f_0_Im
	
	#quad_strength, skew_strength = nonzero , nonzero
	change_value(change_dict_local,"%twiss_pattern","BPM")
	change_value(change_dict_local,"%quad_strength",change_dict["%quad_strength"])
	change_value(change_dict_local,"%pattern_1",error_name)
	change_value(change_dict_local,"%pattern_2",error_name)
	change_value(change_dict_local,"%error_strength",str(error_strength))
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	
	J = error_strength
	B_0 = 1 / (4 * (1 - np.exp(2*np.pi*1j*(Qx - Qy)))) * np.exp(2*np.pi*1j*(deltaPhi(mux0_BPM,mux0_error,Qx) - deltaPhi(muy0_BPM,muy0_error,Qy)))
	B_0.real *=-1
	sqrt_beta0 = np.sqrt(betx0_error * bety0_error)
	sqrt_beta1 = np.sqrt(betx1_error * bety1_error)
	N_BPM = float(len(f_0))
	f_model = f_0 + B_0 * sqrt_beta0 * J
	del_sqrt_beta = sqrt_beta1 - sqrt_beta0
	del_sqrt_beta_est = 1 / N_BPM * np.sum(np.real(1 / (B_0 * J) * (f_res - f_model)))
	
	
	#testing tobias formula
	a = 4 * abs(1 - np.exp(2*np.pi*1j*(Qx - Qy)))
	del_beta_tob = a / J * abs(f_res - f_0) - sqrt_beta0 * np.ones(len(f_res))
	
	
	print((1 / (B_0 * J) * (f_res - f_model)))
	print("del_sqrt_beta",del_sqrt_beta[0])
	print("del_sqrt_beta_relative",del_sqrt_beta[0]/sqrt_beta0)
	print("del_sqrt_beta_est",del_sqrt_beta_est)
	print("del_sqrt_beta_est_tob",np.mean(del_beta_tob))
	print("max estimate",max(np.real(1 / (B_0 * J) * (f_res - f_model))))
	print("error_position",S_all[error_index][0])	
	
	
	 
	#test
	fig1 , ax0 = plt.subplots(1)
	ax0.plot(S, np.real(1 / (B_0 * J) * (f_res - f_model)),label = r"$Re\{\frac{f_{res} - f_{model}}{B_0J}\}$")
	ax0.plot(S,del_beta_tob,label = "del_beta_tob")
	ax0.plot(S, del_sqrt_beta * np.ones(len(S)),label = r"$\Delta\sqrt{\beta_{x,s'}\beta_{y,s'}}$")
	ax0.set_xlabel("S")
	ax0.set_ylabel(r"$\sqrt{\beta_{x,s'}\beta_{y,s'}}$")
	ax0.legend()
	fig1.tight_layout()
	
	fig2 , (ax1,ax2,ax3) = plt.subplots(3,1)
	ax1.plot(S,abs(f_0),label = "|f_0|")
	ax1.plot(S,abs(f_model),label = "|f_model|")
	ax1.plot(S,abs(f_res),label ="|f_res|")
	ax1.set_xlabel("S")
	ax1.set_ylabel("|f_1001|")
	ax1.legend()
	
	ax2.plot(S,np.real(f_0),label = "Re{f_0}")
	ax2.plot(S,np.real(f_model),label = "Re{f_model}")
	ax2.plot(S,np.real(f_res),label ="Re{f_res}")
	ax2.set_xlabel("S")
	ax2.set_ylabel("Re{f_1001}")
	ax2.legend()
	
	ax3.plot(S,np.imag(f_0),label = "Im{f_0}")
	ax3.plot(S,np.imag(f_model),label = "Im{f_model}")
	ax3.plot(S,np.imag(f_res),label = "Im{f_res}")
	ax3.set_xlabel("S")
	ax3.set_ylabel("Im{f_1001}")
	ax3.legend()
	fig2.tight_layout()
	
	fig3 , ax4 = plt.subplots(1)
	ax4.plot(S_all,(betx1_all - betx0_all) / betx0_all,label = "betx")
	ax4.plot(S_all,(bety1_all - bety0_all) / bety0_all,label = "bety")
	ax4.plot(S_all,(np.sqrt(betx1_all*bety1_all) - np.sqrt(betx0_all * bety0_all)) / np.sqrt(betx0_all * bety0_all),label ="delta_sqrt/sqrt_0")
	ax4.legend()
	fig3.tight_layout()
	
	plt.savefig("plots/" + savepath)
	plt.show()
	
	
def beta_beating_estimate_oneerror_abs(change_dict,savepath,error_name,error_strength):
	change_dict_local = copy.deepcopy(change_dict)
	
	#quad_strength, skew_strength = zero , zero
	change_value(change_dict_local,"%twiss_pattern",".")
	change_value(change_dict_local,"%quad_strength","0.")
	change_value(change_dict_local,"%error_strength","0.")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	name_l = np.array(tw40cm.NAME)
	func_BPM = np.frompyfunc(lambda x: 'BPM' in x, 1,1)
	func_name = np.frompyfunc(lambda x: error_name in x, 1,1)
	error_index = func_name(name_l).astype(bool)
	BPM_index = func_BPM(name_l).astype(bool)
	
	Qx , Qy = tw40cm.Q1 , tw40cm.Q2
	mux0_all , muy0_all = np.array(tw40cm.MUX) , np.array(tw40cm.MUY)
	betx0_all , bety0_all = np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	
	mux0_error , muy0_error = mux0_all[error_index] , muy0_all[error_index]
	mux0_BPM , muy0_BPM = mux0_all[BPM_index] , muy0_all[BPM_index]
	betx0_error , bety0_error = betx0_all[error_index] , bety0_all[error_index]
	
	
	#quad_strength, skew_strength = nozero , zero
	change_value(change_dict_local,"%twiss_pattern",".")
	change_value(change_dict_local,"%quad_strength",change_dict["%quad_strength"])
	change_value(change_dict_local,"%error_strength","0.")
	
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	betx1_all , bety1_all = np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	S_all = np.array(tw40cm.S)
	
	betx1_error , bety1_error = betx1_all[error_index] , bety1_all[error_index]
	S = S_all[BPM_index]
	f_0_Re = np.array(tw40cm.F1001R)[BPM_index]
	f_0_Im = np.array(tw40cm.F1001I)[BPM_index]
	f_0 = f_0_Re + 1j * f_0_Im
	
	#quad_strength, skew_strength = nonzero , nonzero
	change_value(change_dict_local,"%twiss_pattern","BPM")
	change_value(change_dict_local,"%quad_strength",change_dict["%quad_strength"])
	change_value(change_dict_local,"%pattern_1",error_name)
	change_value(change_dict_local,"%pattern_2",error_name)
	change_value(change_dict_local,"%error_strength",str(error_strength))
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	
	
	
	J = error_strength
	sqrt_beta0 = np.sqrt(betx0_error * bety0_error)
	sqrt_beta1 = np.sqrt(betx1_error * bety1_error)
	del_sqrt_beta = sqrt_beta1 - sqrt_beta0

	a = 4 * abs(1 - np.exp(2*np.pi*1j*(Qx - Qy)))	
	sqrt_beta_est = abs(f_res) * a / J 
	del_sqrt_beta_est = sqrt_beta_est - sqrt_beta0 * np.ones(len(f_res))
	
	print("sqrt_beta0",sqrt_beta0)
	print("sqrt_beta1",sqrt_beta1)
	
	fig , ax = plt.subplots(1)
	ax.plot(S,del_sqrt_beta_est,label = "del_sqrt_beta_est")
	ax.plot(S,del_sqrt_beta*np.ones(len(f_res)),label = "del_sqrt_beta")
	#ax.plot(S,abs(f_res - f_0) * a,label = "a*abs(f_res - f_0)")
	ax.legend()
	plt.show()	
	
def plot_taylor_expansion(change_dict,savepath,error_name,error_strength):
	change_dict_local = copy.deepcopy(change_dict)
	
	#quad_strength, skew_strength = zero , zero
	change_value(change_dict_local,"%twiss_pattern",".")
	change_value(change_dict_local,"%quad_strength","0.")
	change_value(change_dict_local,"%error_strength","0.")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	name_l = np.array(tw40cm.NAME)
	f_BPM = np.frompyfunc(lambda x: 'BPM' in x, 1,1)
	f_name = np.frompyfunc(lambda x: error_name in x, 1,1)
	error_index = f_name(name_l).astype(bool)
	BPM_index = f_BPM(name_l).astype(bool)
	
	Qx , Qy = tw40cm.Q1 , tw40cm.Q2
	mux0_all , muy0_all = np.array(tw40cm.MUX) , np.array(tw40cm.MUY)
	betx0_all , bety0_all = np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	
	mux0_error , muy0_error = mux0_all[error_index] , muy0_all[error_index]
	mux0_BPM , muy0_BPM = mux0_all[BPM_index] , muy0_all[BPM_index]
	betx0_error , bety0_error = betx0_all[error_index] , bety0_all[error_index]
	
	
	#quad_strength, skew_strength = nozero , zero
	change_value(change_dict_local,"%twiss_pattern",".")
	change_value(change_dict_local,"%quad_strength",change_dict["%quad_strength"])
	change_value(change_dict_local,"%error_strength","0.")
	
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	betx1_all , bety1_all = np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	mux1_all , muy1_all = np.array(tw40cm.MUX) , np.array(tw40cm.MUY)
	S_all = np.array(tw40cm.S)
	
	betx1_error , bety1_error = betx1_all[error_index] , bety1_all[error_index]
	mux1_error , muy1_error = mux1_all[error_index] , muy1_all[error_index]
	mux1_BPM , muy1_BPM = mux1_all[BPM_index] , muy1_all[BPM_index]
	S = S_all[BPM_index]
	f_0_Re = np.array(tw40cm.F1001R)[BPM_index]
	f_0_Im = np.array(tw40cm.F1001I)[BPM_index]
	f_0 = f_0_Re + 1j * f_0_Im
	
	#quad_strength, skew_strength = nonzero , nonzero
	change_value(change_dict_local,"%twiss_pattern","BPM")
	change_value(change_dict_local,"%quad_strength",change_dict["%quad_strength"])
	change_value(change_dict_local,"%pattern_1",error_name)
	change_value(change_dict_local,"%pattern_2",error_name)
	change_value(change_dict_local,"%error_strength",str(error_strength))
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	
	J = error_strength
	B_0 = 1 / (4 * (1 - np.exp(2*np.pi*1j*(Qx - Qy)))) * np.exp(2*np.pi*1j*(deltaPhi(mux0_BPM,mux0_error,Qx) - deltaPhi(muy0_BPM,muy0_error,Qy)))
	B_0.real *=-1
	sqrt_beta0 = np.sqrt(betx0_error * bety0_error)
	sqrt_beta1 = np.sqrt(betx1_error * bety1_error)
	del_sqrt_beta = sqrt_beta1 - sqrt_beta0
	N_BPM = float(len(f_0))
	f_model = f_0 + B_0 * sqrt_beta0 * J
	
	delta_f_beta = B_0 * J * del_sqrt_beta
	delta_mux_BPM , delta_muy_BPM =  mux1_BPM - mux0_BPM , muy1_BPM - muy0_BPM
	delta_f_mu_BPM = 1j * B_0 * J * sqrt_beta0 * (delta_mux_BPM - delta_muy_BPM)
	delta_mux_error , delta_muy_error = mux1_error - mux0_error , muy1_error - muy0_error 
	delta_f_mu_error = 1j * B_0 * J * sqrt_beta0 * (-delta_mux_error + delta_muy_error) * np.ones(len(S))
	f_taylor = f_model + delta_f_beta - delta_f_mu_BPM - delta_f_mu_error
	f_taylor_beta = f_model + delta_f_beta
	
	fig0 , ax0 = plt.subplots(1)
	ax0.plot(S,abs(f_0),label = "f_0")
	ax0.plot(S,abs(f_model),label = "f_model")
	ax0.plot(S,abs(f_taylor_beta),label = "f_taylor_beta")
	ax0.plot(S,abs(f_taylor),label = "f_taylor")
	ax0.plot(S,abs(f_res),label = "f_res")
	ax0.set_title("|f1001|")
	ax0.set_xlabel("S")
	ax0.set_ylabel("|f1001|")
	ax0.legend()
	fig0.tight_layout()
	
	
	fig1 , (ax1,ax2) = plt.subplots(2,1)
	ax1.plot(S,np.real(f_0),label = "f_0")
	ax1.plot(S,np.real(f_model),label = "f_model")
	ax1.plot(S,np.real(f_taylor),label = "f_taylor")
	ax1.plot(S,np.real(f_res),label = "f_res")
	ax1.set_title("Re{f1001}")
	ax1.set_xlabel("S")
	ax1.set_ylabel("Re{f1001}")
	ax1.legend()
	
	
	ax2.plot(S,np.imag(f_0),label = "f_0")
	ax2.plot(S,np.imag(f_model),label = "f_model")
	ax2.plot(S,np.imag(f_taylor),label = "f_taylor")
	ax2.plot(S,np.imag(f_res),label = "f_res")
	ax2.set_title("Re{f1001}")
	ax2.set_xlabel("S")
	ax2.set_ylabel("Re{f1001}")
	ax2.legend()
	fig1.tight_layout()
	
	
	delta_f_mu = delta_f_mu_BPM + delta_f_mu_error
	delta_f_mu_exp = f_res - f_model - delta_f_beta
	fig2 , (ax3,ax4) = plt.subplots(2,1)
	ax3.plot(S,np.real(delta_f_beta),label = "df_beta")	
	ax3.plot(S,np.real(delta_f_mu),label = "df_mu")
	ax3.plot(S,np.real(delta_f_mu_exp),label = "df_mu_exp")
	ax3.set_xlabel("S")
	ax3.set_ylabel("Re{f1001}")
	ax3.legend()
	
	ax4.plot(S,np.imag(delta_f_beta),label = "df_beta")	
	ax4.plot(S,np.imag(delta_f_mu),label = "df_mu")
	ax4.plot(S,np.imag(delta_f_mu_exp),label = "df_mu_exp")
	ax4.set_xlabel("S")
	ax4.set_ylabel("Im{f1001}")
	ax4.legend()
	fig2.tight_layout()
	
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



change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%error_strength","0.")
change_value(change_dict_local,"%quad_strength","0.0004")
#beta_beating_estimate_oneerror(change_dict_local,"systematic_beating_impact_MQS.23R1.B1_strength1e-4_R5.pdf","MQS.23R1.B1",0.004)
#plot_taylor_expansion(change_dict_local,"plot_taylor_expansion_MQS.23R1.B1_strength1e-4_R5.pdf","MQS.23R1.B1",0.0001)


change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%error_strength","0.")
change_value(change_dict_local,"%quad_strength","0.0001")
change_value(change_dict_local,"%Qx","62.28")
change_value(change_dict_local,"%Qy","60.31")
#beta_beating_estimate_oneerror_abs(change_dict_local,"systematic_beating_impact_MQS.23R2.B1_strength1e-4_R5.pdf","MQS.23R2.B1",0.0001)
#beta_beating_estimate_oneerror(change_dict_local,"systematic_beating_impact_MQS.23R2.B1_strength1e-4_R5.pdf","MQS.23R2.B1",0.0001)
#plot_taylor_expansion(change_dict_local,"plot_taylor_expansion_MQS.23R2.B1_strength1e-4_R5.pdf","MQS.23R2.B1",0.0001)


change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%error_strength","0.")
change_value(change_dict_local,"%quad_strength","0.0001")
change_value(change_dict_local,"%quad_pattern_1","R8")
change_value(change_dict_local,"%quad_pattern_2","R8")
change_value(change_dict_local,"%Qx","62.28")
change_value(change_dict_local,"%Qy","60.31")
beta_beating_estimate_oneerror(change_dict_local,"systematic_beating_impact_MQS.23R2.B1_strength1e-4_R8.pdf","MQS.23R2.B1",0.00001)
#plot_taylor_expansion(change_dict_local,"plot_taylor_expansion_MQS.23R2.B1_strength1e-4_R5.pdf","MQS.23R2.B1",0.0001)

change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%error_strength","0.")
change_value(change_dict_local,"%quad_strength","0.0005")
change_value(change_dict_local,"%Qx","62.31")
change_value(change_dict_local,"%Qy","60.32")
#beta_beating_estimate_oneerror(change_dict_local,"systematic_beating_impact_MQS.23R1.B1_strength1e-4_R5.pdf","MQS.23R1.B1",0.0001)
#plot_taylor_expansion(change_dict_local,"plot_taylor_expansion_MQS.23R1.B1_strength1e-4_R5_normal_tunes.pdf","MQS.23R1.B1",0.0001)

"""
R_Re , R_Im , B_c , C , B_1_c , B_2_c , KQS_matrix , KQS_index_l , KQS_name_l , tw40cm = get_responsematrix_components(change_dict_local)
"""
"""
/afs/cern.ch/work/t/tpersson/public/ 
example does not constrain the tunes
Q 10
changeparameter_uter
"""










