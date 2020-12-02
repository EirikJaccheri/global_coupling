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
	change_value(change_dict_local,"%error_strength","0.")
	tw40cm_after_BPM = get_twiss(response_path,"twiss.original",change_dict_local)
	S = np.array(tw40cm_after_BPM.S)
	f_0 = np.array(tw40cm_after_BPM.F1001R) + 1j * np.array(tw40cm_after_BPM.F1001I)
	mux_after_BPM , muy_after_BPM = np.array(tw40cm_after_BPM.MUX) , np.array(tw40cm_after_BPM.MUY)
	
	change_value(change_dict_local,"%twiss_pattern",".")
	#getting lattice parameters with betabeating
	tw40cm_after = get_twiss(response_path,"twiss.original",change_dict_local)
	S_all = np.array(tw40cm_after.S)
	betx_after , bety_after = np.array(tw40cm_after.BETX) , np.array(tw40cm_after.BETY)
	mux_after , muy_after = np.array(tw40cm_after.MUX) , np.array(tw40cm_after.MUY)
	
	change_value(change_dict_local,"%quad_strength","0.")
	#getting lattice parameters without betabeating
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	Qx , Qy = tw40cm.Q1 , tw40cm.Q2
	betx , bety = np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	mux , muy = np.array(tw40cm.MUX) , np.array(tw40cm.MUY)
	name_l = tw40cm.NAME
	
	error_index = name_l.index(error_name)
	S_error = S_all[error_index]
	betx_error , bety_error = betx[error_index] , bety[error_index]
	mux_error , muy_error = mux[error_index] , muy[error_index]
	betx_after_error , bety_after_error = betx_after[error_index] , bety_after[error_index]
	mux_after_error , muy_after_error = mux_after[error_index] , muy_after[error_index]
	
	del_sqrt_beta = np.sqrt(betx_after_error * bety_after_error) - np.sqrt(betx_error * bety_error)
	change_value(change_dict_local,"%twiss_pattern","BPM")
	
	tw40cm_BPM = get_twiss(response_path,"twiss.original",change_dict_local)
	mux_BPM , muy_BPM = np.array(tw40cm_BPM.MUX) , np.array(tw40cm_BPM.MUY)
	
	#f_exp with only model values
	B = 1 / (4*(1 - np.exp(2*np.pi*1j*(Qx - Qy)))) * np.exp(2*np.pi*1j*(deltaPhi(mux_BPM,mux_error,Qx) - deltaPhi(muy_BPM,muy_error,Qy))) 
	B.real *= -1 #becouse of sign connvention
	del_f = np.sqrt(betx_error * bety_error) * error_strength * B
	f_exp =  f_0 + del_f
	
	#f_exp with correct beta at corrector
	del_f_beta = np.sqrt(betx_after_error * bety_after_error) * error_strength * B
	f_exp_beta = f_0 + del_f_beta
	
	#f_exp with correct beta at corrector and correct phaseadvance at corrector
	B_mu = 1 / (4*(1 - np.exp(2*np.pi*1j*(Qx - Qy)))) * np.exp(2*np.pi*1j*(deltaPhi(mux_BPM,mux_after_error,Qx) - deltaPhi(muy_BPM,muy_after_error,Qy))) #husk fortegnkonvensjon
	B_mu.real *= -1 #because of sign_convention
	del_f_beta_mu =   np.sqrt(betx_after_error * bety_after_error) * error_strength * B_mu
	f_exp_beta_mu = f_0 + del_f_beta_mu
	
	#f_exp with correct beta beat at corrector and correct phaseadvance at all BPMs
	B_mu_all = 1 / (4*(1 - np.exp(2*np.pi*1j*(Qx - Qy)))) * np.exp(2*np.pi*1j*(deltaPhi(mux_after_BPM,mux_after_error,Qx) - deltaPhi(muy_after_BPM,muy_after_error,Qy)))
	B_mu_all.real *= -1 #because of sign_convention
	del_f_beta_mu_all =  np.sqrt(betx_after_error * bety_after_error) * error_strength * B_mu_all
	f_exp_beta_mu_all = f_0 + del_f_beta_mu_all
	
	change_value(change_dict_local,"%pattern_1",error_name)
	change_value(change_dict_local,"%pattern_2",error_name)
	change_value(change_dict_local,"%error_strength",str(error_strength))
	change_value(change_dict_local,"%quad_strength",change_dict["%quad_strength"])
	
	
	f_after = get_f(response_path,"twiss.original",change_dict_local)
	
	B_inverse = (1 / B).T / len(B)
	del_sqrt_beta_est_c = 1 / error_strength * np.dot(B_inverse,f_after - f_exp)
	
	#calculating del_sqrt with (f_Re,f_Im) format
	B_l = np.concatenate((B.real,B.imag))
	B_l_inverse = (1 / B_l).T / len(B_l)
	f_after_l = np.concatenate((f_after.real,f_after.imag))
	f_exp_l = np.concatenate((f_exp.real,f_exp.imag))
	del_sqrt_beta_est = 1 / error_strength * np.dot(B_l_inverse,f_after_l - f_exp_l)
	
	
	print(((f_after - f_exp) / B / error_strength).real ) 
	print(np.mean(((f_after - f_exp) / B / error_strength).real))
	
	print("del_sqrt_beta_est")
	print(del_sqrt_beta_est)
	print("del_sqrt_beta_est_c")
	print(del_sqrt_beta_est_c)
	print("del_sqrt_beta")
	print(del_sqrt_beta)
	print("del_sqrt_beta /sqrt_beta_0") 
	print(del_sqrt_beta / np.sqrt(betx_error * bety_error))
	print("mu difference")
	print((mux_after_error - mux_error)/mux_error)
	print("max mu difference")
	print(max((mux_after[np.nonzero(mux_after)] - mux[np.nonzero(mux)])/mux[np.nonzero(mux)]))
	print("max {sqrt(betxbety) * error_strength * B * (delPhi_x - delPhi_y)}")
	print(max(np.sqrt(betx_error * bety_error) * error_strength * B * (deltaPhi(mux_BPM,mux_error,Qx) - deltaPhi(muy_BPM,muy_error,Qy))))
	print("max {B * error_strength}")
	print(max(B*error_strength))
	
	#test comparing complex and concatenate
	del_f_test = 12 * error_strength * B
	f_test = f_after + del_f_test
	f_test_l = np.concatenate((f_test.real,f_test.imag))
	print("complex")
	print(1 / error_strength * np.dot(B_inverse,f_after - f_test))
	print("concatenate")
	print(1 / error_strength * np.dot(B_l_inverse,f_after_l - f_test_l))
	
	
	fig, ax1 = plt.subplots(1)
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_exp),label = "f_exp")
	ax1.plot(S,abs(f_exp_beta),label = "f_exp_beta")
	ax1.plot(S,abs(f_exp_beta_mu),label = "f_exp_beta_mu")
	ax1.plot(S,abs(f_exp_beta_mu_all),label = "f_exp_beta_mu_all")
	ax1.plot(S,abs(f_after),label = "f_after")
	ax1.set_xlabel("S")
	ax1.set_ylabel("|f_{1001}|")
	ax1.legend()
	plt.savefig("plots/" + savepath)
	plt.show()
	
	plt.plot(S,(mux_after_BPM - mux_BPM) / mux_BPM)
	plt.plot(S,(muy_after_BPM - muy_BPM) / muy_BPM)
	plt.show()
	
	
	#test
	fig , (ax1,ax2) = plt.subplots(2,1)
	ax1.plot(S,(f_exp - f_after).real,label = "f_exp - f_after")
	ax1.plot(S,(f_exp_beta - f_after).real,label = "f_exp_beta - f_after")
	ax1.plot(S,(f_exp_beta_mu - f_after).real,label = "f_exp_beta_mu - f_after")
	ax1.plot(S,(f_exp_beta_mu_all - f_after).real,label = "f_exp_beta_mu_all - f_after")
	ax1.plot(S,del_f_test.real,label = "del_f_test")
	ax1.plot(S_error,0.,'x',label = "skew location")
	ax1.legend()
	
	ax2.plot(S,(f_exp - f_after).imag,label = "f_exp - f_after")
	ax2.plot(S,(f_exp_beta - f_after).imag,label = "f_exp_beta - f_after")
	ax2.plot(S,(f_exp_beta_mu - f_after).imag,label = "f_exp_beta_mu - f_after")
	ax2.plot(S,(f_exp_beta_mu_all - f_after).imag,label = "f_exp_beta_mu_all - f_after")
	ax2.plot(S,del_f_test.imag,label = "del_f_test")
	ax2.legend()
	plt.show()


def beta_beating_estimate(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)	
	
	S , f_0 = get_S_f(response_path,"twiss.original",change_dict_local)
	C_min_0 = get_C_min(change_dict_local)
	
	
	change_dict_local = copy.deepcopy(change_dict)
	R_Re , R_Im , B_c , C , B_1_c , B_2_c , KQS_matrix , KQS_index_l , KQS_name_l , tw40cm = get_responsematrix_components(change_dict_local)
	nonzero_index = np.nonzero(KQS_matrix.T[0])
	zero_index = np.nonzero(KQS_matrix.T[0] == 0)
	KQS_matrix_nonzero = np.delete(KQS_matrix,zero_index,axis = 0)
	KQS_index_l_nonzero = KQS_index_l[nonzero_index]
	
	B_c_nonzero = np.delete(B_c,zero_index,axis = 1)
	B_1_c_nonzero = np.delete(B_1_c,zero_index,axis = 1)
	B_2_c_nonzero = np.delete(B_2_c,zero_index,axis = 1)
	B_bar_c = -np.conjugate(C * B_2_c_nonzero) #-* is done to go from madx to franchi sign convention
	
	B_bar = np.block([[np.real(B_bar_c)],[np.imag(B_bar_c)]])
	B_bar_inverse = np.linalg.pinv(B_bar)
	R_z = R_Re + 1j * R_Im
	R = np.block([[R_Re],[R_Im]])
	R_inverse = np.linalg.pinv(R)
	

	#getting betafunctions at nonzero corrector magnets
	betx0 , bety0 = np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	betx0_error , bety0_error = np.take(betx0,KQS_index_l_nonzero) , np.take(bety0,KQS_index_l_nonzero)
	
	
	#calculating response
	knob_Re_res , knob_Im_res = get_response_knobs(R_inverse,change_dict_local)
	knobs_res = np.array([knob_Re_res,knob_Im_res])
	set_knobs(change_dict_local,knob_Re_res,knob_Im_res)
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	
	#finding betabeating estimate
	l = 0.32
	f_predicted = f_0 + np.dot(R_z,knobs_res)
	f_res_l = np.concatenate((np.real(f_res),np.imag(f_res)))
	f_predicted_l = np.concatenate((np.real(f_predicted),np.imag(f_predicted)))
	J = np.dot(KQS_matrix_nonzero,knobs_res) * l
	delta_sqrt_beta_est = 1 / J * np.dot(B_bar_inverse,(f_res_l - f_predicted_l))
		
	#finding true betabeating
	change_value(change_dict_local,"%twiss_pattern",".")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	S_all = np.array(tw40cm.S)
	betx1 , bety1 = np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	betx1_error , bety1_error = np.take(betx1,KQS_index_l_nonzero) , np.take(bety1,KQS_index_l_nonzero)
	delta_sqrt_beta = np.sqrt(betx1_error * bety1_error) - np.sqrt(betx0_error * bety0_error)
	
	delta_sqrt_beta_error = (delta_sqrt_beta_est - delta_sqrt_beta) / delta_sqrt_beta
	print(delta_sqrt_beta_est)
	print(delta_sqrt_beta)
	print(delta_sqrt_beta_error)
		
	fig, ax1 = plt.subplots(1)
	fig.suptitle("$C^-_0 =$" + "{:.3e}".format(C_min_0))		
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_res),label = "f_res")
	ax1.plot(S,abs(f_predicted), label = "f_predicted")
	ax1.legend()
	plt.savefig("plots/" + savepath)
	plt.show()
	
	plt.plot(S_all,np.sqrt(betx1*bety1) - np.sqrt(betx0 * bety0))
	plt.show()
	
	
def plot_taylor_expansion(change_dict,savepath,error_name,error_strength):
	change_dict_local = copy.deepcopy(change_dict)
	change_value(change_dict_local,"%error_strength","0.")
	tw40cm_after_BPM = get_twiss(response_path,"twiss.original",change_dict_local)
	S = np.array(tw40cm_after_BPM.S)
	f_0 = np.array(tw40cm_after_BPM.F1001R) + 1j * np.array(tw40cm_after_BPM.F1001I)
	mux_after_BPM , muy_after_BPM = np.array(tw40cm_after_BPM.MUX) , np.array(tw40cm_after_BPM.MUY)
	
	change_value(change_dict_local,"%twiss_pattern",".")
	#getting lattice parameters with betabeating
	tw40cm_after = get_twiss(response_path,"twiss.original",change_dict_local)
	S_all = np.array(tw40cm_after.S)
	betx_after , bety_after = np.array(tw40cm_after.BETX) , np.array(tw40cm_after.BETY)
	mux_after , muy_after = np.array(tw40cm_after.MUX) , np.array(tw40cm_after.MUY)
	
	change_value(change_dict_local,"%quad_strength","0.")
	#getting lattice parameters without betabeating
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	Qx , Qy = tw40cm.Q1 , tw40cm.Q2
	betx , bety = np.array(tw40cm.BETX) , np.array(tw40cm.BETY)
	mux , muy = np.array(tw40cm.MUX) , np.array(tw40cm.MUY)
	name_l = tw40cm.NAME
	
	error_index = name_l.index(error_name)
	S_error = S_all[error_index]
	betx_error , bety_error = betx[error_index] , bety[error_index]
	mux_error , muy_error = mux[error_index] , muy[error_index]
	betx_after_error , bety_after_error = betx_after[error_index] , bety_after[error_index]
	mux_after_error , muy_after_error = mux_after[error_index] , muy_after[error_index]
	
	del_sqrt_beta = np.sqrt(betx_after_error * bety_after_error) - np.sqrt(betx_error * bety_error)
	change_value(change_dict_local,"%twiss_pattern","BPM")
	
	tw40cm_BPM = get_twiss(response_path,"twiss.original",change_dict_local)
	mux_BPM , muy_BPM = np.array(tw40cm_BPM.MUX) , np.array(tw40cm_BPM.MUY)
	
	#f_exp with only model values
	B = 1 / (4*(1 - np.exp(2*np.pi*1j*(Qx - Qy)))) * np.exp(2*np.pi*1j*(deltaPhi(mux_BPM,mux_error,Qx) - deltaPhi(muy_BPM,muy_error,Qy))) 
	B.real *= -1 #becouse of sign connvention
	del_f = np.sqrt(betx_error * bety_error) * error_strength * B
	f_exp =  f_0 + del_f
	
	#f_taylor
	df_beta = B * del_sqrt_beta * error_strength
	dphi_sx , dphi_sy  = mux_after_BPM - mux_BPM , muy_after_BPM - muy_BPM 
	dphi_sx_error , dphi_sy_error = (mux_after_error - mux_error) * np.ones(len(S)) , (muy_after_error - muy_error) * np.ones(len(S))
	df_phi = 2 * np.pi * 1j * np.sqrt(betx_error * bety_error) * error_strength * B * (deltaPhi(mux_BPM,mux_error,Qx) - deltaPhi(muy_BPM,muy_error,Qy)) * (dphi_sx - dphi_sy - dphi_sx_error + dphi_sy_error)
	f_taylor = f_exp + df_beta - df_phi

	
	change_value(change_dict_local,"%pattern_1",error_name)
	change_value(change_dict_local,"%pattern_2",error_name)
	change_value(change_dict_local,"%error_strength",str(error_strength))
	change_value(change_dict_local,"%quad_strength",change_dict["%quad_strength"])
	
	
	f_after = get_f(response_path,"twiss.original",change_dict_local)
	
	#estimating df_phi
	df_phi_est = f_after - f_exp - df_beta
	
	print(np.max(dphi_sx))
	print(np.mean(dphi_sx))
	print(np.max(dphi_sy))
	print(np.mean(dphi_sy))
	

	fig , ax1 = plt.subplots(1)
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_exp),label = "f_exp")
	ax1.plot(S,abs(f_taylor),label = "f_taylor")
	ax1.plot(S,abs(f_after),label = "f_after")
	ax1.legend()
	plt.show()
	
	fig , ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
	ax1.set_title("Re{f}")
	ax1.plot(S,f_0.real,label = "f_0")
	ax1.plot(S,f_exp.real,label = "f_exp")
	ax1.plot(S,f_taylor.real,label = "f_taylor")
	ax1.plot(S,f_after.real,label = "f_after")
	ax1.legend()
	
	ax2.set_title("Im{f}")
	ax2.plot(S,f_0.imag,label = "f_0")
	ax2.plot(S,f_exp.imag,label = "f_exp")
	ax2.plot(S,f_taylor.imag,label = "f_taylor")
	ax2.plot(S,f_after.imag,label = "f_after")
	ax2.legend()
	
	#fig , (ax1,ax2) = plt.subplots(2,1)
	ax3.set_title("Re{f}")
	ax3.plot(S,df_beta.real,label = "df_beta")
	ax3.plot(S,-df_phi.real,label = "df_phi")
	ax3.plot(S,df_phi_est.real,label = "df_phi_est")
	ax3.legend()
	
	ax4.set_title("Im{f}")
	ax4.plot(S,df_beta.imag,label = "df_beta")
	ax4.plot(S,-df_phi.imag,label = "df_phi")
	ax4.plot(S,df_phi_est.imag,label = "df_phi_est")
	ax4.legend()
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
change_value(change_dict_local,"%error_strength","0.00003*gauss()")
#change_value(change_dict_local,"%quad_strength","0.00018")
#beta_beating_estimate(change_dict_local,"beta_beating_estimate.pdf")

change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%error_strength","0.")
change_value(change_dict_local,"%quad_strength","0.00018")
#beta_beating_estimate_oneerror(change_dict_local,"systematic_beating_impact_MQS.23R1.B1_strength1e-4_R5.pdf","MQS.23R1.B1",0.0001)
#plot_taylor_expansion(change_dict_local,"plot_taylor_expansion_MQS.23R1.B1_strength1e-4_R5.pdf","MQS.23R1.B1",0.0001)


change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%error_strength","0.")
change_value(change_dict_local,"%quad_strength","0.00015")
change_value(change_dict_local,"%quad_pattern_1","R3")
change_value(change_dict_local,"%quad_pattern_2","R3")
beta_beating_estimate_oneerror(change_dict_local,"systematic_beating_impact_MQS.23R1.B1_strength1e-4_R8.pdf","MQS.23R2.B1",0.0001)
#plot_taylor_expansion(change_dict_local,"plot_taylor_expansion_MQS.23R1.B1_strength1e-4_R5.pdf","MQS.23R1.B1",0.0001)
