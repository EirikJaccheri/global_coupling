import sys
sys.path.append("../global_coupling/")
from global_coupling import *
      
modulename = "global_coupling"
if modulename in sys.modules:
    print("You have imported the {} module".format(modulename))

N_error = 12

#home
madx_path = "/home/eirik/madx "
folder_path = "/home/eirik/CERN/global_coupling_correction/colknob_correction/"
lhc_path = "/home/eirik/CERN/lhc2018/2018"

#work
#madx_path = "/home/ehoydals/madx "
#folder_path = "/home/ehoydals/global_coupling_correction/analytical_test/"
#lhc_path =  "/afs/cern.ch/eng/lhc/optics/runII/2018"

response_path = folder_path + "response_colknob_correction.madx"
C_min_path = folder_path + "exact_C_min_colknob_correction.madx"
FineTuneCoupling_path = folder_path + "FineTuneCoupling_colknob_correction.madx"

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


reset_beta_error_dict_new = {
	"%error_strength" : "0.", 
	"%quad_strength" : "0.",
	"%colknob1" : "0.",
	"%colknob5" : "0."
}


set_reset_dict(reset_response_dict,reset_response_dict_new)
set_reset_dict(reset_beta_error_dict,reset_beta_error_dict_new)


def plot_colknob_correction(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	S , f_0 = get_S_f(response_path,"twiss.original",change_dict_local)
	betax_rms , betay_rms , beta_beatx, beta_beaty = get_beta_error_rms(change_dict_local)
	C_min_0 = get_C_min(change_dict_local)
	
	R_inverse = get_responsematrix(change_dict_local)
	knob_Re_res,knob_Im_res = get_response_knobs(R_inverse,change_dict_local)
	knob_Re , knob_Im = float(change_dict_local["%knob_Re_value"]) + knob_Re_res , float(change_dict_local["%knob_Im_value"])  + knob_Im_res
	set_knobs(change_dict_local,knob_Re,knob_Im)
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res = get_C_min(change_dict_local)
	
	fig = plt.figure(figsize = plt.figaspect(0.5))
	fig.suptitle(r"$C_-^0$ = " + "{:.3e}".format(C_min_0) + r"	$C_-^{res}$ = " + "{:.3e}".format(C_min_res) + "	colknob1 = " + change_dict["%colknob1"] + "	colknob5 = "+ change_dict_local["%colknob5"] +'\n' + r"$\beta_{x rms}$ = " + str(round(betax_rms,2))+ r"	$\beta_{y rms}$ = " + str(round(betay_rms,2)))
	ax1 = fig.add_subplot(1,2,1)
	ax1.plot(S,abs(f_0),label = "|f_0|")
	ax1.plot(S,abs(f_res),label = "|f_res|")
	ax1.set_xlabel("S")
	ax1.set_ylabel("|f1001|")
	ax1.legend()
	
	ax2 = fig.add_subplot(1,2,2)
	ax2.plot(S,beta_beatx,label = "beta_beatx")
	ax2.plot(S,beta_beaty,label = "beta_beaty")
	ax2.set_xlabel("S")
	ax2.set_ylabel(r"$\frac{\Delta\beta}{\beta}$")
	ax2.legend()
	
	fig.tight_layout(rect=[0, 0.03, 1, 0.9])
	plt.savefig("plots/" + savepath)
	plt.show()
	

change_dict = {}
change_dict["%lhc_path"] = lhc_path
change_dict["%opticsfile"] = "opticsfile.19"
change_dict["%knob_Re_value"] = "0."
change_dict["%knob_Im_value"] = "0."
change_dict["%knob_Re_type"] = "CMRS.b1_sq"
change_dict["%knob_Im_type"] = "CMIS.b1_sq"
change_dict["%error_component"] = "quadrupole"
change_dict["%error_strength"] = "0.00004*gauss()"
change_dict["%pattern_1"] = "R3"
change_dict["%pattern_2"] = "R3"
change_dict["%quad_component"] = "quadrupole"
change_dict["%quad_pattern_1"] = "R7"
change_dict["%quad_pattern_2"] = "R7"
change_dict["%quad_strength"] = "0.0001"
change_dict["%twiss_pattern"] = "BPM"
change_dict["%colknob1"] = "0."
change_dict["%colknob5"] = "0."


change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%colknob1","5.")
#plot_colknob_correction(change_dict_local,"colknobCorrection_colknob1_randomSQ_localQ_gauss3_R5.pdf")

change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%colknob5","5.")
#plot_colknob_correction(change_dict_local,"colknobCorrection_colknob5_randomSQ_localQ_gauss3_R5.pdf")

change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%quad_pattern_1","R7")
change_value(change_dict_local,"%quad_pattern_2","R7")
change_value(change_dict_local,"%pattern_1",".")
change_value(change_dict_local,"%pattern_2",".")
change_value(change_dict_local,"%quad_strength","0.0001")
change_value(change_dict_local,"%error_strength","0.000001*gauss()")
change_value(change_dict_local,"%colknob1","5")
#plot_colknob_correction(change_dict_local,"colknobCorrection_randomSQ_LocalQR7_colknob1.pdf")

change_dict_local = copy.deepcopy(change_dict)
change_value(change_dict_local,"%quad_pattern_1","R4")
change_value(change_dict_local,"%quad_pattern_2","R4")
change_value(change_dict_local,"%pattern_1",".")
change_value(change_dict_local,"%pattern_2",".")
change_value(change_dict_local,"%quad_strength","0.0002")
change_value(change_dict_local,"%error_strength","0.0000015*gauss()")
change_value(change_dict_local,"%colknob1","5")
plot_colknob_correction(change_dict_local,"colknobCorrection_randomSQ_LocalQR4_colknob1_5.pdf")

"""
#bootleg double correction
knob_Re, knob_Im = 0.0057630052860769025, -1.9595073787698426e-05

set_knobs(change_dict_local,knob_Re,knob_Im)
plot_colknob_correction(change_dict_local,"colknobCorrection_randomSQ_LocalQR4_colknob1_bootlegDoubleCorrection.pdf")
"""




