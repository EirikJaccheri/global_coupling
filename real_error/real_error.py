import sys
sys.path.append("../global_coupling/")
from global_coupling import *
      
modulename = "global_coupling"
if modulename in sys.modules:
    print("You have imported the {} module".format(modulename))

madx_path = "/home/ehoydals/madx "
response_path = "/home/ehoydals/global_coupling_correction/real_error/real_error_response.madx"
C_min_path = "/home/ehoydals/global_coupling_correction/real_error/real_error_C_min_matching.madx"
FineTuneCoupling_path = "/home/ehoydals/global_coupling_correction/real_error/real_error_FineTuneCoupling.madx"

path_dict = {}
path_dict["madx_path"] = madx_path
path_dict["response_path"] = response_path
path_dict["C_min_path"] = C_min_path
path_dict["FineTuneCoupling_path"] = FineTuneCoupling_path
set_global_paths(path_dict)

lhc_path =  "/afs/cern.ch/eng/lhc/optics/runII/2018"

reset_response_dict_new = {
"%error_strength" : "0.",
"%measured_C_min_comment" : "!"
}
set_reset_dict(reset_response_dict,reset_response_dict_new)

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

change_dict = {}
change_dict["%lhc_path"] = lhc_path
change_dict["%knob_Re_value"] = "0."
change_dict["%knob_Im_value"] = "0."
change_dict["%knob_Re_type"] = "CMRS.b1_sq"
change_dict["%knob_Im_type"] = "CMIS.b1_sq"
change_dict["%error_strength"] = "0."
change_dict["%measured_C_min_comment"] = ""
change_dict["%measured_C_min_Re_value"] = "-0.013198519341363707"
change_dict["%measured_C_min_Im_value"] = "0.007487633112887618"
change_dict["%twiss_pattern"] = "BPM"
plot_real_error_correction(change_dict,5,"real_error_correction_madx.pdf")
