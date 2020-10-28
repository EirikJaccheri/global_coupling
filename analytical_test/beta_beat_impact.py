import sys
sys.path.append("../global_coupling/")
from global_coupling import *
      
modulename = "global_coupling"
if modulename in sys.modules:
    print("You have imported the {} module".format(modulename))

N_error = 12

madx_path = "/home/eirik/madx "

folder_path = "/home/eirik/CERN/global_coupling_correction/analytical_test/"
response_path = folder_path + "response_beta_beat_impact.madx"
C_min_path = folder_path + "exact_C_min_analytic_test.madx"
#FineTuneCoupling_path = folder_path + "FineTuneCoupling_squeeze.madx"

path_dict = {}
path_dict["madx_path"] = madx_path
path_dict["response_path"] = response_path
path_dict["C_min_path"] = C_min_path
#path_dict["FineTuneCoupling_path"] = FineTuneCoupling_path
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

def B_matrix(mux_BPM,muy_BPM,Qx,Qy,betax_error,betay_error,mux_error,muy_error,KS):
	A = np.exp(1j*2*np.pi*(muy_BPM - mux_BPM)) / 4 / (1 - np.exp(2*np.pi*1j*(Qx - Qy)))
	B = np.sqrt(betax_error*betay_error)*np.exp(1j*2*np.pi*(mux_error - muy_error))
	return np.tensordot(A,B,axes = 0)


def f_1001(mux_BPM,muy_BPM,Qx,Qy,betax_error,betay_error,mux_error,muy_error,KS):
	B = B_matrix(mux_BPM,muy_BPM,Qx,Qy,betax_error,betay_error,mux_error,muy_error,KS)
	return np.dot(B,KS)

	
def get_knob_matrix(change_dict):
	change_dict_local = copy.deepcopy(change_dict)
	change_value(change_dict_local,"%twiss_pattern",".")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	name_l = tw40cm.NAME
	

	N_KQS = 12
	knob_Re_name = change_dict_local["%knob_Re_type"]
	knob_Im_name = change_dict_local["%knob_Im_type"]
	
	with open("lhc/lhc_as-built.seq",'r') as read_obj:
		data_lhc = read_obj.readlines()
	
	
	KQS_matrix = []
	KQS_index_l = []
	KQS_name_l = []
	with open("knob_matrix.txt",'r') as read_obj:
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
	return np.array(KQS_matrix), np.array(KQS_index_l),KQS_name_l

def plot_f1001_comparison(change_dict,KS_names,KS):
	change_dict_local = copy.deepcopy(change_dict)
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
	
	S = np.array(tw40cm.S)
	f1001_anal = f_1001(mux_BPM,muy_BPM,Qx,Qy,betx_error,bety_error,mux_error,muy_error,KS)
	
	change_value(change_dict_local,"%pattern_1",str(KS_names[0]))
	change_value(change_dict_local,"%pattern_2",str(KS_names[1]))
	change_value(change_dict_local,"%error_strength",str(KS[0]))
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
	plt.show()

def plot_f1001_knob_comparison(change_dict):
	change_dict_local = copy.deepcopy(change_dict)
	knob_Re = float(change_dict_local["%knob_Re_value"])
	knob_Im = float(change_dict_local["%knob_Im_value"])
	knob_arr = np.array([knob_Re,knob_Im])
	
	
	KQS_matrix, KQS_index_l, KQS_name_l = get_knob_matrix(change_dict)
	
	
	change_value(change_dict_local,"%twiss_pattern","BPM")
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm.Cmatrix()
	S = np.array(tw40cm.S)
	f1001_madx = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	
	
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
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(S,abs(f1001_0),label = "madx_knobs0")
	ax1.plot(S,abs(f1001_anal),label = "analytical")
	ax1.plot(S,abs(f1001_madx),label = "madx")
	ax1.set_ylim(bottom=0)
	ax1.set_ylim(top = 2*max(np.max(abs(f1001_madx)),np.max(abs(f1001_anal))))
	ax1.legend()
	plt.show()
	
def get_analytic_response_matrix(change_dict):
	change_dict_local = copy.deepcopy(change_dict)
	change_value(change_dict_local,"%knob_Re_value","0.") #should knobs be zero?
	change_value(change_dict_local,"%knob_Im_value","0.")
	change_value(change_dict_local,"%error_strength","0.")

	KQS_matrix, KQS_index_l, KQS_name_l = get_knob_matrix(change_dict)
	
	
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
	B = B_matrix(mux_BPM,muy_BPM,Qx,Qy,betx_error,bety_error,mux_error,muy_error,KS)
	Z = np.dot(B,KQS_matrix) * l
	R = np.block([[np.real(Z)],[np.imag(Z)]])
	R_inverse = np.linalg.pinv(R)
	return R_inverse
	
	

def plot_response_matrix_comparison(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	f_0 = get_f(response_path,"twiss.original",change_dict_local)
	C_min_0 = get_C_min(change_dict_local)
	
	R_inverse_analytic =  get_analytic_response_matrix(change_dict_local)
	knob_Re_analytic_res, knob_Im_analytic_res = get_response_knobs(R_inverse_analytic,change_dict_local)
	change_value(change_dict_local,"%knob_Re_value",str(knob_Re_analytic_res))
	change_value(change_dict_local,"%knob_Im_value",str(knob_Im_analytic_res))
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm.Cmatrix()
	S = np.array(tw40cm.S)
	f_analytic_res = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	C_min_analytic_res = get_C_min(change_dict_local)
	
	change_value(change_dict_local,"%knob_Re_value","0.")
	change_value(change_dict_local,"%knob_Im_value","0.")
	R_inverse = get_responsematrix(change_dict_local)
	tw40cm0 = get_twiss(response_path,"twiss.original",change_dict_local)
	tw40cm0.Cmatrix()
	knob_Re_res, knob_Im_res = get_response_knobs(R_inverse,change_dict_local)
	set_knobs(change_dict_local,knob_Re_res, knob_Im_res) 
	f_res = get_f(response_path,"twiss.original",change_dict_local)
	C_min_res = get_C_min(change_dict_local)
	
	fig = plt.figure()
	fig.suptitle(r"$C_0$ = " + "{:.3e}".format(C_min_0) + r"	$C_{analytic-res}$ = " + "{:.3e}".format(C_min_analytic_res) + r"	$C_{res}$ = " + "{:.3e}".format(C_min_res))
	
	print(np.linalg.norm(R_inverse_analytic - R_inverse) / np.linalg.norm(R_inverse))
	
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(S,abs(f_0),label = "f_0")
	ax1.plot(S,abs(f_analytic_res),label = "f_analytic_res")
	ax1.plot(S,abs(f_res),label = "f_res")
	ax1.legend()
	plt.savefig('plots/' + savepath)
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


change_dict["%error_strength"] = "0."
change_dict["%quad_strength"] = "0.00018"
KS_names = np.array(["MQS.23R1.B1","MQS.27L2.B1"])
KS = np.array([0.0001,0.0001])
#plot_f1001_comparison(change_dict,KS_names,KS)


change_dict["%knob_Re_type"] = "CMRS.b1_sq"
change_dict["%knob_Im_type"] = "CMIS.b1_sq"
change_dict["%error_strength"] = "0."
change_dict["%quad_strength"] = "0."
change_dict["%knob_Re_value"] = "0."
change_dict["%knob_Im_value"] = "0."
#plot_f1001_knob_comparison(change_dict)



change_dict["%knob_Re_type"] = "CMRS.b1_sq"
change_dict["%knob_Im_type"] = "CMIS.b1_sq"
change_dict["%error_strength"] = "0.00003*gauss()"
change_dict["%quad_strength"] = "0."
change_dict["%knob_Re_value"] = "0."
change_dict["%knob_Im_value"] = "0."
plot_response_matrix_comparison(change_dict,"response_matrix_comparison2.pdf")


