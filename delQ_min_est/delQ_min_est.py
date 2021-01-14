import sys
sys.path.append("../global_coupling/")
from global_coupling import *
      
modulename = "global_coupling"
if modulename in sys.modules:
    print("You have imported the {} module".format(modulename))

N_error = 12

madx_path = "/home/eirik/madx "

folder_path = "/home/eirik/CERN/global_coupling_correction/delQ_min_est/"
response_path = folder_path + "response_delQ_min_est.madx"
C_min_path = folder_path + "exact_C_min_delQ_min_est.madx"
FineTuneCoupling_path = folder_path + "FineTuneCoupling_delQ_min_est.madx"

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

import tfs
import pandas as pd

def delQ_min_new(Q1,Q2,f,S,R,MUX,MUY):
	N = len(f)
	#should this be abs
	delta = Q1 - Q2 - 2

	deltaS_list = [np.abs(S[i +1] - S[i]) for i in range(N-1)]
	deltaS_list.append(S[-1] - S[-2])
	deltaS = np.array(deltaS_list)
	return abs(4 * delta / (2*np.pi*R) * np.sum(deltaS * f * np.exp(-1j * 2*np.pi*(MUY - MUX -  0 * S * delta / R))))


def delQ_min_old(Q1,Q2,f,S,MUX,MUY):
	N = len(f)
	delta = np.abs(Q1%1 - Q2%1)
	return 4 * delta/N * np.sum(np.abs(f))
	
def delQ_min_Cmatrix(Q1,Q2,C):
	det_C = np.linalg.det(C)
	gamma = np.sqrt(1 - det_C)
	delQ_min_l =  2 * gamma * (np.cos(2*np.pi*Q1) - np.cos(2*np.pi*Q2)) / (np.pi*(np.sin(2*np.pi*Q1) + np.sin(2*np.pi*Q2))) * np.sqrt(det_C)
	return np.mean(delQ_min_l)	


def delQ_min_combinedRDT(Q1,Q2,f1001,f1010):
	delQ_min_l = (np.cos(2*np.pi*Q1) - np.cos(2*np.pi*Q2)) / (np.pi*(np.sin(2*np.pi*Q1) + np.sin(2*np.pi*Q2))) * 4 *np.sqrt(abs(f1001)**2 - abs(f1010)**2) / (1 + 4*(abs(f1001)**2 - abs(f1010)**2))
	return np.mean(delQ_min_l)

def delQ_min_teapot(Q1,Q2,C_madx,B_madx_bar):
	return np.mean(np.sqrt((np.linalg.det(C_madx + B_madx_bar))) / (np.pi*(np.sin(2*np.pi*Q1) + np.sin(2*np.pi*Q2))))
	
def delQ_min_Guignard(tw40cm):
	os.system("python3 tfs_to_pandas.py")
	error_dict = pd.read_csv("output_files/K1SL.out")
	
	name_l = tw40cm.NAME
	
	error_names = error_dict["NAME"].to_list()
	error_index = np.array([name_l.index(error_name) for error_name in error_names])
	
	k = error_dict["K1SL"]
	
	Q1 , Q2 = tw40cm.Q1 , tw40cm.Q2
	mux , muy = np.array(tw40cm.MUX)[error_index] , np.array(tw40cm.MUY)[error_index]
	betx , bety = np.array(tw40cm.BETX)[error_index] , np.array(tw40cm.BETY)[error_index]
	S = np.array(tw40cm.S)[error_index]
	delta = Q1 - Q2 - 2
	R = tw40cm.LENGTH / (2*np.pi)
	
	return np.abs(np.sum(np.sqrt(betx*bety) * k * np.exp(-2*np.pi*1j*(mux - muy - 0 * delta * S / R )))/(2*np.pi)) #control the sign 

def plot_delQ_min_comparison(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	
	Q1 , Q2 = tw40cm.Q1 , tw40cm.Q2
	f = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	f1010 = np.array(tw40cm.F1010R) + 1j * np.array(tw40cm.F1010I)
	S = np.array(tw40cm.S)
	R = tw40cm.LENGTH / 2 / np.pi
	MUX , MUY = np.array(tw40cm.MUX) , np.array(tw40cm.MUY)
	C = np.array(tw40cm.C).reshape((-1,2,2))
	
	re13, re14, re23, re24 = np.array(tw40cm.RE13), np.array(tw40cm.RE14), np.array(tw40cm.RE23), np.array(tw40cm.RE24)
	re31, re32, re41, re42 = np.array(tw40cm.RE31), np.array(tw40cm.RE32), np.array(tw40cm.RE41), np.array(tw40cm.RE42)
	B_madx = np.reshape(numpy.array([re13, re14, re23, re24]).T, (-1,2, 2))
	C_madx = np.reshape(numpy.array([re31, re32, re41, re42]).T, (-1,2, 2))
	
	S_matrix = np.array([[0,-1],[1,0]])
	B_madx_bar = -np.matmul(S_matrix,np.matmul(np.transpose(B_madx,(0,2,1)),S_matrix))
	
	C_min_Guignard = delQ_min_Guignard(tw40cm)
	C_min_madx = get_C_min(change_dict_local)
	C_min_Cmatrix = delQ_min_Cmatrix(Q1,Q2,C)
	C_min_combinedRDT = delQ_min_combinedRDT(Q1,Q2,f,f1010)
	C_min_teapot  = delQ_min_teapot(Q1,Q2,C_madx,B_madx_bar)
	C_min_new = delQ_min_new(Q1,Q2,f,S,R,MUX,MUY)
	C_min_old = delQ_min_old(Q1,Q2,f,S,MUX,MUY)
	
	print(C_min_teapot - C_min_madx)
	
	fig , ax = plt.subplots(1)
	fig.suptitle(r"$\Delta Q_{min}$ = " + "{:.3e}".format(C_min_madx) + r"	$\Delta Q_{min}^{Cmatrix}$ = " + "{:.3e}".format(C_min_Cmatrix) + "\n" +  r"$\Delta Q_{min}^{teapot}$ = "  + "{:.3e}".format(C_min_teapot) + r"	$\Delta Q_{min}^{Guignard}$ = "  + "{:.3e}".format(C_min_Guignard) +"\n" + r"$|C^-_{old}|$ = " + "{:.3e}".format(C_min_old) + r"	$|C^-_{new}|$ = " + "{:.3e}".format(C_min_new))
	
	ax.plot(S,abs(f),label = r"$|f_{1001}|$")
	ax.set_xlabel("Position [m]")
	ax.set_ylabel(r"$|f_{1001}|$")
	ax.set_ylim(0,1.2*max(abs(f)))
	fig.tight_layout(rect=[0, 0.03, 1, 0.8])
	plt.savefig("plots/" + savepath)
	plt.show()
	
def plot_delQ_min_comparison_bug(change_dict,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	
	tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	
	f = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
	S = np.array(tw40cm.S)
	
	#calculating with square root
	C = np.array(tw40cm.C)
	gamma = tw40cm.gamma
	print(np.shape(C[0]))
	print(np.shape(gamma))
	f_sqrt = ((C[:,0] + C[:,3]) * 1j + (C[:,1] - C[:,2])) / 4 / gamma
	
	C_min_madx = get_C_min(change_dict_local)
	
	fig , ax = plt.subplots(1)
	fig.suptitle(r"$\Delta Q_{min}^{matching}$ = " + "{:.3e}".format(C_min_madx))
	ax.plot(S,abs(f),label = r"$\gamma = 1 - |\bar{C}|$")
	ax.plot(S,abs(f_sqrt),label = r"$\gamma = \sqrt{1 - |\bar{C}|}$")
	ax.set_xlabel("Position [m]")
	ax.set_ylabel(r"$|f_{1001}|$")
	ax.set_ylim(0,1.2*max(abs(f)))
	ax.legend()
	
	plt.savefig("plots/" + savepath)
	plt.show()

def plot_delQ_min_table(change_dict_dict,savepath):
	column_labels = ["Matching","Cmatrix","CombRDT","Guignard","TobRog","Franchi"]
	row_labels = []
	N_rows = len(change_dict_dict)
	N_columns = len(column_labels)
	Data = np.zeros((N_rows,N_columns)).astype(float)
	for i,key in enumerate(change_dict_dict):
		row_labels.append(key)
		
		change_dict_local = copy.deepcopy(change_dict_dict[key])
		tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	
		Q1 , Q2 = tw40cm.Q1 , tw40cm.Q2
		f = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
		f1010 = np.array(tw40cm.F1010R) + 1j * np.array(tw40cm.F1010I)
		S = np.array(tw40cm.S)
		R = tw40cm.LENGTH / 2 / np.pi
		MUX , MUY = np.array(tw40cm.MUX) , np.array(tw40cm.MUY)
		C = np.array(tw40cm.C).reshape((-1,2,2))
		
		Data[i][0] = get_C_min(change_dict_local)
		Data[i][1] = delQ_min_Cmatrix(Q1,Q2,C)
		Data[i][2] = delQ_min_combinedRDT(Q1,Q2,f,f1010)
		Data[i][3] = delQ_min_Guignard(tw40cm)
		Data[i][4] = delQ_min_new(Q1,Q2,f,S,R,MUX,MUY)
		Data[i][5] = delQ_min_old(Q1,Q2,f,S,MUX,MUY)
		
	df = pd.DataFrame(Data, columns = column_labels, index = row_labels)
	buf_str = "plots/" + savepath
	df.to_latex(buf=buf_str)
	
	fig, ax = plt.subplots()

	# hide axes
	fig.patch.set_visible(False)
	ax.axis('off')
	ax.axis('tight')
	
	rcolors = plt.cm.BuPu(np.full(len(row_labels), 0.1))
	ccolors = plt.cm.BuPu(np.full(len(column_labels), 0.1))
	
	ax.table(cellText=Data, rowLabels = row_labels, rowColours = rcolors, colLabels=column_labels, colColours = ccolors, loc='center')
	#plt.savefig("plots/" + savepath)
	#plt.show()
	
def plot_bar_diagram(change_dict_dict,savepath):
	column_labels = ["Cmatrix","CombRDT","Guignard","TobRog","Franchi"]
	row_labels = []
	N_rows = len(change_dict_dict)
	N_columns = len(column_labels)
	Data = np.zeros((N_rows,N_columns)).astype(float)
	for i,key in enumerate(change_dict_dict):
		row_labels.append(key)
		
		change_dict_local = copy.deepcopy(change_dict_dict[key])
		tw40cm = get_twiss(response_path,"twiss.original",change_dict_local)
	
		Q1 , Q2 = tw40cm.Q1 , tw40cm.Q2
		f = np.array(tw40cm.F1001R) + 1j * np.array(tw40cm.F1001I)
		f1010 = np.array(tw40cm.F1010R) + 1j * np.array(tw40cm.F1010I)
		S = np.array(tw40cm.S)
		R = tw40cm.LENGTH / 2 / np.pi
		MUX , MUY = np.array(tw40cm.MUX) , np.array(tw40cm.MUY)
		C = np.array(tw40cm.C).reshape((-1,2,2))
		
		delta_Q_min_matching = get_C_min(change_dict_local)
		Data[i][0] = delQ_min_Cmatrix(Q1,Q2,C) / delta_Q_min_matching
		Data[i][1] = delQ_min_combinedRDT(Q1,Q2,f,f1010) / delta_Q_min_matching
		Data[i][2] = delQ_min_Guignard(tw40cm) / delta_Q_min_matching
		Data[i][3] = delQ_min_new(Q1,Q2,f,S,R,MUX,MUY) / delta_Q_min_matching
		Data[i][4] = delQ_min_old(Q1,Q2,f,S,MUX,MUY) / delta_Q_min_matching
		
	x = np.arange(N_rows)  # the label locations
	width = 0.15  # the width of the bars
	
	labels  = []
	for i in range(1,N_rows + 1):
		labels.append("E"+str(i))
	
	fig , ax = plt.subplots()
	
	rects_Cmatrix = ax.bar(x - 2 * width,Data[:,0],width,label = "Cmatrix")
	rects_combinedRDT = ax.bar(x - width,Data[:,1],width,label = "CombRDT")
	rects_Guignard = ax.bar(x,Data[:,2],width,label = "Guignard")
	rects_TobRog = ax.bar(x + width, Data[:,3],width,label = "TobRog")
	rects_Franchi = ax.bar(x + 2 * width, Data[:,4],width,label = "Franchi")
	
	ax.set_ylabel('$\Delta Q_{min}$/$\Delta Q_{min}^{matching}$')
	ax.set_title('$Q_x$ = ' + str(round(Q1,3)) + '	$Q_y$ = ' + str(round(Q2,3)))
	
	ax.set_xticks(x)
	ax.set_xticklabels(labels)
	ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=False, ncol=5)
	plt.tight_layout()
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
change_dict["%quad_strength"] = "0."
change_dict["%twiss_pattern"] = "."
change_dict["%colknob1"] = "0."
change_dict["%colknob5"] = "0."


a = np.array([1,1,1])
b = np.array([2,2,2])
c = np.array([3,3,3])
d = np.array([4,4,4])
A = np.array([a,b,c,d])
B = np.reshape(A.T,(-1,2,2))
print(B)

change_dict_local1 = copy.deepcopy(change_dict)
change_value(change_dict_local1,"%error_strength","0.00003*gauss()")
#plot_delQ_min_comparison(change_dict_local1,"delQ_min_comparison_gaussError_uniform.pdf")
#plot_delQ_min_comparison_bug(change_dict_local1,"delQ_min_comparison_bug.pdf")

change_dict_local2 = copy.deepcopy(change_dict)
change_value(change_dict_local2,"%error_strength","0.00001*gauss()")
change_value(change_dict_local2,"%pattern_1","R5")
change_value(change_dict_local2,"%pattern_2","R7")
#plot_delQ_min_comparison(change_dict_local2,"delQ_min_comparison_gaussError_R5R7.pdf")

change_dict_local3 = copy.deepcopy(change_dict)
change_value(change_dict_local3,"%error_strength","0.00001")
change_value(change_dict_local3,"%pattern_1","R5")
change_value(change_dict_local3,"%pattern_2","R5")
#plot_delQ_min_comparison(change_dict_local3,"delQ_min_comparison_uniformError_R5.pdf")

change_dict_local4 = copy.deepcopy(change_dict)
change_value(change_dict_local4,"%error_strength","0.00001")
change_value(change_dict_local4,"%pattern_1","R4")
change_value(change_dict_local4,"%pattern_2","R8")
#plot_delQ_min_comparison(change_dict_local3,"delQ_min_comparison_uniformError_R4R8.pdf")

change_dict_dict = {"1":change_dict_local1,"2":change_dict_local2,"3":change_dict_local3,"4":change_dict_local4}
#change_dict_dict = {"1":change_dict_local1,"2":change_dict_local2,"3":change_dict_local3}
#plot_delQ_min_table(change_dict_dict,"table_woBUG_311319.txt")
plot_bar_diagram(change_dict_dict,"bar_diagram_woBUG_3033.pdf")

#select, flag=twiss, column=name, s, betx, bety, mux, muy, r11, r12, r21, r22, alfx, alfy, re31, r32, re41, r42; twiss, file="minusI.twiss", rmatrix; 

	
	
	
