import sys
sys.path.append("../global_coupling/")
from global_coupling import *
     
modulename = 'global_coupling'
if modulename in sys.modules:
    print 'You have imported the {} module'.format(modulename)

path_dict = {}
madx_path = "/home/ehoydals/madx "
path_dict["madx_path"] = madx_path
set_global_paths(path_dict)

lhc_path =  "/afs/cern.ch/eng/lhc/optics/runII/2018"




def plot_C_min(chagne_dict,n_colknob,colknob_step,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	change_value(change_dict_local,"%MUX","62.315")
	change_value(change_dict_local,"%MUY","60.315")
	change_value(change_dict_local,"%colknob1","0.")
	change_value(change_dict_local,"%colknob5","0.")
	
	colknob_l = colknob_step * np.arange(n_colknob)
	C_min_1_l = np.zeros(n_colknob)
	C_min_5_l = np.zeros(n_colknob)
	for i in range(n_colknob):
		change_value(change_dict_local,"%colknob1",str(colknob_l[i]))
		tw40cm = get_twiss("colknob_test.madx","twiss.test",change_dict_local)
		Q1, Q2 = tw40cm.Q1, tw40cm.Q2
		C_min_1_l[i] = abs(Q1 - Q2 - 2)
		
	change_value(change_dict_local,"%colknob1","0.")
	change_value(change_dict_local,"%colknob5","0.")
	for i in range(n_colknob):
		change_value(change_dict_local,"%colknob5",str(colknob_l[i]))
		tw40cm = get_twiss("colknob_test.madx","twiss.test",change_dict_local)
		Q1, Q2 = tw40cm.Q1, tw40cm.Q2
		Q1, Q2 = tw40cm.Q1, tw40cm.Q2
		C_min_5_l[i] = abs(Q1 - Q2 - 2)
	
	return colknob_l,C_min_1_l,C_min_5_l



def local_coupling(K,betax,betay,mux,muy,l,Q1,Q2,theta):
	p = 2
	print(mux,muy,mux-muy)
	return 1 / (2*np.pi) * K * l * np.sqrt(betax*betay)*np.exp(1j*(2*np.pi*(mux - muy) - (Q1 - Q2 -p)*theta))


def plot_C_min_predicted(change_dict,n_colknob,colknob_step,savepath):
	change_dict_local = copy.deepcopy(change_dict)
	
	change_value(change_dict_local,"%MUX","62.310")
	change_value(change_dict_local,"%MUY","60.320")
	change_value(change_dict_local,"%colknob1","0.")
	change_value(change_dict_local,"%colknob5","0.")
	colknob_l = colknob_step * np.arange(n_colknob)
	C_min_1_l = np.zeros(n_colknob)
	C_min_5_l = np.zeros(n_colknob)
	l = 0.223	
	
	tw40cm = get_twiss("colknob_test.madx","twiss.test",change_dict_local)
	name_l = tw40cm.NAME
	Q1, Q2 = tw40cm.Q1, tw40cm.Q2 
	R = tw40cm.LENGTH / 2 / np.pi
	S = tw40cm.S
	betx_l = np.array(tw40cm.BETX)
	bety_l = np.array(tw40cm.BETY)
	mux_l = np.array(tw40cm.MUX)
	muy_l = np.array(tw40cm.MUY)

	index_L1 = name_l.index("MQSX.3L1")
	index_R1 = name_l.index("MQSX.3R1")
	index_L5 = name_l.index("MQSX.3L5")
	index_R5 = name_l.index("MQSX.3R5")

	theta_L1 = S[index_L1] / R
	theta_R1 = S[index_R1] / R
	theta_L5 = S[index_L5] / R
	theta_R5 = S[index_R5] / R


	#calculating C_- at IP1
	for i in range(n_colknob):
		K1 = -colknob_l[i] * 1e-4
		K2 = float(change_dict_local["%percentage"]) * colknob_l[i] * 1e-4	
		C_min_1 = local_coupling(K1,betx_l[index_L1],bety_l[index_L1],mux_l[index_L1],muy_l[index_L1],l,Q1,Q2,theta_L1) +	local_coupling(K2,betx_l[index_R1],bety_l[index_R1],mux_l[index_R1],muy_l[index_R1],l,Q1,Q2,theta_R1)		  
		C_min_1_l[i] = abs(C_min_1)

		C_min_5 = local_coupling(K1,betx_l[index_L5],bety_l[index_L5],mux_l[index_L5],muy_l[index_L5],l,Q1,Q2,theta_L5) + local_coupling(K2,betx_l[index_R5],bety_l[index_R5],mux_l[index_R5],muy_l[index_R5],l,Q1,Q2,theta_R5)
		C_min_5_l[i] = abs(C_min_5)
	
	return colknob_l,C_min_1_l,C_min_5_l
		

def plot_colknob_percentage_comparison(change_dict,savepath,colknob_value,start,steps):
	change_dict_local = copy.deepcopy(change_dict)
	
	change_value(change_dict_local,"%MUX","62.310")
	change_value(change_dict_local,"%MUY","60.320")
	change_value(change_dict_local,"%colknob1","0.")
	change_value(change_dict_local,"%colknob5","0.")	
	
	change_value(change_dict_local,"%percentage","1.")
	l = 0.223	
	
	tw40cm = get_twiss("colknob_test.madx","twiss.test",change_dict_local)
	name_l = tw40cm.NAME
	Q1, Q2 = tw40cm.Q1, tw40cm.Q2 
	R = tw40cm.LENGTH / 2 / np.pi
	S = tw40cm.S
	betx_l = np.array(tw40cm.BETX)
	bety_l = np.array(tw40cm.BETY)
	mux_l = np.array(tw40cm.MUX)
	muy_l = np.array(tw40cm.MUY)

	index_L1 = name_l.index("MQSX.3L1")
	index_R1 = name_l.index("MQSX.3R1")
	index_L5 = name_l.index("MQSX.3L5")
	index_R5 = name_l.index("MQSX.3R5")

	theta_L1 = S[index_L1] / R
	theta_R1 = S[index_R1] / R
	theta_L5 = S[index_L5] / R
	theta_R5 = S[index_R5] / R

	Q10, Q20 = tw40cm.Q1, tw40cm.Q2	

	change_value(change_dict_local,"%MUX","62.315")
	change_value(change_dict_local,"%MUY","60.315")
	change_value(change_dict_local,"%colknob1",colknob_value)
	change_value(change_dict_local,"%colknob5","0.")	
	
	percentage_l = np.linspace(start,1,steps)
	
	C_min1_l = np.zeros(steps)
	C_min1_l_perturbation = np.zeros(steps)
	for i in range(steps):
		change_value(change_dict_local,"%percentage",str(percentage_l[i]))
		tw40cm = get_twiss("colknob_test.madx","twiss.test",change_dict_local)
		Q1, Q2 = tw40cm.Q1, tw40cm.Q2
		C_min1_l[i] = abs(Q1 - Q2 - 2)
		
		K1 = -float(colknob_value)*1e-4
		K2 = percentage_l[i] * float(colknob_value)*1e-4
		C_min1 = local_coupling(K1,betx_l[index_L1],bety_l[index_L1],mux_l[index_L1],muy_l[index_L1],l,Q10,Q20,theta_L1) +	local_coupling(K2,betx_l[index_R1],bety_l[index_R1],mux_l[index_R1],muy_l[index_R1],l,Q10,Q20,theta_R1)		  
		C_min1_l_perturbation[i] = abs(C_min1)
	
	change_value(change_dict_local,"%colknob1","0.")
	change_value(change_dict_local,"%colknob5",colknob_value)

	C_min5_l = np.zeros(steps)
	C_min5_l_perturbation = np.zeros(steps)
	for i in range(steps):
		change_value(change_dict_local,"%percentage",str(percentage_l[i]))
		tw40cm = get_twiss("colknob_test.madx","twiss.test",change_dict_local)
		Q1, Q2 = tw40cm.Q1, tw40cm.Q2
		C_min5_l[i] = abs(Q1 - Q2 - 2)	
		
		K1 = -float(colknob_value)*1e-4
		K2 = percentage_l[i] * float(colknob_value)*1e-4
		C_min5= local_coupling(K1,betx_l[index_L5],bety_l[index_L5],mux_l[index_L5],muy_l[index_L5],l,Q10,Q20,theta_L5) + local_coupling(K2,betx_l[index_R5],bety_l[index_R5],mux_l[index_R5],muy_l[index_R5],l,Q10,Q20,theta_R5)
		C_min5_l_perturbation[i] = abs(C_min5)
		
	fig = plt.figure()
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(percentage_l,C_min1_l,label = "colknob1")
	ax1.plot(percentage_l,C_min1_l_perturbation,label = "colknob1 perturbation")
	ax1.plot(percentage_l,C_min5_l,label = "colknob5")
	ax1.plot(percentage_l,C_min5_l_perturbation,label = "colknob5 perturbation")
	ax1.legend()
	ax1.set_xlabel("percentage correction")
	ax1.set_ylabel("$C_-$")
	plt.savefig("plots/" + savepath)
	plt.show()

def plot_colknob_perturbation_comparison(change_dict,savepath):
	colknob_l_perturbation , C_min_1_l_perturbation , C_min_5_l_perturbation = 			plot_C_min_predicted(change_dict,5,0.2,"C_min_predicted_temp1.pdf")
	colknob_l_madx , C_min_1_l_madx , C_min_5_l_madx = plot_C_min(change_dict,5,0.2,"C_min_measured_temp1.pdf")
	fig = plt.figure(figsize = plt.figaspect(0.6))

	ax1 = fig.add_subplot(1,2,1)
	ax1.plot(colknob_l_perturbation, C_min_1_l_perturbation, label = "colknob1")
	ax1.plot(colknob_l_perturbation, C_min_5_l_perturbation, label = "colknob5")
	ax1.set_xlabel("colknob")
	ax1.set_ylabel("$C_-$")
	ax1.set_title("perturbation")
	ax1.legend()

	ax2 = fig.add_subplot(1,2,2)
	ax2.plot(colknob_l_madx, C_min_1_l_madx, label = "colknob1")
	ax2.plot(colknob_l_madx, C_min_5_l_madx, label = "colknob5")
	ax2.set_xlabel("colknob")
	ax2.set_ylabel("$C_-$")
	ax2.set_title("madx matching")
	ax2.legend()

	plt.tight_layout()
	plt.savefig("plots/" + savepath)
	plt.show()
	

change_dict = {}
change_dict["%colknob1"] = "0."
change_dict["%colknob5"] = "0."
change_dict["%MUX"] = "62.310"
change_dict["%MUY"] = "60.320"
change_dict["%percentage"] = "1"
change_dict["%lhc_path"] = lhc_path
#plot_C_min(change_dict,2,0.4,"plot_C_min_TEST.pdf")
#plot_C_min_predicted(change_dict,2,0.4,"plot_C_min_predicted_TEST.pdf")
#plot_colknob_perturbation_comparison(change_dict,"pertubation_comparison99.pdf")
plot_colknob_percentage_comparison(change_dict,"knobcomparison_percentage.pdf","5",0.99,5)
"""
colknob1 = "0"
colknob5 = "0"
matching = False
percentage = "1"
error = [colknob1,colknob5,matching,percentage]
#plot_colknob_percentage_comparison(error,"knobcomparison_percentage.pdf","5",0.99,5)
	
colknob1 = "0"
colknob5 = "0"
matching = False
percentage = "0.99"
error = [colknob1,colknob5,matching,percentage]
plot_colknob_perturbation_comparison(error,"pertubation_comparison99.pdf")
"""
