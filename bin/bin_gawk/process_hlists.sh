# preprocessing for emerge using awk
# id pid snapnum mvir rvir rs mpeak mpeak_scale Acc_Rate_1Tdyn Time_to_future_merger Future_merger_MMP_ID t_dynamical rho_nfw_factor_at_rvir

gravG = 4.499753324353495e-24;
# kpc**3/(solMass * yr**2)

rho_crit * delta_vir *

rho_at_rvir = rho_crit * delta_vir * ($12/$13)*($12/$13) / ((1+($12/$13))*((1+($12/$13))*n.log(1+($12/$13))-($12/$13)))

prefactor = str(rho_crit)+" * "+str(delta_vir)+" * "

gawk_command = """gawk 'NR>63 {print $2, $6, $32, $11, $12, $13, $61, $70, $67, $78, $79, sqrt($12*$12*$12/(4.499753324353495e-24*$11)), """+prefactor+"""4 * ($12/$13)*($12/$13) / ((1+($12/$13))*((1+($12/$13))*log(1+($12/$13))-($12/$13)))}' /home/comparat/data/MultiDark/head10k_hlist_0.80130.list > /home/comparat/data/MultiDark/tmp1"""
print gawk_command
os.system(gawk_command)

gawk 'NR>63 {if ( $11 >= 1.51e+11 ) print $2, $6, $32, $11, $12, $13, $61, $70, $67, $78, $79, sqrt($12*$12*$12/(4.49975332435e-24*$11)), 195849652.248 * 4 *  $12 * $12 * ($12/$13)*($12/$13) , ((1+($12/$13))*((1+($12/$13))*log(1+($12/$13))-($12/$13)))}' /home/comparat/data/MultiDark/head10k_hlist_0.80130.list > /home/comparat/data/MultiDark/tmp1

gawk 'NR>63 {if ( $11 >= 1.51e+11 ) print $2, $6, $32, $11, $12, $13, $12/$13, $61, $70, $67, $78, $79, sqrt($12*$12*$12/(4.499753324353495e-24*$11)), 416.689951731 * 201.207074128 * 4.0 * ($12/$13)*($12/$13) , ((1.0+($12/$13))*((1.0+($12/$13))*log(1.0+($12/$13))-($12/$13)))}' /home/comparat/data/MultiDark/head10k_hlist_0.80130.list > /home/comparat/data/MultiDark/tmp1

if ($3 =="" || $4 == "" || $5 == "")
	print "Some score for the student",$1,"is missing";'
}' student-marks

def tau_quenching(tdyn, tau_0, tau_s, m_star):
	if m_star < 1e10 :
		return tdyn * tau_0
	else :
		return tdyn * tau_0 * (m_star * 10.**(-10.))**(tau_s)

f_loss = lambda t : 0.05*n.log( 1 + t / (1.4*10**6))

