SELECT           
     	FOF.Group_M_Crit2500 as m2500c,           
     	SH.BlackHoleMass as bhm           
     FROM           
     	RefL0100N1504_Subhalo as SH,           
     	RefL0100N1504_FOF as FOF     
     WHERE           
     	SH.SnapNum = 27           
     	and SH.GroupID = FOF.GroupID  
		and SH.BlackHoleMass > 0
		and FOF.RandomNumber <0.1
	ORDER BY 
		SH.BlackHoleMass asc