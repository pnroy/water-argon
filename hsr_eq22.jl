using WignerSymbols
eHtokHz=1.0/1.51982850071586e-13

function SR_value(jp,kp,mp,jo,ko,mo,MI)
	# justin
	# C1, l = 1, k = abs(1) = 1/2(C1,ab - C1,ba) = 1/2(49,79- 21,33) = 14,23 et C2, l = 1, k = abs(1) = -14,23
	# C1, l = 2, k = 1 = -1/2(49,79+21,33) = -35,56 ; C1, l = 2, k = -1 = 35,56 et *(-1) pour les C2,...
	# les valeurs sont les mêmes pour +/- 1
	Cl1_i1_plus=14.23/eHtokHz
	Cl1_i2_plus=-14.23/eHtokHz
	Cl1_i1_minus=14.23/eHtokHz
	Cl1_i2_minus=-14.23/eHtokHz

	Cl2_i1_plus=-35.56/eHtokHz
	Cl2_i2_plus=35.56/eHtokHz
	Cl2_i1_minus=35.56/eHtokHz
	Cl2_i2_minus=-35.56/eHtokHz
	value=0.0
	factor=1.
	if abs(jo-jp) <= 1 && abs(ko-kp)==1
		#l=1
		c1=0
		if (kp-ko)>0
			c1=Cl1_i1_plus-Cl1_i2_plus
		else
			c1=Cl1_i1_minus-Cl1_i2_minus
		end
		value+=sqrt(3.)*c1*wigner3j(jo,1,jp,ko,kp-ko,-kp)*
		(sqrt((jp*(jp+1)*(2*jp+1)))*(-1.0)*wigner6j(1,1,1,jp,jp,jo)
		+sqrt((jo*(jo+1)*(2*jo+1)))*wigner6j(1,1,1,jo,jo,jp))
		#l=2
		c2=0
		if (kp-ko)>0 
			c2=Cl2_i1_plus-Cl2_i2_plus
		else
			c2=Cl2_i1_minus-Cl2_i2_minus
		end  
		value+=sqrt(5.)*c2*wigner3j(jo,2,jp,ko,kp-ko,-kp)*
		(sqrt((jp*(jp+1)*(2*jp+1)))*(-1.0)*wigner6j(1,2,1,jp,jp,jo)
		+sqrt((jo*(jo+1)*(2*jo+1)))*wigner6j(1,2,1,jo,jo,jp))

		factor=.25*(-1.)^(mp-ko)*sqrt((2*jo+1)*(2*jp+1))*wigner3j(jp,1,jo,-mp,MI,mo)
		# typo Eq. 35
		#factor=.25*(-1.)^(-mo-ko)*sqrt((2*jo+1)*(2*jp+1))*wigner3j(jp,1,jo,-mp,MI,mo) 
	end
	return factor*value
end

function SR_valueEq22(jp,kp,mp,jo,ko,mo,MI)
	# justin
	# C1, l = 1, k = abs(1) = 1/2(C1,ab - C1,ba) = 1/2(49,79- 21,33) = 14,23 et C2, l = 1, k = abs(1) = -14,23
	# C1, l = 2, k = 1 = -1/2(49,79+21,33) = -35,56 ; C1, l = 2, k = -1 = 35,56 et *(-1) pour les C2,...
	# les valeurs sont les mêmes pour +/- 1
	Cl1_i1_plus=14.23/eHtokHz
	Cl1_i2_plus=-14.23/eHtokHz
	Cl1_i1_minus=14.23/eHtokHz
	Cl1_i2_minus=-14.23/eHtokHz

	Cl2_i1_plus=-35.56/eHtokHz
	Cl2_i2_plus=35.56/eHtokHz
	Cl2_i1_minus=35.56/eHtokHz
	Cl2_i2_minus=-35.56/eHtokHz
	value=0.0
	factor=1.
	if abs(jo-jp) <= 1 && abs(ko-kp)==1

		for i=1:2
			#l=1
			c1=0
			if (kp-ko)>0
				c1=Cl1_i1_plus-Cl1_i2_plus
			else
				c1=Cl1_i1_minus-Cl1_i2_minus
			end
			value+=sqrt(3.)*c1*wigner3j(jo,1,jp,ko,kp-ko,-kp)*
			(sqrt((jp*(jp+1)*(2*jp+1)))*(-1.0)*wigner6j(1,1,1,jp,jp,jo)
			+sqrt((jo*(jo+1)*(2*jo+1)))*wigner6j(1,1,1,jo,jo,jp))
			#l=2
			c2=0
			if (kp-ko)>0 
				c2=Cl2_i1_plus-Cl2_i2_plus
			else
				c2=Cl2_i1_minus-Cl2_i2_minus
			end  
			value+=sqrt(5.)*c2*wigner3j(jo,2,jp,ko,kp-ko,-kp)*
			(sqrt((jp*(jp+1)*(2*jp+1)))*(-1.0)*wigner6j(1,2,1,jp,jp,jo)
			+sqrt((jo*(jo+1)*(2*jo+1)))*wigner6j(1,2,1,jo,jo,jp))
		end
		factor=.25*(-1.)^(1-MI+mo-ko)*sqrt((2*jo+1)*(2*jp+1))*wigner3j(1,1,0,-MI.MI,0)*wigner3j(jp,1,jo,-mp,MI,mo)
	end
	return factor*value
end


jmax=1
for MI=-1:1
for jo=0:jmax
	for mo=-jo:jo
		for ko=-jo:jo
			if mod(ko,2) != 0
				for jp=0:jmax
					if abs(jo-jp) <= 1
						for mp=-jp:jp
							for kp=-jp:jp
								if mod(kp,2) == 0
									value=SR_value(jp,kp,mp,jo,ko,mo,MI)*eHtokHz
									if abs(value) >0
										println("MI=",MI," jo=",jo," ko=",ko," mo=",mo," jp=",jp," kp=",kp," mp=",mp," ",value)
									end
								end
							end
						end		
					end
				end
			end
		end
	end
end
end