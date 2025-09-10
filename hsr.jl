using WignerSymbols
using LinearAlgebra
using LinearMaps
using Arpack
push!(LOAD_PATH,pwd())
using Printf


eHtocm1=219474.631363	#Eh to cm-1
eHtokHz=1.0/1.51982850071586e-13
Ae=     27.38549        #cm-1
Be=     14.58056
Ce=     9.51474
Ae=Ae/eHtocm1
Be=Be/eHtocm1
Ce=Ce/eHtocm1

function spin_isomer(isomer,jmax,matrix)
#Calculate dimension of para/ortho Hilbert space#
	para = []
	ortho = []
	i=0
	for j=0:jmax
		for m=-j:j
			for k=-j:j
				i+=1
				if mod(k,2) == 0
					append!(para,i)
				elseif mod(k,2) != 0
					append!(ortho,i)
				end
			end
		end
	end

	if isomer == "para"
		indices=sort(para)
	elseif isomer == "ortho"
		indices=sort(ortho)
	end

	Nspec=length(indices)

	#################################################
	#Construct new kinetic and potential matrices in#
	#Hilbert space associated with the respective   #
	#spin isomer                                    #
	#################################################

	matrix_out=zeros(Float64,(Nspec,Nspec))

	i=0
	for i1 in indices
		i+=1
		j=0
		for j1 in indices
			j+=1
			matrix_out[i,j]=matrix[i1,j1]
		end
	end

	return matrix_out
end

#####################################################################################################################
function kinetic_rotation(jmax,mmax,kmax,Ae,Be,Ce)

N=0
for j=0:jmax
	mj=min(mmax,j)
	kj=min(kmax,j)
	N=N+(2*mj+1)*(2*kj+1)
end
#Rotational constants#
#Ir-convention#
#D1=0.5*(Be+Ce)
#D2=Ae-0.5*(Be+Ce)
#D3=0.25*(Be-Ce)

#IIr-convention#
#D1=0.5*(Ce+Ae)
#D2=Be-0.5*(Ce+Ae)
#D3=0.25*(Ce-Ae)

#IIl-convention#
D1=0.5*(Ce+Ae)
D2=Be-0.5*(Ce+Ae)
D3=0.25*(Ae-Ce)

#IIIr-convention#
#D1=0.5*(Ae+Be)
#D2=Ce-0.5*(Ae+Be)
#D3=0.25*(Ae-Be)

Trot=zeros(Float64,(N,N))

i1=0
i2=0
i2_old=0
for j=0:jmax
	mj=min(mmax,j)
	kj=min(kmax,j)
	for m=-mj:mj
		i2_old=i2
		for k1=-kj:kj
			i1+=1
			i2=i2_old
        	for k2=-kj:kj
                i2+=1
				if k1 == k2-2
                    Trot[i1,i2]=D3*sqrt(1.0*j*(j+1)-1.0*k2*(k2-1))*sqrt(1.0*j*(j+1)-1.0*(k2-1)*(k2-2))
                end
                if k1 == k2+2
                    Trot[i1,i2]=D3*sqrt(1.0*j*(j+1)-1.0*k2*(k2+1))*sqrt(1.0*j*(j+1)-1.0*(k2+1)*(k2+2))
                end
        	end
            Trot[i1,i1]=j*(j+1)*D1+D2*k1^2
        end
	end
end

return Trot
end

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
function HSR_elements(MI,jmax,Northo,Npara)
	HSRmat=zeros(Float64,(Northo,Npara))
	value=0.0
	index_o=0
	for jo=0:jmax
		for mo=-jo:jo
			for ko=-jo:jo
				if mod(ko,2) != 0
					index_o+=1
					index_p=0
					for jp=0:jmax
						for mp=-jp:jp
							for kp=-jp:jp
								if mod(kp,2) == 0
									index_p+=1
									HSRmat[index_o,index_p]=SR_value(jp,kp,mp,jo,ko,mo,MI)*eHtokHz
									#println(MI," ",jo," ",ko," ",mo," ",jp," ",kp," ",mp," ",HSRmat[index_o,index_p]*eHtokHz)
								end
							end
						end		
					end
				end
			end
		end
	end
	return HSRmat
end


#################################################################
#Calculate kinetic matrices for translational and rotational DOF#
#################################################################
#Ttrans = kinetic_trans_analytic(nmax,omega)
let
	mmax=1
	jmax=1
	kmax=1

	Trot = kinetic_rotation(jmax,mmax,kmax,Ae,Be,Ce)



	#Project kinetic matrix to spin sub space#
	Tpara = spin_isomer("para",jmax,Trot)
	Tortho = spin_isomer("ortho",jmax,Trot)
	Npara = size(Tpara,1)
	Northo = size(Tortho,1)

	e_o,w_o=eigen(Tortho)
	e_p,w_p=eigen(Tpara)

	println()

	for s=1:Npara
		print("Para state ",s,": ")
		@printf("%0.2f",e_p[s]*eHtocm1)
		print(" 1/cm")
		println()
		ii=0
		for jp=0:jmax
			for mp=-jp:jp
				for kp=-jp:jp
					if mod(kp,2) == 0
						ii+=1
						print(w_p[ii,s],"|jp=",jp,",kp=",kp,",mp=",mp,"> ")
						if ii<Npara
							print("+ ")
						end
					end
				end
			end
		end
		println()
	end

	println()

	for s=1:Northo
		print("Ortho state ",s,": ")
		@printf("%0.2f",e_o[s]*eHtocm1)
		print(" 1/cm")
		println()
		ii=0
		for jo=0:jmax
			for mo=-jo:jo
				for ko=-jo:jo
					if mod(ko,2) != 0
						ii+=1
						@printf("%0.2f",w_o[ii,s])
						print("|jo=",jo,",ko=",ko,",mo=",mo,"> ")
						if ii<Northo
							print("+ ")
						end
					end
				end
			end
		end
		println()
	end

	println()

	HSR0_wig=HSR_elements(0,jmax,Northo,Npara)
	HSR1_wig=HSR_elements(1,jmax,Northo,Npara)
	HSRminus1_wig=HSR_elements(-1,jmax,Northo,Npara)

	tmp=BLAS.gemm('C','N', w_o,HSR0_wig)
	HSR0_asym=BLAS.gemm('N','N', tmp,w_p)
	tmp=BLAS.gemm('C','N', w_o,HSR1_wig)
	HSR1_asym=BLAS.gemm('N','N', tmp,w_p)
	tmp=BLAS.gemm('C','N', w_o,HSRminus1_wig)
	HSRminus1_asym=BLAS.gemm('N','N', tmp,w_p)

	println()

	for so=1:Northo
		for sp=1:Npara
			println("ortho state ",so," para state ",sp," HSR0 MI=0 ",HSR0_asym[so,sp])
			println("ortho state ",so," para state ",sp," HSR MI=1 ",HSR1_asym[so,sp])
			println("ortho state ",so," para state ",sp," HSR MI=-1 ",HSRminus1_asym[so,sp])
		end
	end


	# for MI=-1:1
	# for jo=0:jmax
	# 	for mo=-jo:jo
	# 		for ko=-jo:jo
	# 			if mod(ko,2) != 0
	# 				for jp=0:jmax
	# 					if abs(jo-jp) <= 1
	# 						for mp=-jp:jp
	# 							for kp=-jp:jp
	# 								if mod(kp,2) == 0
	# 									value=SR_value(jp,kp,mp,jo,ko,mo,MI)*eHtokHz
	# 									if abs(value) >0
	# 										println("MI=",MI," jo=",jo," ko=",ko," mo=",mo," jp=",jp," kp=",kp," mp=",mp," ",value)
	# 									end
	# 								end
	# 							end
	# 						end		
	# 					end
	# 				end
	# 			end
	# 		end
	# 	end
	# end
	# end
end