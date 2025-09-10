module dipolemoment

using WignerSymbols
push!(LOAD_PATH,pwd())
using trans_real_complex
using spin

export dipole

#######################################################################################################
function dipole(isomer,Utrans,Urot,nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,Nrot2,svd_error)

mmax=jmax
kmax=jmax

Tx,Ty,Tz,RotX,RotY,RotZ = dipole_components(jmax,Ntrans,Nrot) 

RotX_tmp = zeros(Nrot,Nrot)
RotY_tmp = zeros(Nrot,Nrot)
RotZ_tmp = zeros(Nrot,Nrot)

Rx = zeros(Nrot2,Nrot2)
Ry = zeros(Nrot2,Nrot2)
Rz = zeros(Nrot2,Nrot2)

#Transform to spin sub space#
tmpX = spin_isomer(isomer,jmax,RotX)
tmpY = spin_isomer(isomer,jmax,RotY)
tmpZ = spin_isomer(isomer,jmax,RotZ)

#Transform complex Wigner basis to real rotational eigenbasis#
Rx = transform_realbasis(Urot,tmpX)
Ry = transform_realbasis(Urot,tmpY)
Rz = transform_realbasis(Urot,tmpZ)

Dx_matrix = zeros((Ntrans*Nrot2,Ntrans*Nrot2))
Dy_matrix = zeros((Ntrans*Nrot2,Ntrans*Nrot2))
Dz_matrix = zeros((Ntrans*Nrot2,Ntrans*Nrot2))

a1=0
for i1=1:Ntrans
for j1=1:Nrot2
	a1+=1
	a2=0
	for i2=1:Ntrans
	for j2=1:Nrot2
		a2+=1
		Dx_matrix[a1,a2] = Tx[i1,i2]*Rx[j1,j2]
		Dy_matrix[a1,a2] = Ty[i1,i2]*Ry[j1,j2]
		Dz_matrix[a1,a2] = Tz[i1,i2]*Rz[j1,j2]
	end
	end
end
end


return Dx_matrix,Dy_matrix,Dz_matrix 
end
#######################################################################################################
function dipole_components(jmax,Ntrans,Nrot) 

Rot1 = zeros(ComplexF64,(Nrot,Nrot))
Rot2 = zeros(ComplexF64,(Nrot,Nrot))
Rot3 = zeros(ComplexF64,(Nrot,Nrot))


i1=0
for j1=0:jmax
for m1=-j1:j1
for k1=-j1:j1
	i1+=1
	i2=0
	for j2=0:jmax
	for m2=-j2:j2
	for k2=-j2:j2
		i2+=1
		if k1 == k2
			#q=-1 component#
			Rot1[i1,i2] =((-1.0)^(m1+k1))*sqrt((2*j1+1)*(2*j2+1))*wigner3j(BigFloat,j1,1,j2,-k1,0,k1)*wigner3j(BigFloat,j1,1,j2,-m1,-1,m2)		
			#q=0 component#
			Rot2[i1,i2] =((-1.0)^(m1+k1))*sqrt((2*j1+1)*(2*j2+1))*wigner3j(BigFloat,j1,1,j2,-k1,0,k1)*wigner3j(BigFloat,j1,1,j2,-m1,0,m2)		
			#q=1 component#
			Rot3[i1,i2] =((-1.0)^(m1+k1))*sqrt((2*j1+1)*(2*j2+1))*wigner3j(BigFloat,j1,1,j2,-k1,0,k1)*wigner3j(BigFloat,j1,1,j2,-m1,1,m2)		
		end
	end
	end
	end
end
end
end

#Recombine spherical components to Cartesian#
TransX = zeros(ComplexF64,(Ntrans,Ntrans))
TransY = zeros(ComplexF64,(Ntrans,Ntrans))
TransZ = zeros(ComplexF64,(Ntrans,Ntrans))

for i=1:Ntrans
	TransX[i,i] = 1.0
	TransY[i,i] = 1.0
	TransZ[i,i] = 1.0
end

RotX = zeros(ComplexF64,(Nrot,Nrot))
RotY = zeros(ComplexF64,(Nrot,Nrot))
RotZ = zeros(ComplexF64,(Nrot,Nrot))

RotX = -(1.0/sqrt(2.0))*(Rot3-Rot1)
RotY = (1.0im/sqrt(2.0))*(Rot3+Rot1)
RotZ = Rot2

return TransX,TransY,TransZ,RotX,RotY,RotZ
end
#######################################################################################################
end
