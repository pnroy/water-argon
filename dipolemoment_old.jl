module dipolemoment

using LinearAlgebra
using GaussQuadrature: laguerre,legendre 
push!(LOAD_PATH,pwd())
using potential
using trans_real_complex
using spin

export dipole

#######################################################################################################
function dipole(isomer,Utrans,Urot,nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,Nrot2,svd_error)

mmax=jmax
kmax=jmax

Cx,Cy,Cz,TransX,TransY,TransZ,RotX,RotY,RotZ = integration(Property,nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,eq_struct,error)

Nprodx = length(Cx)
Nprody = length(Cy)
Nprodz = length(Cz)

Tx = zeros(Ntrans,Ntrans,Nprodx)
Ty = zeros(Ntrans,Ntrans,Nprody)
Tz = zeros(Ntrans,Ntrans,Nprodz)

RotX_tmp = zeros(Nrot,Nrot,Nprodx)
RotY_tmp = zeros(Nrot,Nrot,Nprody)
RotZ_tmp = zeros(Nrot,Nrot,Nprodz)

Rx = zeros(Nrot2,Nrot2,Nprodx)
Ry = zeros(Nrot2,Nrot2,Nprody)
Rz = zeros(Nrot2,Nrot2,Nprodz)

for ip=1:Nprodx
	#Transform 3D-HO to real basis#
	tmp1 = transform_basis(Utrans,TransX[:,:,ip])
	Tx[:,:,ip] .= real(tmp1)
	#Transform complex Wigner basis to real rotational eigenbasis#
	tmp2 = transform_basis(Urot,RotX[:,:,ip])
	RotX_tmp[:,:,ip] .= real(tmp2)
	#Transform to spin sub space#
	Rx[:,:,ip] .= spin_isomer_matrix(isomer,jmax,mmax,kmax,Urot,RotX_tmp[:,:,ip])
end

for ip=1:Nprody
	#Transform 3D-HO to real basis#
	tmp1 = transform_basis(Utrans,TransY[:,:,ip])
	Ty[:,:,ip] .= real(tmp1)
	#Transform complex Wigner basis to real rotational eigenbasis#
	tmp2 = transform_basis(Urot,RotY[:,:,ip])
	RotY_tmp[:,:,ip] .= real(tmp2)
	#Transform to spin sub space#
	Ry[:,:,ip] .= spin_isomer_matrix(isomer,jmax,mmax,kmax,Urot,RotY_tmp[:,:,ip])
end

for ip=1:Nprodz
	#Transform 3D-HO to real basis#
	tmp1 = transform_basis(Utrans,TransZ[:,:,ip])
	Tz[:,:,ip] .= real(tmp1)
	#Transform complex Wigner basis to real rotational eigenbasis#
	tmp2 = transform_basis(Urot,RotZ[:,:,ip])
	RotZ_tmp[:,:,ip] .= real(tmp2)
	#Transform to spin sub space#
	Rz[:,:,ip] .= spin_isomer_matrix(isomer,jmax,mmax,kmax,Urot,RotZ_tmp[:,:,ip])
end

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
		dum=0.0
		for ip=1:Nprodx
			dum = Cx[ip]*Tx[i1,i2,ip]*Rx[j1,j2,ip]
		end
		Dx_matrix[a1,a2] = dum
	
		dum=0.0
		for ip=1:Nprody
			dum = Cy[ip]*Ty[i1,i2,ip]*Ry[j1,j2,ip]
		end
		Dy_matrix[a1,a2] = dum

		dum=0.0
		for ip=1:Nprodz
			dum = Cz[ip]*Tz[i1,i2,ip]*Rz[j1,j2,ip]
		end
		Dz_matrix[a1,a2] = dum
	end
	end
end
end


return Dx_matrix,Dy_matrix,Dz_matrix 
end
#######################################################################################################
function dipole_integration(nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,svd_error)

mmax=jmax
kmax=jmax
nu=mass*omega*0.5

#Grids#

Ngrid1 = NR*Ntheta*Nphi
Ngrid2 = Nalpha*Nm*Nk
Ngrid = Ngrid1*Ngrid2

R,Rweight = laguerre(NR,0.5)
ctheta1, ct_weight1 = legendre(Ntheta)
ctheta2, ct_weight2 = legendre(Nalpha)
phi1=zeros(Nphi)
for i=1:Nphi
	phi1[i]=(i-1)*2.0*pi/Nphi
end
phi2=zeros(Nm)
for i=1:Nm
	phi2[i]=(i-1)*2.0*pi/Nm
end
chi=zeros(Nk)
for i=1:Nk
	chi[i]=(i-1)*2.0*pi/Nk
end

#Write potential on grid#
dmom1 = zeros(Ngrid1,Ngrid2)
dmom2 = zeros(Ngrid1,Ngrid2)
dmom3 = zeros(Ngrid1,Ngrid2)

i1=0
for ir=1:NR
for it1=1:Ntheta
for ip1=1:Nphi
	i1+=1
	i2=0
	for it2=1:Nalpha
	for ip2=1:Nm
	for ic=1:Nk 
		i2 += 1
		dmom1[i1,i2],dmom2[i1,i2],dmom3[i1,i2] = mu(sqrt(0.5*R[ir]/nu),acos(ctheta1[it1]),phi1[ip1],acos(ctheta2[it2]),phi2[ip2],chi[ic]) 
	end
	end
	end
end
end
end

######################################
#############SVD Approach#############
######################################

f=open("log","a")
println(f,"##Decompose dipole moment on grid by SVD (x component)##")
close(f)
Cx,dx1,dx2 = decomposition(dmom1,svd_error)
f=open("log","a")
println(f,"##Decompose dipole moment on grid by SVD (y component)##")
close(f)
Cy,dy1,dy2 = decomposition(dmom2,svd_error)
f=open("log","a")
println(f,"##Decompose dipole moment on grid by SVD (z component)##")
close(f)
Cz,dz1,dz2 = decomposition(dmom3,svd_error)
Nprodx = length(Cx)
Nprody = length(Cy)
Nprodz = length(Cz)
####################################################################################################
coeff_start1=zeros(ComplexF64,Ntrans)
coeff_start2=zeros(ComplexF64,Nrot)
coeff_grid1=zeros(ComplexF64,Ngrid1)
coeff_grid2=zeros(ComplexF64,Ngrid2)
coeff_out1=zeros(ComplexF64,Ngrid1)
coeff_out2=zeros(ComplexF64,Ngrid2)
coeff_end1=zeros(ComplexF64,(Ntrans))
coeff_end2=zeros(ComplexF64,(Nrot))
transX_matrix=zeros(ComplexF64,(Ntrans,Ntrans,Nprodx))
transY_matrix=zeros(ComplexF64,(Ntrans,Ntrans,Nprody))
transZ_matrix=zeros(ComplexF64,(Ntrans,Ntrans,Nprodz))
rotX_matrix=zeros(ComplexF64,(Nrot,Nrot,Nprodx))
rotY_matrix=zeros(ComplexF64,(Nrot,Nrot,Nprody))
rotZ_matrix=zeros(ComplexF64,(Nrot,Nrot,Nprodz))

#########################################################
##Calculation of dipole matrices in iso. 3D-HO basis##
#########################################################

for ispec1=1:Ntrans
	coeff_start1.=0.0
	coeff_start1[ispec1]=1.0	

	###################################################
	##Transformation from NLM basis to spherical grid##
	###################################################

	coeff_grid1 = transformation_forward_HO(coeff_start1,nmax,NR,Ntheta,Nphi,nu,mass)

	#######################################
	##Act with dipole operator on grid##
	#######################################
	for iprod=1:Nprodx
		for i=1:Ngrid1
			coeff_out1[i]=coeff_grid1[i]*dx1[i,iprod]
		end
		###################################################
		##Transformation from spherical grid to NLM basis##
		###################################################
		coeff_end1 = transformation_backward_HO(coeff_out1,nmax,NR,Ntheta,Nphi,nu,mass,Ntrans)
		for ispec2=1:Ntrans
			transX_matrix[ispec2,ispec1,iprod]=coeff_end1[ispec2]
		end
	end
	for iprod=1:Nprody
		for i=1:Ngrid1
			coeff_out1[i]=coeff_grid1[i]*dy1[i,iprod]
		end
		###################################################
		##Transformation from spherical grid to NLM basis##
		###################################################
		coeff_end1 = transformation_backward_HO(coeff_out1,nmax,NR,Ntheta,Nphi,nu,mass,Ntrans)
		for ispec2=1:Ntrans
			transY_matrix[ispec2,ispec1,iprod]=coeff_end1[ispec2]
		end
	end
	for iprod=1:Nprodz
		for i=1:Ngrid1
			coeff_out1[i]=coeff_grid1[i]*dz1[i,iprod]
		end
		###################################################
		##Transformation from spherical grid to NLM basis##
		###################################################
		coeff_end1 = transformation_backward_HO(coeff_out1,nmax,NR,Ntheta,Nphi,nu,mass,Ntrans)
		for ispec2=1:Ntrans
			transZ_matrix[ispec2,ispec1,iprod]=coeff_end1[ispec2]
		end
	end
end

#####################################################
##Calculation of dipole matrices in Wigner basis##
#####################################################

for ispec1=1:Nrot
	coeff_start2.=0.0
	coeff_start2[ispec1]=1.0	

	#################################################
	##Transformation from JKM basis to angular grid##
	#################################################

	coeff_grid2 = transformation_forward_Wigner(coeff_start2,jmax,mmax,kmax,Nalpha,Nm,Nk)

	#######################################
	##Act with dipole operator on grid##
	#######################################
	for iprod=1:Nprodx
		for i=1:Ngrid2
			coeff_out2[i]=coeff_grid2[i]*dx2[i,iprod]
		end
		#################################################
		##Transformation from angular grid to JMK basis##
		#################################################
		coeff_end2 .= transformation_backward_Wigner(coeff_out2,jmax,mmax,kmax,Nalpha,Nm,Nk,Nrot)
		for ispec2=1:Nrot
			rotX_matrix[ispec2,ispec1,iprod]=coeff_end2[ispec2]
		end
	end
	for iprod=1:Nprody
		for i=1:Ngrid2
			coeff_out2[i]=coeff_grid2[i]*dy2[i,iprod]
		end
		#################################################
		##Transformation from angular grid to JMK basis##
		#################################################
		coeff_end2 .= transformation_backward_Wigner(coeff_out2,jmax,mmax,kmax,Nalpha,Nm,Nk,Nrot)
		for ispec2=1:Nrot
			rotY_matrix[ispec2,ispec1,iprod]=coeff_end2[ispec2]
		end
	end
	for iprod=1:Nprodz
		for i=1:Ngrid2
			coeff_out2[i]=coeff_grid2[i]*dz2[i,iprod]
		end
		#################################################
		##Transformation from angular grid to JMK basis##
		#################################################
		coeff_end2 .= transformation_backward_Wigner(coeff_out2,jmax,mmax,kmax,Nalpha,Nm,Nk,Nrot)
		for ispec2=1:Nrot
			rotZ_matrix[ispec2,ispec1,iprod]=coeff_end2[ispec2]
		end
	end
end

return Cx,Cy,Cz,transX_matrix,transY_matrix,transZ_matrix,rotX_matrix,rotY_matrix,rotZ_matrix
end
#######################################################################################################
function mu(R,theta1,phi1,theta2,phi2,chi)

dm = zeros(3)

dm[1] = sin(theta2)*cos(phi2)+R*sin(theta1)*cos(phi1)
dm[2] = sin(theta2)*sin(phi2)+R*sin(theta1)*sin(phi1)
dm[3] = cos(theta2)+R*cos(theta1)

#L = 1+R^2+2*R*(sin(theta1)*cos(phi1)*sin(theta2)*cos(phi2)+sin(theta1)*sin(phi1)*sin(theta2)*sin(phi2)+cos(theta1)*cos(theta2))
L = dm[1]^2+dm[2]^2+dm[3]^2


dm.=dm/sqrt(L)

return dm[1],dm[2],dm[3]
end
#######################################################################################################
function decomposition(matrix,error)

Ngrid1 = size(matrix,1)
Ngrid2 = size(matrix,2)

#Decompose potential by SVD#
F=svd(matrix)

#Determine the number of important singular values#
values=F.S/sum(F.S)

e=0.0
Nprod=0
while e <1.0-error
	Nprod+=1
	e=e+values[Nprod]
end
log_file=open("log","a")
println(log_file,"Number of products: ",Nprod)
println(log_file,"Chosen SV error: ",error)
println(log_file,"Sum of first Nprod SV: ",e)

#Determine the Frobenius norm of the difference between truncated SVD and initial potential#
mprod=zeros(ComplexF64,(Ngrid1,Ngrid2))
for i=1:Ngrid1
for j=1:Ngrid2
	prod=0.0
	for ip=1:Nprod
		prod=prod+F.U[i,ip]*F.S[ip]*F.Vt[ip,j]
	end
	mprod[i,j]=prod
end
end

diff=norm(abs.(mprod.-matrix), 2)
println(log_file,"Frobenius norm = ",diff)
println(log_file,"Fnorm divided by matrix size= ",diff/(Ngrid1*Ngrid2))
println(log_file)
close(log_file)

trans=zeros((Ngrid1,Nprod))
rot=zeros((Ngrid2,Nprod))
C=zeros(Nprod)	

for i=1:Ngrid1
	trans[i,1:Nprod].=F.U[i,1:Nprod]
end
for i=1:Ngrid2
	rot[i,1:Nprod].=F.Vt[1:Nprod,i]
end
C[:].=F.S[1:Nprod]

return C,trans,rot
end
#######################################################################################################
end
