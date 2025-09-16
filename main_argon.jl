using LinearAlgebra
using LinearMaps
using Arpack
push!(LOAD_PATH,pwd())
using kinetic
using potential
using trans_real_complex
using spin
using thermal_properties
using dipolemoment
using sym_break
using WignerSymbols
using DelimitedFiles


#Conversion factors#
a0toangst=0.529177210903 #Bohr to angstrom
eHtocm1=219474.631363	#Eh to cm-1
utome=1.0/(5.48579909065*1e-4)
eHtoJ=4.3597447222071e-18
eHtokHz=1.0/1.51982850071586e-13
a0tom=0.529177210903e-10 #Bohr to angstrom
a0topm=52.9177210903


#######################
#Matrix-vector product#
#######################

function Hv!(y::AbstractVector,x::AbstractVector)

	k=0
	for it=1:Ntrans
	for ir=1:Nrot
		k+=1
		wf_R[ir,it]=real(x[k])
		wf_I[ir,it]=imag(x[k])
		wf[ir,it]=x[k]
	end
	end

	one=1.0
	for i=1:Nrot
		trans_R = BLAS.symv('U',one, Ttrans_real,wf_R[i,:])
		trans_I = BLAS.symv('U',one, Ttrans_real,wf_I[i,:])
		Mv1[i,:] = trans_R + im*trans_I
	end
	for i=1:Ntrans
		rot_R = BLAS.symv('U',one, Trot_spin,wf_R[:,i])
		rot_I = BLAS.symv('U',one, Trot_spin,wf_I[:,i])
		Mv2[:,i] = rot_R + im*rot_I
	end

	Mv3 .= 0.0
	for ip=1:Nprod
		for i=1:Nrot
			trans_R = BLAS.symv('U',one, Vtrans[:,:,ip],wf_R[i,:])
			trans_I = BLAS.symv('U',one, Vtrans[:,:,ip],wf_I[i,:])
			Vv[i,:] = trans_R + im*trans_I
		end
		for i=1:Ntrans
			rot_R = BLAS.symv('U',Cpot[ip], Vrot[:,:,ip],real(Vv[:,i]))
			rot_I = BLAS.symv('U',Cpot[ip], Vrot[:,:,ip],imag(Vv[:,i]))
			Mv3[:,i] .+= rot_R + im*rot_I
		end
	end


	#Transform from NLM-basis to (R,theta,phi)-grid#
	# for i=1:Nrot
	# tempv_R = transformation_forward_HO(wf_R[i,:],nmax,NR,Ntheta,Nphi,nu,mass)
	# tempv_I = transformation_forward_HO(wf_I[i,:],nmax,NR,Ntheta,Nphi,nu,mass)
	# end
	# for i=1:Ntrans
	# 	wf_grid_R = transformation_forward_Wigner(tempv_R[:,i],jmax,mmax,kmax,Nalpha,Nm,Nk)
	# 	wf_grid_I = transformation_forward_Wigner(tempv_I[:,i],jmax,mmax,kmax,Nalpha,Nm,Nk)
	# end
	# 	#Multipliy transformed vectors with potential on grid#
	# 	for ig1=1:Ngrid1
	# 		for ig2=1:Ngrid2
	# 			coeff_grid[ig1,ig2] = Vpot[ig1,ig2]*coeff1_out[ig1]*coeff2_out[ig2]
	# 		end
	# 		#Transform from (Theta,phi,chi)-grid to JKM-basis#
	# 		coeff_tmp[ig1,:] .= transformation_backward_Wigner(coeff_grid[ig1,:],jmax,mmax,kmax,Nalpha,Nm,Nk,Nrot)
	# 	end	
	# 	#Transform from (R,theta,phi)-grid to NLM-basis#
	# 	for j=1:Nrot
	# 		coeff_end[:,j] = transformation_backward_HO(coeff_tmp[:,j],nmax,NR,Ntheta,Nphi,nu,mass,Ntrans)
	# 	end
	
	# 	for i=1:Ntrans
	# 	for j=1:Nrot
	# 		matrix[i,ispec1,j,ispec2] = coeff_end[i,j]
	# 	end
	# 	end
	# end

	

	k=0
	for it=1:Ntrans
	for ir=1:Nrot
		k+=1
		y[k]=Mv1[ir,it] + Mv2[ir,it] + Mv3[ir,it] # no SB

	end
	end

	return y
end 
function Hv_pn(x)

	y=zeros(ComplexF64,Ntrans*Nrot)

	k=0
	for it=1:Ntrans
	for ir=1:Nrot
		k+=1
		wf_R[ir,it]=real(x[k])
		wf_I[ir,it]=imag(x[k])
		wf[ir,it]=x[k]
	end
	end

	one=1.0
	for i=1:Nrot
		trans_R = BLAS.symv('U',one, Ttrans_real,wf_R[i,:])
		trans_I = BLAS.symv('U',one, Ttrans_real,wf_I[i,:])
		Mv1[i,:] = trans_R + im*trans_I
	end
	for i=1:Ntrans
		rot_R = BLAS.symv('U',one, Trot_spin,wf_R[:,i])
		rot_I = BLAS.symv('U',one, Trot_spin,wf_I[:,i])
		Mv2[:,i] = rot_R + im*rot_I
	end

	Mv3 .= 0.0
	for ip=1:Nprod
		for i=1:Nrot
			trans_R = BLAS.symv('U',one, Vtrans[:,:,ip],wf_R[i,:])
			trans_I = BLAS.symv('U',one, Vtrans[:,:,ip],wf_I[i,:])
			Vv[i,:] = trans_R + im*trans_I
		end
		for i=1:Ntrans
			rot_R = BLAS.symv('U',Cpot[ip], Vrot[:,:,ip],real(Vv[:,i]))
			rot_I = BLAS.symv('U',Cpot[ip], Vrot[:,:,ip],imag(Vv[:,i]))
			Mv3[:,i] .+= rot_R + im*rot_I
		end
	end
	
	for i=1:Ntrans
		rot_R = BLAS.symv('U',one,Vq1,wf_R[:,i])	
		rot_I = BLAS.symv('U',one,Vq1,wf_I[:,i])
		Mv5[:,i] = rot_R+im*rot_I
	end

	one2 = 1.0 + 0.0im
	Mv6 .= 0.0
	for m=-2:2
		for m2=-1:1
		for m1=-1:1
			for i=1:Nrot
				Vv2[i,:] = BLAS.gemv('N',one2, Y[:,:,2+m1],wf[i,:])
			end
			for i=1:Ntrans
				Mv6[:,i] += BLAS.gemv('N',wig[2+m1,2+m2,3+m], D1[:,:,2+m2],Vv2[:,i])
			end
		end
		end
	end

	k=0
	for it=1:Ntrans
	for ir=1:Nrot
		k+=1
		y[k]=Mv1[ir,it] + Mv2[ir,it] + Mv3[ir,it] # no SB
	end
	end

	return y
end

###################################################################
function matrix_prod_ABC(U,A,B)
	tmp = BLAS.gemm('C','N', U,A)
	tmp2 = BLAS.gemm('N','N', tmp,B)
	C = BLAS.gemm('N','N', tmp2,U)
	return C
end

###################################################################
function matrix_prod_AB(A,B)
	C = BLAS.gemm('N','N', A,B)
	return C
end

####################################
##Read information from input file##
####################################
f=open("input.txt") 
lines=readlines(f)
lsplit=split(lines[2])
jmax=parse(Int64, lsplit[2])
kmax=jmax
mmax=jmax
lsplit=split(lines[5])
Nalpha=parse(Int64, lsplit[2])
lsplit=split(lines[6])
Nm=parse(Int64, lsplit[2])
lsplit=split(lines[7])
Nk=parse(Int64, lsplit[2])
lsplit=split(lines[10])
Ae=parse(Float64, lsplit[2])
lsplit=split(lines[11])
Be=parse(Float64, lsplit[2])
lsplit=split(lines[12])
Ce=parse(Float64, lsplit[2])
Ae=Ae/eHtocm1
Be=Be/eHtocm1
Ce=Ce/eHtocm1
lsplit=split(lines[15])
nmax=parse(Int64, lsplit[2])
lsplit=split(lines[16])
mH=parse(Float64, lsplit[2])
lsplit=split(lines[17])
mO=parse(Float64, lsplit[2])
mass=(mO+2*mH)*utome
lsplit=split(lines[18])
kconst=parse(Float64, lsplit[2])
kconst=(kconst/eHtoJ)*(a0tom^2)
omega=sqrt(kconst/mass)
nu=mass*omega*0.5
lsplit=split(lines[19])
dCI=parse(Float64, lsplit[2])
dCI=dCI/a0topm


#lsplit=split(lines[22])
#mass=parse(Float64, lsplit[2])
#mass=mass*utome
#lsplit=split(lines[23])
#omega=parse(Float64, lsplit[2])
#omega=24.38/mass
#println("Denk' an die korrekte Einheit von omega")
#nu=mass*omega*0.5

lsplit=split(lines[22])
NR=parse(Int64, lsplit[2])
lsplit=split(lines[23])
Ntheta=parse(Int64, lsplit[2])
lsplit=split(lines[24])
Nphi=parse(Int64, lsplit[2])

lsplit=split(lines[27])
eq_struct=zeros(2)
eq_struct[1]=parse(Float64,lsplit[2])
eq_struct[1]=eq_struct[1]/a0toangst
lsplit=split(lines[28])
eq_struct[2]=parse(Float64,lsplit[2])
eq_struct[2]=eq_struct[2]*pi/180.0

lsplit=split(lines[31])
isomer=lsplit[2]

lsplit=split(lines[34])
Nstates=parse(Int64, lsplit[2])

lsplit=split(lines[37])
Tstart=parse(Float64, lsplit[2])
lsplit=split(lines[38])
dT=parse(Float64, lsplit[2])
lsplit=split(lines[39])
Ntemp=parse(Int64, lsplit[2])

lsplit=split(lines[42])
svd_err=parse(Float64, lsplit[2])

lsplit=split(lines[45])
vib_state=parse(Int64, lsplit[2])

lsplit=split(lines[48])
model_flag=parse(Int64, lsplit[2])
close(f)
model="AI"
if model_flag==0
	model="LO"
end
if model_flag==1
	model="AI"
end

#Check if Nm and Nk are odd#
if mod(Nm,2) == 0
	println("Even number of phi points (Euler angle)")
	exit()
elseif mod(Nk,2) == 0
	println("Even number of chi points (Euler angle)")
	exit()
elseif mod(Nphi,2) == 0
	println("Even number of phi points (3D-HO)")
	exit()
end

#Check if Nm/Nk is at least 2*mmax+1/2*kmax+1#
if Nm < 2*mmax+1
	println("Nphi must be at least 2*mmax+1 (Euler angles)")
	exit()
elseif Nk < 2*kmax+1
	println("Nchi must be at least 2*kmax+1 (Euler angles)")
	exit()
end

#Check if phi-grid is equal or larger than the spectral m-basis#
if Nphi < 2*nmax+1
	println("Increase Nphi, too small for chosen basis (3D-HO)")
	println("Try Nphi >= ", 2*nmax+1)
	exit()
end

#Basis size of the 3D-HO eigenbasis#
Ntrans=0
for n=0:nmax
	global Ntrans += Int(((n+1)*(n+2))/2)
end
#Basis size of the Wigner basis#
Nrot=0
for j=0:jmax
	mj=min(mmax,j)
	kj=min(kmax,j)
	global Nrot=Nrot+(2*mj+1)*(2*kj+1)
end

# test representtaion transformation
d_theta=0.0
d_chi=0.0
d_phi=0.0

#mat=change_rep(Nrot,jmax,Nalpha,Nm,Nk,d_theta,d_chi,d_phi)
#println(mat)

#exit()

#################################################################
#Calculate kinetic matrices for translational and rotational DOF#
#################################################################
#Ttrans = kinetic_trans_analytic(nmax,omega)
Ttrans = kinetic_translation(nmax,omega,mass)
H3D=H3D_translation(nmax,omega,mass)
Trot,Jx,Jz = kinetic_rotation(jmax,mmax,kmax,Ae,Be,Ce)

#Project kinetic matrix to spin sub space#
Tpara = spin_isomer("para",jmax,Trot)
Tortho = spin_isomer("ortho",jmax,Trot)
Npara = size(Tpara,1)
Northo = size(Tortho,1)

Identity=0.0*Trot

Nhalf=0
for j=0:jmax
	mj=min(mmax,j)
	kj=min(kmax,j)
	global Nhalf+=(2*mj+1)*(2*kj+1)
end

for i=1:Nhalf
	Identity[i,i]+=1.0
end
T_half=spin_isomer_half("para",jmax,Identity)

T_half_t=transpose(T_half)

Tparatest_tmp=T_half_t*Trot

Tparatest=Tparatest_tmp*T_half

println()

#show(Tparatest-Tpara)
#exit()

##Rotations##
#Transform kinetic matrix to real basis#
Upara = trans_realwigner(jmax,"para")
Tpara_spin = real(transform_realbasis(Upara,Tpara))
Uortho = trans_realwigner(jmax,"ortho")
Tortho_spin = real(transform_realbasis(Uortho,Tortho))

#Diagonalize free rotor
e_Tpara,ev_Tpara=eigen(Tpara)
e_Tpara_spin,ev_Tpara_spin=eigen(Tpara_spin)

e_Tortho,ev_Tortho=eigen(Tortho)
e_Tortho_spin,ev_Tortho_spin=eigen(Tortho_spin)

##Translations##
Utrans_trans = trans_realSPH(nmax,Ntrans)
Ttrans_real = real(transform_realbasis(Utrans_trans,Ttrans))
H3D_real = real(transform_realbasis(Utrans_trans,H3D))

f=open("log","w")
println(f,"#######################################")
println(f,"###########Basis information###########")
println(f,"#######################################")
println(f,"#Rotation#")
println(f,"jmax= ",jmax,", mmax= ",mmax,", kmax= ",kmax)
println(f,"Ntheta= ",Nalpha,", Nphi= ",Nm,", Nchi= ",Nk)
println(f,"Dimension of local Hilbert space: ",Npara+Northo)
println(f,"Dimension of para Hilbert space: ",Npara)
println(f,"Dimension of ortho Hilbert space: ",Northo)
println(f,"#Translation#")
println(f,"nmax= ",nmax,", lmax= ",nmax,", mmax= ",nmax)
println(f,"NR= ",NR,", Ntheta= ",Ntheta,", Nphi= ",Nphi)
println(f,"mass= ",mass/utome," omega= ",omega)
println(f,"Dimension of local Hilbert space: ",Ntrans)
println(f,"#######################################")
println(f,"#####Rotational constants in cm^-1#####")
println(f,"#######################################")
println(f,"Ae= ",Ae*eHtocm1 ,", Be= ",Be*eHtocm1 ,", Ce= ",Ce*eHtocm1)
println(f,"#######################################")
println(f)
println(f,"Calculate cage potential matrix")
close(f)

#Calculate cage potential in product basis#
Cpot,Vtrans_matrix,Vrot_matrix = potential_matrix(nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,eq_struct,svd_err,dCI,kconst,vib_state,model)
#Cpot,Vtrans_matrix,Vrot_matrix = potential_matrix(nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,eq_struct,svd_err,dCI,kconst,vib_state,model)
trig_matrix=trigo_euler_matrix(jmax,Nalpha,Nm,Nk,Nrot)
# 1 cϕsθ Jx
# 2 sϕsθ Jx
# 3 cθ Jx
# 4 cϕcθcχ−sϕsχ Jz
# 5 sϕcθcχ+ cϕsχ Jz
# 6 −sθcχ Jz

hsr1=matrix_prod_AB(trig_matrix[1,:,:],Jx)
hsr2=matrix_prod_AB(trig_matrix[2,:,:],Jx)
hsr3=matrix_prod_AB(trig_matrix[3,:,:],Jx)
hsr4=matrix_prod_AB(trig_matrix[4,:,:],Jz)
hsr5=matrix_prod_AB(trig_matrix[5,:,:],Jz)
hsr6=matrix_prod_AB(trig_matrix[6,:,:],Jz)

Nprod=size(Vtrans_matrix,3)
f=open("log","a")
println(f,"Cage potential matrix done")
close(f)

#Transform cage potential to real 3D-HO basis#
Vtrans = zeros(Float64,(Ntrans,Ntrans,Nprod))
Vpara = zeros(Float64,(Npara,Npara,Nprod))
Vortho = zeros(Float64,(Northo,Northo,Nprod))

for ip=1:Nprod
	#Transform 3D-HO to real basis#
	Vtrans[:,:,ip] .= real(transform_realbasis(Utrans_trans,Vtrans_matrix[:,:,ip]))
	
	#Transform to spin sub space and to real rotational eigenbasis#
	tmp = spin_isomer("para",jmax,Vrot_matrix[:,:,ip])
	Vpara[:,:,ip] = real(transform_realbasis(Upara,tmp))
	tmp = spin_isomer("ortho",jmax,Vrot_matrix[:,:,ip])
	Vortho[:,:,ip] = real(transform_realbasis(Uortho,tmp))
end

###############################################
#Calculate symmetry-breaking (quadrupole) term#
###############################################
mu=-0.737196
QBF=[1.53843,0,-0.09973,0,1.53843]
A=6.809e-6

#Vsb = Vquad(nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,A,QBF,mu)
term1,Y1m,wig,D1m0 = Vquad(nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,A,QBF,mu)

Vq1_p = zeros(Float64,(Npara,Npara))
Vq1_o = zeros(Float64,(Northo,Northo))

tmp_p = spin_isomer("para",jmax,term1)
tmp_o = spin_isomer("ortho",jmax,term1)
Vq1_p = real(transform_realbasis(Upara,tmp_p))
Vq1_o = real(transform_realbasis(Uortho,tmp_o))

D1_p = zeros(ComplexF64,(Npara,Npara,3))
D1_o = zeros(ComplexF64,(Northo,Northo,3))
for m=-1:1
	global tmp_p = spin_isomer("para",jmax,D1m0[:,:,2+m])
	global tmp_o = spin_isomer("ortho",jmax,D1m0[:,:,2+m])
	D1_p[:,:,2+m] = transform_realbasis(Upara,tmp_p)
	D1_o[:,:,2+m] = transform_realbasis(Uortho,tmp_o)
end

Y = zeros(ComplexF64,(Ntrans,Ntrans,3))
for m=-1:1
	Y[:,:,2+m] = transform_realbasis(Utrans_trans,Y1m[:,:,2+m])
end

f=open("energies_para.txt","w")
println(f,"#nmax= ",nmax," NR= ",NR," Ntheta= ",Ntheta," Nphi= ",Nphi)
println(f,"#jmax= ",jmax," Nalpha= ",Nalpha," Nphi= ",Nm," Nchi= ",Nk)
println(f,"#Ntrans= ",Ntrans)
println(f,"#Nrot= ",Npara)
println(f)
close(f)

f=open("energies_ortho.txt","w")
println(f,"#nmax= ",nmax," NR= ",NR," Ntheta= ",Ntheta," Nphi= ",Nphi)
println(f,"#jmax= ",jmax," Nalpha= ",Nalpha," Nphi= ",Nm," Nchi= ",Nk)
println(f,"#Ntrans= ",Ntrans)
println(f,"#Nrot= ",Northo)
println(f)
close(f)

f=open("heat_capacity.txt","w")
println(f,"#nmax= ",nmax," NR= ",NR," Ntheta= ",Ntheta," Nphi= ",Nphi)
println(f,"#jmax= ",jmax," Nalpha= ",Nalpha," Nphi= ",Nm," Nchi= ",Nk)
println(f,"#Ntrans= ",Ntrans)
println(f,"#Nrot= ",Nrot)
println(f,"#Npara= ",Npara)
println(f,"#Nportho= ",Northo)
println(f)
close(f)
#
#f=open("polarizability.txt","w")
#println(f,"nmax= ",nmax," NR= ",NR," Ntheta= ",Ntheta," Nphi= ",Nphi)
#println(f,"jmax= ",jmax," Nalpha= ",Nalpha," Nphi= ",Nm," Nchi= ",Nk)
#println(f,"Ntrans= ",Ntrans)
#println(f,"Nrot= ",Nrot2)
#println(f,"Spin isomer: ",isomer)
#println(f)
#close(f)


I2m=A*[im,im-1.0,0,im+1,-im]
let

	f=open("log","a")
	println(f,"Start diagonalization")

	println(f,"Para isomer:")
	close(f) 
	
	global Trot_spin = Tpara_spin
	global Tpara_spin_SR = Tpara_spin
	global Tortho_spin_SR = Tortho_spin
	global Vrot = Vpara
	global Nrot=Npara
	global wf = zeros(ComplexF64,(Npara,Ntrans))
	global wf_R = zeros(Float64,(Npara,Ntrans))
	global wf_I = zeros(Float64,(Npara,Ntrans))
	global Mv1 = zeros(ComplexF64,(Npara,Ntrans))
	global Mv2 = zeros(ComplexF64,(Npara,Ntrans))
	global Mv3 = zeros(ComplexF64,(Npara,Ntrans))
	global Mv5 = zeros(ComplexF64,(Npara,Ntrans))
	global Mv6 = zeros(ComplexF64,(Npara,Ntrans))
	global Vv = zeros(ComplexF64,(Npara,Ntrans))
	global Vv2 = zeros(ComplexF64,(Npara,Ntrans))

	global trans_R = zeros(Float64,(Ntrans))
	global trans_I = zeros(Float64,(Ntrans))
	global rot_R = zeros(Float64,(Npara))
	global rot_I = zeros(Float64,(Npara))
	global Vq1 = zeros(Float64,(Npara,Npara))
	global D1 = zeros(ComplexF64,(Npara,Npara,3))
	
	global Vq1 = Vq1_p
	global D1 = D1_p

	Nsize_para=Nstates

	D = LinearMap{ComplexF64}(Hv!,Ntrans*Nrot,ishermitian=true,ismutating=true)

	e_para,Wpara = Arpack.eigs(D,nev=Nstates,which=:SR)

	#rotational analysis
	println("rotational analysis")
	for state=1:Nstates
		println(real(e_para[state]-e_para[1])*eHtocm1)

		for n=1:Nrot
			ov=0.0	
			k=0
			for it=1:Ntrans
			for ir=1:Nrot
				k+=1
				ov+=real(ev_Tpara_spin[ir,n]*Wpara[k,state]*conj(ev_Tpara_spin[ir,n]*Wpara[k,state]))
			end
			end
			if ov>0.01
				println(n," ",e_Tpara_spin[n]*eHtocm1," ",ov)
			end
		end
		println()
	end
	#translational analysis
	e_t,ev_t=eigen(H3D_real)
	println("translational analysis")

	for state=1:Nstates
		println(real(e_para[state]-e_para[1])*eHtocm1)

		for n=1:Ntrans
			ov=0.0	
			k=0
			for it=1:Ntrans
			for ir=1:Nrot
				k+=1
				ov+=real(ev_t[it,n]*Wpara[k,state]*conj(ev_t[it,n]*Wpara[k,state]))
			end
			end
			if ov>0.01
				println(n," ",e_t[n]*eHtocm1," ",ov)
			end
		end
		println()
	end

	# explicit diag
	#Nsize_para=Ntrans*Nrot
	# Hexp= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))
	# v=zeros(ComplexF64,Ntrans*Nrot)
	# u=zeros(ComplexF64,Ntrans*Nrot)
	# for eye=1:Ntrans*Nrot
	# 	u[eye]=1.0+0.0im
	# 	for jay=1:Ntrans*Nrot
	# 		v[jay]=1.0+0.0im

	# 		Hexp[eye,jay]=dot(u,Hv_pn(v))
	# 		v[jay]=0.0+0.0im
	# 	end
	# 	u[eye]=0.0+0.0im
	# end

	#e_para,Wpara=eigen(Hexp)
	# for eye in e_para
	# 	println(real(eye*eHtocm1))
	# end
	
	f=open("log","a")
	println(f,"Ortho isomer:")
	close(f) 

	global Trot_spin = Tortho_spin
	global Vrot = Vortho
	global Nrot=Northo
	global wf = zeros(ComplexF64,(Northo,Ntrans))
	global wf_R = zeros(Float64,(Northo,Ntrans))
	global wf_I = zeros(Float64,(Northo,Ntrans))
	global Mv1 = zeros(ComplexF64,(Northo,Ntrans))
	global Mv2 = zeros(ComplexF64,(Northo,Ntrans))
	global Mv3 = zeros(ComplexF64,(Northo,Ntrans))
	global Mv5 = zeros(ComplexF64,(Northo,Ntrans))
	global Mv6 = zeros(ComplexF64,(Northo,Ntrans))
	global Vv = zeros(ComplexF64,(Northo,Ntrans))
	global Vv2 = zeros(ComplexF64,(Northo,Ntrans))

	global Vq1 = zeros(Float64,(Northo,Northo))
	global D1 = zeros(ComplexF64,(Northo,Northo,3))
	global rot_R = zeros(Float64,(Northo))
	global rot_I = zeros(Float64,(Northo))
	global Vq1 = Vq1_o
	global D1 =  D1_o

	D = LinearMap{ComplexF64}(Hv!,Ntrans*Nrot,ishermitian=true,ismutating=true)

	e_ortho,Wortho = Arpack.eigs(D,nev=Nstates,which=:SR)

	Nsize_ortho	= Nstates

	#rotational analysis
	println("rotational analysis")
	for state=1:Nstates
		println(real(e_ortho[state]-e_ortho[1])*eHtocm1)

		for n=1:Nrot
			ov=0.0	
			k=0
			for it=1:Ntrans
			for ir=1:Nrot
				k+=1
				ov+=real(ev_Tortho_spin[ir,n]*Wortho[k,state]*conj(ev_Tortho_spin[ir,n]*Wortho[k,state]))
			end
			end
			if ov>0.01
				println(n," ",e_Tortho_spin[n]*eHtocm1," ",ov)
			end
		end
		println()
	end
	#translational analysis
	e_t,ev_t=eigen(H3D_real)
	println("translational analysis")

	for state=1:Nstates
		println(real(e_ortho[state]-e_ortho[1])*eHtocm1)

		for n=1:Ntrans
			ov=0.0	
			k=0
			for it=1:Ntrans
			for ir=1:Nrot
				k+=1
				ov+=real(ev_t[it,n]*Wortho[k,state]*conj(ev_t[it,n]*Wortho[k,state]))
			end
			end
			if ov>0.01
				println(n," ",e_t[n]*eHtocm1," ",ov)
			end
		end
		println()
	end



	f=open("log","a")
	println(f,"Write eigenvalues")
	println(f)
	close(f)
	
	fp=open("energies_para.txt","w")
	fo=open("energies_ortho.txt","w")
	for istates=1:Nsize_para
		println(fp,round(real(e_para[istates])*eHtocm1,digits=18)," ",round(real(e_para[istates]-e_para[1])*eHtocm1,digits=16))
	end
	for istates=1:Nsize_ortho
		println(fo,round(real(e_ortho[istates])*eHtocm1,digits=18)," ",round(real(e_ortho[istates]-e_para[1])*eHtocm1,digits=16))
	end

	close(fp)
	close(fo)
	
	f=open("log","a")
	println(f,"Calculate heat capacity")
	f3=open("heat_capacity.txt","a")
	for it=1:Ntemp
		cv_p,cv_o,cv_full=heat_capacity(real(e_para),real(e_ortho),Tstart+(it-1)*dT)
		println(f3,round(Tstart+(it-1)*dT,digits=3),"  ",cv_p*eHtocm1,"  ",cv_o*eHtocm1,"  ",cv_full*eHtocm1)
	end
	close(f3)
	println(f,"Heat capacity done")
	println(f)
	close(f)

end