using LinearAlgebra
using LinearMaps
using Arpack
using Printf
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
using Random


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
		#trans_R = BLAS.symv('U',one, Ttrans_real,wf_R[i,:])
		#trans_I = BLAS.symv('U',one, Ttrans_real,wf_I[i,:])
		#Mv1[i,:] = trans_R + im*trans_I
		#Mv1[i,:] = BLAS.gemv('N', Ttrans,wf[i,:]) #mod Ttrans
		Mv1[i,:] = BLAS.hemv('U', Ttrans,wf[i,:]) #mod Ttrans
	end
	for i=1:Ntrans
		#rot_R = BLAS.symv('U',one, Trot_spin,wf_R[:,i])
		#rot_I = BLAS.symv('U',one, Trot_spin,wf_I[:,i])
		#Mv2[:,i] = rot_R + im*rot_I
		#Mv2[:,i] = BLAS.gemv('N', Trot_C,wf[:,i])
		Mv2[:,i] = BLAS.hemv('U', Trot_C,wf[:,i])
	end

	Mv3 .= 0.0
	for ip=1:Nprod
		for i=1:Nrot
			#trans_R = BLAS.symv('U',one, Vtrans[:,:,ip],wf_R[i,:])
			#trans_I = BLAS.symv('U',one, Vtrans[:,:,ip],wf_I[i,:])
			#Vv[i,:] = trans_R + im*trans_I
			Vv[i,:] = BLAS.hemv('U', Vtrans_matrix[:,:,ip],wf[i,:]) #mod
		end
		for i=1:Ntrans
			#rot_R = BLAS.symv('U',Cpot[ip], Vrot[:,:,ip],real(Vv[:,i]))
			#rot_I = BLAS.symv('U',Cpot[ip], Vrot[:,:,ip],imag(Vv[:,i]))
			#Mv3[:,i] .+= rot_R + im*rot_I
			Mv3[:,i] .+= BLAS.hemv('U',Cpot[ip]+0im, Vrot_C[:,:,ip],(Vv[:,i]))
		end
	end	



	k=0
	for it=1:Ntrans
	for ir=1:Nrot
		k+=1
		y[k]=Mv1[ir,it] + Mv2[ir,it] + Mv3[ir,it] 
		#y[k]=Mv1[ir,it] + Mv2[ir,it]  # no pot
		#y[k]=Mv2[ir,it]   # rot kin

		#noise=(1.0+0.01*(randn(rng)-0.5))
		# noise_mag=abs(noise)
		# norm_noise=noise/noise_mag
		#y[k]*=noise
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

function commutator(A,B)
	C= matrix_prod_AB(A,B)-matrix_prod_AB(B,A)

	n=size(A)[1]

	sum=0.0
	for eye=1:n
		for jay=1:n
			mag=abs(C[eye,jay])
			sum+=mag
		end
	end
	return sum
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
	Ae = 27.8806/eHtocm1;Be = 14.5216/eHtocm1;Ce=9.2778/eHtocm1
	kconst=(2.0/eHtoJ)*(a0tom^2)
	omega=sqrt(kconst/mass)
	nu=mass*omega*0.5
	if vib_state==4
		Ae = 26.6303/eHtocm1;Be = 14.4225/eHtocm1;Ce=9.1418/eHtocm1
	end
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
dummy,LX,LY,LZ=kinetic_trans_analytic(nmax,omega)
H3D=H3D_translation(nmax,omega,mass)
Trot,Jx,Jz,JX,JY,JZ = kinetic_rotation(jmax,mmax,kmax,Ae,Be,Ce)

# eLZ,evLZ=eigen(LZ)
# println("LZ eigenvalues: ",eLZ)
#Project kinetic matrix to spin sub space#
Tpara = spin_isomer("para",jmax,Trot)
JXpara=spin_isomer("para",jmax,JX)
JYpara=spin_isomer("para",jmax,JY)
JZpara=spin_isomer("para",jmax,JZ)
Tortho = spin_isomer("ortho",jmax,Trot)
JXortho = spin_isomer("ortho",jmax,JX)
JYortho = spin_isomer("ortho",jmax,JY)
JZortho = spin_isomer("ortho",jmax,JZ)

# irot=0
# irotNS=0
# for j=0:jmax
# 	for m=-j:j
# 		for k=-j:j
# 			global irotNS+=1
# 			if mod(k,2) == 0
# 				global irot+=1
# 				irotp=0
# 				irotNSp=0
# 				for jp=0:jmax
# 					for mp=-jp:jp
# 						for kp=-jp:jp
# 							irotNSp+=1
# 							if mod(kp,2) == 0
# 								irotp+=1	
# 								println(Trot[irotNS,irotNSp]," ",Tpara[irot,irotp])
# 							end
# 						end
# 					end
# 				end
# 			end
# 		end
# 	end
# end

# exit()




#Jxpara = spin_isomer("para",jmax,Jx)
#Jxortho = spin_isomer("ortho",jmax,Jx)
#Jzpara = spin_isomer("para",jmax,Jz)
#Jzortho = spin_isomer("ortho",jmax,Jz)
Npara = size(Tpara,1)
Northo = size(Tortho,1)
println()
# eJZpara,evLZ=eigen(JZpara)
# println("JZ eigenvalues: ",eJZpara)

# Tpara_mat=0.0*Tpara
# v=zeros(ComplexF64,Npara)
# u=zeros(ComplexF64,Npara)
# for eye=1:Npara
# 	u[eye]=1.0+0.0im
# 	for jay=1:Npara
# 		v[jay]=1.0+0.0im
# 		Tpara_mat[eye,jay]=dot(conj(u),BLAS.hemv('U', Tpara,v))
# 		v[jay]=0.0+0.0im
# 	end
# 	u[eye]=0.0+0.0im
# end
# AB= zeros(ComplexF64,(Npara,Npara))
# BA= zeros(ComplexF64,(Npara,Npara))

# mul!(AB,Tpara,Tpara_mat)
# mul!(BA,Tpara_mat,Tpara)
# C=AB-BA
# sum=0.0
# for eye=1:Npara
# 	for jay=1:Npara
# 		mag=abs(C[eye,jay])
# 		global sum+=mag
# 		if mag>1e-10
# 			println(mag)
# 		end
# 	end
# end
# println("para [Tpara_mat,Tpara] = ",sum)

# Tortho_mat=0.0*Tortho
# v=zeros(ComplexF64,Northo)
# u=zeros(ComplexF64,Northo)
# for eye=1:Northo
# 	u[eye]=1.0+0.0im
# 	for jay=1:Northo
# 		v[jay]=1.0+0.0im
# 		Tortho_mat[eye,jay]=dot(conj(u),BLAS.hemv('U', Tortho,v))
# 		v[jay]=0.0+0.0im
# 	end
# 	u[eye]=0.0+0.0im
# end

# AB= zeros(ComplexF64,(Northo,Northo))
# BA= zeros(ComplexF64,(Northo,Northo))

# mul!(AB,Tortho,Tortho_mat)
# mul!(BA,Tortho_mat,Tortho)
# C=AB-BA
# sum=0.0
# for eye=1:Northo
# 	for jay=1:Northo
# 		mag=abs(C[eye,jay])
# 		global sum+=mag
# 		if mag>1e-10
# 			println(mag)
# 		end
# 	end
# end
# println("ortho [Tortho_mat,Tortho] = ",sum)


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
T_half_para=spin_isomer_half("para",jmax,Identity)
T_half_ortho=spin_isomer_half("ortho",jmax,Identity)

T_half_para_t=transpose(T_half_para)

Tparatest_tmp=T_half_para_t*Trot

Tparatest=Tparatest_tmp*T_half_para



#show(Tparatest-Tpara)
#exit()

##Rotations##
#Transform kinetic matrix to real basis#
Upara = trans_realwigner(jmax,"para")
Tpara_spin = real(transform_realbasis(Upara,Tpara))
#Jxpara_spin = real(transform_realbasis(Upara,Jxpara))
#Jzpara_spin = real(transform_realbasis(Upara,Jzpara))
Uortho = trans_realwigner(jmax,"ortho")
Tortho_spin = real(transform_realbasis(Uortho,Tortho))
#Jxortho_spin = real(transform_realbasis(Uortho,Jxortho))
#Jzortho_spin = real(transform_realbasis(Uortho,Jzortho))

#Diagonalize free rotor
e_Tpara,ev_Tpara=eigen(Tpara)
e_Tpara_spin,ev_Tpara_spin=eigen(Tpara_spin)

e_Tortho,ev_Tortho=eigen(Tortho)
e_Tortho_spin,ev_Tortho_spin=eigen(Tortho_spin)

##Translations##
Utrans_trans = trans_realSPH(nmax,Ntrans)
Ttrans_real = real(transform_realbasis(Utrans_trans,Ttrans))
H3D_real = real(transform_realbasis(Utrans_trans,H3D))
e_t,ev_t=eigen(H3D_real)


path = "./JMax_$jmax"*"_NMax_$nmax/"
mkpath("./JMax_$jmax"*"_NMax_$nmax")

f=open(path*"log","w")
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
Cpot,Vtrans_matrix,Vrot_matrix = potential_matrix(nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,eq_struct,svd_err,dCI,kconst,vib_state,model, path)
#Cpot,Vtrans_matrix,Vrot_matrix = potential_matrix(nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,eq_struct,svd_err,dCI,kconst,vib_state,model)
trig_matrix=trigo_euler_matrix(jmax,Nalpha,Nm,Nk,Nrot)
#dipole_matrix=dipole_basis(jmax,Nalpha,Nm,Nk,Nrot)
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
f=open(path*"log","a")
println(f,"Cage potential matrix done")
close(f)

#Transform cage potential to real 3D-HO basis#
Vtrans = zeros(Float64,(Ntrans,Ntrans,Nprod))
Vpara = zeros(Float64,(Npara,Npara,Nprod))
Vortho = zeros(Float64,(Northo,Northo,Nprod))
Vpara_C = zeros(ComplexF64,(Npara,Npara,Nprod))
Vortho_C = zeros(ComplexF64,(Northo,Northo,Nprod))
mu_para = zeros(Float64,(Npara,Npara,3))
mu_ortho = zeros(Float64,(Northo,Northo,3))

for ip=1:Nprod
	#Transform 3D-HO to real basis#
	Vtrans[:,:,ip] .= real(transform_realbasis(Utrans_trans,Vtrans_matrix[:,:,ip]))
	
	#Transform to spin sub space and to real rotational eigenbasis#
	tmp = spin_isomer("para",jmax,Vrot_matrix[:,:,ip])
	Vpara_C[:,:,ip] = spin_isomer("para",jmax,Vrot_matrix[:,:,ip])
	Vpara[:,:,ip] = real(transform_realbasis(Upara,tmp))
	tmp = spin_isomer("ortho",jmax,Vrot_matrix[:,:,ip])
	Vortho_C[:,:,ip] = spin_isomer("ortho",jmax,Vrot_matrix[:,:,ip])
	Vortho[:,:,ip] = real(transform_realbasis(Uortho,tmp))
end

for ip=1:3
	tmp = spin_isomer("para",jmax,trig_matrix[ip,:,:])
	mu_para[:,:,ip] = real(transform_realbasis(Upara,tmp))
	tmp = spin_isomer("ortho",jmax,trig_matrix[ip,:,:])
	mu_ortho[:,:,ip] = real(transform_realbasis(Uortho,tmp))
end

f=open(path*"energies_para.txt","w")
println(f,"#nmax= ",nmax," NR= ",NR," Ntheta= ",Ntheta," Nphi= ",Nphi)
println(f,"#jmax= ",jmax," Nalpha= ",Nalpha," Nphi= ",Nm," Nchi= ",Nk)
println(f,"#Ntrans= ",Ntrans)
println(f,"#Nrot= ",Npara)
println(f)
close(f)

f=open(path*"energies_ortho.txt","w")
println(f,"#nmax= ",nmax," NR= ",NR," Ntheta= ",Ntheta," Nphi= ",Nphi)
println(f,"#jmax= ",jmax," Nalpha= ",Nalpha," Nphi= ",Nm," Nchi= ",Nk)
println(f,"#Ntrans= ",Ntrans)
println(f,"#Nrot= ",Northo)
println(f)
close(f)

f=open(path*"heat_capacity.txt","w")
println(f,"#nmax= ",nmax," NR= ",NR," Ntheta= ",Ntheta," Nphi= ",Nphi)
println(f,"#jmax= ",jmax," Nalpha= ",Nalpha," Nphi= ",Nm," Nchi= ",Nk)
println(f,"#Ntrans= ",Ntrans)
println(f,"#Nrot= ",Nrot)
println(f,"#Npara= ",Npara)
println(f,"#Nportho= ",Northo)
println(f)
close(f)


let

	f=open(path*"log","a")
	println(f,"Start diagonalization")

	println(f,"Para isomer:")
	close(f) 
	
	global Trot_spin = Tpara_spin
	global Trot_C = Tpara
	global Vrot = Vpara
	global Vrot_C = Vpara_C
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
	global rng = Xoshiro(1234);

	D = LinearMap{ComplexF64}(Hv!,Ntrans*Nrot,ishermitian=true,ismutating=true)

	e_para,Wpara = Arpack.eigs(D,nev=Nstates,which=:SR)

	global Nstates=size(e_para)[1]
	#JZ+LZ expectation

	JX_plus_LX=zeros(ComplexF64,(Npara*Ntrans,Npara*Ntrans))
	JY_plus_LY=zeros(ComplexF64,(Npara*Ntrans,Npara*Ntrans))
	JZ_plus_LZ=zeros(ComplexF64,(Npara*Ntrans,Npara*Ntrans))
	JX_mat=zeros(ComplexF64,(Npara*Ntrans,Npara*Ntrans))
	JY_mat=zeros(ComplexF64,(Npara*Ntrans,Npara*Ntrans))
	JZ_mat=zeros(ComplexF64,(Npara*Ntrans,Npara*Ntrans))
	LZ_mat=zeros(ComplexF64,(Npara*Ntrans,Npara*Ntrans))
	LX_mat=zeros(ComplexF64,(Npara*Ntrans,Npara*Ntrans))
	LY_mat=zeros(ComplexF64,(Npara*Ntrans,Npara*Ntrans))
	Trot_test=zeros(ComplexF64,(Npara*Ntrans,Npara*Ntrans))
	Ttrans_test=zeros(ComplexF64,(Npara*Ntrans,Npara*Ntrans))

	#for i=1:Npara
	#	println(real(JZpara[i,:]))
	#end

	for it=1:Ntrans
		for ir=1:Npara
			for itp=1:Ntrans
				Ttrans_test[(it-1)*Npara+ir,(itp-1)*Npara+ir]=Ttrans[it,itp]
			end
		end
	end

	i=0
	for n1=0:nmax
		for l1=n1:-2:0
			for m1=-l1:l1
				i+=1
				for irot=1:Npara
					for irotp=1:Npara
						if irot==irotp
							JZ_plus_LZ[(i-1)*(Npara)+irot,(i-1)*(Npara)+irotp]=m1+JZpara[irot,irot]
							JZ_mat[(i-1)*(Npara)+irot,(i-1)*(Npara)+irotp]=JZpara[irot,irot]
							LZ_mat[(i-1)*(Npara)+irot,(i-1)*(Npara)+irotp]=m1
						end
						Trot_test[(i-1)*(Npara)+irot,(i-1)*(Npara)+irotp]=Trot_C[irot,irotp]
					end
				end
			end
		end
	end

	i=0
	for n1=0:nmax
		for l1=n1:-2:0
			for m1=-l1:l1
				i+=1
				for irot=1:Npara
					for irotp=1:Npara
						JX_plus_LX[(i-1)*(Npara)+irot,(i-1)*(Npara)+irotp]=JXpara[irot,irotp]
						JX_mat[(i-1)*(Npara)+irot,(i-1)*(Npara)+irotp]=JXpara[irot,irotp]
						JY_plus_LY[(i-1)*(Npara)+irot,(i-1)*(Npara)+irotp]=JYpara[irot,irotp]
						JY_mat[(i-1)*(Npara)+irot,(i-1)*(Npara)+irotp]=JYpara[irot,irotp]
					end
				end
			end
		end
	end

	i1=0
	for n1=0:nmax
		for l1=n1:-2:0
			for m1=-l1:l1
				i1+=1
				i2=0
				for n2=0:nmax
					for l2=n2:-2:0
						for m2=-l2:l2
							i2+=1
							if n1 == n2 && l1 == l2 
								for irot=1:Npara
									JX_plus_LX[(i1-1)*(Npara)+irot,(i2-1)*(Npara)+irot]+=LX[i1,i2]
									JY_plus_LY[(i1-1)*(Npara)+irot,(i2-1)*(Npara)+irot]+=LY[i1,i2]
									LX_mat[(i1-1)*(Npara)+irot,(i2-1)*(Npara)+irot]=LX[i1,i2]
									LY_mat[(i1-1)*(Npara)+irot,(i2-1)*(Npara)+irot]=LY[i1,i2]
								end
							end
						end
					end
				end
			end
		end
	end


	#println(sum," ",norm)
	#println()

	sum=zeros(Float64,(3,3))
	for s=2:4
		for sp=2:4
			for i=1:Ntrans*Nrot
				sum[s-1,sp-1]+=real(conj(Wpara[i,s])*JZ_plus_LZ[i,i]*Wpara[i,sp])*(1.0+0.000001*(rand()-0.5))
			end
		end
	end

	println("<i|LZ+JZ|j>")
	for i=1:3
		for j=1:3
			print(round(real(sum[i,j]),digits=5),"   ")
		end
		println()
	end
	e_sum,ev_sum=eigen(sum)
	println(e_sum)

	
	#sum+=(Wpara[(itrans-1)*Npara+irot,state]*conj(Wpara[(itrans-1)*Npara+irot,state]))*(m1+(JZpara[irot,irot]))
	# 					norm+=(Wpara[(itrans-1)*Npara+irot,state]*conj(Wpara[(itrans-1)*Npara+irot,state]))
	# 				end
	# 			end
	# 		end
	# 	end
	# 	println("<JZ+LZ> ",sum," ",norm)
	# 	println()
	# end


	#rotational analysis
	f = open(path*"rotational_wavefunction_para.txt", "w")
	f_rt = open(path*"rot_trans_wavefunction_para.txt", "w")
	println(f, "rotational analysis para")
	println(f,"#nmax= ",nmax," NR= ",NR," Ntheta= ",Ntheta," Nphi= ",Nphi)
	println(f,"#jmax= ",jmax," Nalpha= ",Nalpha," Nphi= ",Nm," Nchi= ",Nk)

	for state=1:Nstates
		println(f, real(e_para[state]-e_para[1])*eHtocm1)

		for n=1:Nrot
			ov=0.0	
			k=0
			for it=1:Ntrans
				for ir=1:Nrot
					k+=1
					ov+=ev_Tpara_spin[ir,n]*conj(Wpara[k,state])
				end
			end
			ov2=real(conj(ov)*ov)

			if ov2>0.01
				println(f, n," ",e_Tpara_spin[n]*eHtocm1," ",ov2)
			end
			# joint analysis
			for nt=1:Ntrans
				ov_rt=0.0	
				k=0
				for it=1:Ntrans
					for ir=1:Nrot
						k+=1
						ov_rt+=ev_Tpara_spin[ir,n]*ev_t[it,nt]*conj(Wpara[k,state])
					end
				end
				ov_rt2=real(conj(ov_rt)*ov_rt)
				if ov_rt2>0.01
					s=String[]
					sf=@sprintf("%-5d %-7d %-10g %-10g %-10g",n,nt,round(e_Tpara_spin[n]*eHtocm1,digits=2),round(e_t[nt]*eHtocm1,digits=2),round(ov_rt2,digits=2))
					push!(s, sf)
					println(f_rt,s[1])
				end
			end
		end
		println(f)
		println(f_rt)
	end
	close(f)
	close(f_rt)
	#translational analysis
	f = open(path*"translational_wavefunction_para.txt", "w")

	println(f, "translational analysis para")
	println(f,"#nmax= ",nmax," NR= ",NR," Ntheta= ",Ntheta," Nphi= ",Nphi)
	println(f,"#jmax= ",jmax," Nalpha= ",Nalpha," Nphi= ",Nm," Nchi= ",Nk)

	for state=1:Nstates
		println(f, real(e_para[state]-e_para[1])*eHtocm1)

		for n=1:Ntrans
			ov=0.0	
			k=0
			for it=1:Ntrans
				for ir=1:Nrot
					k+=1
					ov+=ev_t[it,n]*conj(Wpara[k,state])
				end
			end
			ov2=real(conj(ov)*ov)
			if ov2>0.01
				println(f, n," ",e_t[n]*eHtocm1," ",ov2)
			end
		end
		println(f)
	end
	close(f)
	testQN=0
	if testQN==0
		# explicit diag
		Nsize_para=Ntrans*Nrot
		Hexp= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))
		v=zeros(ComplexF64,Ntrans*Nrot)
		u=zeros(ComplexF64,Ntrans*Nrot)
		w=zeros(ComplexF64,Ntrans*Nrot)
		for eye=1:Ntrans*Nrot
			u[eye]=1.0+0.0im
			for jay=1:Ntrans*Nrot
				v[jay]=1.0+0.0im
				Hexp[eye,jay]=dot(conj(u),Hv!(w,v))*(1.0+0.000001*(rand()-0.5))
				#Hexp[eye,jay]=dot(conj(u),Hv!(w,v))
				v[jay]=0.0+0.0im
			end
			u[eye]=0.0+0.0im
		end
		e_para,ev_para=eigen(Hexp)
		temp_mat1= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))
		temp_mat2= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))
		mul!(temp_mat1,JX_plus_LX,ev_para)
		mul!(temp_mat2,conj(transpose(ev_para)),temp_mat1)
		println("<i|LX+JX|j>")
		for i=1:10
			for j=1:10
				print(round(real(temp_mat2[i,j]),digits=2),"   ")
			end
			println()
		end

		sum2=zeros(Float64,(3,3))
		for i=2:4
			for j=2:4
				sum2[i-1,j-1]=real(temp_mat2[i,j])
			end
		end

		e_2,ev_2=eigen(sum2)
		println(e_2)

		temp_mat1= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))
		temp_mat2= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))
		mul!(temp_mat1,JY_plus_LY,ev_para)
		mul!(temp_mat2,conj(transpose(ev_para)),temp_mat1)
		println("<i|LY+JY|j>")
		for i=1:10
			for j=1:10
				print(round(real(temp_mat2[i,j]),digits=2),"   ")
			end
			println()
		end

		sum2=zeros(Float64,(3,3))
		for i=2:4
			for j=2:4
				sum2[i-1,j-1]=real(temp_mat2[i,j])
			end
		end

		e_2,ev_2=eigen(sum2)
		println(e_2)

		temp_mat1= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))
		temp_mat2= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))
		mul!(temp_mat1,JZ_plus_LZ,ev_para)
		mul!(temp_mat2,conj(transpose(ev_para)),temp_mat1)
		println("<i|LZ+JZ|j>")
		for i=1:10
			for j=1:10
				print(round(real(temp_mat2[i,j]),digits=2),"   ")
			end
			println()
		end

		sum2=zeros(Float64,(3,3))
		for i=2:4
			for j=2:4
				sum2[i-1,j-1]=real(temp_mat2[i,j])
			end
		end

		e_2,ev_2=eigen(sum2)
		println(e_2)

	
		println("para [JX+LX,H] = "," ",commutator(JX_plus_LX,Hexp))
		println("para [JY+LY,H] = "," ",commutator(JY_plus_LY,Hexp))
		println("para [JZ+LZ,H] = "," ",commutator(JZ_plus_LZ,Hexp))
		println()
		println("para [JX,H] = "," ",commutator(JX_mat,Hexp))
		println("para [JY,H] = "," ",commutator(JY_mat,Hexp))
		println("para [JZ,H] = "," ",commutator(JZ_mat,Hexp))
		println()
		println("para [LX,H] = "," ",commutator(LX_mat,Hexp))
		println("para [LY,H] = "," ",commutator(LY_mat,Hexp))
		println("para [LZ,H] = "," ",commutator(LZ_mat,Hexp))
		println()


		println("para [JX+LX,Trot+Ttrans] = "," ",commutator(JX_plus_LX,Trot_test+Ttrans_test))
		println("para [JY+LY,Trot+Ttrans] = "," ",commutator(JY_plus_LY,Trot_test+Ttrans_test))
		println("para [JZ+LZ,Trot+Ttrans] = "," ",commutator(JZ_plus_LZ,Trot_test+Ttrans_test))
		println()

		println("para [JX,Trot+Ttrans] = "," ",commutator(JX_mat,Trot_test+Ttrans_test))
		println("para [JY,Trot+Ttrans] = "," ",commutator(JY_mat,Trot_test+Ttrans_test))
		println("para [JZ,Trot+Ttrans] = "," ",commutator(JZ_mat,Trot_test+Ttrans_test))
		println()

		println("para [LX,Trot+Ttrans] = "," ",commutator(LX_mat,Trot_test+Ttrans_test))
		println("para [LY,Trot+Ttrans] = "," ",commutator(LY_mat,Trot_test+Ttrans_test))
		println("para [LZ,Trot+Ttrans] = "," ",commutator(LZ_mat,Trot_test+Ttrans_test))
	
	end

	#e_para,Wpara=eigen(Hexp)
	# for eye in e_para
	# 	println(real(eye*eHtocm1))
	# end
	
	f=open(path*"log","a")
	println(f,"Ortho isomer:")
	close(f) 

	global Trot_spin = Tortho_spin
	global Trot_C = Tortho
	global Vrot_C = Vortho_C
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

	global rot_R = zeros(Float64,(Northo))
	global rot_I = zeros(Float64,(Northo))

# test JZ+LZ
	JZ_plus_LZ=zeros(ComplexF64,(Northo*Ntrans,Northo*Ntrans))
	Trot_test=zeros(ComplexF64,(Northo*Ntrans,Northo*Ntrans))

	i=0
	for n1=0:nmax
		for l1=n1:-2:0
			for m1=-l1:l1
				i+=1
				for irot=1:Northo
					for irotp=1:Northo
						if irot==irotp
							JZ_plus_LZ[(i-1)*(Northo)+irot,(i-1)*(Northo)+irotp]=m1
						end
						JZ_plus_LZ[(i-1)*(Northo)+irot,(i-1)*(Northo)+irotp]+=JZortho[irot,irotp]

						Trot_test[(i-1)*(Northo)+irot,(i-1)*(Northo)+irotp]=Tortho[irot,irotp]
					end
				end
			end
		end
	end



	Nsize_ortho=Ntrans*Nrot

	# Hexp= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))
	# v=zeros(ComplexF64,Ntrans*Nrot)
	# u=zeros(ComplexF64,Ntrans*Nrot)
	# w=zeros(ComplexF64,Ntrans*Nrot)
	# for eye=1:Ntrans*Nrot
	# 	u[eye]=1.0+0.0im
	# 	for jay=1:Ntrans*Nrot
	# 		v[jay]=1.0+0.0im

	# 		Hexp[eye,jay]=dot(conj(u),Hv!(w,v))
	# 		v[jay]=0.0+0.0im
	# 	end
	# 	u[eye]=0.0+0.0im
	# end

	# AB= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))
	# BA= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))

	# mul!(AB,JZ_plus_LZ,Hexp)
	# mul!(BA,Hexp,JZ_plus_LZ)
	# C=AB-BA

	# sum=0.0
	# for eye=1:Ntrans*Nrot
	# 	for jay=1:Ntrans*Nrot
	# 		mag=abs(C[eye,jay])
	# 		sum+=mag
	# 		# if mag>1e-10
	# 		# 	println(mag)
	# 		# end
	# 	end
	# end
	# println("ortho [JZ+LZ,H] = ",sum)

	# exit()


	D = LinearMap{ComplexF64}(Hv!,Ntrans*Nrot,ishermitian=true,ismutating=true)

	e_ortho,Wortho = Arpack.eigs(D,nev=Nstates,which=:SR)

	Nstates=size(e_ortho)[1]

	#rotational analysis
	f = open(path*"rotational_wavefunction_ortho.txt", "w")
	f_rt = open(path*"rot_trans_wavefunction_ortho.txt", "w")
	println(f, "rotational analysis ortho")
	println(f,"#nmax= ",nmax," NR= ",NR," Ntheta= ",Ntheta," Nphi= ",Nphi)
	println(f,"#jmax= ",jmax," Nalpha= ",Nalpha," Nphi= ",Nm," Nchi= ",Nk)
	for state=1:Nstates
		println(f, real(e_ortho[state]-e_para[1])*eHtocm1)
		println(f_rt, "state: ",state," E_state =",real(e_ortho[state]-e_para[1])*eHtocm1)
		println(f_rt,"n_rot n_trans E_rot      E_trans    |<state|n_rot n_trans>|^2")

		for n=1:Nrot
			ov=0.0	
			k=0
			for it=1:Ntrans
				for ir=1:Nrot
					k+=1
					ov+=ev_Tortho_spin[ir,n]*conj(Wortho[k,state])
				end
			end
			ov2=real(conj(ov)*ov)
			if ov2>0.01
				println(f, n," ",e_Tortho_spin[n]*eHtocm1," ",ov2)
			end

			# joint analysis
			for nt=1:Ntrans
				ov_rt=0.0	
				k=0
				for it=1:Ntrans
					for ir=1:Nrot
						k+=1
						#ov_rt+=real(ev_Tortho_spin[ir,n]*ev_t[it,nt]*conj(Wortho[k,state])*conj(ev_t[it,nt]*ev_Tortho_spin[ir,n]*conj(Wortho[k,state])))
						ov_rt+=ev_Tortho_spin[ir,n]*ev_t[it,nt]*conj(Wortho[k,state])
					end
				end
				ov_rt2=real(conj(ov_rt)*ov_rt)
				if ov_rt2>0.01
					s=String[]
					sf=@sprintf("%-5d %-7d %-10g %-10g %-10g",n,nt,round(e_Tortho_spin[n]*eHtocm1,digits=2),round(e_t[nt]*eHtocm1,digits=2),round(ov_rt2,digits=2))
					push!(s, sf)
					println(f_rt,s[1])
				end
			end
		end	
		println(f_rt)
		println(f)
	end
	close(f)
	close(f_rt)
	#translational analysis
	f = open(path*"translational_wavefunction_ortho.txt", "w")
	println(f, "translational analysis ortho")
	println(f,"#nmax= ",nmax," NR= ",NR," Ntheta= ",Ntheta," Nphi= ",Nphi)
	println(f,"#jmax= ",jmax," Nalpha= ",Nalpha," Nphi= ",Nm," Nchi= ",Nk)

	for state=1:Nstates
		println(f, real(e_ortho[state]-e_para[1])*eHtocm1)

		for n=1:Ntrans
			ov=0.0	
			k=0
			for it=1:Ntrans
				for ir=1:Nrot
					k+=1
					ov+=ev_t[it,n]*conj(Wortho[k,state])
				end
			end
			ov2=real(conj(ov)*ov)
			if ov2>0.01
				println(f, n," ",e_t[n]*eHtocm1," ",ov2)
			end
		end
		println(f)
	end
	close(f)




	f=open(path*"log","a")
	println(f,"Write eigenvalues")
	println(f)
	close(f)
	
	fp=open(path*"energies_para.txt","a")
	fo=open(path*"energies_ortho.txt","a")
	for istates=1:Nstates
		println(fp,round(real(e_para[istates])*eHtocm1,digits=18)," ",round(real(e_para[istates]-e_para[1])*eHtocm1,digits=16)," ",round(real(e_para[istates]-0.5*kconst*dCI*dCI)*eHtocm1,digits=18))
	end
	for istates=1:Nstates
		println(fo,round(real(e_ortho[istates])*eHtocm1,digits=5)," ",round(real(e_ortho[istates]-e_para[1])*eHtocm1,digits=5)," ",round(real(e_ortho[istates]-0.5*kconst*dCI*dCI)*eHtocm1,digits=5))
	end

	close(fo)
	close(fp)

	spectrum=0
	Emax=100.0/eHtocm1
	DeltaEmax=220.0/eHtocm1
	

	if spectrum==0

		ftp=open(path*"transitions_para.txt","w")
		fto=open(path*"transitions_ortho.txt","w")
		for istates=1:Nstates
			for jstates=(istates):Nstates
				deltaE=real(e_para[jstates]-e_para[istates])
				if real(e_para[istates]-e_para[1])<Emax && deltaE<DeltaEmax
					mu_element = zeros(ComplexF64,(3))
					k=0
					for it=1:Ntrans
						for ir=1:Npara
							k+=1
							kp=0
							for itp=1:Ntrans
								for irp=1:Npara
									kp+=1
									if it==itp
										for a=1:3
											mu_element[a]+=conj(Wpara[k,jstates])*mu_para[ir,irp,a]*Wpara[kp,istates]
										end
									 end
								end
							end
						end
					end
					dipole_mag=0.0
					for a=1:3
						dipole_mag+=real(conj(mu_element[a])*mu_element[a])
					end
					println(ftp,istates," ",round(real(e_para[istates]-e_para[1])*eHtocm1,digits=2)," ",jstates,"  ",round(deltaE*eHtocm1,digits=2)," ",dipole_mag)
				end
			end
		end
		close(ftp)

		for istates=1:Nstates
			for jstates=(istates):Nstates
				deltaE=real(e_ortho[jstates]-e_ortho[istates])
				if real(e_ortho[istates]-e_para[1])<Emax && deltaE<DeltaEmax
				# dipole elementsk=0
					mu_element = zeros(ComplexF64,(3))
					k=0
					for it=1:Ntrans
						for ir=1:Northo
							k+=1
							kp=0
							for itp=1:Ntrans
								for irp=1:Northo
									kp+=1
									if it==itp
										for a=1:3
											mu_element[a]+=conj(Wortho[k,istates])*mu_ortho[ir,irp,a]*Wortho[kp,jstates]
										end
									end
								end
							end
						end
					end
					dipole_mag=0.0
					for a=1:3
						dipole_mag+=real(conj(mu_element[a])*mu_element[a])
					end
					println(fto,istates," ",round(real(e_ortho[istates]-e_para[1])*eHtocm1,digits=2)," ",jstates," ",round(deltaE*eHtocm1,digits=2)," ",dipole_mag)		
				end
			end
		end
		close(fto)
	end
	
	f=open(path*"log","a")
	println(f,"Calculate heat capacity")
	f3=open(path*"heat_capacity.txt","a")
	for it=1:Ntemp
		cv_p,cv_o,cv_full=heat_capacity(real(e_para),real(e_ortho),Tstart+(it-1)*dT)
		println(f3,round(Tstart+(it-1)*dT,digits=3),"  ",cv_p*eHtocm1,"  ",cv_o*eHtocm1,"  ",cv_full*eHtocm1)
	end
	close(f3)
	println(f,"Heat capacity done")
	println(f)
	close(f)

end