using LinearAlgebra
using LinearMaps
using Arpack
push!(LOAD_PATH,pwd())
using kinetic
using potential_nosvd
#using potential
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

	#
	


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

	k=0
	for it=1:Ntrans
	for ir=1:Nrot
		k+=1
		y[k]=Mv1[ir,it] + Mv2[ir,it]  # no V
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

#################################################################
#Calculate kinetic matrices for translational and rotational DOF#
#################################################################
#Ttrans = kinetic_trans_analytic(nmax,omega)
Ttrans = kinetic_translation(nmax,omega,mass)
H3D=H3D_translation(nmax,omega,mass)
Trot,Jx,Jz = kinetic_rotation(jmax,mmax,kmax,Ae,Be,Ce)

f=open("log","w")
println(f,"#######################################")
println(f,"###########Basis information###########")
println(f,"#######################################")
println(f,"#Rotation#")
println(f,"jmax= ",jmax,", mmax= ",mmax,", kmax= ",kmax)
println(f,"Ntheta= ",Nalpha,", Nphi= ",Nm,", Nchi= ",Nk)
println(f,"#Translation#")
println(f,"nmax= ",nmax,", lmax= ",nmax,", mmax= ",nmax)
println(f,"NR= ",NR,", Ntheta= ",Ntheta,", Nphi= ",Nphi)
println(f,"mass= ",mass/utome," omega= ",omega)
println(f,"Dimension of local Hilbert space: ",Ntrans*Nrot)
println(f,"#######################################")
println(f,"#####Rotational constants in cm^-1#####")
println(f,"#######################################")
println(f,"Ae= ",Ae*eHtocm1 ,", Be= ",Be*eHtocm1 ,", Ce= ",Ce*eHtocm1)
println(f,"#######################################")
println(f)
println(f,"Calculate cage potential matrix")
close(f)

#Calculate cage potential in product basis#
Vmatrix = potential_matrix(nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,eq_struct,dCI,kconst,vib_state,model)


f=open("log","a")
println(f,"Cage potential matrix done")
close(f)


let
	
	#explicit diag
	Hexp= zeros(ComplexF64,(Ntrans*Nrot,Ntrans*Nrot))
	v=zeros(ComplexF64,Ntrans*Nrot)
	u=zeros(ComplexF64,Ntrans*Nrot)
	for eye_t=1:Ntrans
		for eye_r=1:Nrot
			eye=(eye_t-1)*Nrot+eye_r
			u[eye]=1.0+0.0im
			for jay_t=1:Ntrans
				for jay_r=1:Nrot
					jay=(jay_t-1)*Nrot+jay_r
					v[jay]=1.0+0.0im
					if eye_r==jay_r
						Hexp[eye,jay]+=Ttrans[eye_t,jay_t]
					end
					if eye_t==jay_t
						Hexp[eye,jay]+=Trot[eye_r,jay_r]
					end
					Hexp[eye,jay]+=Vmatrix[eye_t,jay_t,eye_r,jay_r]
					v[jay]=0.0+0.0im
				end
			end
			u[eye]=0.0+0.0im
		end
	end

	e,W=eigen(Hexp)
	for eye in e
		println(real((eye-e[1])*eHtocm1))
	end
	
	
	# f=open("log","a")
	# println(f,"Write eigenvalues")
	# println(f)
	# close(f)
	
	# f=open("energies.txt","w")
	# for istates=1:Nsize_para
	# 	println(fp,round(real(e_para[istates])*eHtocm1,digits=18)," ",round(real(e_para[istates]-e_para[1])*eHtocm1,digits=16))
	# end
	# for istates=1:Nsize_ortho
	# 	println(fo,round(real(e_ortho[istates])*eHtocm1,digits=18)," ",round(real(e_ortho[istates]-e_para[1])*eHtocm1,digits=16))
	# end

	# close(fp)
	# close(fo)
	
	# f=open("log","a")
	# println(f,"Calculate heat capacity")
	# f3=open("heat_capacity.txt","a")
	# for it=1:Ntemp
	# 	cv_p,cv_o,cv_full=heat_capacity(real(e_para),real(e_ortho),Tstart+(it-1)*dT)
	# 	println(f3,round(Tstart+(it-1)*dT,digits=3),"  ",cv_p*eHtocm1,"  ",cv_o*eHtocm1,"  ",cv_full*eHtocm1)
	# end
	# close(f3)
	# println(f,"Heat capacity done")
	# println(f)
	# close(f)

end