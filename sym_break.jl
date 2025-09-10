module sym_break

using LinearAlgebra
using WignerSymbols
using SpecialFunctions
using GaussQuadrature: laguerre,legendre 
push!(LOAD_PATH,pwd())
using potential

export Vquad
#############################################################################################
function Vquad(nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,A,QBF,mu)

mmax=jmax
kmax=jmax
nu=mass*omega*0.5

################
##Define grids##
################

Ngrid_trans = NR*Ntheta*Nphi
Ngrid_rot = Nalpha*Nm*Nk

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

##############################
##Calculate rotational terms##
##############################
D2mq = zeros(Float64,(Nrot,Nrot,5,5))
D1m0 = zeros(ComplexF64,(Nrot,Nrot,3))

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
		for m=-2:2
		for q=-2:2
			#if m1-m-m2 == 0 && k1-q-k2 == 0
				dum=wigner3j(BigFloat,j1,2,j2,m1,-m,-m2)*wigner3j(BigFloat,j1,2,j2,k1,-q,-k2)
				D2mq[i1,i2,2+1+m,2+1+q]=dum*sqrt((2*j1+1)*(2*j2+1))*(-1.0)^(m-q+m2-k2)
			#end
		end
		end
		for m=-1:1
			#if m1-m-m2 == 0 && k1-k2 == 0
				dum=wigner3j(BigFloat,j1,1,j2,m1,-m,-m2)*wigner3j(BigFloat,j1,1,j2,k1,0,-k2)
				D1m0[i1,i2,1+1+m]=dum*sqrt((2*j1+1)*(2*j2+1))*(-1.0)^(m+m2-k2)
			#end
		end
	end
	end
	end
end
end
end

#################################
##Calculate translational terms##
#################################
U = Urn(nu,nmax,NR)
lag = zeros(Float64,(nmax+1,nmax+1,nmax+1,nmax+1))


for n1=0:nmax
for l1=n1:-2:0
for n2=0:nmax
for l2=n2:-2:0
	dum = 0.0
	for ir=1:NR
		dum+=U[ir,n1+1,l1+1]*sqrt(0.5*R[ir]/nu)*U[ir,n2+1,l2+1]
	end
		lag[n1+1,l1+1,n2+1,l2+1]=dum
end
end
end
end

Y1m=zeros(ComplexF64,(Ntrans,Ntrans,3))
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
		for q=-1:1
			Y=sqrt((2*l1+1)*3*(2*l2+1)/(4*pi))*wigner3j(BigFloat,l1,1,l2,0,0,0)*wigner3j(BigFloat,l1,1,l2,-m1,q,m2)*(-1)^m1
			Y1m[i1,i2,1+1+q]=lag[n1+1,l1+1,n2+1,l2+1]*Y
		end
	end
	end
	end
end
end
end

##############
##Define Q2m##
##############
I2m=A*[im,im-1.0,0,im+1,-im]
delta = zeros(Float64,(Ntrans,Ntrans))
for it=1:Ntrans
	delta[it,it] = 1.0
end

wig = zeros(ComplexF64,(3,3,5))

for m=-2:2
for m2=-1:1
for m1=-1:1	
	wig[2+m1,2+m2,3+m]=sqrt(40*pi)*mu*wigner3j(BigFloat,1,1,2,m1,m2,-m)*I2m[3-m]
end
end
end

term1 = zeros(ComplexF64,(Nrot,Nrot))

for ir2=1:Nrot
for ir1=1:Nrot 
	dum1=0.0
	for q=1:5
	for m=-2:2
        	dum1+=D2mq[ir1,ir2,3+m,q]*QBF[q]*I2m[3-m]*(-1)^m
	end
	end
	term1[ir1,ir2] = dum1
end
end

return term1,Y1m,wig,D1m0
end
#############################################################################################
end
