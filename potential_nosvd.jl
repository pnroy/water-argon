module potential

using LinearAlgebra
using SpecialFunctions
using Jacobi
using AssociatedLegendrePolynomials
using ClassicalOrthogonalPolynomials: laguerrel
using GaussQuadrature: laguerre,legendre 
using FFTW
push!(LOAD_PATH,pwd())
using kinetic
using h2o_c60

export potential_matrix

###############################################################################################################
function potential_matrix(nmax,NR,Ntheta,Nphi,jmax,Nalpha,Nm,Nk,omega,mass,Ntrans,Nrot,eq_struct)

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
Vpot = zeros(Ngrid1,Ngrid2)

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
		Vpot[i1,i2] = cage_potential(mass,omega,sqrt(0.5*R[ir]/nu),acos(ctheta1[it1]),phi1[ip1],acos(ctheta2[it2]),phi2[ip2],chi[ic],eq_struct) 
		#Vpot[i1,i2]=0.5*mass*(omega^2)*(0.5*R[ir]/nu)
	end
	end
	end
end
end
end

coeff1 = zeros(ComplexF64,Ntrans)
coeff2 = zeros(ComplexF64,Nrot)
coeff_grid = zeros(ComplexF64,(Ngrid1,Ngrid2))
coeff_tmp = zeros(ComplexF64,(Ngrid1,Nrot))
coeff_end = zeros(ComplexF64,(Ntrans,Nrot))
matrix = zeros(ComplexF64,(Ntrans,Ntrans,Nrot,Nrot))


for ispec1=1:Ntrans
	coeff1 .= 0.0
	coeff1[ispec1] = 1.0
	#Transform from NLM-basis to (R,theta,phi)-grid#
	coeff1_out = transformation_forward_HO(coeff1,nmax,NR,Ntheta,Nphi,nu,mass)
	for ispec2=1:Nrot
		coeff2 .= 0.0	
		coeff2[ispec2] = 1.0
		#Transform from JKM-basis to (Theta,phi,chi)-grid#
		coeff2_out = transformation_forward_Wigner(coeff2,jmax,mmax,kmax,Nalpha,Nm,Nk)
		#Multipliy transformed vectors with potential on grid#
		for ig1=1:Ngrid1
			for ig2=1:Ngrid2
				coeff_grid[ig1,ig2] = Vpot[ig1,ig2]*coeff1_out[ig1]*coeff2_out[ig2]
			end
			#Transform from (Theta,phi,chi)-grid to JKM-basis#
			coeff_tmp[ig1,:] .= transformation_backward_Wigner(coeff_grid[ig1,:],jmax,mmax,kmax,Nalpha,Nm,Nk,Nrot)
		end	
		#Transform from (R,theta,phi)-grid to NLM-basis#
		for j=1:Nrot
			coeff_end[:,j] = transformation_backward_HO(coeff_tmp[:,j],nmax,NR,Ntheta,Nphi,nu,mass,Ntrans)
		end
	
		for i=1:Ntrans
		for j=1:Nrot
			matrix[i,ispec1,j,ispec2] = coeff_end[i,j]
		end
		end
	end
end

return matrix
end
####################################################################################################################
function transformation_forward_HO(coeff_in,nmax,NR,Ntheta,Nphi,nu,mass)
####################################################################################################################
##Function to perform a transformation from the NLM (isotropic 3D-HO eigenfunctions) basis to the R,theta,phi-grid##
####################################################################################################################

coeff_start = zeros(ComplexF64,(nmax+1,nmax+1,2*nmax+1))

ii=0
for n=0:nmax
for l=n:-2:0
for m=-l:l
	ii += 1
	coeff_start[n+1,l+1,nmax+1+m] = coeff_in[ii]
end
end
end


#Transform from the N basis to the R (Gauss-Laguerre) grid#
coeff_out = RN_transformation("forward",coeff_start,nmax+1,NR,nmax,nu,mass)

#Transform from the L basis to the theta (Gauss-Legendre) grid#
coeff_out2 = ThetaL_transformation("forward",coeff_out,nmax+1,Ntheta,nmax,NR)

#Transform from the M basis to the phi grid#
coeff=zeros(ComplexF64,Nphi)
coeff_grid=zeros(ComplexF64,(NR*Ntheta*Nphi))
tmp=zeros(ComplexF64,Nphi)

Pm=plan_ifft(tmp)
ig=1
for ir=1:NR
for it=1:Ntheta
        for m=1:2*nmax+1
		coeff[m]=coeff_out2[ir,it,m]
        end
	tmp=Pm*coeff
	for ip=1:Nphi
		coeff_grid[ig]=tmp[ip]
		ig+=1
	end
end
end

return coeff_grid 
end
###############################################################################################################
function transformation_backward_HO(coeff_grid,nmax,NR,Ntheta,Nphi,nu,mass,Nspec)
####################################################################################################################
##Function to perform a transformation from the R,theta,phi-grid to the NLM (isotropic 3D-HO eigenfunctions) basis##
####################################################################################################################

##Conversion from R,theta,phi indices to composite index##
coeff_start=zeros(ComplexF64,(NR,Ntheta,Nphi))

ii=0
for ir=1:NR
for it=1:Ntheta
for ip=1:Nphi
    ii+=1
    coeff_start[ir,it,ip]=coeff_grid[ii]
end
end
end

##Transform from phi grid to M basis##
coeff_tmp=zeros(ComplexF64,Nphi)
coeff_out=zeros(ComplexF64,(Ntheta,NR,2*nmax+1))
tmp=zeros(ComplexF64,Nphi)

Pm=plan_fft(coeff_tmp)
for ir=1:NR
for it=1:Ntheta
	for ip=1:Nphi
		tmp[ip]=coeff_start[ir,it,ip]
	end
        coeff_tmp=Pm*tmp
	for m=1:2*nmax+1
		coeff_out[it,ir,m]=coeff_tmp[m]
	end
end
end

#Transform from the theta (Gauss-Legendre) grid to the L basis#

coeff_out2 = ThetaL_transformation("backward",coeff_out,Ntheta,nmax+1,nmax,NR)

#Transform from the R (Gauss-Laguerre) grid to the N basis#
coeff_out3 = RN_transformation("backward",coeff_out2,NR,nmax+1,nmax,nu,mass)

coeff_out4 = zeros(ComplexF64,Nspec)
ii=0
for n=0:nmax
for l=n:-2:0
for m=-l:l
	ii += 1
	coeff_out4[ii] = coeff_out3[l+1,n+1,nmax+1+m]
end
end
end

return coeff_out4
end
###############################################################################################################
function RN_transformation(direction,coeff,N1,N2,nmax,nu,mass)

if direction == "forward"
        trans='N'
	Utrans = Urn(nu,nmax,N2)
elseif direction == "backward"
        trans='T'
	Utrans = Urn(nu,nmax,N1)
end

#Transform between R- and n-grid#
coeff2_R = zeros((N2,nmax+1,2*nmax+1))
coeff2_I = zeros((N2,nmax+1,2*nmax+1))
for l=0:nmax
for m=-l:l
	coeff2_R[:,l+1,nmax+1+m].=BLAS.gemv(trans, 1.0, Utrans[:,:,l+1], real(coeff[:,l+1,nmax+1+m]))
	coeff2_I[:,l+1,nmax+1+m].=BLAS.gemv(trans, 1.0, Utrans[:,:,l+1], imag(coeff[:,l+1,nmax+1+m]))
end
end

coeff_out = zeros(ComplexF64,(nmax+1,N2,2*nmax+1))
for n=1:N2
for l=0:nmax
for m=-l:l
	coeff_out[l+1,n,nmax+1+m] = coeff2_R[n,l+1,nmax+1+m] + im*coeff2_I[n,l+1,nmax+1+m]
end
end
end

return coeff_out
end
###############################################################################################################
function ThetaL_transformation(direction,coeff,N1,N2,nmax,NR)

if direction == "forward"
        trans='N'
	Utrans = Utl(nmax,N2)
elseif direction == "backward"
        trans='T'
	Utrans = Utl(nmax,N1)
end

#Transform between L- and theta-grid#
coeff2_R = zeros((N2,NR,2*nmax+1))
coeff2_I = zeros((N2,NR,2*nmax+1))

for ir=1:NR
for m=1:2*nmax+1
	coeff2_R[:,ir,m].=BLAS.gemv(trans, 1.0, Utrans[:,:,m], real(coeff[:,ir,m]))
	coeff2_I[:,ir,m].=BLAS.gemv(trans, 1.0, Utrans[:,:,m], imag(coeff[:,ir,m]))
end
end

coeff_out = zeros(ComplexF64,(NR,N2,2*nmax+1))
for ir=1:NR
for it=1:N2
for m=1:2*nmax+1
	coeff_out[ir,it,m] = coeff2_R[it,ir,m] + im*coeff2_I[it,ir,m]
end
end
end

return coeff_out
end
###############################################################################################################
function Urn(nu,nmax,NR)

Rgrid,weight = laguerre(NR,0.5)
Utrans = zeros((NR,nmax+1,nmax+1))

for n=0:nmax
	for l=n:-2:0
	k=Int((n-l)/2)
	prefactor = laguerre_normalization(k,l,nu)/(2.0*sqrt(nu*(2.0*nu)^(l+0.5)))
	for ir=1:NR
		lag = laguerrel(k,l+0.5,Rgrid[ir])
		Utrans[ir,n+1,l+1] = prefactor*sqrt(weight[ir]*Rgrid[ir]^l)*lag
	end
	end
end


return Utrans
end
###############################################################################################################
function Utl(nmax,Ntheta)

Utrans = zeros((Ntheta,nmax+1,2*nmax+1))

ctheta, ct_weight = legendre(Ntheta)

for l=0:nmax
	for m=-l:-1
		norm = legendre_normalization(l,m)
		for it=1:Ntheta
			Utrans[it,l+1,nmax+1+m] = (-1.0)^(abs(m))*(factorial(l-abs(m))/factorial(l+abs(m)))*Plm(l,abs(m),ctheta[it])*sqrt(ct_weight[it])*norm
		end
	end
	for m=0:l
		norm = legendre_normalization(l,m)
		for it=1:Ntheta
			Utrans[it,l+1,nmax+1+m] = Plm(l,m,ctheta[it])*sqrt(ct_weight[it])*norm
		end
	end


end

return Utrans
end
###############################################################################################################
function laguerre_normalization(k,l,nu)

norm = sqrt(sqrt((2*nu^3)/pi)*(2^(k+2*l+3)*factorial(k)*nu^(l)/dfactorial(2*k+2*l+1)))

return norm
end
###############################################################################################################
function legendre_normalization(l,m)

norm = sqrt(0.5*(2*l+1)*factorial(l-m)/factorial(l+m))
norm = norm*(-1)^m	#Condon-Shortley phase factor

return norm
end
###############################################################################################################
function dfactorial(n)

if mod(n,2) == 0
	Nfac = Int(n/2)
elseif mod(n,2) == 1
	Nfac = Int((n+1)/2)
end

fac = 1.0
for ii=0:Nfac-1
	fac*=(n-2*ii)
end

return fac
end
###############################################################################################################
###############################################################################################################
function transformation_forward_Wigner(coeff_in,jmax,mmax,kmax,Nalpha,Nm,Nk)

coeff_start=zeros(ComplexF64,(jmax+1,2*mmax+1,2*kmax+1))
i=0
for j=0:jmax
	mj=min(j,mmax)
	kj=min(j,kmax)
	for m=-mj:mj
	for k=-kj:kj
		i+=1
		coeff_start[j+1,mmax+1+m,kmax+1+k]=coeff_in[i]
	end
	end
end


##Transform from J basis to theta grid##
coeff_out=ThetaJ_transformation("forward" ,jmax+1,Nalpha,jmax,mmax,kmax,coeff_start)

##Transform from M basis to phi grid##
coeff=zeros(ComplexF64,Nm)
coeff_out2=zeros(ComplexF64,(Nalpha,2*kmax+1,Nm))
tmp=zeros(ComplexF64,Nm)

Pm=plan_ifft(tmp)
for i2=1:2*kmax+1
for alpha=1:Nalpha
        for m=1:2*mmax+1
		coeff[m]=coeff_out[alpha,m,i2]
        end
	tmp=Pm*coeff
	for m=1:Nm
		coeff_out2[alpha,i2,m]=tmp[m]
	end
end
end

##Transform from K basis to chi grid##
coeff2=zeros(ComplexF64,Nk)
coeff_out3=zeros(ComplexF64,(Nk,Nm,Nalpha))
tmp2=zeros(ComplexF64,Nk)

Pk=plan_ifft(tmp2)
for i2=1:Nm
for alpha=1:Nalpha
	for k=1:2*kmax+1
		coeff2[k]=coeff_out2[alpha,k,i2]
	end
	tmp2=Pk*coeff2
	for k=1:Nk
		coeff_out3[k,i2,alpha]=tmp2[k]
	end
end
end

##Conversion from theta,phi,chi indices to composite index##
coeff_grid=zeros(ComplexF64,Nalpha*Nm*Nk)

i=0
for it=1:Nalpha
for ip=1:Nm
for ic=1:Nk
    i+=1
    coeff_grid[i]=coeff_out3[ic,ip,it]
end
end
end

return coeff_grid

end
###############################################################################################################
function transformation_backward_Wigner(coeff_grid,jmax,mmax,kmax,Nalpha,Nm,Nk,Nspec)

##Conversion from theta,phi,chi indices to composite index##
coeff_start=zeros(ComplexF64,(Nalpha,Nm,Nk))

i=0
for alpha=1:Nalpha
for ip=1:Nm
for ic=1:Nk
    i+=1
    coeff_start[alpha,ip,ic]=coeff_grid[i]
end
end
end


##Transform from chi grid to K basis##
coeff_tmp=zeros(ComplexF64,Nk)
coeff_out=zeros(ComplexF64,(Nalpha,2*kmax+1,Nm))
tmp=zeros(ComplexF64,Nk)

Pk=plan_fft(coeff_tmp)
for i2=1:Nm
for alpha=1:Nalpha
	for k=1:Nk
		tmp[k]=coeff_start[alpha,i2,k]
	end
        coeff_tmp=Pk*tmp
	for k=1:2*kmax+1
		coeff_out[alpha,k,i2]=coeff_tmp[k]
	end
end
end

##Transform from phi grid to M basis##
coeff_tmp2=zeros(ComplexF64,Nm)
coeff_out2=zeros(ComplexF64,(Nalpha,2*mmax+1,2*kmax+1))
tmp2=zeros(ComplexF64,Nm)

Pm=plan_fft(coeff_tmp2)
for i2=1:2*kmax+1
for alpha=1:Nalpha
	for m=1:Nm
		tmp2[m]=coeff_out[alpha,i2,m]
	end
	coeff_tmp2=Pm*tmp2
	for m=1:2*mmax+1
        	coeff_out2[alpha,m,i2]=coeff_tmp2[m]
	end
end
end

#Transform from theta grid to J basis#

coeff_out3= ThetaJ_transformation("backward",Nalpha,jmax+1,jmax,mmax,kmax,coeff_out2)

coeff_out4=zeros(ComplexF64,Nspec)
i=0
for j=0:jmax
	mj=min(j,mmax)
	kj=min(j,kmax)
	for m=-mj:mj
	for k=-kj:kj
		i+=1
		coeff_out4[i]=coeff_out3[j+1,mmax+1+m,kmax+1+k]
	end
	end
end

return coeff_out4
end
###############################################################################################################
function ThetaJ_transformation(direction,N1,N2,jmax,mmax,kmax,coeff_in)

coeff_R=zeros((N1,2*mmax+1,2*kmax+1))
coeff_I=zeros((N1,2*mmax+1,2*kmax+1))

for k=1:2*kmax+1
for m=1:2*mmax+1
for i=1:N1
	coeff_R[i,m,k]=real(coeff_in[i,m,k])
	coeff_I[i,m,k]=imag(coeff_in[i,m,k])
end
end
end

if direction == "forward"
        trans='N'
        matrix=Utj(N2,jmax,mmax,kmax)
elseif direction == "backward"
        trans='T'
        matrix=Utj(N1,jmax,mmax,kmax)
end

coeff2_R=zeros((N2,2*mmax+1,2*kmax+1))
coeff2_I=zeros((N2,2*mmax+1,2*kmax+1))


#Act with collocation matrix on vector#
for k=1:2*kmax+1
for m=1:2*mmax+1
	coeff2_R[:,m,k].=BLAS.gemv(trans, 1.0, matrix[:,:,m,k], coeff_R[:,m,k])
	coeff2_I[:,m,k].=BLAS.gemv(trans, 1.0, matrix[:,:,m,k], coeff_I[:,m,k])
end
end
coeff_out=zeros(ComplexF64,(N2,2*mmax+1,2*kmax+1))

for k=1:2*kmax+1
for m=1:2*mmax+1
for i=1:N2
	coeff_out[i,m,k]=coeff2_R[i,m,k]+im*coeff2_I[i,m,k]
end
end
end

return coeff_out
end 
###############################################################################################################
function Utj(Nalpha,jmax,mmax,kmax)

#Gauss-Legendre points#
ctheta, ct_weight = legendre(Nalpha)

matrix=zeros((Nalpha,jmax+1,2*mmax+1,2*kmax+1))
for m=-mmax:mmax
for k=-kmax:kmax
for j=0:jmax
	if j < abs(m) || j < abs(k) 
        	matrix[:,j+1,m+mmax+1,k+kmax+1].=0.0
	else
        	wig=wigner_dmatrix(j,m,k,ctheta,Nalpha)
        	for alpha=1:Nalpha
            		matrix[alpha,j+1,m+mmax+1,k+kmax+1]=sqrt((2.0*j+1.0)/2.0)*wig[alpha]*sqrt(ct_weight[alpha])
        	end
    	end
end
end
end

return matrix
end
###############################################################################################################
function wigner_dmatrix(j,m,k,u,Nu)

mone=-1.0
two=2.0
one=1.0

alpha=(m-k)
beta=m+k
n=j-m

factor=factorial(big(n+alpha))*factorial(big(n+beta))/(factorial(big(n))*factorial(big(n+alpha+beta)))

poly=zeros(Nu)
P=zeros(Nu)
if alpha < 0.0 && beta < 0.0
	for i=1:Nu
		P[i]=Jacobi.jacobi(u[i], n+alpha+beta, -alpha, -beta)
        end
	poly[:].=P[:].*((0.5*(u[:].-1.0)).^(-alpha)).*((0.5*(u[:].+1.0)).^(-beta))
end
if alpha < 0.0 && beta >= 0.0
 	for i=1:Nu
		P[i]=Jacobi.jacobi(u[i], n+alpha, -alpha, beta)
        end
	poly[:].=P[:]*factor.*((0.5*(u[:].-1.0)).^(-alpha))
end
if alpha >= 0.0 && beta < 0.0 
	for i=1:Nu        
		P[i]=Jacobi.jacobi(u[i], n+beta, alpha, -beta)
        end
	poly[:].=P[:]*factor.*((0.5*(u[:].+1.0)).^(-beta))
end
if alpha >= 0.0 && beta >= 0.0
        for i=1:Nu
		P[i]=Jacobi.jacobi(u[i], n, alpha, beta)
        end
	poly[:].=P[:]
end

djmk=zeros(Nu)
for i=1:Nu
        djmk[i]=((mone^alpha)/(two^m))*sqrt(factorial(big(j+m))*factorial(big(j-m))/(factorial(big(j+k))*factorial(big(j-k))))
        djmk[i]=djmk[i]*((1.0-u[i])^(alpha*0.5))*((1.0+u[i])^(beta*0.5))*poly[i]
end

return djmk
end
###############################################################################################################

end
