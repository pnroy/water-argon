module kinetic

using SpecialFunctions

export kinetic_rotation,kinetic_translation,kinetic_trans_analytic,harmonic_integral,H3D_translation

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

Trot=zeros(ComplexF64,(N,N))
Jx=zeros(ComplexF64,(N,N))
Jz=zeros(ComplexF64,(N,N))

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
				if k1 == k2-1
                    Jx[i1,i2]=0.5*sqrt(1.0*j*(j+1)-1.0*k2*(k2-1))
                end
                if k1 == k2+1
                    Jx[i1,i2]=0.5*sqrt(1.0*j*(j+1)-1.0*k2*(k2+1))
                end
        	end
            Trot[i1,i1]=j*(j+1)*D1+D2*k1^2
            Jz[i1,i1]=k1
        end
	end
end

return Trot,Jx,Jz
end

####################################################################################################################
function kinetic_trans_analytic(nmax,omega)
####################################################################################################################
##Analytic kinetic matrix elements for translational motion in isotropic 3D-HO basis, from Felker and Bacic, JCP,### 
#########################################152, 2020 (Supplementary material)#########################################
####################################################################################################################

Nspec=0
for n=0:nmax
	Nspec += Int(((n+1)*(n+2))/2)
end

Ttrans = zeros(ComplexF64,(Nspec,Nspec))

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
		if n1 == n2+2 && l1 == l2 && m1 == m2
			Ttrans[i1,i2] = 0.25*omega*sqrt((n2+2)*(n2+3)-l1*(l1+1))
		end
		if n1+2 == n2 && l1 == l2 && m1 == m2
			Ttrans[i1,i2] = 0.25*omega*sqrt((n1+2)*(n1+3)-l1*(l1+1))
		end
	end
	end
	end
end
end
end
	
i=0
for n=0:nmax
for l=n:-2:0
for m=-l:l
	i+=1
	Ttrans[i,i] = 0.5*omega*(n+1.5)
end
end
end

return Ttrans
end
####################################################################################################################
function kinetic_translation(nmax,omega,mass)
####################################################################################################################
##Calculate kinetic matrix elements for translational motion in isotropic 3D-HO basis, from Kalugina and Roy, JCP,##
#####################################################147, 2017######################################################
####################################################################################################################

Nspec=0
for n=0:nmax
	Nspec += Int(((n+1)*(n+2))/2)
end

#Calculate Hamiltonian of isotropic 3D HO#
H3D = zeros((Nspec,Nspec))
i1=0
for n1=0:nmax
for l1=n1:-2:0
for m1=-l1:l1
	i1+=1
	i2=0
	k1 = Int((n1-l1)/2)
	for n2=0:nmax
	for l2=n2:-2:0
	for m2=-l2:l2	
		i2+=1
		k2 = Int((n2-l2)/2)
		if k1 == k2 && l1 == l2 && m1 == m2
			H3D[i1,i2] = omega*(n1+1.5) 
		end		
	end
	end
	end
end
end
end

#Calculate integral over harmonic potential in isotropic 3D-HO basis#
I3 = harmonic_integral(nmax,Nspec,omega,mass)
	
Ttrans = zeros(ComplexF64,(Nspec,Nspec)) 

Ttrans = H3D - I3

return Ttrans
end

function H3D_translation(nmax,omega,mass)
####################################################################################################################
##Calculate kinetic matrix elements for translational motion in isotropic 3D-HO basis, from Kalugina and Roy, JCP,##
#####################################################147, 2017######################################################
####################################################################################################################

Nspec=0
for n=0:nmax
	Nspec += Int(((n+1)*(n+2))/2)
end

#Calculate Hamiltonian of isotropic 3D HO#
H3D = zeros(ComplexF64,(Nspec,Nspec))
i1=0
for n1=0:nmax
for l1=n1:-2:0
for m1=-l1:l1
	i1+=1
	i2=0
	k1 = Int((n1-l1)/2)
	for n2=0:nmax
	for l2=n2:-2:0
	for m2=-l2:l2	
		i2+=1
		k2 = Int((n2-l2)/2)
		if k1 == k2 && l1 == l2 && m1 == m2
			H3D[i1,i2] = omega*(n1+1.5) 
		end		
	end
	end
	end
end
end
end


return H3D
end

####################################################################################################################
function harmonic_integral(nmax,Nspec,omega,mass)

lmax = nmax
if mod(nmax,2) == 0
	kmax = Int(nmax/2)
elseif mod(nmax,2) == 1
	kmax = Int((nmax-1)/2)
end

nu=mass*omega*0.5

integral=zeros(kmax+1,kmax+1,lmax+1)

for k1=0:kmax
for k2=0:kmax
	for l=0:lmax
		tau = k1-k2
		m = 0.5*(2*l+3)
		alpha = l+0.5
		dum = 0.0
		if tau >= 0
			for sigma=0:k2
				dum+=binomial(Int(m-alpha),tau+sigma)*binomial(Int(m-alpha),sigma)*gamma(m+k2-sigma+1)/factorial(k2-sigma)
			end	
			dum=dum*mass*(omega^2)*((-1)^(k1+k2))/(8*nu*(2*nu)^m)
		else
			for sigma=0:k1
				dum+=binomial(Int(m-alpha),abs(tau)+sigma)*binomial(Int(m-alpha),sigma)*gamma(m+k1-sigma+1)/factorial(k1-sigma)
			end	
			dum=dum*mass*(omega^2)*((-1)^(k1+k2))/(8*nu*(2*nu)^m)
		end
		integral[k1+1,k2+1,l+1]=dum*laguerre_normalization(k1,l,nu)*laguerre_normalization(k2,l,nu)
	end
end
end

int_out = zeros(ComplexF64,(Nspec,Nspec))
i1=0
for n1=0:nmax
for l1=n1:-2:0
for m1=-l1:l1
	i1+=1
	i2=0
	k1 = Int((n1-l1)/2)
	for n2=0:nmax
	for l2=n2:-2:0
	for m2=-l2:l2	
		i2+=1
		k2 = Int((n2-l2)/2)
		if l1 == l2 && m1 == m2
			int_out[i1,i2] = integral[k1+1,k2+1,l1+1]
		end	
	end	
	end	
	end	
end	
end	
end	

return int_out
end
####################################################################################################################
function laguerre_normalization(k,l,nu)

norm = sqrt(sqrt((2*nu^3)/pi)*(2^(k+2*l+3)*factorial(k)*nu^(l)/dfactorial(2*k+2*l+1)))

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
end
