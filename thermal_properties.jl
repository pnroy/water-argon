module thermal_properties

export heat_capacity,polarizability

#####################################################################################################################
function heat_capacity(e_para,e_ortho,T)

#Constants, conversion factors#
kB=1.380649e-23
EhtoJ = 4.3597447222071e-18

Nstates=length(e_para)+length(e_ortho)
Npara = length(e_para)
Northo = length(e_ortho)
beta=EhtoJ/(kB*T)

energy = zeros(Float64,(Nstates))
k=1
for ii=1:Npara
	energy[ii] = e_para[ii]
	k+=1
end
for ii=1:Northo
	energy[k] = e_ortho[ii]
end
sort!(energy)

#Cv in full Hilbert space#
E=0.0
E2=0.0
Z=0.0
for istates=1:Nstates
	E+=energy[istates]*exp(-(energy[istates]-energy[1])*beta)
	E2+=(energy[istates]^2)*exp(-(energy[istates]-energy[1])*beta)
	Z+=exp(-(energy[istates]-energy[1])*beta)
end
E=E/Z
E2=E2/Z

Cv_full=(beta/T)*(E2-E^2)

#Cv in para space#
E=0.0
E2=0.0
Z=0.0
for istates=1:Npara
	E+=e_para[istates]*exp(-(e_para[istates]-e_para[1])*beta)
	E2+=(e_para[istates]^2)*exp(-(e_para[istates]-e_para[1])*beta)
	Z+=exp(-(e_para[istates]-e_para[1])*beta)
end
E=E/Z
E2=E2/Z

Cv_para=(beta/T)*(E2-E^2)

#Cv in ortho space#
E=0.0
E2=0.0
Z=0.0
for istates=1:Northo
	E+=e_ortho[istates]*exp(-(e_ortho[istates]-e_ortho[1])*beta)
	E2+=(e_ortho[istates]^2)*exp(-(e_ortho[istates]-e_ortho[1])*beta)
	Z+=exp(-(e_ortho[istates]-e_ortho[1])*beta)
end
E=E/Z
E2=E2/Z

Cv_ortho=(beta/T)*(E2-E^2)

return Cv_para,Cv_ortho,Cv_full
end
#####################################################################################################################
function polarizability(energy,muX,muY,muZ,T,omega_ext)

Nstates = length(energy)

mu2 = zeros(Nstates,Nstates)
for i1=1:Nstates
for i2=1:Nstates
	mu2[i1,i2] = muX[i1,i2]^2+muY[i1,i2]^2+muZ[i1,i2]^2
end
end

#Constants, conversion factors#
kB=1.380649e-23
EhtoJ = 4.3597447222071e-18
eHtocm1=219474.631363	#Eh to cm-1

#Energy threshold to distinguish degenerate states#
Ethresh = 0.1	#cm-1
Ethresh = Ethresh/eHtocm1

beta=EhtoJ/(kB*T)

Z=0.0
alpha=0.0
for i1=1:Nstates-1
	for i2=(i1+1):Nstates
		if abs(energy[i1]-energy[i2]) > Ethresh
			alpha+=exp(-(energy[i1]-energy[1])*beta)*(energy[i2]-energy[i1])*mu2[i1,i2]/((energy[i2]-energy[i1])^2-omega_ext^2)	
		end
	end
	Z+=exp(-(energy[i1]-energy[1])*beta)
end

alpha = 2*alpha/(3.0*Z)

return alpha
end
#####################################################################################################################
end
