module h2o_Argon

using LinearAlgebra

export cage_potential

###########################################################################################
function cage_potential(R,theta1,phi1,theta2,phi2,chi,eq_struct,iv)

#Conversion factor#
Ehtokcalmol=627.509474	#Hartree to kcal/mol
a0toangst=0.529177210903
eHtocm1=219474.631363	#Eh to cm-1

#Convert R from bohr to A#
R=R*a0toangst

X=R*sin(theta1)*cos(phi1)
Y=R*sin(theta1)*sin(phi1)
Z=R*cos(theta1)

# COM = zeros(3)
# COM[1] = R*sin(theta1)*cos(phi1)
# COM[2] = R*sin(theta1)*sin(phi1)
# COM[3] = R*cos(theta1)

#Read coordinates of cage atoms#
f = open("Ar.xyz") 
lines = readlines(f)
close(f)
lsplit = split(lines[1])
Ncage = parse(Int64, lsplit[1])

Aratoms = zeros((Ncage,3))
for i=1:Ncage
	lsplit = split(lines[i+2])
 	Aratoms[i,1] = parse(Float64,lsplit[2])
 	Aratoms[i,2] = parse(Float64,lsplit[3])
 	Aratoms[i,3] = parse(Float64,lsplit[4])
end

V=0.0 

rAr_COM=zeros(3)

for i=1:Ncage
	rAr_COM[1]=Aratoms[i,1]-X
	rAr_COM[2]=Aratoms[i,2]-Y
	rAr_COM[3]=Aratoms[i,3]-Z
	dist=norm(rAr_COM)
	# get theta and phi from rotation matrix
	#Rotate molecule by Euler angles around COM, matrix is for BFF to SFF, use transpose
	Rot=rotation_matrix(theta2,phi2,chi)
	rAr_BFF=BLAS.gemv('T' , 1.0, Rot, rAr_COM)
	theta_BFF=acos(rAr_BFF[3]/dist)
	# atan2 below  
	phi_BFF=atan(rAr_BFF[2],rAr_BFF[1])

#  iv=1 for (v1,v2,v3)=(0,0,0); 
#  iv=2 for (v1,v2,v3)=(1,0,0); 
#  iv=3 for (v1,v2,v3)=(0,1,0); 
#  iv=4 for (v1,v2,v3)=(0,0,1). 

	result_ref = Ref{Float64}()
	ccall((:h2oar_3dpes_, "./potlibarh20.dylib"), Nothing, (Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},Ref{Int}), phi_BFF,theta_BFF,dist,result_ref,iv)
	value=result_ref[]

	V+=value

end
V = V/eHtocm1

return V
end
###########################################################################################


###########################################################################################
function rotation_matrix(theta,phi,chi)

#body fixed to space fixed
cp=cos(phi)
sp=sin(phi)
ct=cos(theta)
st=sin(theta)
ck=cos(chi)
sk=sin(chi)

rotmat=zeros(Float64,(3,3))

rotmat[1,1]=ct*cp*ck-sp*sk
rotmat[1,2]=-ct*cp*sk-sp*ck
rotmat[1,3]=st*cp
rotmat[2,1]=ct*sp*ck+cp*sk
rotmat[2,2]=-ct*sp*sk+cp*ck
rotmat[2,3]=st*sp
rotmat[3,1]=-st*ck
rotmat[3,2]=st*sk
rotmat[3,3]=ct


return rotmat
end
###########################################################################################
end
