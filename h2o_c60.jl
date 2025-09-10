module h2o_c60

using LinearAlgebra

export cage_potential

###########################################################################################
function cage_potential(R,theta1,phi1,theta2,phi2,chi,eq_struct)

#Conversion factor#
Ehtokcalmol=627.509474	#Hartree to kcal/mol
a0toangst=0.529177210903
#Lennard-Jones parameter#
eps_H = 0.0256	#kcal/mol
sig_H = 2.64 #A
eps_O = 0.1039
sig_O = 3.372

#Convert R from bohr to A#
R=R*a0toangst

#Read coordinates of cage atoms#
f = open("c60.txt") 
lines = readlines(f)
close(f)
lsplit = split(lines[1])
Ncage = parse(Int64, lsplit[1])

Catoms = zeros((Ncage,3))
for i=1:Ncage
	lsplit = split(lines[i+2])
 	Catoms[i,1] = parse(Float64,lsplit[2])
 	Catoms[i,2] = parse(Float64,lsplit[3])
 	Catoms[i,3] = parse(Float64,lsplit[4])
end
Catoms *= a0toangst

rO,rH1,rH2 = get_cart_coord(R,theta1,phi1,theta2,phi2,chi,eq_struct)

V=0.0
for i=1:Ncage
	dist = sqrt((rO[1]-Catoms[i,1])^2+(rO[2]-Catoms[i,2])^2+(rO[3]-Catoms[i,3])^2)
	V += lennard_jones(dist,eps_O,sig_O)
	dist = sqrt((rH1[1]-Catoms[i,1])^2+(rH1[2]-Catoms[i,2])^2+(rH1[3]-Catoms[i,3])^2)
	V += lennard_jones(dist,eps_H,sig_H)
	dist = sqrt((rH2[1]-Catoms[i,1])^2+(rH2[2]-Catoms[i,2])^2+(rH2[3]-Catoms[i,3])^2)
	V += lennard_jones(dist,eps_H,sig_H)
end
V = V/Ehtokcalmol

return V
end
###########################################################################################
function get_cart_coord(R,theta1,phi1,theta2,phi2,chi,eq_struct)

#Convert R distance to angstrom#
a0toangst=0.529177210903
rOH=eq_struct[1]*a0toangst
hh_angle=eq_struct[2]
hh_angle=hh_angle*0.5

mH = 1.00782503207
mO = 15.99491461956

O_init=zeros(3)
H1_init=zeros(3)
H2_init=zeros(3)

O_init[3] = 0.125534*a0toangst

H1_init[1] = 1.453650*a0toangst
H1_init[3] = -0.996156*a0toangst

H2_init[1] = -1.453650*a0toangst
H2_init[3] = -0.996156*a0toangst

#H1_init[1]=-rOH*sin(hh_angle)
#H1_init[3]=-rOH*cos(hh_angle)
#
#H2_init[1]=rOH*sin(hh_angle)
#H2_init[3]=-rOH*cos(hh_angle)

#Shift center-of mass to coordinate center#
#COM=zeros(3)
#for i=1:3
#	COM[i]=(mH*(H1_init[i]+H2_init[i])+mO*O_init[i])/(2*mH+mO)
#	O_init[i]=O_init[i]-COM[i]
#	H1_init[i]=H1_init[i]-COM[i]
#	H2_init[i]=H2_init[i]-COM[i]
#end

rO=zeros(3)
rH1=zeros(3)
rH2=zeros(3)

#Rotate molecule by Euler angles around COM#
Rot=rotation_matrix(theta2,phi2,chi)

rO=BLAS.gemv('N' , 1.0, Rot, O_init)
rH1=BLAS.gemv('N' , 1.0, Rot, H1_init)
rH2=BLAS.gemv('N' , 1.0, Rot, H2_init)

#Shift center of mass#
shift = zeros(3)
shift[1] = R*sin(theta1)*cos(phi1)
shift[2] = R*sin(theta1)*sin(phi1)
shift[3] = R*cos(theta1)

for i=1:3
	rO[i] = rO[i] + shift[i]
	rH1[i] = rH1[i] + shift[i]
	rH2[i] = rH2[i] + shift[i]
end

return rO,rH1,rH2
end
###########################################################################################
function lennard_jones(R,epsilon,sigma)

q=sigma/R
V = 4*epsilon*(q^12-q^6)

return V
end
###########################################################################################
function rotation_matrix(theta,phi,chi)

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
