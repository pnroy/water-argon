module h2o_lo

using LinearAlgebra

export lo_potential

###########################################################################################
function lo_potential(R,theta1,phi1,theta2,phi2,chi,dCI,kconst)

rCI=zeros(3)
rCI[3] = dCI

Rot=rotation_matrix(theta2,phi2,chi)

rCI_rot=BLAS.gemv('N' , 1.0, Rot, rCI)

COM = zeros(3)
COM[1] = R*sin(theta1)*cos(phi1)
COM[2] = R*sin(theta1)*sin(phi1)
COM[3] = R*cos(theta1)

V=0.5*kconst*(norm(COM)^2+norm(rCI_rot)^2+2*dot(COM,rCI_rot))

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
