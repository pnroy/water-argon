using LinearAlgebra
using Optim
using ReverseDiff

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
function lo_potential(R,theta1,phi1,theta2,phi2,chi,dCI,kconst)
eHtocm1=219474.631363	#Eh to cm-1
a0toangst=0.529177210903 #Bohr to angstrom

rCI=zeros(3)
rCI[3] = dCI

R=R/a0toangst
Rot=rotation_matrix(theta2,phi2,chi)

rCI_rot=BLAS.gemv('N' , 1.0, Rot, rCI)

COM = zeros(3)
COM[1] = R*sin(theta1)*cos(phi1)
COM[2] = R*sin(theta1)*sin(phi1)
COM[3] = R*cos(theta1)

V=0.5*kconst*(norm(COM)^2+norm(rCI_rot)^2+2*dot(COM,rCI_rot))

V=V*eHtocm1
return V
end

###########################################################################################
function potential(R,theta,phi,chi)

	#Read coordinates of cage atoms#
	f = open("Ar.xyz") 
	#f = open("lattice.xyz") 
	lines = readlines(f)
	close(f)
	lsplit = split(lines[1])
	Ncage = parse(Int64, lsplit[1])

	Aratoms = zeros(Float64,(Ncage,3))
	for i=1:Ncage
		lsplit = split(lines[i+2])
		Aratoms[i,1] = parse(Float64,lsplit[2])
		Aratoms[i,2] = parse(Float64,lsplit[3])
		Aratoms[i,3] = parse(Float64,lsplit[4])
	end


	rAr_COM=zeros(Float64,3)
	V=0.0 
	for i=1:Ncage
#	Threads.@threads for i=1:Ncage
		rAr_COM[1]=Aratoms[i,1]-R[1]
		rAr_COM[2]=Aratoms[i,2]-R[2]
		rAr_COM[3]=Aratoms[i,3]-R[3]
		dist=norm(rAr_COM)
		# get theta and phi from rotation matrix
		#Rotate molecule by Euler angles around COM, matrix is for BFF to SFF, use transpose
		Rot=rotation_matrix(theta,phi,chi)
		rAr_BFF=BLAS.gemv('T' , 1.0, Rot, rAr_COM)
		theta_BFF=acos(rAr_BFF[3]/dist)
		# atan2 below
		phi_BFF=atan(rAr_BFF[2],rAr_BFF[1])

		# call fortran code
		#	cmd_string=string("`./a.out ",dist," ",theta_BFF," ",phi_BFF,"`")
		#	expr = Meta.parse(cmd_string)
		#	cmd_object = eval(expr)
		#	value=parse(Float64,readchomp(cmd_object))
			result_ref = Ref{Float64}()
			ccall((:h2oar_3dpes_, "./potlibarh20.dylib"), Nothing, (Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},Ref{Int}), phi_BFF,theta_BFF ,dist,result_ref,1)
			value=result_ref[]
		V+=value
	end
	return V
end




#######################
let
	#Conversion factors#
	a0toangst=0.529177210903 #Bohr to angstrom
	eHtocm1=219474.631363	#Eh to cm-1
	utome=1.0/(5.48579909065*1e-4)
	eHtoJ=4.3597447222071e-18
	eHtokHz=1.0/1.51982850071586e-13
	a0tom=0.529177210903e-10 
	a0topm=52.9177210903

	# # generate cage coordinates
	# l=[-1/sqrt(2),1/sqrt(2)]
	# l=3.72*l
	# f=open("Ar.xyz","w")
	# println(f,"12")
	# println(f," ")
	# atom_index="Ar"
	# for i=1:2
	# 	for j=1:2
	# 		x=l[i]
	# 		y=l[j]
	# 		z=0
	# 		#global atom_index+=1
	# 		#println("Ar ",x," ",y," ",z)
	# 		println(f,atom_index," ",x," ",y," ",z)
	# 		x=l[i]
	# 		z=l[j]
	# 		y=0
	# 		#global atom_index+=1
	# 		println(f,atom_index," ",x," ",y," ",z)
	# 		y=l[i]
	# 		z=l[j]
	# 		x=0
	# 		#global atom_index+=1
	# 		println(f,atom_index," ",x," ",y," ",z)
	# 	end
	# end
	# close(f)

	kconst=(4.0/eHtoJ)*(a0tom^2)
	dCI=10.0/a0topm	

	N=50

	R=zeros(Float64,3)
	Theta=pi/3
	Phi=pi/7
	Rmag=1.0

	theta=pi/3
	phi=0.0
	chi=pi/7
	dphi=2*pi/N
	f=open("phi.dat","w")
	for i=0:N
		phi=i*dphi
		value=potential(R,theta,phi,chi)
		value_lo=lo_potential(Rmag,Theta,Phi,theta,phi,chi,dCI,kconst)
		println(f,phi," ",value," ",value_lo)
	end
	close(f)

	theta=0.0
	phi=0.0
	chi=0.0
	dchi=2*pi/N
	f=open("chi.dat","w")
	for i=0:N
		chi=i*dchi
		value=potential(R,theta,phi,chi)
		value_lo=lo_potential(Rmag,Theta,Phi,theta,phi,chi,dCI,kconst)
		println(f,chi," ",value," ",value_lo)
	end
	close(f)

	theta=0.0
	phi=0.0
	chi=0.0
	dtheta=pi/N
	f=open("theta.dat","w")
	for i=0:N
		theta=i*dtheta
		value=potential(R,theta,phi,chi)
		value_lo=lo_potential(Rmag,Theta,Phi,theta,phi,chi,dCI,kconst)
		println(f,theta," ",value," ",value_lo)
	end
	close(f)

	theta=0.0
	phi=0.0
	chi=0.0
	Rmax=.3
	dR=Rmax/N
	f=open("R.dat","w")

	f2=open("R.xyz","w")

	println(f2,N)

	println(f2," ")


	Theta=pi/4.0
	#Theta=pi/2.0
	Theta=0.0
	Phi=pi
	for i=0:N
		Rmag=i*dR
		R[1]=Rmag*sin(Theta)*cos(Phi)
		R[2]=Rmag*sin(Theta)*sin(Phi)
		R[3]=Rmag*cos(Theta)
		value=potential(R,theta,phi,chi)
		value_lo=lo_potential(Rmag,Theta,Phi,theta,phi,chi,dCI,kconst)
		println(f,Rmag," ",value," ",value_lo)
		println(f2,"H2O ",R[1]," ",R[2]," ",R[3])

	end
	close(f)
	close(f2)

end
