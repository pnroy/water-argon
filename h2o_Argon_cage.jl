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
function potential(R,theta,phi,chi,vmat,Ngrid)
	#Conversion factor#
	a0toangst=0.529177210903
	eHtocm1=219474.631363	#Eh to cm-1

	#Read coordinates of cage atoms#
	f = open("Ar.xyz") 
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


	rmax=4.0
	rmin=2.5
	dr=(rmax-rmin)/Ngrid
	dt=pi/Ngrid
	dp=2*pi/Ngrid
	
	gridV=0

	rAr_COM=zeros(Float64,3)
	V=0.0 
	for i=1:Ncage
#	Threads.@threads for i=1:Ncage
		rAr_COM=Aratoms[i,:]-R
		dist=norm(rAr_COM)
		# get theta and phi from rotation matrix
		#Rotate molecule by Euler angles around COM, matrix is for BFF to SFF, use transpose
		Rot=rotation_matrix(theta,phi,chi)
		rAr_BFF=BLAS.gemv('T' , 1.0, Rot, rAr_COM)
		theta_BFF=acos(rAr_BFF[3]/dist)
		# atan2 below
		phi_BFF=atan(rAr_BFF[2],rAr_BFF[1])
		if gridV==1
		# call grid potential
			ir=floor(Int,1+(dist-rmin)/dr)
			if ir<1
				ir=1
			end
			if ir>Ngrid
				ir=Ngrid
			end
			it=floor(Int,1+(theta_BFF)/dt)
			if it<1
				it=1
			end
			if it>Ngrid
				it=Ngrid
			end
			ip=floor(Int,1+(phi_BFF)/dp)
			if ip<1
				ip=1
			end
			if ip>Ngrid
				ip=Ngrid
			end
			value=vmat[ir,it,ip]
		else
			# call fortran code
			#	cmd_string=string("`./a.out ",dist," ",theta_BFF," ",phi_BFF,"`")
			#	expr = Meta.parse(cmd_string)
			#	cmd_object = eval(expr)
			#	value=parse(Float64,readchomp(cmd_object))
				result_ref = Ref{Float64}()
				ccall((:h2oar_3dpes_, "./potlibarh20.dylib"), Nothing, (Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},Ref{Float64}), theta_BFF, phi_BFF,dist,result_ref,1)
				value=result_ref[]
		end
		V+=value
	end
	return V
end


#######################
let

	Ngrid=3
	rmax=7.0
	rmin=2.5
	dr=(rmax-rmin)/Ngrid
	dt=pi/Ngrid
	dp=2*pi/Ngrid

	#pre-compute potential on grid
	vmat=zeros(Float64,(Ngrid,Ngrid,Ngrid))
	# for ir=1:Ngrid
	# 	r=rmin+(ir-1)*dr
	# 	for it=1:Ngrid
	# 		theta=(it-1)*dt
	# 		for ip=1:Ngrid
	# 			phi=(ip-1)*dp
	# 			cmd_string=string("`./a.out ",r," ",theta," ",phi,"`")
	# 			expr = Meta.parse(cmd_string)
	# 			cmd_object = eval(expr)
	# 			vmat[ir,it,ip]=parse(Float64,readchomp(cmd_object))
	# 		end
	# 	end
	# end


	N=50

	R=zeros(Float64,3)

	theta=0.0
	phi=0.0
	chi=0.0
	dphi=2*pi/N
	f=open("phi.dat","w")
	for i=0:N
		phi=i*dphi
		value=potential(R,theta,phi,chi,vmat,Ngrid)
		println(f,phi," ",value)
	end
	close(f)

	theta=0.0
	phi=0.0
	chi=0.0
	dchi=2*pi/N
	f=open("chi.dat","w")
	for i=0:N
		chi=i*dchi
		value=potential(R,theta,phi,chi,vmat,Ngrid)
		println(f,chi," ",value)
	end
	close(f)

	theta=0.0
	phi=0.0
	chi=0.0
	dtheta=pi/N
	f=open("theta.dat","w")
	for i=0:N
		theta=i*dtheta
		value=potential(R,theta,phi,chi,vmat,Ngrid)
		println(f,theta," ",value)
	end
	close(f)

	theta=0.0
	phi=0.0
	chi=0.0
	Rmax=8.0
	dR=Rmax/N
	f=open("R.dat","w")

	f2=open("R.xyz","w")

	println(f2,N)

	println(f2," ")

	Theta=pi/4.0
	Theta=0.0
	Phi=0
	for i=0:N
		Rmag=i*dR
		R[1]=Rmag*sin(Theta)*cos(Phi)
		R[2]=Rmag*sin(Theta)*sin(Phi)
		R[3]=Rmag*cos(Theta)
		value=potential(R,theta,phi,chi,vmat,Ngrid)
		println(f,Rmag," ",value)
		println(f2,"H2O ",R[1]," ",R[2]," ",R[3])

	end
	close(f)
	close(f2)

end
