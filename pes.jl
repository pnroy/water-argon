using LinearAlgebra


let

	N=29

	rmax=4.0
	rmin=2.5
	dr=(rmax-rmin)/N
	dt=pi/N
	dp=2*pi/N

	f=open("pes.dat","w")
	for ir=0:N
		r=rmin+ir*dr
		for it=0:N
			theta=it*dt
			for ip=0:N
				phi=ip*dp
					cmd_string=string("`./a.out ",r," ",theta," ",phi,"`")
					expr = Meta.parse(cmd_string)
					cmd_object = eval(expr)
					v=parse(Float64,readchomp(cmd_object))
					println(f,r," ",theta," ",phi," ",v)
			end
		end
	end
	close(f)

end
