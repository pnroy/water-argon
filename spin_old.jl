module spin

using LinearAlgebra

export spin_isomer,spin_isomer_matrix,spin_isomer_matrix_complex,spin_dimension
####################################################################
function spin_isomer(isomer,jmax,mmax,kmax,Utrans,V1,V2,Nprod)

#Get inidices of eigenbasis sorted by the j(j+1) values#
jlist=jstates(Utrans,jmax,mmax,kmax)

#############################################
#Determine the indices for the allowed basis#
#function for the chosen spin isomer        #
#############################################
para=[]
ortho=[]
iso=0
i=0
for j=0:jmax
	for i1=1:(2*j+1)
		for i2=1:(2*j+1)
			i+=1
			if iso == 0
				append!(para,jlist[i])
			elseif iso == 1
				append!(ortho,jlist[i])
			end
		end
		if iso == 0
			iso=1
		elseif iso == 1
			iso=0
		end
	end
end

indices=[]

#Choose para or ortho and sort the indices back#
if isomer == "para"
	indices=sort(para)
elseif isomer == "ortho"
	indices=sort(ortho)
end

Nspec=length(indices)

#################################################
#Construct new kinetic and potential matrices in#
#Hilbert space associated with the respective   #
#spin isomer                                    #
#################################################

V1new=zeros(Float64,(Nspec,Nspec,Nprod))
V2new=zeros(Float64,(Nspec,Nspec,Nprod))

i=0
for i1 in indices
    	i+=1
	j=0
	for j1 in indices
    		j+=1
		V1new[i,j,:].=real(V1[i1,j1,:])
		V2new[i,j,:].=real(V2[i1,j1,:])
    	end
end

return V1new,V2new
end
##########################################################################
function spin_dimension(isomer,jmax,mmax,kmax,Utrans)

#Get inidices of eigenbasis sorted by the j(j+1) values#
jlist=jstates(Utrans,jmax,mmax,kmax)

#############################################
#Determine the indices for the allowed basis#
#function for the chosen spin isomer        #
#############################################
para=[]
ortho=[]
iso=0
i=0
for j=0:jmax
	for i1=1:(2*j+1)
		for i2=1:(2*j+1)
			i+=1
			if iso == 0
				append!(para,jlist[i])
			elseif iso == 1
				append!(ortho,jlist[i])
			end
		end
		if iso == 0
			iso=1
		elseif iso == 1
			iso=0
		end
	end
end

indices=[]

#Choose para or ortho and sort the indices back#
if isomer == "para"
	indices=sort(para)
elseif isomer == "ortho"
	indices=sort(ortho)
end

Nspec=length(indices)

return Nspec
end
##########################################################################
function spin_isomer_matrix(isomer,jmax,mmax,kmax,Utrans,matrix)

#Get inidices of eigenbasis sorted by the j(j+1) values#
jlist=jstates(Utrans,jmax,mmax,kmax)

#############################################
#Determine the indices for the allowed basis#
#function for the chosen spin isomer        #
#############################################
para=[]
ortho=[]
iso=0
i=0
for j=0:jmax
	for i1=1:(2*j+1)
		for i2=1:(2*j+1)
			i+=1
			if iso == 0
				append!(para,jlist[i])
			elseif iso == 1
				append!(ortho,jlist[i])
			end
		end
		if iso == 0
			iso=1
		elseif iso == 1
			iso=0
		end
	end
end

indices=[]

#Choose para or ortho and sort the indices back#
if isomer == "para"
	indices=sort(para)
elseif isomer == "ortho"
	indices=sort(ortho)
end

Nspec=length(indices)

#################################################
#Construct new kinetic and potential matrices in#
#Hilbert space associated with the respective   #
#spin isomer                                    #
#################################################

matrix_new=zeros(Float64,(Nspec,Nspec))

i=0
for i1 in indices
    	i+=1
	j=0
	for j1 in indices
    		j+=1
		matrix_new[i,j]=real(matrix[i1,j1])
    	end
end

return matrix_new
end
##########################################################################
function spin_isomer_matrix_complex(isomer,jmax,mmax,kmax,Utrans,matrix)

#Get inidices of eigenbasis sorted by the j(j+1) values#
jlist=jstates(Utrans,jmax,mmax,kmax)

#############################################
#Determine the indices for the allowed basis#
#function for the chosen spin isomer        #
#############################################
para=[]
ortho=[]
iso=0
i=0
for j=0:jmax
	for i1=1:(2*j+1)
		for i2=1:(2*j+1)
			i+=1
			if iso == 0
				append!(para,jlist[i])
			elseif iso == 1
				append!(ortho,jlist[i])
			end
		end
		if iso == 0
			iso=1
		elseif iso == 1
			iso=0
		end
	end
end

indices=[]

#Choose para or ortho and sort the indices back#
if isomer == "para"
	indices=sort(para)
elseif isomer == "ortho"
	indices=sort(ortho)
end

Nspec=length(indices)

#################################################
#Construct new kinetic and potential matrices in#
#Hilbert space associated with the respective   #
#spin isomer                                    #
#################################################

matrix_new=zeros(ComplexF64,(Nspec,Nspec))

i=0
for i1 in indices
    	i+=1
	j=0
	for j1 in indices
    		j+=1
		matrix_new[i,j]=matrix[i1,j1]
    	end
end

return matrix_new
end
##############################################################
function jstates(Utrans,jmax,mmax,kmax)

Nspec=size(Utrans,1)
#Define J^2 operator in Wigner basis#
J2=zeros(ComplexF64,(Nspec,Nspec))
i=0
for j=0:jmax
	mj=min(j,mmax)
	kj=min(j,kmax)
	for m=-mj:mj
	for k=-kj:kj
		i+=1
		J2[i,i]=j*(j+1)
	end
	end
end

#Transform J^2 to real Wigner and eigenstate basis#
tmp=BLAS.gemm('C', 'N', Utrans,J2)
J2new=BLAS.gemm('N','N', tmp,Utrans)

jvalue=zeros(Float64,Nspec)

for i=1:Nspec
	jvalue[i]=round(real(J2new[i,i]))
end

#Create list#
jlist=[(i,jvalue[i]) for i=1:Nspec]
#Sort list by j(j+1)#
jlist_ord=sort(jlist, by = x -> x[2])
list_final=zeros(Int64,Nspec)
list_final.=first.(jlist_ord)


return list_final

end
##############################################################
end
