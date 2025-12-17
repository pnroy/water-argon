module spin

using LinearAlgebra

export spin_isomer, spin_isomer_half, spin_isomerPN
####################################################################
function spin_isomer(isomer,jmax,matrix)

#Calculate dimension of para/ortho Hilbert space#

para = []
ortho = []
i=0
for j=0:jmax
for m=-j:j
for k=-j:j
	i+=1
	if mod(k,2) == 0
		append!(para,i)
	elseif mod(k,2) != 0
		append!(ortho,i)
	end
end
end
end

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

matrix_out=zeros(ComplexF64,(Nspec,Nspec))

i=0
for i1 in indices
    i+=1
	j=0
	for j1 in indices
    	j+=1
		matrix_out[i,j]=matrix[i1,j1]
    end
end

return matrix_out
end
##########################################################################


function spin_isomerPN(isomer,jmax,matrix)

#Calculate dimension of para/ortho Hilbert space#

para = []
ortho = []
i=0
for j=0:jmax
for m=-j:j
for k=-j:j
	i+=1
	if mod(k,2) == 0
		append!(para,i)
	elseif mod(k,2) != 0
		append!(ortho,i)
	end
end
end
end

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

matrix_out=zeros(ComplexF64,(Nspec,Nspec))

# i=0
# for j=0:jmax
# for m=-j:j
# for k=-j:j
# 	if isomer == "para"
# 		if mod(k,2) == 0
# 			i+=1
# 	i=0
# for j=0:jmax
# for m=-j:j
# for k=-j:j
# 	if mod(k,2) == 0
# 		i+=1
		




# 		matrix_out[i,ip]=matrix[i1,j1]
#     end
#end

return matrix_out
end
##########################################################################


function spin_isomer_half(isomer,jmax,matrix)

#Calculate dimension of para/ortho Hilbert space#

para = []
ortho = []
i=0
for j=0:jmax
for m=-j:j
for k=-j:j
	i+=1
	if mod(k,2) == 0
		append!(para,i)
	elseif mod(k,2) != 0
		append!(ortho,i)
	end
end
end
end

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

mmax=2*jmax+1
kmax=2*jmax+1
N=0
for j=0:jmax
	mj=min(mmax,j)
	kj=min(kmax,j)
	N=N+(2*mj+1)*(2*kj+1)
end


matrix_out=zeros(ComplexF64,(N,Nspec))

i=0
for i=1:N
	j=0
	for j1 in indices
		j+=1
		matrix_out[i,j]=matrix[i,j1]
	end
end

return matrix_out
end
##########################################################################
end