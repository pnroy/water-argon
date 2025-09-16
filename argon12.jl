using LinearAlgebra
using LinearMaps
using Arpack
push!(LOAD_PATH,pwd())
l=[-1/sqrt(2),1/sqrt(2)]

#distance from lattice constant
dist=5.311/sqrt(2)

#distance from PES
#dist=3.63
l=dist*l

f=open("Ar.xyz","w")
println(f,"12")
println(f," ")
atom_index="Ar"
for i=1:2
    for j=1:2
        x=l[i]
        y=l[j]
        z=0
        #global atom_index+=1
        #println("Ar ",x," ",y," ",z)
        println(f,atom_index," ",x," ",y," ",z)
        x=l[i]
        z=l[j]
        y=0
        #global atom_index+=1
        println(f,atom_index," ",x," ",y," ",z)
        y=l[i]
        z=l[j]
        x=0
        #global atom_index+=1
        println(f,atom_index," ",x," ",y," ",z)
    end
end
close(f)