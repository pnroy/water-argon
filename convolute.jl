f_para = open("transitions_para.txt") 
f_ortho = open("transitions_ortho.txt") 
lines_para = readlines(f_para)
lines_ortho = readlines(f_ortho)
close(f_para)
close(f_ortho)



Nw=1000
wmax=135
wmin=1
dw=(wmax-wmin)/Nw

fwhmmin=0.01
fwhmmax=2.5
dfwhm=(fwhmmax-fwhmmin)/Nw

Iw=zeros(Float64,Nw)

Npara=length(lines_para)
Northo=length(lines_ortho)

kB_cm=0.695034800
Temperature=20.0
kBT=kB_cm*Temperature

for line in lines_para
    lsplit=split(line)
    E0=parse(Float64,lsplit[2])
    w0=parse(Float64,lsplit[4])
    intensity=parse(Float64,lsplit[5])
    #println(deltaw," ",intensity)
    for i=1:Nw
        w=wmin+(i-1)*dw
        #boltz=exp(-(E0)/kBT)
        boltz=exp(-(E0)/kBT)-exp(-(E0+w0)/kBT)
        #Iw[i]+=intensity*1.0/(1.0+x*x)*boltz
        fwhm=fwhmmin+i*dfwhm
        x=(w-w0)/fwhm/2.0
        Iw[i]+=w*intensity*1.0/(1.0+x*x)*boltz
    end
end

count=0
gs_ortho=0.
for line in lines_ortho
    global count+=1
    lsplit=split(line)
    E0=parse(Float64,lsplit[2])
    if count ==1
        global gs_ortho=E0
    end
    w0=parse(Float64,lsplit[4])
    intensity=parse(Float64,lsplit[5])
    #println(w0," ",intensity)
    for i=1:Nw
        w=wmin+(i-1)*dw
        #boltz=exp(-(E0-gs_ortho)/kBT)
        #boltz=exp(-(E0)/kBT)
        boltz=exp(-(E0)/kBT)-exp(-(E0+w0)/kBT)
        #Iw[i]+=intensity*1.0/(1.0+x*x)*boltz*3
        fwhm=fwhmmin+i*dfwhm
        x=(w-w0)/fwhm/2.0
        Iw[i]+=w*intensity*1.0/(1.0+x*x)*boltz*3
    end
end
for i=1:Nw
    w=wmin+(i-1)*dw
    println(w," ",Iw[i])
end

