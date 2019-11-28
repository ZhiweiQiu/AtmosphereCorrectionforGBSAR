function Phuidu=bilinc(Baseimg,Bx,By)


Bxz=fix(Bx);
Byz=fix(By);

detx=Bx-Bxz;
dety=By-Byz;

if (detx==0 && dety==0)
    Phuidu=Baseimg(Bx,By);
else
    Phuidu=Baseimg(Bxz,Byz)*(1-detx)*(1-dety)+Baseimg(Bxz,Byz+1)*(1-detx)*dety+detx*(1-dety)*Baseimg(Bxz+1,Byz)+detx*dety*Baseimg(Bxz+1,Byz+1);
end