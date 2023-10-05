% This function computes the vortex horseshoe axial and tangential 
% interaction influence functions UAHIF_int and UTHIF_int respectively 
% UAHIF_int(n,m)=influence of mth horseshoe vortex shed from one propulsor 
% component (Mact panels) on nth control point of the other component  
% (Mpas panels).  
  
function [UAHIF_int,UTHIF_int]=Horseshoe_int(Mpas,Mact,Zact,TANBIVact,... 
                                        RCpas,RVact,Xf,Hub_Flag,Rhub_oRact) 
UAHIF_int=zeros(Mpas,Mact); 
UTHIF_int=zeros(Mpas,Mact); 
for n=1:Mpas 
    for m=1:Mact+1 
        [UAHough(m)]=Hough(Zact,Xf,TANBIVact(m),RCpas(n),RVact(m)); 
        if Hub_Flag == 1  
            RCW    = RCpas(n); 
            RVW    = Rhub_oRact^2/RVact(m);             
            TANBIW = TANBIVact(m)*RVact(m)/RVW;    
  
            [UAHough_h] = Hough(Zact,Xf,TANBIW,RCW,RVW); 
  
            UAHough(m) = UAHough(m)+UAHough_h; 
        end 
    end 
     
    for m=1:Mact 
        UAHIF_int(n,m)=UAHough(m+1)-UAHough(m); 
        S=(RVact(m)-RCpas(n))*(RVact(m+1)-RCpas(n)); 
        if S<0 && Xf>0 
            UTHIF_int(n,m)=Zact/RCpas(n); 
        else 
            UTHIF_int(n,m)=0; 
        end 
    end 
end 