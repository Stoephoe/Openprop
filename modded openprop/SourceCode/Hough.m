% Hough function 
% Returns circumferrential mean axial 'induction factor'  
% Xf:axial distance between propulsors in terms of R (Xf=distance/R) 
% Xf:positive for downstream, negative for upstream. 
 
function[UA_Hough]=Hough(Z,Xf,tanbi,rc,rv) 
q=1+(Xf^2+(rc-rv)^2)/(2*rc*rv); 
s=asin(Xf/sqrt(Xf^2+(rc-rv)^2));        
                                    %amplitude wrt elliptical integrals 
t=sqrt(4*rc*rv/(Xf^2+(rc+rv)^2));   
                                    %t=k (modulus wrt elliptical integrals) 
if rc>rv          
    C1=   Xf/(2*sqrt(rc*rv))*Q2Mhalf(q)-pi/2*Heuman(s,asin(t)); 
else 
    C1=pi+Xf/(2*sqrt(rc*rv))*Q2Mhalf(q)+pi/2*Heuman(s,asin(t)); 
end 

UA_Hough=Z*C1/(2*pi*rv*tanbi); 