%Q2Mhalf: Legendre fuction of the second kind and minus half order 
        %Ref: Handbook of Math Functions, Abramowitz and Stegun, 1972 
        %section 8.13.3, p.337 
        %uses modulus k for elliptic integrals (m=k^2) 
        %ellipke uses parameter m.  
         
function [Q2M] = Q2Mhalf(q)        
  
k=sqrt(2/(q+1)); 
[K,E]=ellipke(k^2); 
Q2M=k*K; 
end
  
%checked with ref p.340 example 
%Validated with the National Bureau of Standards Tables of  
%Associated Legendre Functions  
%(Columbia University Press, New York, 1945), p.264. 
%%From the tables: Q2Mhalf(1.5)=2.01891, Q2Mhalf(2.7)=1.38958, 
%Q2Mhalf(6)=0.911696, Q2Mhalf(8.4)=0.768523, Q2Mhalf(10)=0.703806 
  