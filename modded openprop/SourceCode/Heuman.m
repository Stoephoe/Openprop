function [H] = Heuman(phi,alpha)   
[K,E]=ellipke(sin(alpha)^2); 

F=ellipticF(sin(phi),sin(pi/2-alpha));

EE=ellipticE(sin(phi),sin(pi/2-alpha));

%F=mfun('EllipticF',sin(phi),sin(pi/2-alpha));   %Incomplete elliptic  
%                                                integral, 1st kind 
%EE=mfun('EllipticE',sin(phi),sin(pi/2-alpha));  %Incomplete elliptic  
%                                                Integral, 2nd kind 
H=2/pi*(K*EE-(K-E)*F); 
end