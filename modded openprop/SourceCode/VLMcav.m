% Last modified: MAY/01/10 by Dimitrios Laskos  
% Original codes by Hsin-Lung Chung 
% 2D Vortex/Source Lattice with Lighthill Correction Program (VLM) 
% This file contains the algorithms for VLM. 
function [xt, CPU, CPL, CLNum] = VLMcav(N,CL,Alpha,TOC); 

 
% global N CL Alpha TOC; 
%========================================================================== 
U = 1;          c = 1;          % Free Stream Velocity and Chord Length 
for i = 1:N 
   xv(i) = c/2 * (1-cos((i-1/2)*pi/N));     % Vortex Position 
   xc(i) = c/2 * (1-cos(i*pi/N));           % CP Position 
   dx(i) = pi * sqrt(xv(i)*(c-xv(i))) / N;  % Interval between vrotices 
end 
for i = 1:N         % Influence Matrix A (i:CP; j:vortex) 
    for j = 1:N 
        A(i,j) = 1/(2*pi*(xv(j)-xc(i)));         
    end 
end 
% ============================================================= Camber Term 
[B,F,Gexact] = MeanLine(xv,xc);     % Function for NACA a = 0.8 Mean Line 
for i = 1:N 
   B(i) = CL*B(i) - Alpha*pi/180; 
   F(i) = CL*F(i);                  % Camber F 
end 
Gamma = (B/A');                     % Point Vortex Strength 
G = Gamma./dx;                      % Vortex Sheet Strength 
CLNum = 2*sum(Gamma);               % Numerical Lift Coeff 
% ========================================================== Thickness Term 
xt(1) = 0;                          % Thickness at the leading edge 
for i = 1:length(xc) 
   xt(i+1) = xc(i);  
end 
  
% =================================== 
thick_toggle='NACA65A'; % or NACA66TMB 
% =================================== 
  
[RLE,yt,dydx] = Thickness(TOC, xt, thick_toggle); 
for i = 1:N                         % i for CP; j for Vortices 
    for j = 1:N 
        ut(i,j) = (yt(j+1)-yt(j))/(xc(i)-xv(j))/(2*pi); 
    end 
    UT(i) = sum(ut(i,:));           % UT @ Control Points 
end 
UTVP = spline(xc,UT,xv);            % UT @ Vortex Points 
% =========================================== Leading Edge Surface Velocity 
QU = Alpha*pi/180*sqrt(2*c/RLE);    % Surface Velocity 
CPU(1) = QU^2-1;                    % Minus Cp on the upper surface at LE 
CPL(1) = CPU(1);                    % Minus Cp on the lower surface at LE 
% ======================================================== Surface Velocity 
for i = 1:N 
    if dydx(i)>0 
        FLH(i) = 1/sqrt(1+dydx(i)^2); 
    else 
        FLH(i) = 1; 
    end 
    QU(i) = (1+UT(i)+1/2*G(i))*FLH(i);  % Velocity on Upper Surface 
    CPU(i+1) = QU(i)^2-1;               % -Cp     
    QL(i) = (1+UT(i)-1/2*G(i))*FLH(i);  % Velocity on Lower Surface 
    CPL(i+1) = QL(i)^2-1;               % -Cp 
end 

 
% ================================================================ Plotting 
plot_toggle='no'; %or 'yes' 
if strcmp(plot_toggle,'yes') 
    figure; 
    plot(xt,CPU,'-r','LineWidth',2);        hold on; 
    plot(xt,CPL,':b','LineWidth',2);        hold off; 
    grid on;            xlim([0 1]);        xlabel('X/C');   ylabel('-Cp'); 
    legend('Upper Surface','Lower Surface'); 
    % title(strcat('thickness form : ',thick_toggle)) 
    title(['thickness form ',thick_toggle]) 
end 
  
  
  
% Former FORTRAN Subroutine "AEIGHT" ============== APR/27/07 by H.L. Chung 
function [B,F,Gexact] = MeanLine(xv,xc) 
a = 0.8;            % For NACA a=0.8 
MC = length(xv); 
g = -1/(1-a) * (a^2*(log(a)/2-1/4)+1/4); 
h = 1/(1-a) * ((1-a)^2*log(1-a)/2 - (1-a)^2/4) + g; 
AlphaIdeal = -h / (2*pi*(a+1)); 
for i = 1:MC 
   C1 = max(1- xv(i),1e-6); 
   CA = a - xv(i); 
   if (abs(CA)<1e-6) 
       CA = CA+1e-5; 
   end 
      P = 1/2*CA^2*log(abs(CA))-1/2*C1^2*log(C1)+1/4*(C1^2-CA^2); 
      F(i)=(P/(1-a)-xv(i)*log(xv(i))+g-h*xv(i))/(2*pi*(a+1))+C1*AlphaIdeal; 
   if (xv(i)<=a) 
       Gexact(i) = 1/(a+1); 
   else 
       Gexact(i) = 1/(a+1) * (1-xv(i))/(1-a); 
   end 
end 
for j = 1:MC 
   C1 = max(1-xc(j),1e-6); 
   CA = a - xc(j); 
   if (abs(CA)<1e-6) 
       CA = CA+1e-5; 
   end 
   R = -(a-xc(j))*log(abs(CA))-1/2*CA+C1*log(C1)+1/2*C1; 
   S = -1/2*C1+1/2*CA; 
   T = -log(xc(j))-1-h; 
   B(j) = ((R+S)/(1-a)+T)/(2*pi*(a+1)) - AlphaIdeal;  
end 
  
  
  
% Function for thickness ========================  
    function[RLE,YT,DYDX] = Thickness(thk, xt, thick_toggle) 
if strcmp(thick_toggle,'NACA66TMB') 
    PC=[0.000, 0.010, 0.025, 0.050, 0.100, 0.200, 0.300, 0.400, 0.450,... 
      0.500, 0.600, 0.700, 0.800, 0.900, 0.950, 0.975, 0.990, 1.000]; 
  THICK = [0.0000, 0.1870, 0.2932, 0.4132, 0.5814, 0.8000, 0.9274,... 
       0.9904, 1.0000, 0.9917, 0.9256, 0.7934, 0.5950, 0.3306,... 
 
       0.1736, 0.0888, 0.0360, 0.0000]; 
  RLE_CONST = 0.448; 
elseif strcmp(thick_toggle,'NACA65A') 
    PC = [0.000, 0.005, 0.0075, 0.0125, 0.0250, 0.05, 0.075, 0.1, 0.15,... 
      0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,... 
      0.75, 0.8, 0.85, 0.9, 0.95, 1]; 
  THICK = [0, 0.1556, 0.1879, 0.2387, 0.3265, 0.4379, 0.5311, 0.6094,... 
       0.7331, 0.8268, 0.8978, 0.9495, 0.9835, 1.0000, 0.9976, 0.9736,... 
       0.9272, 0.8612, 0.7803, 0.6868, 0.5828, 0.4709, 0.3548, 0.2382,... 
       0.1208, 0.0031]; 
   RLE_CONST = 0.654; 
end 
NT = length(xt); 
RLE = RLE_CONST*thk^2; 
PSQ = sqrt(PC); 
TRLE = 2*sqrt(2*RLE_CONST); 
XSQ = sqrt(xt); 
YSPLN = spline(PSQ,THICK,XSQ); 
YT = thk.*YSPLN; 
for N=1:NT-1 
   DYDX(N) = (YT(N+1)-YT(N))/(xt(N+1)-xt(N))/2; 
end 
 