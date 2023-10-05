% ------------------------------------------------------- Cavitation inputs
rho  = 1025;        % kg/m^3
H    = 3.048;       % m
g    = 9.81;        % m/s^2
Patm = 101325;      % Pa
Pv   = 2500;        % Pa
Meanline  = 'NACA a=0.8 (modified)';
Thickness = 'NACA 65A010';
t0oc01       = [.1551 .1181 .0902 .0694 .0541 .0419 .0332 .0324... 
                .0204 .005];
t0oc02       = t0oc01;
XR1          = [0.2,0.3,0.4,0.5,0.6,...
                   0.7,0.8,0.9,0.95,1];
XR2          = [0.2,0.3,0.4,0.5,0.6,...
                   0.7,0.8,0.9,0.99,1];

SIGMAs  = (Patm + rho*g*H - Pv)/(0.5*rho*Vs^2);
%SIGMAs  = SIGMAs/10;

[f0octilde, CLItilde, junk1, junk2, junk3, junk4, junk5] = GeometryFoil2D(Meanline,Thickness);
    
    Xf0octilde = f0octilde * ones(size(XR1));
     XCLItilde =  CLItilde * ones(size(XR1));

f0octilde = pchip(XR1,Xf0octilde,RC1);
 CLItilde = pchip(XR1, XCLItilde,RC1);

Xt0oD1 = t0oc01 .* XCoD1;
Xt0oD2 = t0oc02 .* XCoD2;
t0oD1 = pchip(XR1,Xt0oD1 ,RC1);   % section thickness / propeller dia. at ctrl pts
t0oD2 = pchip(XR2,Xt0oD2 ,RC2);   % section thickness / propeller dia. at ctrl pts

% -----------------------------------------------------------------
% Method: (Coney, 1989) cavitation method -- ASSUMES GIVEN THICKNESS DISTRIBUTION t0oD      
% -----------------------------------------------------------------   
SIGMA1 = SIGMAs./VSTAR1.^2;    % local cavitation number
SIGMA2 = SIGMAs./VSTAR2.^2;    % local cavitation number

f0oD1 = (2*pi*G1'./ VSTAR1) .* f0octilde ./ CLItilde;
f0oD2 = (2*pi*G2'./ VSTAR2) .* f0octilde ./ CLItilde;

CoD1  = (8.09*f0oD1+3.033*t0oD1)./(2*SIGMA1) + sqrt((8.09*f0oD1+3.033*t0oD1).^2 + 4* SIGMA1 .* (26.67*f0oD1.^2 + 10*f0oD1.*t0oD1) )./(2*SIGMA1);
CoD2  = (8.09*f0oD2+3.033*t0oD2)./(2*SIGMA2) + sqrt((8.09*f0oD2+3.033*t0oD2).^2 + 4* SIGMA2 .* (26.67*f0oD2.^2 + 10*f0oD2.*t0oD2) )./(2*SIGMA2);
CoD1(M1) = 0;
CoD2(M2) = 0;
% -----------------------------------------------------------------   