% This function computes the thrust - torque coefficients, and it computes 
% the total efficiency of the CRP set, Kerwin eqns 161-162, p.138, and 
% eqns 196-197, p. 152, Coney eq. 2-65, p.45 with the inclusion of hub drag 
% 
% ------------------------------------------------------------------------- 
% indices 1,2 refer to the forward and the aft propellers respectively 
% Input Variables: 
    % CD        [ ],    section drag coefficient 
    % RV        [ ],    radius of vortex point / propeller radius 
    % VAC       [ ],    axial inflow velocity at c. points / ship velocity 
    % TANBC     [ ],    tangent of beta at the control points 
    % UASTAR    [ ],    axial      induced velocity / ship velocity 
    % UTSTAR    [ ],    tangential induced velocity / ship velocity 
    % CoD       [ ],    section chord length / propeller diameter 
    % G         [ ],    circulation / (2*pi * prop radius * ship velocity) 
    % RC        [ ],    radius of control point / propeller radius 
    % Fh        [N],    Hub drag 
    % Z         [ ],    number of blades  
    % Js        [ ],    advance coefficient 
    % VMIV      [ ],    Volumetric Mean Inflow Velocity / ship velocity 
    % N         [RPM],  Propeller speed 
    % Vs        [m/s],  Ship speed 
    % R         [m],    Propeller radius 
  
% 
% Output variables: 
    % CT        [ ],    thrust coefficient, eqn (161) p.138 
    % CQ        [ ],    torque coefficient, eqn (161) p.138 
    % CP        [ ],    power coefficient based on torque 
    % KT        [ ],    thrust coefficient, eqn (162) p.138 
    % KQ        [ ],    torque coefficient, eqn (162) p.138 

 
    % EFFY      [ ],    total efficiency of the CRP set 
    % VSTAR     [ ],    total inflow velocity / ship velocity 
% 
% ------------------------------------------------------------------------- 
  
function [CT1,CQ1,KT1,KQ1,CT2,CQ2,KT2,KQ2,EFFY,VSTAR1,VSTAR2] =... 
   Forces(CD1,CD2,DR1,DR2,VAC1,VAC2,TANBC1,TANBC2,... 
          UASTAR1,UASTAR2,UTSTAR1,UTSTAR2,CoD1,CoD2,G1,G2,M1,M2,RC1,RC2,... 
          Fh,Z1,Z2,Js1,Js2,VMIV1,VMIV2,N1,N2,Vs,R1,R2) 
  
VASTAR1    = VAC1         + UASTAR1;         % total axial vel. / ship vel. 
VASTAR2    = VAC2         + UASTAR2;           
VTSTAR1    = VAC1./TANBC1 + UTSTAR1;        % total tang.  vel. / ship vel.  
VTSTAR2    = VAC2./TANBC2 + UTSTAR2;         
VSTAR1     = sqrt(VTSTAR1.^2 + VASTAR1.^2); % total inflow vel. / ship vel. 
VSTAR2     = sqrt(VTSTAR2.^2 + VASTAR2.^2);   
  
sin_BetaI1 = VASTAR1./VSTAR1; 
sin_BetaI2 = VASTAR2./VSTAR2; 
cos_BetaI1 = VTSTAR1./VSTAR1; 
cos_BetaI2 = VTSTAR2./VSTAR2; 
if CD1 < 1 
    DVISC1 = VSTAR1.^2.*CoD1.*CD1/(2*pi);   % normalized viscous drag force 
    DVISC2 = VSTAR2.^2.*CoD2.*CD2/(2*pi);    
else                          % CD > 1 means the input is L/D (legacy code) 
    DVISC1 = VSTAR1.*G1./CD1; 
    DVISC2 = VSTAR2.*G2./CD2; 
end 
  
% ----------------------- Compute CT and CQ, Kerwin eqns. (196-197), p. 152   
CT1 = 0; 
CQ1 = 0; 
CT2 = 0; 
CQ2 = 0; 
  
for m=1:M1 
    CT1 = CT1 + ... 
            (VSTAR1(m)*G1(m)*cos_BetaI1(m)-DVISC1(m)*sin_BetaI1(m))*DR1(m); 
    CQ1 = CQ1 + ... 
     (VSTAR1(m)*G1(m)*sin_BetaI1(m)+DVISC1(m)*cos_BetaI1(m))*RC1(m)*DR1(m); 
end 
for m=1:M2 
    CT2 = CT2 + ... 
            (VSTAR2(m)*G2(m)*cos_BetaI2(m)-DVISC2(m)*sin_BetaI2(m))*DR2(m); 
    CQ2 = CQ2 + ... 
     (VSTAR2(m)*G2(m)*sin_BetaI2(m)+DVISC2(m)*cos_BetaI2(m))*RC2(m)*DR2(m); 
end 
  
CT1   = CT1*4*Z1;             % eqn 196, p.152 (w/ addition for CTD) 
CQ1   = CQ1*4*Z1;             % eqn 197, p.152 
%CP1   = CQ1*pi/Js1;          % power coefficient based on torque 
KT1   = CT1*Js1^2*pi/8;       % eqn 167, p.139 
KQ1   = CQ1*Js1^2*pi/16;      % eqn 167, p.139 
%EFFY1 = CT1*VMIV1/CP1;       % efficiency 
  
 
CT2   = CT2*4*Z2;             % eqn 196, p.152 (w/ addition for CTD) 
CQ2   = CQ2*4*Z2;             % eqn 197, p.152 
%CP2   = CQ2*pi/Js2;          % power coefficient based on torque 
KT2   = CT2*Js2^2*pi/8;       % eqn 167, p.139 
KQ2   = CQ2*Js2^2*pi/16;      % eqn 167, p.139 
%EFFY2 = CT2*VMIV2/CP2;       % efficiency 
  
%--------------- Include effect of hub drag on efficiency ----------------- 
T1= CT1*(.5*1025*Vs^2*pi*R1^2); 
T2= CT2*(.5*1025*Vs^2*pi*R2^2); 
Q1= CQ1*(.5*1025*Vs^2*pi*R1^3); 
Q2= CQ2*(.5*1025*Vs^2*pi*R2^3); 
EFFY=Vs*(VMIV1*T1+VMIV2*T2-Fh)/((2*pi/60)*(N1*Q1+N2*Q2)); %total efficiency 
% Coney, eq. 2-65, p.45 
  
% 
% ===================================================== END Forces Function  
% ========================================================================= 
 
