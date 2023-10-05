% ========================================================================= 
% Contra-Rotating Propeller Design Code based on lifting line theory 
% Numerical Implementation of the Variational Optimization Method 
% for Two-Component Propulsors developed by Kerwin, et al. (1986) 
% 
% -------------------- Copyright 2010 Dimitrios Laskos -------------------- 
% This program is free software.  You can redistribute it and/or modify it 
% under the terms of the GNU General Public License version 2, as published 
% by the Free Software Foundation.  This program is distributed in the hope  
% that it will be useful, but WITHOUT ANY WARRANTY; without even the  
% implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License for more details. 
% 
% -------------------------------------------------------------- References 
% 1) J.S. Carlton, "Marine Propellers & Propulsion", chapter 3, 1994. 
% 2) J. Kerwin, W. Coney & C. Hsin, "Optimum Circulation Distributions for 
%           Single and Multi-Component Propulsors", 21st American Towing  
%           Tank Conference (ATTC), 1986. 
% 3) W. Coney, "A Method for the Design of a Class of Optimum Marine 
%           Propulsors", Ph.D. thesis, MIT, 1989. 
% 4) J. Kerwin, "Hydrofoils and Propellers", MIT Course 2.23 notes, 2007. 
% 5) J.W. Wrench, "The Calculation of Propeller Induction Factors", David  
%           Taylor Model Basin (Technical Report 1116), Feb. 1957.  
% 6) M. Wang, "Hub Effects in Propeller Design and Analysis",Ph.D. thesis, 
%           MIT, 1985. 
% 7) G. Hough & D. Ordway, "Generalized Actuator Disk", Developements in 
%           Theoretical and Applied Mechanics, Vol.2, pp.23-31, 1965.  
% 8) B. Epps et al., "OpenProp: An open-source parametric design tool for 
%           propellers", Grand Challenges in Modeling & Simulation 

 
%           Conference (GCMC '09), 2009. 
% ========================================================================= 
% Last modified: 05/03/2010 by Dimitrios Laskos  
% ------------------------------------------------------------------------- 
% indices 1,2 refer to the forward and the aft propellers respectively 
% Input Variables: 
    % 
    % Rhub         [m],    Hub radius (common for both propellers) 
    % R1,R2        [m],    Propeller radius 
    % M1,M2        [ ],    Number of vortex panels over the radius 
    % Z1,Z2        [ ],    Number of blades 
    % Tr           [N],    Required total thrust  
    % q            [ ],    torque ratio Q2/Q1 
    % N1,N2        [RPM],  Propeller speed 
    % XR1,XR2      [ ],    Radial locations for defining inflow velocities   
    %                      and geometric properties 
    % XCoD1,XCoD2  [ ],    chord / propeller diameter 
    % CD           [ ],    section drag coefficient 
    % XVA1,XVA2    [ ],    Va/Vs, axial inflow vel. / ship vel. 
    % XVT1,XVT2    [ ],    Vt/Vs, tangential inflow vel. / ship vel. 
    % Vs           [m/s],  Ship speed 
    % Xf           [m],    Axial separation between propellers 
    % ITER         [ ],    Max. Iterations for Circulation Convergence and  
    %                      Wake Alignment  
    % spacing              Type of radial spacing ('cosine' or 'constant') 
    % Hub_Flag             Inclusion of hub effects (1=YES, 0=NO) 
    % Rhv                  Hub Vortex Radius/Hub Radius 
       
% Output variables: 
    % EFFY         [ ],    total efficiency 
    % CT           [ ],    thrust coefficient, eqn (161) p.138 
    % CQ           [ ],    torque coefficient, eqn (161) p.138 
    % KT           [ ],    thrust coefficient, eqn (162) p.138 
    % KQ           [ ],    torque coefficient, eqn (162) p.138 
    % G            [ ],    non-dimensional Circulation 
    % UA_SELF      [ ],    Axial self-induced velocity vector / Vs 
    % UT_SELF      [ ],    Tangential self-induced velocity vector / Vs 
    % UA_INT1_2    [ ],    Axial interaction velocity vector on component 1 
    %                      induced by component 2 / Vs 
    % UT_INT1_2    [ ],    Tangential interaction velocity vector on  
    %                      component 1 induced by component 2 / Vs 
    % TANBIC       [ ],    Tangent of hydrodynamic pitch angle 
    % VSTAR        [ ],    Total inflow velocity / Vs 
    % Cl           [ ],    Required lift coefficient 
     
% 
% ------------------------------------------------------------------------- 
  
function[EFFY,CT1,CT2,CQ1,CQ2,KT1,KT2,KQ1,KQ2,RC1,RC2,G1,... 
             G2,UA_SELF1,UT_SELF1,UA_INT1_2,UT_INT1_2,UA_SELF2,UT_SELF2,... 
             UA_INT2_1,UT_INT2_1,TANBIC1,TANBIC2,VSTAR1,VSTAR2,Cl1,Cl2,...
             Np,X3D,Y3D,Z3D,X3D_aft,Y3D_aft,Z3D_aft,CoD1,CoD2,A]=... 
        CoupledCRP(Rhub,R1,R2,M1,M2,Z1,Z2,Tr,q,N1,N2,XR1,XR2,XCoD1,... 
              XCoD2,XVA1,XVA2,XVT1,XVT2,Vs,Xf,ITER,spacing,Hub_Flag,Rhv...
              ,geometry_flag,CD)  
  
% ------------------------ Apply spacing ------------------------- 
RV1=zeros(1,M1+1);RC1=zeros(1,M1);                 % initialize RC1 and RV1 
 
 
RV2=zeros(1,M2+1);RC2=zeros(1,M2);                 % initialize RC2 and RV2 
RoR=1; 
Rhub_oR1=Rhub/R1; 
Rhub_oR2=Rhub/R2; 
if strcmp(spacing,'constant')==1             %Constant spacing 
    if  Hub_Flag==0 
        DRR1 = (RoR-Rhub_oR1)/(M1+.5);       % panel size 
        DRR2 = (RoR-Rhub_oR2)/(M2+.5); 
        RV1(M1+1)=RoR-.25*DRR1;              % 25% tip inset 
        RV2(M2+1)=RoR-.25*DRR2; 
        RV1(1)=Rhub_oR1+.25*DRR1;            % 25% hub inset (NO IMAGE HUB) 
        RV2(1)=Rhub_oR2+.25*DRR2;  
    elseif Hub_Flag==1 
        DRR1 = (RoR-Rhub_oR1)/(M1+.25);      % panel size 
        DRR2 = (RoR-Rhub_oR2)/(M2+.25); 
        RV1(M1+1)=RoR-.25*DRR1;              % 25% tip inset 
        RV2(M2+1)=RoR-.25*DRR2; 
        RV1(1)=Rhub_oR1;                     % 25% hub inset (NO IMAGE HUB) 
        RV2(1)=Rhub_oR2; 
    end 
     
    RC1(1)=RV1(1)+.5*DRR1;                   % ctrl pt at mid-panel 
    for m=2:M1 
        RV1(m)=RV1(m-1)+DRR1; 
        RC1(m)=RC1(m-1)+DRR1; 
    end 
    RC2(1)=RV2(1)+.5*DRR2;                   % ctrl pt at mid-panel 
    for m=2:M2 
        RV2(m)=RV2(m-1)+DRR2; 
        RC2(m)=RC2(m-1)+DRR2; 
    end 
     
elseif strcmp(spacing,'cosine')==1           %Cosine spacing 
     
    DEL1 = pi/(2*M1);             
    Rdif1  = 0.5*(RoR - Rhub_oR1); 
    for m = 1:M1+1 
        RV1(m) = Rhub_oR1 + Rdif1*(1-cos(2*(m-1)*DEL1));  
    end 
    for n = 1:M1 
        RC1(n) = Rhub_oR1 + Rdif1*(1-cos((2*n-1)*DEL1)); 
    end 
    DEL2 = pi/(2*M2);             
    Rdif2  = 0.5*(RoR - Rhub_oR2); 
    for m = 1:M2+1 
        RV2(m) = Rhub_oR2 + Rdif2*(1-cos(2*(m-1)*DEL2));  
    end 
    for n = 1:M2 
        RC2(n) = Rhub_oR2 + Rdif2*(1-cos((2*n-1)*DEL2)); 
    end 
end 
% ------------------------------------------------------------------------- 
DR1=diff(RV1);DR2=diff(RV2); 
% ---------Interpolate Va,Vt and CoD at vortex and control points---------- 
VAC1 = pchip(XR1,XVA1,RC1);     % axial inflow vel. / ship vel. at ctrl pts 
VTC1 = pchip(XR1,XVT1,RC1);     % tang. inflow vel. / ship vel. at ctrl pts 

 
CoD1 = pchip(XR1,XCoD1,RC1);    % chord / propeller diameter at ctrl pts 
VAC2 = pchip(XR2,XVA2,RC2);      
VTC2 = pchip(XR2,XVT2,RC2);      
CoD2 = pchip(XR2,XCoD2,RC2);     
% ------------------------------------------------------------------------- 
Js1=Vs/((N1/60)*2*R1); %n=N/60 [rev/sec] 
Js2=Vs/((N2/60)*2*R2); 
om1=pi/Js1;            %tip speed ratio 
om2=pi/Js2; 
  
%--------------- Initialize induced velocity vectors ---------------------- 
UTSTAR1(1:M1)=0; 
UTSTAR2(1:M2)=0; 
UASTAR1(1:M1)=0; 
UASTAR2(1:M2)=0; 
  
% ------------------------- Assign Initial Values ------------------------- 
G1_last=0; 
G2_last=0; 
LT_last=-1; 
LQ_last=0; 
A=zeros(M1+M2+2); 
B=zeros(M1+M2+2,1); 
% ------------ Initial estimates for hydrodynamic pitch angles ------------ 
[TANBIC1,TANBIV1] = find_tan_BetaI(VAC1,VTC1,UASTAR1,UTSTAR1,RC1,RV1,Js1);   
[TANBIC2,TANBIV2] = find_tan_BetaI(VAC2,VTC2,UASTAR2,UTSTAR2,RC2,RV2,Js2); 
TANBC1=TANBIC1; 
TANBC2=TANBIC2; 
% ------------------------------------------------------------------------- 
% Iteration for betaI's. BetaI's are fixed. 
B_iter=1; 
B1_res=1; 
B2_res=1; 
B_res=[B1_res B2_res]; 
TANBIC1_last=TANBIC1; 
TANBIC2_last=TANBIC2; 
% ------------ Compute Horseshoe Influence Functions ---------------------- 
% UAHIF1_2 is the horseshoe influence matrix for the axial interaction 
% velocities induced by component 2 (aft) on component 1 (forward) 
[UAHIF1,UTHIF1]=Horseshoe_self(M1,Z1,TANBIV1,RC1,RV1,Hub_Flag,Rhub_oR1); 
[UAHIF1_2,UTHIF1_2]=Horseshoe_int(M1,M2,Z2,TANBIV2,RC1,RV2,-Xf,Hub_Flag,... 
                                                                 Rhub_oR2); 
[UAHIF2,UTHIF2]=Horseshoe_self(M2,Z2,TANBIV2,RC2,RV2,Hub_Flag,Rhub_oR2); 
[UAHIF2_1,UTHIF2_1]=Horseshoe_int(M2,M1,Z1,TANBIV1,RC2,RV1,Xf,Hub_Flag,... 
                                                                 Rhub_oR1); 
% ------------------------------------------------------------------------- 
figure; 
hold on 
G1(1)=0; 
G2(1)=0; 
while B_iter<ITER & any(B_res)==1  %(WHILE LOOP B1)   
    G_iter=1; 
    G1_res=1; 
    G2_res=1; 
    LT_res=1; 
    LQ_res=1; 

 
    rho=1025; 
    while G_iter<ITER &(G1_res>1e-5 | G2_res>1e-5 | LT_res>1e-5 | ... 
                                               LQ_res>1e-5)%(WHILE LOOP G1) 
    % Solve simultaneous equations for G1,G2,LT and LQ. 
    % Setting up the linear system of M1+M2+2 equations 
    % There are omissions in equations 2.61 and 2.62 which affect the  
    % coefficients of LQ and LT, A(:,M1+M2+1) and A(:,M1+M2+2), as well as  
    % the constant values in matrix B.      
    % First eq. 2.61 (Coney, p.42) 
        for i=1:M1 
            for m=1:M1 
                A(i,m)=(om1+q*LQ_last)*Z1*(UAHIF1(i,m)*RC1(i)*DR1(i)+... 
                                           UAHIF1(m,i)*RC1(m)*DR1(m))+... 
                              LT_last*Z1* (UTHIF1(i,m)*      DR1(i)+... 
                                           UTHIF1(m,i)*      DR1(m)); 
            end 
            for m=1:M2 
             A(i,m+M1)=(om1+q*LQ_last)*Z1*(UAHIF1_2(i,m)*RC1(i)*DR1(i))+... 
                       (om2-LQ_last)*  Z2*(UAHIF2_1(m,i)*RC2(m)*DR2(m))+... 
                        LT_last*       Z1* UTHIF1_2(i,m)*       DR1(i)+... 
                        LT_last*       Z2* UTHIF2_1(m,i)*       DR2(m); 
            end 
            A(i,M1+M2+1)=Z1*(VTC1(i)+om1*RC1(i))*DR1(i); 
            A(i,M1+M2+2)=Z1*q*VAC1(i)*RC1(i)*DR1(i); 
             
%    The circulation coefficients in the thrust constrain equation must be 
%    multiplied by (2*rho*Vs^2*pi*R^2) for dimensional consistency since Tr 
%    has dimensions [N] 
            A(M1+M2+1,i)=(2*rho*Vs^2*pi*R1^2)*... 
                    Z1*(VTC1(i)+om1*RC1(i)+UTSTAR1(i))*DR1(i);%thrust terms 
            A(M1+M2+2,i)=q*Z1*(VAC1(i)+UASTAR1(i))*RC1(i)*DR1(i);%torque t. 
            B(i)=-Z1*om1*VAC1(i)*RC1(i)*DR1(i); 
        end 
         
    % Then eq. 2.62 (Coney, p.43) 
        for i=1:M2 
            for m=1:M1 
               A(i+M1,m)=(om1+q*LQ_last)*Z1*UAHIF1_2(m,i)*RC1(m)*DR1(m)+... 
                         (om2-LQ_last)*Z2*  UAHIF2_1(i,m)*RC2(i)*DR2(i)+... 
                         LT_last*       (Z1*UTHIF1_2(m,i)*       DR1(m)+... 
                                         Z2*UTHIF2_1(i,m)*       DR2(i)); 
            end
            for m=1:M2 
              A(i+M1,m+M1)=(om2-LQ_last)*Z2*(UAHIF2(i,m)*RC2(i)*DR2(i)+... 
                                             UAHIF2(m,i)*RC2(m)*DR2(m))+... 
                           LT_last*Z2*      (UTHIF2(i,m)*       DR2(i)+... 
                                             UTHIF2(m,i)*       DR2(m)); 
            end 
            A(i+M1,M1+M2+1)=Z2*(VTC2(i)+om2*RC2(i))*DR2(i); 
            A(i+M1,M1+M2+2)=-Z2*VAC2(i)*RC2(i)*DR2(i);  
     
%    The circulation coefficients in the thrust constrain equation must be 
%    multiplied by (2*rho*Vs^2*pi*R^2) for dimensional consistency since Tr 
%    has dimensions [N] 
            A(M1+M2+1,i+M1)=(2*rho*Vs^2*pi*R2^2)*... 
                   Z2*(VTC2(i)+om2*RC2(i)+UTSTAR2(i))*DR2(i); %thrust terms 
 
            A(M1+M2+2,i+M1)=-Z2*(VAC2(i)+UASTAR2(i))*RC2(i)*DR2(i);%torque  
            B(M1+i)=-Z2*om2*VAC2(i)*RC2(i)*DR2(i); 
        end 
%         Modify terms related to circulation at the hub (innermost radial 
%         distance). The difference in the results is very small though. 
        if Hub_Flag==1 
            a=LT_last/(16*pi)*(log(1/Rhv)+3); 
            A(1,1)=A(1,1)-a*2*Z1^2; 
            A(1,M1+1)=A(1,M1+1)+a*2*Z1*Z2; 
            A(M1+1,1)=A(M1+1,1)+a*2*Z1*Z2; 
            A(M1+1,M1+1)=A(M1+1,M1+1)-a*2*Z2^2; 
        end 
        
        %---Compute total velocities used in viscous force calculations---- 
        VASTAR1=VAC1+UASTAR1;VASTAR2=VAC2+UASTAR2; 
        VTSTAR1=VTC1+om1*RC1+UTSTAR1;VTSTAR2=VTC2+om2*RC2+UTSTAR2; 
        VSTAR1=sqrt(VASTAR1.^2+VTSTAR1.^2); 
        VSTAR2=sqrt(VASTAR2.^2+VTSTAR2.^2); 

        Tv1=0; 
        Qv1=0; 
         
        for i=1:M1 
            Tv1=Tv1-(1/2)*Z1*VSTAR1(i)*(VAC1(i)+UASTAR1(i))... 
                                       *CoD1(i)*CD(i)*DR1(i); %viscous thrust 
                            
            Qv1=Qv1+(1/2)*Z1*VSTAR1(i)*(VTC1(i)+om1*RC1(i)+UTSTAR1(i))... 
                                       *RC1(i)*CoD1(i)*CD(i)*DR1(i);%v.torque 
        end 
%       or Tv1=Z1*sum(VSTAR1.*(VAC1+UASTAR1).*CoD1.*Cd.*DR1) 
%       and similarly for Qv1. 

        Tv2=0; 
        Qv2=0; 

        for i=1:M2 
            Tv2=Tv2-(1/2)*Z2*VSTAR2(i)*(VAC2(i)+UASTAR2(i))... 
                                       *CoD2(i)*CD(i)*DR2(i); %viscous thrust 
                            
            Qv2=Qv2+(1/2)*Z2*VSTAR2(i)*(VTC2(i)+om2*RC2(i)+UTSTAR2(i))... 
                                       *RC2(i)*CoD2(i)*CD(i)*DR2(i);%v. torque 
        end 
% The viscous thrust terms above must be myltiplied by 2*rho*Vs^2*R^2 in 
% order to represent dimensional values [N]since Tr is dimensional. 
        B(M1+M2+1)=Tr-(2*rho*Vs^2)*(R1^2*Tv1+R2^2*Tv2); 
% ========================================================================= 
% Account for hub drag term if a hub image is present 
% Fh is expressed in dimensional form in order to be consistent  
% with the dimensional value of Tr. 
       if Hub_Flag==1 
            Fh=rho/(16*pi)*(log(1/Rhv)+3)*(Z1*G1(1)*sqrt(0.5*rho*... 
                Vs^2*pi*R1^2)-Z2*G2(1)*sqrt(0.5*rho*Vs^2*pi*R2^2))^2; 
       elseif Hub_Flag==0 
            Fh=0; 
       end 
       B(M1+M2+1)=B(M1+M2+1)+Fh; 
% ========================================================================= 
% Divide viscous torque terms by pi since the non-dimensionalizing parame- 

 
% ters for Qi and Qv differ by the myltiplication parameter pi (Qi=1/pi*Qv) 
        B(M1+M2+2)=(Qv2-q*Qv1)/pi; 
        GL=linsolve(A,B); 
        G1=GL(1:M1); 
        G2=GL(M1+1:M1+M2); 
        LT=GL(M1+M2+1); 
        LQ=GL(M1+M2+2); 
  
% ----------Compute induced velocities------------------------------------- 
        [UA_SELF1,UT_SELF1,UA_INT1_2,UT_INT1_2]=Induced_Velocity(M1,M2,... 
                                    G1,G2,UAHIF1,UTHIF1,UAHIF1_2,UTHIF1_2); 
        UASTAR1=UA_SELF1+UA_INT1_2;  %total axial induced velocity/Vs 
        UTSTAR1=UT_SELF1+UT_INT1_2;  %total tangential induced velocity/Vs 
  
        [UA_SELF2,UT_SELF2,UA_INT2_1,UT_INT2_1]=Induced_Velocity(M2,M1,... 
                                    G2,G1,UAHIF2,UTHIF2,UAHIF2_1,UTHIF2_1); 
        UASTAR2=UA_SELF2+UA_INT2_1;  %total axial induced velocity/Vs 
        UTSTAR2=UT_SELF2+UT_INT2_1;  %total tangential induced velocity/Vs 
%--------------------------------------------------------------------------  
%     
%-----------Update velocities used in viscous force calculations---------- 
%        VAC1 = pchip(XR1,XVA1,RC1);     % axial inflow vel. / ship vel. at ctrl pts 
%        VTC1 = pchip(XR1,XVT1,RC1);     % tang. inflow vel. / ship vel. at ctrl pts
%        VAC2 = pchip(XR2,XVA2,RC2);     % axial inflow vel. / ship vel. at ctrl pts 
%        VTC2 = pchip(XR2,XVT2,RC2);     % tang. inflow vel. / ship vel. at ctrl pts
%        VASTAR1=VAC1+UASTAR1;VASTAR2=VAC2+UASTAR2; 
%        VTSTAR1=VTC1+om1*RC1+UTSTAR1;VTSTAR2=VTC2+om2*RC2+UTSTAR2; 
%        VSTAR1=sqrt(VASTAR1.^2+VTSTAR1.^2); 
%        VSTAR2=sqrt(VASTAR2.^2+VTSTAR2.^2);
%        
%       Added code for chord optimalization.
%       Epps: A Method for Propeller Blade Optimization and Cavitation 
%       Inception Mitigation 3.6 and 3.7
       CLmax1 = 0.5 + (0.2-0.5)*((R1*XR1-Rhub)/(R1-Rhub));
       CLmax2 = 0.5 + (0.2-0.5)*((R2*XR2-Rhub)/(R2-Rhub));
       
       Gamma1=G1*2*pi*R1*Vs;
       Gamma2=G2*2*pi*R2*Vs;
       XCoD1 = ((2*Gamma1')./(VSTAR1*Vs.*CLmax1*2*R1));
       XCoD2 = ((2*Gamma2')./(VSTAR2*Vs.*CLmax2*2*R2));

       CoD1 = pchip(XR1,XCoD1,RC1);
       CoD2 = pchip(XR1,XCoD2,RC1);

       HHcod = plot(RC1,-CoD1,'-',RC1,CoD1,'-',RC2,-CoD2,'*',RC2,CoD2,'*');

       for i=1:M1
           CLmax1(i) = 0.5 + (0.2-0.5)*((XR1(i)*R1-Rhub)/(R1-Rhub));
           CLmax2(i) = 0.5 + (0.2-0.5)*((XR2(i)*R2-Rhub)/(R2-Rhub));
           XCoD1(i) = (2*G1(i))/(Vs*CLmax1(i));
           XCoD2(i) = (2*G2(i))/(Vs*CLmax2(i));
   end
       XCoD1 = XCoD1/R1;
       XCoD2 = XCoD2/R2;
       CoD1 = pchip(XR1,XCoD1,RC1);
       CoD2 = pchip(XR1,XCoD2,RC1);
%-------------------------------------------------------------------------- 

% ------------Prepare for next iteration----------------------------------- 
        G_iter=G_iter+1; 
        G1_res=abs(G1-G1_last); 
        G2_res=abs(G2-G2_last); 
        LT_res=abs(LT-LT_last); 
        LQ_res=abs(LQ-LQ_last); 
        G1_last=G1; 
        G2_last=G2; 
        LT_last=LT; 
        LQ_last=LQ; 

        % WARNING IF LOOP G1 DOESN'T CONVERGE 
        % check for G1,G2 and LM convergence 
        if G_iter > ITER 
            warning('on'), 
            warning('WARNING: While loop G1 did NOT converge.'), 
            warning('off'), 
        end         
    end                                           %(END WHILE LOOP G1)  
  
% -------------Allign wake to new circulation distributions---------------- 
    [UAHIF1,UTHIF1,UAHIF2,UTHIF2,UAHIF1_2,UTHIF1_2,UAHIF2_1,UTHIF2_1,... 
     UASTAR1,UTSTAR1,UASTAR2,UTSTAR2,TANBIC1,TANBIV1,TANBIC2,TANBIV2] = ... 
     Align_wake(TANBIC1,TANBIV1,TANBIC2,TANBIV2,ITER,M1,M2,Z1,Z2,RC1,... 
     RV1,RC2,RV2,G1,G2,VAC1,VTC1,VAC2,VTC2,Js1,Js2,Xf,Hub_Flag,Rhub_oR1,... 
                                                                 Rhub_oR2);     
% -------------------------End of wake alignment--------------------------- 
  
    B_iter=B_iter+1; 
    B1_res=abs(TANBIC1-TANBIC1_last) 
    B2_res=abs(TANBIC2-TANBIC2_last) 
    B_res=[B1_res>1e-2 B2_res>1e-2]; %convergence limit for TANBI 
    TANBIC1_last=TANBIC1; 
    TANBIC2_last=TANBIC2; 
%     -------------Plot Circulation Distributions-------------- 
     
    plot(RC1,G1,RC2,G2,'r') 
 
%    ---------------------------------------------------------- 
    if B_iter > ITER 
        warning('on'), 
        warning('WARNING: While loop B1 did NOT converge.'), 
        warning('off'), 
    end 
end                                               %(END WHILE LOOP B1) 
grid on; 
hold off 
  
% ========== Plotting self-induced and interaction velocities ============= 
  
figure; 
subplot(2,1,1); 
plot(RC1,UA_SELF1,'-*',RC2,UA_SELF2,'-*r');grid on 
title('axial self-induced velocities') 
subplot(2,1,2) 
plot(RC1,UT_SELF1,'-*',RC2,UT_SELF2,'-*r');grid on 
title('tangential self-induced velocities') 
  
figure; 
subplot(2,1,1); 
plot(RC1,UA_INT1_2,'-*',RC2,UA_INT2_1,'-*r');grid on 
title('axial interaction velocities') 
subplot(2,1,2) 
plot(RC1,UT_INT1_2,'-*',RC2,UT_INT2_1,'-*r');grid on 
title('tangential interaction velocities') 
  
% ================= Induced Velocities Far Downstream ===================== 
[UAHIFinf_1,UTHIFinf_1]=Horseshoe_int(M2,M1,Z1,TANBIV1,RC2,RV1,20*Xf,... 
                                                        Hub_Flag,Rhub_oR1); 
[UAHIFinf_2,UTHIFinf_2]=Horseshoe_int(M1,M2,Z2,TANBIV2,RC1,RV2,19*Xf,... 
                                                        Hub_Flag,Rhub_oR2); 
  
[UA_SELFinf2,UT_SELFinf2,UA_INTinf_1,UT_INTinf_1]=Induced_Velocity(M2,... 
                             M1,G2,G1,UAHIF2,UTHIF2,UAHIFinf_1,UTHIFinf_1); 
[UA_SELFinf1,UT_SELFinf1,UA_INTinf_2,UT_INTinf_2]=Induced_Velocity(M1,... 
                             M2,G1,G2,UAHIF1,UTHIF1,UAHIFinf_2,UTHIFinf_2); 
figure; 
plot(RC2,UA_INTinf_1,RC2,UA_INTinf_2,'-r',RC2,... 
                                      UA_INTinf_1+UA_INTinf_2,'--');grid on 
legend('forward propeller','aft propeller','total CRP') 
title('Axial induced velocities far downstream'); 
figure; 
plot(RC2,UT_INTinf_1,RC2,-UT_INTinf_2,'-r',RC2,... 
                                      UT_INTinf_1-UT_INTinf_2,'--');grid on 
legend('forward propeller','aft propeller','total CRP') 
title('Tangential induced velocities far downstream'); 
% ========================= Forces Function================================ 
VMIV1 = 2*trapz(XR1,XR1.*XVA1)/(RoR^2-Rhub_oR1^2);% [ ], VMIV/ship velocity 
VMIV2 = 2*trapz(XR2,XR2.*XVA2)/(RoR^2-Rhub_oR2^2);% [ ], VMIV/ship velocity 
  
[CT1,CQ1,KT1,KQ1,CT2,CQ2,KT2,KQ2,EFFY,VSTAR1,VSTAR2] =... 
   Forces(CD,CD,DR1,DR2,VAC1,VAC2,TANBC1,TANBC2,... 
          UASTAR1,UASTAR2,UTSTAR1,UTSTAR2,CoD1,CoD2,G1,G2,M1,M2,RC1,RC2,... 
          Fh,Z1,Z2,Js1,Js2,VMIV1,VMIV2,N1,N2,Vs,R1,R2); 

Gamma1=G1*2*pi*R1*Vs; 
Gamma2=G2*2*pi*R2*Vs; 

Cl1= 2*Gamma1'./(VSTAR1*Vs.*XCoD1.*2*R1); 
Cl2= 2*Gamma2'./(VSTAR2*Vs.*XCoD2.*2*R2); 
       
% ========================================================================= 
%geometry_flag=1;     %flag for geometry generation 
if geometry_flag==1 
% -------------------------------------- Compute required lift coefficients 
Gamma1=G1*2*pi*R1*Vs; 
Gamma2=G2*2*pi*R2*Vs; 
%Cl1= 2*Gamma1'./(VSTAR1*Vs.*CoD1.*2*R1); 
%Cl2= 2*Gamma2'./(VSTAR2*Vs.*CoD2.*2*R2); 
% ============= Inputs necessary for geometry generation ================== 
skew01       = zeros(1,11);                          % Skew [deg] 
skew02       = zeros(1,11); 
rake01       = zeros(1,11);                          % Xs/D, Rake 
rake02       = zeros(1,11); 
% --------------------- t0/c, thickness / chord --------------------------- 
t0oc01       = [.2056 .1551 .1181 .0902 .0694 .0541 .0419 .0332 .0324... 
                .0204 .005];        
               %[0.0815 0.0771 0.0731 0.0664 0.0608 0.0561 0.0522... 
               %  0.0489 0.0457 0.0457 0.005]; 
t0oc02       = t0oc01;

XR1          = [0.2,0.3,0.4,0.5,0.6,...
                   0.7,0.8,0.9,0.95,0.99,1];
XR2          = [0.2,0.3,0.4,0.5,0.6,...
                   0.7,0.8,0.9,0.95,0.99,1];
% ------------------------------------------------------------------------- 
BetaI_c1=atand(TANBIC1); 
BetaI_c2=atand(TANBIC2); 
Np=40;            % Number of points over the chord 
% ======================= Generate Propeller Geometry ===================== 

[X3D,Y3D,Z3D,X3D_aft,Y3D_aft,Z3D_aft,] = Geometry(XR1,XR2,t0oc01,t0oc02,skew01,skew02,rake01,...
    rake02,RC1,RC2,Cl1,Cl2,BetaI_c1,BetaI_c2,Xf,Z1,Z2,Rhub,CoD1,CoD2,...
    R1,R2,M1,M2,Np); 
                    
% ============ Increase the number of sections and the respective values                    
% (t0oc,f0oc,AlphaI,Cl,Z3D,Vstar,RC,c) such that the computation of the  
% cavitating area is more accurate ======================================== 
% Cosine spacing is used with the end values remaining the same. The new 
% number of sections are given by M1_int, M2_int. 
M1_int=40;M2_int=40; 
DEL1=(RC1(end)-RC1(1))/2; 
DEL2=(RC2(end)-RC2(1))/2; 
for n=1:M1_int 
    RC1_int(n)=RC1(1)+DEL1*(1-cos(n*pi/M1_int)); 
end 
  
for n=1:M2_int 
    RC2_int(n)=RC2(1)+DEL2*(1-cos(n*pi/M2_int)); 
end 
% ======= Now interpolate to find new values at RC1,2_int locations ======= 
Gamma1_int=pchip(RC1,Gamma1,RC1_int); 
Gamma2_int=pchip(RC2,Gamma2,RC2_int); 
VSTAR1_int=pchip(RC1,VSTAR1*Vs,RC1_int);   %dimensional velocity! 
VSTAR2_int=pchip(RC2,VSTAR2*Vs,RC2_int);   %dimensional velocity! 
UASTAR1_int=pchip(RC1,UASTAR1,RC1_int); 
UASTAR2_int=pchip(RC2,UASTAR2,RC2_int); 
CoD1_int=pchip(RC1,CoD1,RC1_int); 
CoD2_int=pchip(RC2,CoD2,RC2_int); 
Cl1_int= 2*Gamma1_int./(VSTAR1_int.*CoD1_int.*2*R1); %new Cl1 
Cl2_int= 2*Gamma2_int./(VSTAR2_int.*CoD2_int.*2*R2); %new Cl2 
BetaI_c1_int=pchip(RC1,BetaI_c1,RC1_int); 
BetaI_c2_int=pchip(RC2,BetaI_c2,RC2_int); 

 
  
% ------------------ Alternatively could extrapolate ----------------------  
% RC1_int=0.9*Rhub_oR1+(1-0.9*Rhub_oR1)*(sin((0:40)*pi/(2*40))); 
% RC2_int=0.9*Rhub_oR1+(1-0.9*Rhub_oR1)*(sin((0:40)*pi/(2*40))); 
% Gamma1_int=interp1(RC1,Gamma1,RC1_int,'pchip','extrap'); 
%... 
% ------------------------------------------------------------------------- 
  
% ===== Run Geometry module again ====== 
[f0oc1,f0oc2,t0oc1,t0oc2,AlphaI1,AlphaI2,X3D,Y3D,Z3D,X3D_aft,Y3D_aft,... 
Z3D_aft,c1,c2,x0_1,x0_2,theta_Z1,theta_Z2] = Geometry(XR1,XR2,t0oc01,... 
t0oc02,skew01,skew02,rake01,rake02,RC1_int,RC2_int,Cl1_int,Cl2_int,... 
BetaI_c1_int,BetaI_c2_int,Xf,Z1,Z2,Rhub,CoD1_int,CoD2_int,R1,R2,M1_int,... 
                                                                M2_int,Np); 
  
% =============== Run Cavitation module for both propellers =============== 
Cp_mode='VLM'; 
H=2;           % shaft centerline depth 
[Color_matrix_upper,Color_matrix_lower,cav_mess1]=Cavitation(Cp_mode,... 
M1_int,t0oc1,f0oc1,AlphaI1,Cl1_int,H,Z3D,VSTAR1_int,Np,Z1,x0_1,RC1_int,... 
                                  c1,R1,theta_Z1,BetaI_c1_int,UASTAR1_int); 
       
[Color_matrix_upper_aft,Color_matrix_lower_aft,cav_mess2]=... 
    Cavitation(Cp_mode,M2_int,t0oc2,f0oc2,AlphaI2,Cl2_int,H,Z3D_aft,... 
    VSTAR2_int,Np,Z2,x0_2,RC2_int,c2,R2,theta_Z2,BetaI_c2_int,UASTAR2_int);  
        
% ============== Plot Propeller Cavitation Image ====================        
figure; 
grid on;   
axis equal; 
axis([-2*Xf*R1 R1 -1.1*R1 1.1*R1 -1.1*R1 1.1*R1]); 
xlabel('X (3D) [m]','FontSize',12);  
ylabel('Y (3D) [m]','FontSize',12);  
zlabel('Z (3D) [m]','FontSize',12);  
title(['3D Cavitation Image using ',Cp_mode],'FontSize',16);  
hold on 
% ================ Plot forward propeller blade surfaces ================== 
for k=1:Z1 
     
    surf(X3D(:,1:Np,1),Y3D(:,1:Np,k),Z3D(:,1:Np,k),... 
        Color_matrix_upper(:,:,k)); 
end 
  
for k=1:Z1 
     
    surf(X3D(:,2*Np:-1:Np+1,1),Y3D(:,2*Np:-1:Np+1,k),... 
        Z3D(:,2*Np:-1:Np+1,k),Color_matrix_lower(:,:,k)); 
end 
  
% ================== Plot aft propeller blade surfaces ==================== 
for k=1:Z2 
     
    surf(X3D_aft(:,1:Np,1),Y3D_aft(:,1:Np,k),Z3D_aft(:,1:Np,k),... 
        Color_matrix_upper_aft(:,:,k)); 
end 


  
for k=1:Z2 
     
    surf(X3D_aft(:,2*Np:-1:Np+1,1),Y3D_aft(:,2*Np:-1:Np+1,k),... 
        Z3D_aft(:,2*Np:-1:Np+1,k),Color_matrix_lower_aft(:,:,k)); 
end 
shading interp; 
  
% Plot the hub using only one color 
% ================================= 
hub_clr=mean(caxis); 
tick = 90:-15:0; 
[yh0,zh0,xh0] = cylinder(Rhub*sind(tick),50);    
xh0 = -0.5*c1(1)*xh0 - 0.75*R1; 
surf(xh0,yh0,zh0,hub_clr*ones(7,51)); 
     
[yh1,zh1,xh1] = cylinder(Rhub,50); 
xh1 = 6*c1(1)*xh1 - 0.75*R1; 
surf(xh1,yh1,zh1,hub_clr*ones(2,51));  
colorbar 
text(0,0,R1,cav_mess1,'FontSize',15,'HorizontalAlignment','center') 
text(0,0,-R1,'-C_p','FontSize',15,'HorizontalAlignment','center') 
text(-Xf*R1,0,R2,cav_mess2,'FontSize',15,'HorizontalAlignment','center') 
else
    Np = 0;
    X3D = 0;
    Y3D = 0;
    Z3D = 0;
end 
