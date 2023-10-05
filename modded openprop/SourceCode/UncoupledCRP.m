% === === Written by Dimitrios Laskos ===
% Contra-Rotating Propeller Design Code based on lifting line theory
% Numerical Implementation of the iterative process for CRP Design
% by coupling Single Propellers Design codes. The Variational Optimization
% Method developed by Kerwin, et al. (1986) is used in order to
% determine optimum circulation distributions for the CRP set.
% - References
% 1) E.B. Caster & T.A. Lafone,"A Computer Program for the Preliminary
% Design of Contrarotating Propellers",DTNSRDC Report SPD-596-01,
% 1975.
% 2) J. Kerwin, W. Coney & C. Hsin, "Optimum Circulation Distributions for
% Single and Multi-Component Propulsors", 21st American Towing
% Tank Conference (ATTC), 1986.
% 3) W. Coney, "A Method for the Design of a Class of Optimum Marine
%Propulsors", Ph.D. thesis, MIT, 1989.
%4) B.D. Cox & A.M. Reed, "Contrarotating Propellers-Design Theory and
%Application", Propellers '88 Symposium, 1988.
% 5) J. Kerwin, "Hydrofoils and Propellers", MIT Course 2.23 notes, 2007.
% 6) J.W. Wrench, "The Calculation of Propeller Induction Factors", David
% Taylor Model Basin (Technical Report 1116), Feb. 1957.
% 7) G. Hough & D. Ordway, "Generalized Actuator Disk", Developements in
% Theoretical and Applied Mechanics, Vol.2, pp.23-31, 1965.
% 8) B. Epps et al., "OpenProp: An open-source parametric design tool for
% propellers", Grand Challenges in Modeling & Simulation
% Conference (GCMC '09), 2009.
% - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% The function 'Coney' is a modified version of the one used in MIT
% OpenProp v2 Version 1.0. For the calculation of the axial and tangential
% interaction velocities a modified version of the 'CMV' function is used.
% indices 1,2 refer to the forward and the aft propellers respectively
% Input Variables:
% q [, torque ratio Q2/Ql
%Js [ ], Advance coefficient (same for both propellers)
% Rhub [m], Hub radius (common for both propellers)
% R [m], Propeller radius
% Zl,Z2 [ ], Number of blades
% Mp [ ], Number of vortex panels over the radius
% CTPDES [ ], Design thrust loading
% Hub Flag Inclusion of hub effects (l=YES, 0=NO)
% XR ], Radial locations for defining inflow velocities
% and geometric properties
% XVA [ ], Va/Vs, axial inflow vel. / ship vel.
% XVT [ ], Vt/Vs, tangential inflow vel. / ship vel.
% XCoD [ ], chord / propeller diameter
% XCD [, section drag coefficient
% spacing Type of radial spacing ('cosine' or 'constant')
.................... ............ . ...... .... ...........
% Xs [m], Axial separation between propellers
% Output variables:
% n total [ ], total efficiency
% RC,RC2 [, Control Points for forward and aft propellers
% G1,G2 [ ], Non-dimensional circulation
% --- -------------------------------------
function [RC,RC2,G1,G2,n_total]=CRPUncoupled(q,Js,Rhub,R,Zl,Z2,Mp,...
CTPDES,HubFlag,XR,XVA,XVT,XCoD,XCD,spacing,Xs)
% = New Inputs for Variational Optimization =
    Rhv = 0.042;
    SCF =1;
    %Xf = Xs/(R);
    ITER=10;
    q_iter = 1;
    q_res=1;
    q_last2=0;
    q_lastl=0;
    CTPDES1MF_last2=0;
    CTPDES1MF_last1=0;
% Application of Newton method for finding the specific thrust ratio which
% yields the required torque ratio for a given thrust loading
while q_iter<ITER & q_res>1e-5 %(WHILE LOOP A)
    if q_iter==1
        CTPDES1MF=1;
    elseif q_iter==2
        CTPDES1MF=1+(q-(Kq2/Kql))/(5*q);
    elseif q_iter>2
        CTPDES1MF=CTPDES1MF_lastl+(CTPDES1MF_lastl-CTPDES1MF_last2)* ...
        (q-qlast1)/(q_last1-q_last2);
    end
CTPDES1=CTPDES1MF*(CTPDES/2);
CTPDES2=CTPDES-CTPDES1; %thrust coefficient required by aft propeller
iter flag=l;
if iter_flag==1 % (IF CONDITION B)
% ====== = = == = = == = = == = = == = = == = = == = = == = = == = = =
% iterative procedure for determining circulation distributions
% for the forward and the aft propellers of the CRP set.
G1_last=O;
G2_last=O;
G_iter=1;
G1_res=1;
G2_res=1;
while G_iter<ITER & (G1_res>1e-5 | G2_res>1e-5) % (WHILE LOOP B)
%solve for Gl,G2 and update respective onset flows
% Variational Optimization for forward prop
if G_iter==1
    [RV,G1,TANBIV,TANBIC,VAC,VTC,UASTAR,UTSTAR,RC,CD,CoD] ...
        = Coney(Rhub,R,Z1,Mp,ITER,Rhv,SCF,Js,CTPDES1,Hub_Flag,...
        XR,XCoD,XCD,XVA,XVT,spacing,'normal',O);
% Xs
elseif G_iter~=1
[RV,G1,TANBIV,TANBIC,VAC1,VTC1,UASTAR,...
UTSTAR,RC,CD,CoD,Kql,Ktl,CT1,CP1]...
= Coney(Rhub,R,Z1,Mp,ITER,Rhv,SCF,Js,CTPDES1,Hub_Flag,...
RC,CoD,CD,VA1,VT1,spacing,'none',RV);
end
% Calculate interaction velocities at aft propeller plane
Vinter2=zeros(2,length(G1));
for i=1:length(Gl)
[Vinter2(:,i)] = CMV(Xf,RC(i),RV,G1,TANBIV,Zl);
end
UA2_INT=Vinter2(1,:);
UT2_INT=Vinter2(2,:);
VA2=VAC-UA2_INT;
VT2=VTC-UT2_INT;
% Variational Optimization for aft prop
[RV2,G2,TANBIV2,TANBIC2,VAC2,VTC2,UASTAR2,UTSTAR2,RC2...
CD,CoD,Kq2,Kt2,CT2,CP2] ...
= Coney(Rhub,R,Z2,Mp,ITER,Rhv,SCF,Js,CTPDES2,HubFlag,...
RC,CoD,CD,VA2,VT2,spacing, 'none',RV);
% calculate interaction velocities at forward propeller plane
Vinterl=zeros(2,length(G2));
for i=1:length(G2)
[Vinterl(:,i)] = CMV(-Xf,RC(i),RV2,G2,TANBIV2,Z2);
end
UA1_INT=Vinterl(1,:);
UT1_INT=Vinterl(2,:);
VA1=VAC-UA1_INT;
VT1=VTC-UT1INT;
G_iter=G_iter+l
G1_res=abs(Gl-Gl_last);
G2_res=abs(G2-G2_last);
G1_last=G1;
G2_last=G2;
end %(END OF WHILE LOOP B)
elseif iter_flag~=1
% Variational Opt for forward prop
[RV,G,TANBIV,TANBIC,VAC,VTC,UASTAR,UTSTAR,RC,CD,CoD] ...
= Coney(Rhub,R,Z1,Mp,ITER,Rhv,SCF,Js,CTPDES1,HubFlag,...
XR,XCoD,XCD,XVA,XVT,spacing,'normal',0);
% First calculate interaction velocities at aft propeller plane
Vinter2=zeros(2,length(G));
for i=1:length(G)
[Vinter2(:,i)] = CMV(Xf,RC(i),RV,G,TANBIV,Zl);
end
VA2=VAC-Vinter2(1,:);
VT2=VTC-Vinter2(2,:);
....... I, , . j& . .. .. .. ..... ...
CTPDES2=CTPDES-CTPDES1; %thrust coefficient required by aft propeller
[RV2,G2,TANBIV2,TANBIC2,VAC2,VTC2,UASTAR2,UTSTAR2,RC2] ...
= Coney(Rhub,R,Z2,Mp,ITER,Rhv,SCF,Js,CTPDES2,Hub_Flag,...
RC,CoD,CD,VA2,VT2,spacing,'none',RV);
% calculate interaction velocities at forward propeller plane
Vinter1=zeros(2,length(G2));
for i=1:length(G2)
[Vinterl(:,i)] = CMV(-Xf,RC(i),RV2,G2,TANBIV2,Z2);
end
VA1=VAC-Vinterl(1,:);
VT1=VTC-Vinterl(2,:);
% again run Coney.m for forward prop
[RVnew,G1_new,TANBIVnew,TANBICnew,VAC1,VTC1,UASTAR,UTSTAR,RCnew,...
CD1,CoDl,Kql,Ktl,CT1,CP1,EFFY1,VMIV1]=Coney(Rhub,R,Z1,Mp,ITER,Rhv,...
SCF,Js,CTPDES1,HubFlag,RC,CoD,CD,VA1,VT1,spacing,'none',RV);
% calculate new interaction velocities at aft propeller plane
Vinter2b=zeros(2,length(G));
for i=1:length(G)
[Vinter2b(:,i)] = CMV(Xf,RC2(i),RVnew,Glnew,TANBIVnew,Zl);
end
VA2=VAC-Vinter2b(1,:);
VT2=VTC-Vinter2b(2,:);
% run Coney.m for aft propeller
[RV2 new,G2, new,TANBIV2_new,TANBIC2_new,VAC2,VTC2,UASTAR2,UTSTAR2,...
RC2 new,CD2,CoD2,Kq2,Kt2,CT2,CP2,EFFY2,VMIV2]=Coney(Rhub,R,Z2,Mp,...
ITER,Rhv,SCF,Js,CTPDES2,HubFlag,RC2,CoD,CD,VA2,VT2,spacing,'none',RV2);
end %(END OF IF CONDITION B)
q_iter=q_iter+l;
q_res=abs((Kq2/Kql)-q);
q_last2=qlast1;
q_lastl=Kq2/Kql;
CTPDES1MF_last2=CTPDES1MF_last1;
CTPDESlMF_lastl=CTPDESlMF;
end % (END OF WHILE LOOP A)
% --------------------- Compute total efficiency --------------------------
% This expression however applies only to the case for which the diameters,
% and the speeds of the two propellers are equal. The effect of the hub
% drag is also neglected
n_total=(CT1+CT2)/(CP1+CP2);
figure;
subplot (3,1,1)
plot(RC,Gl,'-*',RC2,G2,'-*r');grid on;xlabel('Control points radii (RC)');
ylabel('Non-dimensional Circulation (G)');subplot(3,1,2)
plot(RC,atand(TANBIC),RC2,atand(TANBIC2),'-r');
xlabel('Control points radii (RC)');
ylabel('Hydrodynamic Pitch Angle \betai (deg)');grid on
................ .............. M . ... .... -
% Plotting total inflow velocities
VASTAR= VAC + UASTAR; % total axial inflow vel. / ship vel.
VTSTAR= VTC + pi*RC/Js + UTSTAR; % total tangential inflow vel. / ship vel.
VSTAR = sqrt(VTSTAR.^2 + VASTAR.^2);%magnitude of the inflow vel./ship vel.
VASTAR2= VAC2 + UASTAR2; %total axial inflow vel./ship vel.
VTSTAR2= VTC2 + pi*RC2/Js + UTSTAR2;%total tangential inflow vel./ship vel.
VSTAR2 =sqrt(VTSTAR2.^2+VASTAR2.^2);%magnitude of the inflow vel./ship vel.
subplot(3,1,3);plot(RC,VSTAR,'-*',RC2,VSTAR2, '-*r');grid on;
xlabel('Control points radii (RC)')
ylabel('total inflow velocity (V s ta r/Vs)')
% = = = Induced Velocities far downstream
M1=Mp;M2=Mp;Z1=Z;Z2=Z;
[UAHIFinf_1,UTHIFinf_l]=Horseshoe_int(M2,Ml,Z1,TANBIV,RC2,RV,20*Xf,...
HubFlag,Rhub_oR);
[UAHIFinf_2,UTHIFinf_2]=Horseshoe_int(Ml,M2,Z2,TANBIV2,RC,RV2,19*Xf,...
HubFlag,Rhub_oR);
[UAINTinf_1,UTINTinf_1]=Induced_Velocity_int(M2,M1,G1,UAHIFinf_1,...
UTHIFinf_1);
[UAINTinf_2,UTINTinf_2]=Induced_Velocity_int(M1,M2,G2,UAHIFinf_2,...
UTHIFinf_2);
% Induced velocities far downstream for Single propeller having double the
% the same number of blades as the CRP set
Z_SR=Z1+Z2;
[RVs,Gs,TANBIVs,TANBICs,VACs,VTCs,UASTARs,UTSTARs,RCs,CDs,CoDs]= ...
Coney(Rhub,R,Z_SR,Mp,ITER,Rhv,SCF,Js,CTPDES,HubFlag,XR,XCoD,XCD,XVA,...
XVT,spacing,'normal',O);
[UAHIFinf_SR,UTHIFinf_SR]=Horseshoe_int(Mp,Mp,Z_SR,TANBIVs,RCs,RVs,...
20*Xf,HubFlag,Rhub_oR);
[UAINTinf_SR,UTINTinf_SR]=Induced_Velocity_int(Mp,Mp,Gs,UAHIFinfSR,...
UTHIFinf_SR);
% -- -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
figure;
plot(RC2,UA_INTinf_1,RC2,UA_INTinf_2,'-r',...
RC2,UAINTinf_1+UAINTinf_2,'--g*',RCs,UAINTinfSR, '--b*');grid on
legend('forward propeller','aft propeller','total CRP','single propeller')
title('Axial induced velocities far downstream / Uncoupled method');
figure;
plot(RC2,UTINTinf_1,RC2,-UT_INTinf_2,'-r',...
RC2,UTINTinf_1-UTINTinf_2,'--g*',RCs,UTINTinf_SR, '--b*');grid on
legend('forward propeller','aft propeller','total CRP','single propeller')
title('Tangential induced velocities far downstream / Uncoupled method');
% Plot interaction and self induced velocities
figure;
subplot (2,1,1);
plot(UASTAR,RC,'-*',Xs+UASTAR2,RC2,'-*r');grid on
hold on
plot(zeros(1,length(RC)),RC,'--','Linewidth',2)
plot(Xs*ones(1,length(RC2)),RC2,'--r','Linewidth',2)
hold off
title('axial self-induced velocities')
subplot (2,1,2)
plot(UTSTAR,RC,'-*',Xs+UTSTAR2,RC2,'-*r');grid on
........... O W .......... ...... . ....
hold on
plot(zeros(1,length(RC)),RC,'--','Linewidth',2)
plot(Xs*ones(1,length(RC2)),RC2,'--r','Linewidth',2)
title('tangential self-induced velocities')
% --- Interaction Velocities-------------------
figure;
subplot(2,1,1);
plot(-UA1_INT,RC,'-*',Xs-UA2_INT,RC2,'-*r');grid on
hold on
plot(zeros(1,length(RC)),RC,'--','Linewidth',2)
plot(Xs*ones(1,length(RC2)),RC2,'--r','Linewidth',2)
title('axial interaction velocities')
subplot(2,1,2)
plot(-UT1_INT,RC,'-*',Xs-UT2_INT,RC2,'-*r');grid on
hold on
plot(zeros(1,length(RC)),RC,'--','Linewidth',2)
plot(Xs*ones(1,length(RC2)),RC2,'--r','Linewidth',2)
title('tangential interaction velocities')
end