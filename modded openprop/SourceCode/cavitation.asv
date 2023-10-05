%Written by Dimitrios Laskos
% Function Cavitation calculates pressure coefficients on blades' surfaces
% and assigns colors depending on whether the values exceed cavitation
% numbers (indicating cavitation inception) or not.
% ===============Inputs====================================================
% Cp_mode 'VLM' or 'XFOIL' depending on which method is
% implemented for calculating pressure coefficients
% Mp        [ ],    Number of points over the span
% t0oc      [ ],    Maximum thickness / chord at each radius
% f0oc      [ ],    Maximum camber ratio (f0/c=0.0679*Cl)
% AlphaI    [deg],  Ideal angle of attack (AlphaI=1.54*Cl)
% Cl        [ ],    Lift coefficient (Cl=Clideal)
% H         [m],    Shaft centerline depth
% Z3D       [m],    vertical location
% Vstar     [ ],    Total inflow velocity
% Np        [ ],    Number of points over the chord
% Z         [ ],    Blade number
% x0        [ ],    chordwise location [0:1]
% RC        [ ],    Non-dimensional radius for control points
% c         [m],    Chordlengths of blade sections along span
% R         [m],    Propeller radius
% theta_Z   [deg],  angle between blades
% BetaI_c   [deg],  Hydrodynamic pitch angle
% UASTAR    [ ],    Total axial induced velocity (self- and interaction-)
% =========================================================================
function[Color_matrix_upper,Color_matrix_lower,cav_mess]=...
 cavitation(Cp_mode,Mp,t0oc,f0oc,AlphaI,Cl,H,Z3D,Vstar,...
 Np,Z,x0,RC,c,R,theta_Z,BetaI_c,UASTAR)

% =========================== Execution of XFOIL ==========================
if strcmp(Cp_mode,'XFOIL')
 xdir='.\Xfoil\';
 foil_type='LOAD'; % or 'NACA'
 foil_name='foildata';
 for i=1:Mp % for each section along the span
 makefoil(t0oc(i),f0oc(i),'NACAa=08.txt','65A010.txt',foil_name);
 cmd=[xdir,'xfoil.exe',' ',foil_type,...
 ' ',foil_name,' NORM ',' GDES TSET ',num2str(t0oc(i)),' ',...
 num2str(f0oc(i)),' ','GDES EXEC ',' PANE',' OPER ALFA ',...
 num2str(AlphaI(i)),' OPER CPWR ',' ','CParray'];
 system(cmd)


 fid=fopen('CParray');
 datain=textscan(fid,'%f64 %f64','headerlines',1);
 fclose(fid);
 Length=length(datain{1,2})
 if length(datain{1,2})>=160
 cpi{1,i}=datain{1,2};
 xcpi{1,i}=datain{1,1};
 else
 cmd=[xdir,'xfoil.exe',' ',foil_type,...
 ' ',foil_name,'NORM ',' GDES TSET ',num2str(t0oc(i)),' ',...
 num2str(f0oc(i)),' ','GDES EXEC ',' OPER ALFA ',...
 num2str(AlphaI(i)),' OPER CPWR ',' ','CParray'];
 system(cmd)

 end

 fid=fopen('CParray');
 datain=textscan(fid,'%f64 %f64','headerlines',1);
 fclose(fid);
 cpi{1,i}=datain{1,2};
 xcpi{1,i}=datain{1,1};

% Remove double values from xcpi arrays and keep only those appearing
% first such that the interpolation routine doesn't crash
 n=length(xcpi{1,i})-1;
 Bpos=[];
 counter=0;
 for l=1:n
 if xcpi{1,i}(l)==xcpi{1,i}(l+1)
 counter=counter+1;
 Bpos(counter)=l+1;
 end
 end

 if isempty(Bpos)~=1
 counter1=0;
 ind_matrix=[];
 for l=1:n+1

     if (l~=Bpos)==1
 counter1=counter1+1;
 ind_matrix(counter1)=l;
 end
 end
 xcpi{1,i}=xcpi{1,i}(ind_matrix);
 cpi{1,i}=cpi{1,i}(ind_matrix);
 end

% plot Cp distribution for each section
% figure;grid on;
% plot(xcpi{1,i},-(cpi{1,i}));
% title({['section # ',num2str(i)]});

 end
 for i=1:Mp;
 for j=1:length(xcpi{1,i})-1
 xcpi_compare{1,i}(j,1)=xcpi{1,i}(j)-xcpi{1,i}(j+1);
 end
 end
 % Indexing begins from TE (xcpi=1), goes to LE (xcpi=0)along upper side
 % and returns to TE (xcpi=1) again along the lower foil side
 for i=1:Mp
 ind_upper{1,i}=find(xcpi_compare{1,i}>=0);
 ind_upper{1,i}=[ind_upper{1,i};ind_upper{1,i}(end)+1];
 ind_lower{1,i}=find(xcpi_compare{1,i}<0);
 ind_lower{1,i}=ind_lower{1,i}+ones(length(ind_lower{1,i}),1);
 xcpi_upper{1,i}=xcpi{1,i}(ind_upper{1,i});
 xcpi_lower{1,i}=xcpi{1,i}(ind_lower{1,i});
 cpi_upper{1,i}=cpi{1,i}(ind_upper{1,i});
 cpi_lower{1,i}=cpi{1,i}(ind_lower{1,i});

 end
 % Interpolate to find Cp values at Np positions along the chord
 for i=1:Mp
 Cpi_upper(i,:)=pchip(xcpi_upper{1,i},cpi_upper{1,i},x0);
 Cpi_lower(i,:)=pchip(xcpi_lower{1,i},cpi_lower{1,i},x0);
 end
%======================== End of XFOIL Execution ==========================
elseif strcmp(Cp_mode,'VLM')
 unsteady_flag=0;
 % ===========Cp calculation using VLM code ============================
 if unsteady_flag==0
 for i=1:Mp
 [xt, CPU(i,:), CPL(i,:)]=VLMcav(40, Cl(i),0,t0oc(i));
 Cpi_upper(i,:)=pchip(xt,-CPU(i,:),x0);
 Cpi_lower(i,:)=pchip(xt,-CPL(i,:),x0);
 end
 end
 % ================ Unsteady Cavitation Calculation ====================
 if unsteady_flag==1
 load wake_030910 wake_full
 theta=[0:5:360];
 roR_wake=[0.2:0.05:1];
 [THETA,ROR_WAKE]=meshgrid(theta,roR_wake);
 [THETA_Z,R_C]=meshgrid(theta_Z(1:end-1),RC);

 wake_int=interp2(THETA,ROR_WAKE,wake_full',THETA_Z,R_C);
 VAC_wake=wake_int';
% ========= calculate VTSTAR ==============================================
 Vs=5; % change manually (could be added as a function input)
 VTSTAR=(Vstar/Vs).*cosd(BetaI_c);
 beta_wake=atand((ones(Z,1)*UASTAR+VAC_wake)./(ones(Z,1)*VTSTAR));
% ===== Now calculate CP using VLM ========================================
for k=1:Z
 for i=1:Mp
 delta_alpha(k,i)=beta_wake(k,i)-BetaI_c(i);% (a_ideal-alpha)
 end
end
 for k=1:Z
 for i=1:Mp
 [xt, CPU(i,:,k), CPL(i,:,k)]=VLMcav(40, Cl(i),...
 -delta_alpha(k,i),t0oc(i));
 Cpi_upper(i,:,k)=pchip(xt,-CPU(i,:,k),x0);
 Cpi_lower(i,:,k)=pchip(xt,-CPL(i,:,k),x0);
 end
 end
 end
end
% ================End of Cp Calculation===================
%Accurate calculation of sigma for all blades
% variation of Z3D along section is taken into account
rho=1025;
for k=1:Z
 for i=1:Mp
 for j=1:Np
 SIGMA2(i,j,k)=(101000+rho*9.81*(H-Z3D(i,j,k))-2500)./...
 (rho*Vstar(i)^2/2); % cavitation matrix
 end
 end
end
% Check for cavitation on suction side
% ====================================
Cpi_upper_cmp=zeros(Mp,Np,Z);
if unsteady_flag==1
 for k=1:Z
 Cpi_upper_cmp(:,:,k)=Cpi_upper(:,:,k);
 end
elseif unsteady_flag==0
 for k=1:Z
 Cpi_upper_cmp(:,:,k)=Cpi_upper;
 end
end
Cav_matrix_upper=+(SIGMA2<-Cpi_upper_cmp);
% Check for cavitation on pressure side
% ====================================
Cpi_lower_cmp=zeros(Mp,Np,Z);
if unsteady_flag==1
 for k=1:Z
 Cpi_lower_cmp(:,:,k)=Cpi_lower(:,:,k);
 end
elseif unsteady_flag==0
 for k=1:Z
 Cpi_lower_cmp(:,:,k)=Cpi_lower;
 end
end
Cav_matrix_lower=+(SIGMA2<-Cpi_lower_cmp);
% figure;
% grid on;
% axis equal;
% axis([-R/2 R -1.1*R 1.1*R -1.1*R 1.1*R]);
% xlabel('X (3D) [m]','FontSize',12);
% ylabel('Y (3D) [m]','FontSize',12);
% zlabel('Z (3D) [m]','FontSize',12);
% title(['3D Cavitation Image using ',Cav_module],'FontSize',16);
% hold on
% Suction sides
% =====================
for k=1:Z
 for i=1:Mp
 for j=1:Np
 if Cav_matrix_upper(i,j,k)==1
 Color_matrix_upper(i,j,k)=2*Cav_matrix_upper(i,j,k);
 else
 Color_matrix_upper(i,j,k)=Cpi_upper_cmp(i,j,k);
 end
 end
 end
end
% Pressure sides
% =====================
for k=1:Z
 for i=1:Mp
 for j=1:Np
 if Cav_matrix_lower(i,j,k)==1
 Color_matrix_lower(i,j,k)=2*Cav_matrix_lower(i,j,k);
 else
 Color_matrix_lower(i,j,k)=Cpi_lower_cmp(i,j,k);
 end
 end
 end
end
% Print message
% =============
B_lower=+any(any(any(Cav_matrix_lower)));
B_upper=+any(any(any(Cav_matrix_upper)));
% B_cmp=+(B_lower==1)|(B_upper==1);
if (B_lower==1)|(B_upper==1)==1
 message='Cavitation present';
else
 message='No cavitation present'
end
% text(0,0,R,message,'FontSize',15);
Area=0;
for i=1:Mp-1
 RCdif(i)=RC(i+1)-RC(i);
 Rdif(i)=RCdif(i)*R;
 Area=Area+(c(i)+c(i+1))*Rdif(i)/2;
end
% ====Cavitating area=========
% All trapezoids forming the blade area of a specific section
% ((Mp-1) sections in total), have the same area since the length
% of the bases is the same and their height (Rdif) is common
cav_area_lower=0;
cav_area_upper=0;
for i=1:Mp %or from 2:Mp
 num_lower(i)=length(find(Cav_matrix_lower(i,:,1)));
 num_upper(i)=length(find(Cav_matrix_upper(i,:,1)));
 if num_lower(i)~=0
 cav_area_lower=cav_area_lower+num_lower(i)*((c(i)+c(i-1))/...
 (Np-1))*Rdif(i-1)/2;
 end
 if num_upper(i)~=0
 cav_area_upper=cav_area_upper+num_upper(i)*((c(i)+c(i-1))/...
 (Np-1))*Rdif(i-1)/2;
 end
end
perc_lower=100*cav_area_lower/Area;
perc_upper=100*cav_area_upper/Area;
message1=strcat('TDC face cavitation:',num2str(perc_lower),'%');
message2=strcat('TDC back cavitation:',num2str(perc_upper),'%');
cav_mess={message1;message2};