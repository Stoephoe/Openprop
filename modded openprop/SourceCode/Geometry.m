% ========================================================================= 
% =================================== Determine Propeller Geometry Function      
%  
% This function determines the geometry of the CRP set. It outputs the  
% geometry as a 3D image. 
% 
% Reference: J.S. Carlton, "Marine Propellers & Propulsion", ch. 3, 1994.  
% 
% ------------------------------------------------------------------------- 
% Input Variables: 
% 
%   filename            file name prefix for all output files 
%   Date_string         time and date to print on reports 
%   Make2Dplot_flag     flag for whether to make 2D geometry plot 
%   Make3Dplot_flag     flag for whether to make 3D geometry plot 
%   Make_Rhino_flag     flag for whetehr to make a Rhino output file 
%   Meanline            flag for choice of meanline  form 
%   Thickness           flag for choice of thickness form 
% 
%   XR          [ ],    input radii / propeller radius 
%   t0oc0       [ ],    input thickness / chord at each radius 
%   skew0       [deg],  input skew              at each radius 
%   rake0       [ ],    input rake / diameter   at each radius 
% 
%   RC          [ ],    control point radii / propeller radius 
%   Cl          [ ],    section lift coefficients 
%   BetaI_c     [deg],  BetaI at the control points 
%   Xf          [m]     Axial separation between propellers 
%   AlphaI      [deg],  ideal angle of attack 
%   Z           [ ],    number of blades 
%   Rhub        [m],    hub radius  
%   CoD         [ ],    chord / diameter at each control point radius 
%   R           [m],    propeller radius  
%   M           [ ],    number of radial 2D cross-sections 
%   Np          [ ],    number of points in each 2D section 
% 
% ------------------------------------------------------------------------- 
  
function [f0oc1,f0oc2,t0oc1,t0oc2,AlphaI1,AlphaI2,X3D,Y3D,Z3D,... 
          X3D_aft,Y3D_aft,Z3D_aft,c1,c2,x0_1,x0_2,theta_Z1,theta_Z2] = ... 
          Geometry(XR1,XR2,t0oc01,t0oc02,skew01,skew02,rake01,... 
                   rake02,RC1,RC2,Cl1,Cl2,BetaI_c1,BetaI_c2,Xf,Z1,Z2,... 
                   Rhub,CoD1,CoD2,R1,R2,M1,M2,Np) 
  
% ---------------------------- Interpolate input geometry at control points 
f0oc1=0.0679*Cl1;                   %max camber ratio (NACA a=0.8 meanline) 
f0oc2=0.0679*Cl2;                    
t0oc1 = pchip(XR1,t0oc01,RC1);       % [ ],   thickness ratio 
t0oc2 = pchip(XR2,t0oc02,RC2);        
skew1 = pchip(XR1,skew01,RC1);       % [deg], angular translation along  
%                                             mid-chord helix 
skew2 = pchip(XR2,skew02,RC2); 
D1=2*R1; 
D2=2*R2; 
rake1 = pchip(XR1,rake01,RC1)*D1;     % [m],   translation along propeller  
%                                              axis (3D X-axis) 
rake2 = pchip(XR2,rake02,RC2)*D2; 
AlphaI1=1.54*Cl1; 
AlphaI2=1.54*Cl2; 
  
% --------------- Find basic geometry parameters chord, radius, pitch, etc. 
theta_nt1 = BetaI_c1 + AlphaI1;      % Nose-tail pitch angle, [deg] 
theta_nt2 = BetaI_c2 + AlphaI2; 
  
PoD1     = tand(theta_nt1).*pi.*RC1; % Pitch / propeller diameter, [ ] 
PoD2     = tand(theta_nt2).*pi.*RC2; 
c1       = CoD1.*D1;                 % section chord at the c. points [m]                  
c2       = CoD2.*D2; 
r1       = RC1.*R1;                  % radius of the c. points [m] 

 
r2       = RC2.*R2; 
theta_Z1 = 0:360/Z1:360;             % angle between blades [deg] 
theta_Z2 = 360/(2*Z2):360/Z2:360+360/(2*Z2);   % angle between blades [deg] 
%       or 0:360/Z2:360 
% ---------------------------------------- Lay out the 2D coordinate system 
% 
% xN   [ ], x/c coordinate in 2D NACA foil tables 
%               At the Leading  Edge: xN = 0, x1 =  c/2, x0 = 0 
%               At the Trailing Edge: xN = 1, x1 = -c/2, x0 = 1 
% x0   [ ], x/c distance along mid-chord line to interpolate NACA data. 
% x1   [m], x   distance along mid-chord line to evaluate elliptical or  
%               parabolic formulae. By definition, x1 == c/2 - c*x0. 
% x2D  [m], x   position in 2D space on upper and lower foil surfaces 
% y2D  [m], y   position in 2D space on upper and lower foil surfaces 
% x2Dr [m], x   position in 2D space after rotation for pitch angle 
% y2Dr [m], y   position in 2D space after rotation for pitch angle 
%  
  
xN = [0 .5 .75 1.25 2.5 5 7.5 10 15 20 25 30 35 40 45 50 ... 
     55 60 65 70 75 80 85 90 95 100]./100;  
  
for i = 1:M1                     % for each radial section along the span 
    for j = 1:Np                 % for each point          along the chord 
        x0_1(1,j)      =               (j-1)/(Np-1);    % [0   :    1] 
        x1_1(i,j)      = c1(i)/2 - c1(i)*(j-1)/(Np-1);  % [c/2 : -c/2] 
    end 
end 
  
for i = 1:M2                     % for each radial section along the span 
    for j = 1:Np                 % for each point          along the chord 
        x0_2(1,j)      =               (j-1)/(Np-1);    % [0   :    1] 
        x1_2(i,j)      = c2(i)/2 - c2(i)*(j-1)/(Np-1);  % [c/2 : -c/2] 
    end 
end 
  
% ------------------ Find meanline and thickness profiles (at x1 positions) 
% 
% foc    = camber / chord ratio (NACA data at xN positions) 
% dfdxN  = slope of camber line (NACA data at xN positions) 
% fscale = scale to set max camber    ratio to f0oc for each section 
% tscale = scale to set max thickness ratio to t0oc for each section 
% f      = camber               at x1 positions  
% dfdx   = slope of camber line at x1 positions 
% t      = thickness            at x1 positions 
  
         % ------------------------- Use NACA a=0.8 meanline 
    foc = [0 .287 .404 .616 1.077 1.841 2.483 3.043 3.985 4.748 ... 
           5.367 5.863 6.248 6.528 6.709 6.79 6.77 6.644 6.405  ... 
           6.037 5.514 4.771 3.683 2.435 1.163 0]./100; 
  
    dfdxN = [.48535 .44925 .40359 .34104 .27718 .23868 .21050 ... 
             .16892 .13734 .11101 .08775 .06634 .04601 .02613 ... 
             .00620 -.01433 -.03611 -.06010 -.08790 -.12311   ... 
             -.18412 -.23921 -.25583 -.24904 -.20385]; 
  

 
    fscale1 = f0oc1 / max(foc); 
    fscale2 = f0oc2 / max(foc); 
     
    for i = 1:M1 
        for j = 1:Np 
            f1(i,:)    = pchip(xN       ,foc  .*fscale1(i).*c1(i),x0_1);  
            dfdx1(i,:) = pchip(xN(2:end),dfdxN.*fscale1(i)      ,x0_1);  
        end 
    end 
  
    for i = 1:M2 
        for j = 1:Np 
            f2(i,:)    = pchip(xN       ,foc  .*fscale2(i).*c2(i),x0_2);  
            dfdx2(i,:) = pchip(xN(2:end),dfdxN.*fscale2(i)      ,x0_2);  
        end 
    end 
  
    %this is for NACA 66mod with t0/c=0.1 
%    toc_66 = [0 .665 .812 1.044 1.466 2.066 2.525 2.907 3.521 4 ... 
%              4.363 4.637 4.832 4.952 5 4.962 4.846 4.653     ... 
%              4.383 4.035 3.612 3.11 2.532 1.877 1.433 .333]./100; 
    %this is for the NACA 65A010 
     toc_65 = [0 .765 .928 1.183 1.623 2.182 2.65 3.04 3.658 4.127 ... 
               4.483 4.742 4.912 4.995 4.983 4.863 4.632 4.304     ... 
               3.899 3.432 2.912 2.352 1.771 1.188 .604 .021]./100; 
  
    tscale1 = t0oc1 / max(toc_65); 
    tscale2 = t0oc2 / max(toc_65); 
  
    for i = 1:M1    
        for j = 1:Np    
            t1(i,:) = pchip(xN,toc_65.*tscale1(i).*c1(i),x0_1);   
        end 
    end 
  
    for i = 1:M2    
        for j = 1:Np    
            t2(i,:) = pchip(xN,toc_65.*tscale2(i).*c2(i),x0_2);   
        end 
    end 
  
% ------------------------------------- Find 2D unroatated section profiles 
% x2D  [m], x position in 2D space on upper (x2D_u) and lower (x2D_l) surf. 
% y2D  [m], y position in 2D space on upper (y2D_u) and lower (y2D_l) surf. 
for i = 1:M1                             % for each section along the span 
    for j = 1:Np                         % for each point   along the chord 
        x2D_u(i,j) = x1_1(i,j) + (t1(i,j)/2)*sin(atan(dfdx1(i,j)));  
        x2D_l(i,j) = x1_1(i,j) - (t1(i,j)/2)*sin(atan(dfdx1(i,j)));   
        y2D_u(i,j) =  f1(i,j) + (t1(i,j)/2)*cos(atan(dfdx1(i,j)));  
        y2D_l(i,j) =  f1(i,j) - (t1(i,j)/2)*cos(atan(dfdx1(i,j)));  
    end 
end 
  
for i = 1:M2                             % for each section along the span 
    for j = 1:Np                         % for each point   along the chord 

 
        x2D_u_aft(i,j) = x1_2(i,j) + (t2(i,j)/2)*sin(atan(dfdx2(i,j)));  
        x2D_l_aft(i,j) = x1_2(i,j) - (t2(i,j)/2)*sin(atan(dfdx2(i,j)));   
        % ------------ For aft propeller signs are reversed---------------- 
        y2D_u_aft(i,j) =  -f2(i,j) - (t2(i,j)/2)*cos(atan(dfdx2(i,j)));  
        y2D_l_aft(i,j) =  -f2(i,j) + (t2(i,j)/2)*cos(atan(dfdx2(i,j)));  
    end 
end 
  
% -----------------------Put all the numbers in one list------------------- 
% First Np values are the upper surface (suction side),and the second Np  
% values are the lower surface (pressure side). 
x2D(:,   1:Np   ) = x2D_u(:,1:Np);      
x2D(:,1+Np:Np+Np) = x2D_l(:,Np:-1:1);   
y2D(:,   1:Np   ) = y2D_u(:,1:Np); 
y2D(:,1+Np:Np+Np) = y2D_l(:,Np:-1:1); 
  
% ---------------------------- Put all the numbers in one list for aft prop 
x2D_aft(:,   1:Np   ) = x2D_u_aft(:,1:Np);      
x2D_aft(:,1+Np:Np+Np) = x2D_l_aft(:,Np:-1:1);   
y2D_aft(:,   1:Np   ) = y2D_u_aft(:,1:Np); 
y2D_aft(:,1+Np:Np+Np) = y2D_l_aft(:,Np:-1:1); 
  
% --------------------------------------- Find 2D rotated section profiles 
% x2Dr [m], x position in 2D space after rotation for pitch angle 
% y2Dr [m], y position in 2D space after rotation for pitch angle 
for i = 1:M1          % for each section along the span 
    for j = 1:2*Np    % for each point   along the upper and lower surfaces 
        x2Dr(i,j) = x2D(i,j)*cosd(theta_nt1(i))... 
                - y2D(i,j)*sind(theta_nt1(i)); % rotated 2D upper surface x 
        y2Dr(i,j) = x2D(i,j)*sind(theta_nt1(i))... 
                + y2D(i,j)*cosd(theta_nt1(i)); % rotated 2D upper surface y 
    end 
end 
  
% --------------------------- Find 2D rotated section profiles for aft prop 
theta_nt_aft=180-theta_nt2; 
for i = 1:M2          % for each section along the span 
    for j = 1:2*Np    % for each point   along the upper and lower surfaces 
        x2Dr_aft(i,j) = x2D_aft(i,j)*cosd(theta_nt_aft(i))... 
            - y2D_aft(i,j)*sind(theta_nt_aft(i)); % rotated upper surface x 
        y2Dr_aft(i,j) = x2D_aft(i,j)*sind(theta_nt_aft(i))... 
            + y2D_aft(i,j)*cosd(theta_nt_aft(i)); % rotated upper surface y 
    end 
end 
  
% --------------------------- Invoke skew and rake, and find 3D coordinates 
% X3D [m], X position in 3D space (corresponds to y position in 2D space) 
% Y2D [m], Y position in 3D space 
% Z3D [m], Z position in 3D space 
  
for i = 1:M1          % for each section along the span 
    for j = 1:2*Np    % for each point   along the upper and lower surfaces 
        X3D(i,j,1) = - rake1(i) ... 
            - r1(i)*(pi*skew1(i)/180)*tand(theta_nt1(i)) + y2Dr(i,j); 
      
 
 
        for k = 1:Z1   % for each blade 
            Y3D(i,j,k) = r1(i)*sind(skew1(i)... 
                - (180/pi)*x2Dr(i,j)/r1(i) - theta_Z1(k)); 
            Z3D(i,j,k) = r1(i)*cosd(skew1(i)... 
                - (180/pi)*x2Dr(i,j)/r1(i) - theta_Z1(k)); 
        end 
    end 
end 
% -------------- Invoke skew and rake, and find 3D coordinates for aft prop 
  
for i = 1:M2          % for each section along the span 
    for j = 1:2*Np    % for each point   along the upper and lower surfaces 
        X3D_aft(i,j,1) = - rake2(i) - ... 
            r2(i)*(pi*skew2(i)/180)*tand(theta_nt_aft(i)) + y2Dr_aft(i,j); 
        X3D_aft(i,j,1)= X3D_aft(i,j,1)-Xf; %*R1; 
        for k = 1:Z2   % for each blade 
            Y3D_aft(i,j,k) = r2(i)*sind(skew2(i) ... 
                - (180/pi)*x2Dr_aft(i,j)/r2(i) - theta_Z2(k)); 
            Z3D_aft(i,j,k) = r2(i)*cosd(skew2(i)... 
                - (180/pi)*x2Dr_aft(i,j)/r2(i) - theta_Z2(k)); 
        end 
    end 
end 
  
% ----------------------------------------------- Create 3D Propeller Image 
  
    Fig3_S = figure('units','normalized','position',[.61 .06 .4 .3],... 
                    'name','Propeller Image','numbertitle','off'); 
    hold on; 
        
    % ------------------------------------------ Plot the propeller surface    
    for k = 1:Z1 
        surf(X3D(:,:,1),Y3D(:,:,k),Z3D(:,:,k));   
    end 
     
    for k = 1:Z2 
        surf(X3D_aft(:,:,1),Y3D_aft(:,:,k),Z3D_aft(:,:,k)); 
    end 
  
    colormap gray;      
    shading interp;   
%     shading faceted; 
    grid on;         
    axis equal; 
    axis([-0.2 R1 -1.1*R1 1.1*R1 -1.1*R1 1.1*R1]); 
    xlabel('X (3D) [m]','FontSize',12);  
    ylabel('Y (3D) [m]','FontSize',12);  
    zlabel('Z (3D) [m]','FontSize',12);  
    title('3D Propeller Image','FontSize',16); 
  
    % -------------------------------------------------------- Plot the hub 
    tick = 90:-15:0; 
    [yh0,zh0,xh0] = cylinder(Rhub*sind(tick),50);    
    xh0 = -0.5*c1(1)*xh0 - 0.75*R1; 
    surf(xh0,yh0,zh0); 

 
     
    [yh1,zh1,xh1] = cylinder(Rhub,50); 
    xh1 = 6*c1(1)*xh1 - 0.75*R1; 
    surf(xh1,yh1,zh1);     
  
    % ----- Plot the suction side (green) & pressure side (red) of the prop 
    for i = 1:M1          % for each section along the span 
        for k = 1:Z1       % for each blade  
            plot3(X3D(i,1:Np,1),Y3D(i,1:Np,k),Z3D(i,1:Np,k),... 
                'g','Linewidth',1); % suction surface 
            plot3(X3D(i,Np+1:2*Np,1),Y3D(i,Np+1:2*Np,k),... 
                Z3D(i,Np+1:2*Np,k),'r','Linewidth',1); % pressure surface 
        end 
    end 
     
    for i = 1:M2          % for each section along the span 
        for k = 1:Z2       % for each blade  
%               Now for aft prop 
            plot3(X3D_aft(i,1:Np,1),Y3D_aft(i,1:Np,k),Z3D_aft(i,1:Np,k),... 
                'g','Linewidth',1); % suction surface 
            plot3(X3D_aft(i,Np+1:2*Np,1),Y3D_aft(i,Np+1:2*Np,k),... 
                Z3D_aft(i,Np+1:2*Np,k),'r','Linewidth',1);%pressure surface 
        end 
    end 
     
    for j = 1:Np          % for each point along the chord 
        for k = 1:Z1       % for each blade  
            plot3(X3D(:,j,1),Y3D(:,j,k),Z3D(:,j,k),... 
                'g','Linewidth',1); % suction surface 
            plot3(X3D(:,j+Np,1),Y3D(:,j+Np,k),Z3D(:,j+Np,k),... 
                'r','Linewidth',1); % pressure surface 
        end 
    end   
     
    for j = 1:Np          % for each point along the chord 
        for k = 1:Z2       % for each blade  
            %  Now for aft prop 
            plot3(X3D_aft(:,j,1),Y3D_aft(:,j,k),Z3D_aft(:,j,k),... 
                'g','Linewidth',1); % suction surface 
            plot3(X3D_aft(:,j+Np,1),Y3D_aft(:,j+Np,k),... 
                Z3D_aft(:,j+Np,k),'r','Linewidth',1); % pressure surface 
        end 
    end   
     
    % --------------------------------- Plot the leading and trailing edges 
    for k = 1:Z1           % for each blade 
        plot3(X3D(:,1,1), Y3D(:,1,k), Z3D(:,1,k), 'b','Linewidth',2); %L.E. 
        plot3(X3D(:,Np,1),Y3D(:,Np,k),Z3D(:,Np,k),'k','Linewidth',2); %T.E. 
    end 
     
    for k = 1:Z2           % for each blade 
        plot3(X3D_aft(:,1,1),Y3D_aft(:,1,k),Z3D_aft(:,1,k),... 
            'b','Linewidth',2); %L.E. 
        plot3(X3D_aft(:,Np,1),Y3D_aft(:,Np,k),Z3D_aft(:,Np,k),... 
            'k','Linewidth',2); %T.E. 

 
    end 
  
    % ------------------------------------------ Plot the coordinate system 
  
    % Axes 
    plot3([0 R1],[0 0],[0 0],'y','LineWidth',2), 
    plot3([0 0],[0 R1],[0 0],'y','LineWidth',2), 
    plot3([0 0],[0 0],[0 R1],'y','LineWidth',2), 
     
    % Circle at the X = 0 location on the hub 
    phi = 0:0.01:2*pi; 
    Xhc =   zeros(size(phi)); 
    Yhc = - Rhub * sin(phi); 
    Zhc =   Rhub * cos(phi); 
    plot3(Xhc,Yhc,Zhc,'y','LineWidth',2), 
    % Circle at the X = -Xf*R1/2 location on the hub 
    Xhc_mid =   -(Xf/2)*R1*ones(size(phi)); 
    plot3(Xhc_mid,Yhc,Zhc,'k','LineWidth',4) 
     
 