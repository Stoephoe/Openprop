  %%********************************************************************%%    
  %%--------------------------Input Parameters--------------------------%%
  %%********************************************************************%%

    Rhub        = 0.066/2;                     %[m],    Hub radius (common for both propellers) 
    R1          = 0.25/2;                       %[m],    Propeller radius front prop
    R2          = 0.20/2;                       %[m],    Propeller radius back prop
    M1          = 10;                           %[ ],    Number of vortex panels over the radius front prop
    M2          = 10;                           %[ ],    Number of vortex panels over the radius back prop
    Z1          = 3;                            %[ ],    Number of blades front prop
    Z2          = 3;                            %[ ],    Number of blades back prop
    Tr          = 350;                          %[N],    Required total thrust
    q           = 1;                            %[ ],    torque ratio Q2/Q1 
    N1          = 1500;                         %[RPM],  Propeller speed front prop
    N2          = 1500;                         %[RPM],  Propeller speed back prop
    XR1         = [0.284,0.3,0.4,0.5,0.6,...
                   0.7,0.8,0.9,0.95,1];         %[ ],    Radial locations for defining inflow velocities   
                                                %        and geometric properties front prop
    XR2         = [0.284,0.3,0.4,0.5,0.6,...
                   0.7,0.8,0.9,0.95,1];         %[ ],    Radial locations for defining inflow velocities   
                                                %        and geometri properties back prop
    XCoD1       =  [0.1700,0.1912,0.2124,...
                   0.2296,0.2405,0.2411,...
                   0.2273,0.1957,0.1500,0.001]; %[ ],    chord / propeller diameter front prop
    XCoD2       = [0.1700,0.1912,0.2124,...
                   0.2296,0.2405,0.2411,...
                   0.2273,0.1957,0.1500,0.001]; %[ ],    chord / propeller diameter back prop
    CD          = [0.001,0.001,0.001,0.001,...
                   0.001,0.001,0.001,0.001,...
                   0.001,0.001];                %[ ],    Section drag coëfficient
    Xf          = 0.04;                         %[m],    Axial separation between propellers 
    ITER        = 5;                            %[ ],    Max. Iterations for Circulation Convergence and  
                                                %        Wake Alignment  
    spacing     = 'cosine';                     %        Type of radial spacing ('cosine' or 'constant') 
    f_Geometry  = 1;                            %        flag for geometry generation
    Hub_Flag    = 1;                            %        Inclusion of hub effects (1=YES, 0=NO) 
    Rhv         = 0.001;                        %        Hub Vortex Radius/Hub Radius ???
       
    Vs          = 6;                            %[m/s],  Ship speed 
    XVA1        = ones(size(XR1));              %[ ],    Va/Vs, axial inflow vel. / ship vel.
    XVA2        = ones(size(XR2));              %[ ],    Va/Vs, axial inflow vel. / ship vel. 
    XVT1        = zeros(size(XR1));             %[ ],    Vt/Vs, tangential inflow vel. / ship vel. 
    XVT2        = zeros(size(XR2));             %[ ],    Vt/Vs, tangential inflow vel. / ship vel.

  %%********************************************************************%%

addpath  ./SourceCode



[EFFY,CT1,CT2,CQ1,CQ2,KT1,KT2,KQ1,KQ2,RC1,RC2,G1,... 
        G2,UA_SELF1,UT_SELF1,UA_INT1_2,UT_INT1_2,UA_SELF2,UT_SELF2,... 
        UA_INT2_1,UT_INT2_1,TANBIC1,TANBIC2,VSTAR1,VSTAR2,Cl1,Cl2,...
        Np,X3D,Y3D,Z3D,X3D_aft,Y3D_aft,Z3D_aft,CoD1,CoD2,A]=CoupledCRP(Rhub,R1,R2,M1,...
        M2,Z1,Z2,Tr,q,N1,N2,XR1,XR2,XCoD1,XCoD2,XVA1,XVA2,XVT1,XVT2,Vs,...
        Xf,ITER,spacing,Hub_Flag,Rhv,f_Geometry,CD);

        figure(),
        plot(RC1,-CoD1,'-',RC1,CoD1,'-');
        hold on;
%%
if f_Geometry==1
    Export_Solidworks_v19('PropV5_front.csv',Np,size(X3D,1)-1,Z1,X3D,Y3D(:,:,1),Z3D(:,:,1));    
    Export_Solidworks_v19('PropV5_back.csv',Np,size(X3D_aft,1)-1,Z2,X3D_aft,Y3D_aft(:,:,ceil(Z2/2)),Z3D_aft(:,:,ceil(Z2/2)));
end


