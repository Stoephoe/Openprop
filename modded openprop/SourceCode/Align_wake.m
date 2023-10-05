% ========================================================================= 
% ===================================================== Align_wake Function 
% 
% This function aligns the wake to the given circulation distribution by 
% iteratively computing: 
%   UAHIF1,UTHIF1,UAHIF2,UTHIF2         = the horseshoe influence functions 
%                                         for the self-induced velocities 
%   UAHIF1_2,UTHIF1_2,UAHIF2_1,UTHIF2_2 = the horseshoe influence functions 
%                                         for the interaction velocities   
%   UASTAR1,UTSTAR1,UASTAR2,UTSTAR2     = the induced velocities 
%   TANBIC1,TANBIV1,TANBIC2,TANBIV2     = the velocity angles 
% 
% ------------------------------------------------------------------------- 
  
function [UAHIF1,UTHIF1,UAHIF2,UTHIF2,UAHIF1_2,UTHIF1_2,UAHIF2_1,... 
    UTHIF2_1,UASTAR1,UTSTAR1,UASTAR2,UTSTAR2,TANBIC1,TANBIV1,TANBIC2,... 
    TANBIV2] = Align_wake(TANBIC1,TANBIV1,TANBIC2,TANBIV2,ITER,M1,M2,... 
    Z1,Z2,RC1,RV1,RC2,RV2,G1,G2,VAC1,VTC1,VAC2,VTC2,Js1,Js2,Xf,Hub_Flag,... 
    Rhub_oR1,Rhub_oR2) 
         
    % ----------- Iterate to ALIGN WAKE to the new circulation distribution 
    W_iter = 1;                     % iteration in the wake alignment loop 
    W_res1  = 1;                     % residual BetaI between interations 
    W_res2  = 1; 
    TANBIW1_last = TANBIC1;           % the last value of TANBIC 
    TANBIW2_last = TANBIC2; 
  
    while W_iter < ITER & (W_res1 > 1e-5 | W_res2 > 1e-5 )%(WHILE LOOP WA1) 
         
        % --------- Compute the vortex Horseshoe Influence Functions ------ 
                                         

 
        [UAHIF1,UTHIF1]=Horseshoe_self(M1,Z1,TANBIV1,RC1,RV1,Hub_Flag,... 
                                                                 Rhub_oR1); 
        [UAHIF1_2,UTHIF1_2]=Horseshoe_int(M1,M2,Z2,TANBIV2,RC1,RV2,-Xf,... 
                                                        Hub_Flag,Rhub_oR2); 
     
        [UAHIF2,UTHIF2]=Horseshoe_self(M2,Z2,TANBIV2,RC2,RV2,Hub_Flag,... 
                                                                 Rhub_oR2); 
        [UAHIF2_1,UTHIF2_1]=Horseshoe_int(M2,M1,Z1,TANBIV1,RC2,RV1,Xf,... 
                                                        Hub_Flag,Rhub_oR1); 
         
         
         
        % ---- Compute induced velocities at control points. eqn 254, p.179 
        %[UASTAR,UTSTAR] = Induced_Velocity(Mp,G,UAHIF,UTHIF,UADUCT,dCirc); 
         
        [UA_SELF1,UT_SELF1,UA_INT1_2,UT_INT1_2]=Induced_Velocity(M1,M2,... 
                                    G1,G2,UAHIF1,UTHIF1,UAHIF1_2,UTHIF1_2); 
        UASTAR1=UA_SELF1+UA_INT1_2;  %total axial induced velocity/Vs 
        UTSTAR1=UT_SELF1+UT_INT1_2;  %total tangential induced velocity/Vs 
  
        [UA_SELF2,UT_SELF2,UA_INT2_1,UT_INT2_1]=Induced_Velocity(M2,M1,... 
                                    G2,G1,UAHIF2,UTHIF2,UAHIF2_1,UTHIF2_1); 
        UASTAR2=UA_SELF2+UA_INT2_1;  %total axial induced velocity/Vs 
        UTSTAR2=UT_SELF2+UT_INT2_1;  %total tangential induced velocity/Vs 
         
         
         
         
        % --------------- Compute tan(BetaI) for the new induced velocities 
        [TANBIC1,TANBIV1] = find_tan_BetaI(VAC1,VTC1,UASTAR1,UTSTAR1,... 
                                                              RC1,RV1,Js1); 
        [TANBIC2,TANBIV2] = find_tan_BetaI(VAC2,VTC2,UASTAR2,UTSTAR2,... 
                                                              RC2,RV2,Js2); 
         
         
        % ---------------------------------- Prepare for the next iteration 
        W_iter = W_iter + 1                 % iteration in the BetaI loop 
        W_res1  = abs(TANBIC1 - TANBIW1_last); % residual BetaI 
        W_res2  = abs(TANBIC2 - TANBIW2_last); 
        TANBIW1_last = TANBIC1;               % the last value of TANBIC 
        TANBIW2_last = TANBIC2; 
     
        if W_iter > ITER 
            warning('on'), 
            warning('WARNING: While loop WA1 did NOT converge.'), 
            warning('off'), 
        end  
    end                                              % (END WHILE LOOP WA1) 
  
% ================================================= END Align_wake Function  
% ========================================================================= 
 
 