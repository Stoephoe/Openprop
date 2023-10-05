% ========================================================================= 
% =============================================== Induced_Velocity Function 
% 
% This function computes induced velocities at control points, Kerwin  
% eqn 254, p.179, normalized by the ship speed. 
% The self-induced velocities are assumed to be those having an index of 1, 
% while the interaction are those induced by component 2 on component 1 
% having an index of 1_2. 
% ------------------------------------------------------------------------- 
  
function [UA_SELF,UT_SELF,UA_INT,UT_INT] = Induced_Velocity(M1,M2,G1,G2,... 
                                           UAHIF1,UTHIF1,UAHIF1_2,UTHIF1_2) 
  
UA_SELF(1:M1) = 0;                       
UT_SELF(1:M1) = 0; 
UA_INT(1:M1)  = 0; 
UT_INT(1:M1)  = 0; 
  
for n = 1:M1                                    % for each control point, n       
    for m = 1:M1                                % for each vortex  panel, m 
       UA_SELF(n) = UA_SELF(n) + G1(m)*UAHIF1(n,m); % UASTAR / ship speed   
       UT_SELF(n) = UT_SELF(n) + G1(m)*UTHIF1(n,m); % UASTAR / ship speed  
    end 
    for m=1:M2 
        UA_INT(n) = UA_INT(n)  + G2(m)*UAHIF1_2(n,m); 
        UT_INT(n) = UT_INT(n)  + G2(m)*UTHIF1_2(n,m); 
    end 
end  
end 
% =========================================== END Induced_Velocity Function  
% ========================================================================= 