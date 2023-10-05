function [UAHIF,UTHIF] = Horseshoe_self(Mp,Z,TANBIV,RC,RV,Hub_Flag,... 
                                                                  Rhub_oR) 
  
for n = 1:Mp                 % for each control point, n     (FOR LOOP MF2)   
    for m = 1:Mp+1           % for each vortex  point, m     (FOR LOOP MF3) 
  
        % -- Find induction factors for a unit vortex shed at RV(m) 
        % -- Wrench returns 2*pi*R*u_bar 
        [UAW(m),UTW(m)] = Wrench(Z,TANBIV(m),RC(n),RV(m)); 
  
        % ---------------------------- Find hub-image effects, Kerwin p.181 
        if Hub_Flag == 1  
            RCW    = RC(n); 
            RVW    = Rhub_oR^2/RV(m);              
            TANBIW = TANBIV(1)*RV(1)/RVW;      
  
            [UAWh,UTWh] = Wrench(Z,TANBIW,RCW,RVW); 
  
            UAW(m) = UAW(m)-UAWh; 
            UTW(m) = UTW(m)-UTWh;                
        end 
    end                                                % (END FOR LOOP MF3) 
  
    % The Horseshoe Influence Function for vortex panel m is the 
    % effect of the induction by a helical trailing vortex at  
    % vortex point m   with circulation -Gamma(m) and another at  
    % vortex point m+1 with circulation +Gamma(m).  
    % UAHIF(n,m) = u_barA horseshoe influence function in eqn 254. 
    % UAW(m)     = u_barA Wrench velocity given in eqn 202-203. 
    for m = 1:Mp                                % for each vortex  panel, m 
        UAHIF(n,m) = UAW(m+1)-UAW(m);           % 2*pi*R*(HIF)      
        UTHIF(n,m) = UTW(m+1)-UTW(m);           % 2*pi*R*(HIF) 
    end 
     
end  % (END FOR LOOP MF2) 
end 

 
  
% ================================================== END Horseshoe Function   
% ========================================================================= 
  
% ========================================================================= 
% ========================================================= Wrench Function  
% 
% This function evaluates the Wrench u_bar velocity induced on a point on  
% a lifting line due to a helical trailing vortex. These formulae were  
% derived in 1957 by J.W. Wrench. This function returns u_bar given in  
% Kerwin eqns 202-205, p.154. 
% 
% NOTE:  There are TWO ERRORS in Kerwin, as of the Spring 2007 printing. 
%        These have been corrected in the present implementation. 
% 
% 1. Eqn 202, u_bar_a. Should be (y-2*Z*y*y0*F1), not (y-2*Z*rv*F1) to  
%    agree with Wrench, eqn 31. 
% 2. Eqn 204, F2.  Need to kill the leading "-" sign to make F2 agree  
%    with Wrench equation 29.   
% 
% ------------------------------------------------------------------------- 
% Variables: 
    % Z         [ ],    number of blades 
    % tan_betaW [ ],    tangent of the pitch angle of helical wake trail 
    % rc        [ ],    radius of control point  / propeller radius 
    % rv        [ ],    radius of helical vortex / propeller radius 
  
    % u_barA    [ ],    Wrench u_bar velocity in the axial      direction 
    % u_barT    [ ],    Wrench u_bar velocity in the tangential direction 
    % y,y0,U,F1,F2,     auxilary variables.  See Kerwin eqns. 202-205.  
% 
% ------------------------------------------------------------------------- 
  
function [u_barA, u_barT] = Wrench(Z,tan_betaW,rc,rv) 
  
%     % --------------- Enable this to check for infinite bladed propellers 
%     if Z > 20   % Return infinite blade result if Z > 20. 
%         if rc = rv 
%             IF_A = 0; 
%             IF_T = 0;   
%              
%         elseif rc < rv 
%             IF_A = Z/(2*rv*tan_betaW);    % 2*pi*R*(eqn 206) 
%             IF_T = 0;                     % 2*pi*R*(eqn 206) 
%              
%         else % rc > rv             
%             IF_A = 0;                     % 2*pi*R*(eqn 207) 
%             IF_T = Z/(2*rc);              % 2*pi*R*(eqn 207) 
%         end 
%         return; 
%     end 
  
  
y  = rc/(rv*tan_betaW); 
y0 = 1/tan_betaW; 

 
U  = ((y0*(sqrt(1+y ^2)-1))*exp(sqrt(1+y^2)-sqrt(1+y0^2))/... 
      (y *(sqrt(1+y0^2)-1)))^Z; 
  
   
if rc == rv 
    IF_A = 0; 
    IF_T = 0; 
     
elseif rc < rv     
     
    F1     = -(1/(2*Z*y0)) * ((1+y0^2)/(1+y^2))^0.25 * ... 
             ((1/(U^-1-1)) + (1/(24*Z))*... 
             (((9*y0^2+2)/(1+y0^2)^1.5)+((3*y^2-2)/(1+y^2)^1.5))*... 
             log(1+(1/(U^-1-1))) ); 
   
    u_barA = (Z  /(2*rc))*(y-2*Z*y*y0*F1);     % 2*pi*R*(eqn 202) 
    u_barT = (Z^2/rc)    *(y0*F1);             % 2*pi*R*(eqn 202) 
  
else % rc > rv 
     
    F2     = (1/(2*Z*y0)) * ((1+y0^2)/(1+y^2))^0.25 * ... 
             ((1/(U-1))    - (1/(24*Z))*... 
             (((9*y0^2+2)/(1+y0^2)^1.5)+((3*y^2-2)/(1+y^2)^1.5))*... 
             log(1+(1/(U-1))) );     
         
    u_barA = -(Z^2/rc)    *(y*y0*F2);          % 2*pi*R*(eqn 203) 
    u_barT =  (Z  /(2*rc))*(1+2*Z*y0*F2);      % 2*pi*R*(eqn 203) 
    
end 
end 
% ===================================================== END Wrench Function  
% ========================================================================= 