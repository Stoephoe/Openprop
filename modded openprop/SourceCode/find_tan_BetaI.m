% ========================================================================= 
% ================================================= find_tan_BetaI Function 
% 
% This function computes tan(BetaI), Kerwin eqn 193, p. 151. 
% UASTAR, UTSTAR represent the total induced velocities 
% (sum of self-induced and interaction velocities) 
% ------------------------------------------------------------------------- 
  
function [TANBIC,TANBIV] = find_tan_BetaI(VAC,VTC,UASTAR,UTSTAR,RC,RV,Js) 
  
VASTAR = VAC            + UASTAR; % total axial      inflow vel. / ship vel. 
VTSTAR = VTC + pi*RC/Js + UTSTAR; %      / ship vel. 
  
TANBIC = VASTAR./VTSTAR;          % tan(BetaI) at control pts. 
TANBIV = pchip(RC,TANBIC,RV);     % tan(BetaI) at vortex  pts.   
end 
% ============================================= END find_tan_BetaI Function  
% ========================================================================= 