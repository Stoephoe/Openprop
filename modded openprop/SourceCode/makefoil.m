%Code by Chris Peterson. Code will read in specified camber and thickness
% distributions and generate foil geometry file for XFOIL. Thickness and
% camber are scaled to t set and f set.
% Coordinates start at TE, go forward CCW along upper surfact to LE,
% and back to TE along lower surface.

function [] = makefoil(t_set, f_set, mean_type, thick_type, save_as)

% clc; clear all; close all;
% t set 0.1;
% f set = 0.08;
% mean type = 'NACAa=08(Brockett).txt';
% thick type = 'NACA66(Brockett).txt';
% save as = 'brockett';

make_plot   = 'no';     %Generate plot toggle ('yes' or 'no')
N_parab_def = 35;       %Number of points to make nose parabola. Fails at
                        %numbers < -20
N_parab_eval= 6;        %Number of points to include at the nose in
                        %data export;
N_surf_pts  = 80;       %Number of points along body to TE
                        %(not including LE)
                        %N_parab pts + N surf pts must be < 150
fract       = 1-2/N_parab_eval; %Fraction of parabola to use from LE
                                %to 0.005.
                                %Max parabola point must be less than 0.005
                                %to prevent sharp cornder at 0.005

conc_fact = 2;          %Power for exponential disribution at LE. This
                        %concentrates point near tip.

                    
%Get meanline and dy/dx distributions from mean line data base
[x_f fc_o dydx_o]   = getmeanline(mean_type);
[x_t tc_o RLE_o]    = getthickdist(thick_type);

%scale appropriately
t_set   = t_set/2;      %uses 1/2 thickness
if max(fc_o) ~= 0
    f_scale = f_set/max(fc_o);
elseif max(fc_0) == 0
    f_scale = 0;
end

f_c       = fc_o * f_scale;
dydx      = dydx_o * f_scale;
t_scale   = t_set/max(tc_o);
t_c       = tc_o * t_scale;
RLE       = RLE_o * (t_scale)^2;

%Findpoints along RLE nose parabola
x_RLE     = fract*0.005*(0:1/(N_parab_def-1):1).^conc_fact;
t_RLE     = sqrt(2*RLE*(x_RLE));

%Spline parabola and tabulated data for thickness function
x_locs    = [x_RLE x_t(2:end)];   %New combined x/c values

%1e8 sets init slope = ~inf
t_fnct    = csape(x_locs, [1e10 t_RLE t_c(2:end) 1], [1 0]);
    %Make x locations for generating data file
    %Cosine spasing from 0.005 to TE
x_cos_sp  = 0.005 + 0.5*0.995*(1-cos(0:pi/(N_surf_pts-1):pi));
    %Exponetial spacing for nose
x_eval_LE = fract*0.005*(0:1/(N_parab_eval-1):1).^conc_fact;
t_eval_LE = sqrt(2*RLE*(x_eval_LE));
x_eval_mb = [x_cos_sp];                    %Establishes eval points
t_eval_mb = fnval(t_fnct, x_eval_mb);      %Evaluates spline at eval points
x_eval    = [x_eval_LE x_eval_mb];
t_eval    = [t_eval_LE t_eval_mb];

%spline tabulated data for camber at same x/c locations as thickness
f_fnct    = csape(x_f, f_c);
f_eval    = fnval(f_fnct, x_eval);
dydx_eval = fnval(fnder(f_fnct), x_eval);

%Plotting for unrotated parameters
if strcmp(make_plot,'yes')
    figure();
    hold on;
    axis equal;                             %Set X:Y to unity
    title('Camber, Thickness, and LE Graphical Display')
    xlable('X/C');
    xlim([-0.01 0.25]);                     %Set Initial Zoom
        %Plot thickness
        fnplt(t_fnct, 'y'); fnplt(f_fnct,'g')
        plot(x_t, t_c, 'co'); plot(x_f, f_c, 'ro')
        plot(x_RLE, t_RLE,'k.');
        %Plot RLE Circle and parabola for viewing on plot
        plot(RLE - RLE*cos(0:pi/100:pi), RLE*(sin(0:pi/100:pi)), 'b.');
        plot((0:1/10000:0.2), sqrt(2*RLE*(0:1/10000:0.2)), 'r.');
        %Plot camber

    legend('Splined Thickness', 'Splined Camber',...
        'Tabulated Thickness (Scaled)','Tabulated camber (scaled)',...
        'Calcuated Parabola', 'Leading Edge Radius', 'LE Parabola',...
        'Location', 'southeast')
end

%Calculate upper and lower surface ordinates
x_u = x_eval - t_eval.*sin(atan(dydx_eval));
y_u = f_eval + t_eval.*cos(atan(dydx_eval));
x_l = x_eval + t_eval.*sin(atan(dydx_eval));
y_l = f_eval - t_eval.*cos(atan(dydx_eval));

%Solve for most forward point on foil
[x_fwd, min_i] = min(x_u);
y_fwd = y_u(min_i);

%New plot for actual upper and lower surfaces
if strcmp(make_plot,'yes')
    figure();
    hold on;
    axis equal;                             %Set X:Y to unity
    xlim([0 1]);                            %Set Initial Zoom
    plot(x_u, y_u, 'b-', x_u, y_u, 'r.')
    plot(x_l, y_l, 'b-', x_l, y_l, 'r.');
    plot(x_eval, f_eval, 'g-', x_eval, f_eval, 'r.')
    plot(x_fwd, y_fwd, 'kp')
end

%Combine coordinates into a single array of points from TE along upper
%surface around LE back to TE along lower surace
x_comb  = [fliplr(x_u) x_l];
y_comb  = [fliplr(y_u) y_l];

%Rotate and scale such that max forward point is at 0,0, and TE is at 0,1.
%Assumes TE is already at 0,0 (Uses method in Brockett Report)
shift_ang = atan(y_fwd/(1-x_fwd));

%Scaled chord length back to 1 (accounts for portion forward of 0)
x_scaled = (x_comb-x_fwd)./(1-x_fwd);
y_scaled = (y_comb-y_fwd)./(1-x_fwd);
%Rotate so that most forward point is at 0,0
x_rot = (x_scaled.*cos(shift_ang) - y_scaled.*sin(shift_ang))/...
            sqrt(1+(y_fwd/(1-x_fwd))^2);
y_rot = (y_scaled.*cos(shift_ang) + x_scaled.*sin(shift_ang))/...
            sqrt(1+(y_fwd/(1-x_fwd))^2);

%New plot for final upper and lower surfaces
if strcmp(make_plot,'yes')
    figure ();
    hold on;
    title('Final Points exported to Data File.');
    axis equal; %Set X:Y to unity
    xlim([O 1.1]); %Set Initial Zoom
    plot(x_rot, y_rot, x_rot, y_rot, 'r.');
    legend('Connect the dots', 'Actual data points');
end

%Write to text file for use in XFOIL.
cmd = ['del ', 'foildata.txt'];    %save as is file name to be written to
system(cmd);                %Delets previous file
fid = fopen('foildata.txt', 'w');   %permission specifier changed from 'w' to 'wt'
for i = 1:length(x_rot)
    fprintf(fid, '%10.8f %10.8f\n', x_rot(i), y_rot(i));
end
fclose (fid);
