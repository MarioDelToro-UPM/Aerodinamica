%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Vortex Lattice Method Wing Solver                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% Wing deffinition variables

NACA_perf = 2400;
AR = 4; %Aspect Ratio
sup = 1; %Surface area
estr = 0.5;
sweep = deg2rad(45);
tors = deg2rad(-4);

%% Geometry deffinition variables

Nc = 1; %Number of chordwise panels (placeholder, more than 1 WILL break this code)
Ns = 10; %Number of spanwise panels for each semiplane
Dist = 0; %0 for equidist, 1 for more vorteces toward the root, 2 for more vorteces toward the wingtips

%% Wing semiplane geometry calculations

%lets assume the center wing profile is located at x,y = 0,0

%Division of the y axis will be made as a distribution of points from said
%axis, it will NOT be done dividing any of the wing chords

%Position variables are stored for the right semiplane only

semispan = sqrt(AR/sup)/2;
cr = sup/(semispan*(1+estr)); %chord at root
ct = cr * estr; %chord at wingtip

if Dist == 0
    yw = linspace(0, semispan, Ns);
elseif Dist == 1
    yw = linspace(0, pi/2, Ns);
    yw = semispan*(1-cos(yw));
elseif Dist == 2
    yw = linspace(0, pi/2, Ns);
    yw = semispan*sen(yw);
end

xw_4 = zeros(1,Ns); %x at 1/4 chord
chord = zeros(1,Ns);


for i=1:Ns
    chord(i) = cr-(cr-ct)/semispan * yw(i);

   %Position of the 1/4 chord line of the semiplane, each airfoil is defined
   %around these points
    xw_4(i) = cr/4 + yw(i)*tan(sweep);
end


%Definition of the rest of the lattice points in the right semiplane
xw_le = xw_4 - chord/4; %x at trailing edge
xw_te = xw_4 + 3*chord/4; %x at leading edge
xw_ctrl = xw_4 +chord/2; %x at control line

x_ctrl = zeros (1, Ns-1);
y_ctrl = zeros (1, Ns-1);

%Definition of control points for induced velocity calculations
for i=1:(Ns-1)
    y_ctrl(i) = (yw(i)+yw(i+1))/2;
    x_ctrl(i) = cr/4 + y_ctrl(i)*tan(sweep) + 0.5*(cr-(cr-ct)/semispan * y_ctrl(i));
end

c_geo_mean = cr * (1+estr)/2;
c_aero_mean = 2/3 * cr * (1+estr+estr^2)/(1+estr);
y_aero_ctr = semispan * (1+2*estr)/(3*(1+estr));

x_aero_ctr = 0.25*cr;
for i=1:Ns
    x_aero_ctr = x_aero_ctr + 2*tan(sweep)/sup * chord(i)*yw(i);
end








%% Plot temporal para probar cosas

figure
hold on
plot(xw_4, yw)
plot(xw_le, yw)
plot(xw_te, yw)
plot(xw_ctrl, yw)
plot(x_ctrl,y_ctrl, LineStyle="none", Marker="o")

