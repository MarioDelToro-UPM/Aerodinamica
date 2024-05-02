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
Ns = 2; %Number of spanwise panels for each semiplane
Dist = 0; %0 for equidist, 1 for more vorteces toward the root, 2 for more vorteces toward the wingtips

%% Wing semiplane geometry calculations

%lets assume the center wing profile is located at x,y = 0,0

%Division of the y axis will be made as a distribution of points from said
%axis, it will NOT be done dividing any of the wing chords

%Position variables are stored for the right semiplane only

semispan = sqrt(AR/sup)/2;
cr = sup/(semispan*(1+estr)); %chord at root
ct = cr * estr; %chord at wingtip

% if Dist == 0
%     for i = 1:Nc
%         yw(:, i) = linspace(0, semispan, Ns+1);
%     end
% elseif Dist == 1
%     for i=1:Nc
%     yw(:,i) = linspace(0, pi/2, Ns+1);
%     yw(:,i) = semispan*(1-cos(yw));
%     end
% elseif Dist == 2
%     for i=1:Nc
%     yw(:,i) = linspace(0, pi/2, Ns+1);
%     yw(:,i) = semispan*sen(yw);
%     end
% end

if Dist == 0
    yw = linspace(0, semispan, Ns+1);
elseif Dist == 1
    yw = linspace(0, pi/2, Ns+1);
    yw = semispan*(1-cos(yw));
elseif Dist == 2
    yw = linspace(0, pi/2, Ns+1);
    yw = semispan*sen(yw);
end


xw_4 = zeros(1,Ns+1); %x at 1/4 chord
chord = zeros(1,Ns+1);


for i=1:Ns+1
    chord(i) = cr-(cr-ct)/semispan * yw(i);

   %Position of the 1/4 chord line of the semiplane, each airfoil is defined
   %around these points
    xw_4(i) = cr/4 + yw(i)*tan(sweep);
end


%Definition of the rest of the lattice points in the right semiplane
xw_le = xw_4 - chord/4; %x at trailing edge
xw_te = xw_4 + 3*chord/4; %x at leading edge
xw_ctrl = xw_4 +chord/2; %x at control line

x_ctrl = zeros (1, Ns);
y_ctrl = zeros (1, Ns);

%Definition of control points for induced velocity calculations
for i=1:(Ns)
    y_ctrl(i) = (yw(i)+yw(i+1))/2;
    x_ctrl(i) = cr/4 + y_ctrl(i)*tan(sweep) + 0.5*(cr-(cr-ct)/semispan * y_ctrl(i));
end

c_geo_mean = cr * (1+estr)/2;
c_aero_mean = 2/3 * cr * (1+estr+estr^2)/(1+estr);

%Aerodynamic center coordinates
y_aero_ctr = semispan * (1+2*estr)/(3*(1+estr));
x_aero_ctr = 0.25*cr + tan(sweep) * y_aero_ctr;


%% Velocity solving for the whole wing

%v_ind_mat=zeros(Nc, Ns);


% for t=1:Nc
%     for u=1:Ns
%         %first two loops iterate through panel t,u to calculate speed
%         %induced at position i,j in space
%         for i=1:Nc
%             for j=1:Ns
%                 a = x_ctrl(t,u) - xw_4(i,j);
%                 b = y_ctrl(t,u) - yw(i,j);
%                 c = x_ctrl(t,u) - xw_4(i,j+1);
%                 d = y_ctrl(t,u) - yw(i,j+1);
%                 e = sqrt(a^2+b^2);
%                 f = sqrt(c^2+d^2);
% 
%                 g = xw_4(i,j+1) - xw_4(i,j);
%                 h = yw(i,j+1) - yw(i,j);
% 
%                 k = (g*a + h*b)/e - (g*c + h*d)/f;
%                 l = -1/b * (1 + a/e) + 1/d * (1 * c/f);
% 
%                 v_ind_mat(t,u) = v_ind_mat(t,u) + k/((4*pi)*(a*d-c*d)) + l/(4*pi);
%             end
%         end
%     end
% end

v_ind_mat = zeros (Ns, Ns);


for u=1:Ns
    %first loop iterates through panel u to calculate speed
    %induced at position j in space
    for j=1:Ns
        %calculations for right semiplane-induced speed

        a = x_ctrl(j) - xw_4(u);
        b = y_ctrl(j) - yw(u);
        c = x_ctrl(j) - xw_4(u+1);
        d = y_ctrl(j) - yw(u+1);
        e = sqrt(a^2+b^2);
        f = sqrt(c^2+d^2);

        g = xw_4(u+1) - xw_4(u);
        h = yw(u+1) - yw(u);

        k = (g*a + h*b)/e - (g*c + h*d)/f;
        l = -1/b * (1 + a/e) + 1/d * (1 + c/f);

        v_ind_mat(j,u) = v_ind_mat(j,u) + k/((4*pi)*(a*d-c*d)) + l/(4*pi);


        %calculations for left semiplane-induced speed

        a = x_ctrl(j) - xw_4(u+1);
        b = y_ctrl(j) + yw(u+1);
        c = x_ctrl(j) - xw_4(u);
        d = y_ctrl(j) + yw(u);
        e = sqrt(a^2+b^2);
        f = sqrt(c^2+d^2);

        g = xw_4(u+1) - xw_4(u);
        h = yw(u) - yw(u+1);

        k = (g*a + h*b)/e - (g*c + h*d)/f;
        l = -1/b * (1 + a/e) + 1/d * (1 + c/f);
        
        
        v_ind_mat(j,u) = v_ind_mat(j,u) + k/((4*pi)*(a*d-c*d)) + l/(4*pi);
    end
end








%% Plot temporal para probar cosas

figure
hold on
plot(xw_4, yw)
plot(xw_le, yw)
plot(xw_te, yw)
plot(xw_ctrl, yw)
plot(x_ctrl,y_ctrl, LineStyle="none", Marker="o")

