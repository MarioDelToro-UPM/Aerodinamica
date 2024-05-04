%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Vortex Lattice Method Wing Solver                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% Wing deffinition variables

NACA_perf = 2400;
AR = 4; %Aspect Ratio
sup = 1; %Surface area
taper = 0.5;
sweep = deg2rad(45);
tors = deg2rad(-4);


%% Flight deffinition variables
u_inf = 1;
AoA = deg2rad(5);

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
cr = sup/(semispan*(1+taper)); %chord at root
ct = cr * taper; %chord at wingtip

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


xw_4 = zeros(Ns+1,1); %x at 1/4 chord
chord = zeros(Ns+1,1);


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

x_ctrl = zeros(Ns, 1);
y_ctrl = zeros(Ns, 1);
tors_ctrl = zeros(Ns, 1);
chord_ctrl = zeros(Ns, 1);

%Definition of control points for induced velocity calculations
for i=1:(Ns)
    y_ctrl(i) = (yw(i)+yw(i+1))/2;
    x_ctrl(i) = cr/4 + y_ctrl(i)*tan(sweep) + 0.5*(cr-(cr-ct)/semispan * y_ctrl(i));
    tors_ctrl(i) = y_ctrl(i)*tors/semispan;
    chord_ctrl(i) = (chord(i)+chord(i+1))/(2*Nc);
end

c_geo_mean = cr * (1+taper)/2;
c_aero_mean = 2/3 * cr * (1+taper+taper^2)/(1+taper);

%Aerodynamic center coordinates
y_aero_ctr = semispan * (1+2*taper)/(3*(1+taper));
x_aero_ctr = 0.25*cr + tan(sweep) * y_aero_ctr;
%CUIDADO NO CUADRA CON LA VALIDACIÓN


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

v_ind_mat = zeros(Ns, Ns);
w_i_ctrl = zeros(1,Ns);


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

        v_ind_mat(j,u) = v_ind_mat(j,u) + k/((4*pi)*(a*d-c*b)) + l/(4*pi);

        %calculations for left semiplane-induced speed

        a = x_ctrl(j) - xw_4(u+1);
        b = y_ctrl(j) + yw(u+1);
        c = x_ctrl(j) - xw_4(u);
        d = y_ctrl(j) + yw(u);
        e = sqrt(a^2+b^2);
        f = sqrt(c^2+d^2);

        g = xw_4(u) - xw_4(u+1);
        h = yw(u+1) - yw(u);

        k = (g*a + h*b)/e - (g*c + h*d)/f;
        l = -1/b * (1 + a/e) + 1/d * (1 + c/f);
        
        v_ind_mat(j,u) = v_ind_mat(j,u) + k/((4*pi)*(a*d-c*b)) + l/(4*pi);
    end
end

b = zeros(Ns, 1);

%b vector deffinition, angle of each panel measured against flight velocity 

for j=1:Ns
    b(j) = AoA+tors_ctrl(j); 
end

b = -b;


%circulation generated by each lifting line in each panel

gamma = linsolve(v_ind_mat, b);

gamma_real = gamma * u_inf;


%% Compute each panel as a lifting body

cl_panel(:) = 2 * (gamma./chord_ctrl);

%poner cuando se haga para todos los paneles el cl_perf

c_L = 0;

for j=1:Ns
    c_L = c_L + cl_panel(j) * (yw(j+1)-yw(j))*chord_ctrl(j);
end

c_L = 2*c_L/sup;

%% Arodynamic characteristics
cmoy = zeros(Ns, 1);

for j = 1:Ns
    cmoy(j) = -cl_panel(j) * (x_ctrl(j)-chord_ctrl(j)/2)/chord_ctrl(j);
end

CMoy =0;
for j=1:Ns
    CMoy = CMoy + cmoy(j)*((chord(j)+chord(j+1))/2)*(yw(j+1)-yw(j))*chord_ctrl(j);
end 

CMoy = (2/(sup*c_aero_mean)) * CMoy;

CMca = CMoy + c_L*x_aero_ctr/c_aero_mean;
%REVISAR, VALIDACION NO CUADRA
%% RESISTENCIA INDUCIDA DE TREFT
w_drag_d = zeros(Ns,Ns);
w_drag_i = zeros(Ns,Ns);

for u=1:Ns
    for j=1:Ns
        w_drag_d(u,j)= 1/(2*pi*(y_ctrl(j)-yw(u+1)));
        w_drag_i(u,j)= -1/(2*pi*(y_ctrl(j)+yw(u+1)));
    end 
end
w_drag = w_drag_i + w_drag_d;

gamma_neta = zeros(Ns,1);
for j=1:Ns-1
    gamma_neta(j) = gamma(j)-gamma(j+1);
end 
gamma_neta(Ns) = gamma(Ns);

w_inf = gamma_neta'*w_drag;
alpha_inf = w_inf/u_inf;

CdiT = -cl_panel.*alpha_inf/2;

CDi =0;
for j=1:Ns
CDi = CDi + 2/sup * CdiT(j)*(yw(j+1)-yw(j))*chord_ctrl(j);
end
%% Plot temporal para probar cosas

figure
hold on
plot(xw_4, yw)
plot(xw_le, yw)
plot(xw_te, yw)
plot(xw_ctrl, yw)
plot(x_ctrl,y_ctrl, LineStyle="none", Marker="o")

