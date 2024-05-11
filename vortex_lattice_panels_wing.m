%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Vortex Lattice Method Wing Solver                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% Wing deffinition variables

NACA_perf = 3315; %TENEMOS QUE AÑADIR LA PENDIENTE EN 3/4 DE LA CUERDA DONDE TOCA

%Following data is approximately the wing of an Ha.1109, Bf-109G, or
%similar derivatives of the german aircraft
aspect_ratio = 6.1374; %Aspect Ratio
surface = 16.05; %Surface area
Taper = 0.475;
Sweep = deg2rad(0);
torsion = deg2rad(-2);


%% Flight deffinition variables
u_inf = 1;
n_alpha = 5;
alpha = deg2rad(linspace(-8, 8, n_alpha));

%% Geometry deffinition variables

N_chord = 1; %Number of chordwise panels (placeholder, more than 1 WILL break this code)
N_span = 10; %Number of spanwise panels for each semiplane
Distrib = 3; %0 for equidist, 1 for more panels toward the root, 2 for more panels toward the wingtips



[xw_4th, ywing, x_ctr, y_ctr, tors_c, chord_c, x_c, y_c, c_geo_m, c_aero_m, chord]...
    = wing_geo(aspect_ratio, surface, Taper, Sweep, torsion, N_chord, N_span, Distrib);

gamma_real = zeros(n_alpha, N_span);
cL = zeros(n_alpha, 1);
cD = zeros(n_alpha, 1);
cMoy = zeros(n_alpha, 1);
cMca = zeros(n_alpha, 1);
w_ctrl = zeros(n_alpha, N_span);
cl = zeros(n_alpha, N_span);

for i = 1:n_alpha
    [gamma_real(i, :), cL(i), cMoy(i), cMca(i), w_ctrl(i, :), cl(i, :)] = ...
        wing_solve(x_c, y_c, xw_4th, ywing, chord, tors_c, chord_c, ...
        surface, alpha(i), u_inf, c_aero_m, x_ctr);

    cD(i) = induced_drag(surface, y_c, ywing, gamma_real(i), u_inf, cl(i), chord_c);
end

poly_cL = polyfit(alpha, cL, 1);

alpha_0 = -poly_cL(2)/poly_cL(1);
alpha_1 = (1-poly_cL(2))/poly_cL(1);

[~,~,~,~,~,cla] = wing_solve(x_c, y_c, xw_4th, ywing, chord, tors_c, chord_c, ...
        surface, alpha_0, u_inf, c_aero_m, x_ctr);

[~,~,~,~,~,clb] = wing_solve(x_c, y_c, xw_4th, ywing, chord, tors_c, chord_c, ...
        surface, alpha_1, u_inf, c_aero_m, x_ctr);

cla = cla - clb;


%% Plots for the report

colors = [[0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250];...
     [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]];

figure
hold on
for i=1:n_alpha
    alpha_here = num2str(rad2deg(alpha(i)));
    display = strcat('\alpha = ', alpha_here, 'º');
    plot(y_c, cl(i,:), Color=colors(i,:), DisplayName=display)
    plot(-y_c, cl(i,:), Color=colors(i,:), HandleVisibility='off')
end
ylabel('c_l(y)')
legend('show', 'Location', 'northwest')
hold off

clear alpha_here
clear display




figure
hold on
grid on
plot(rad2deg(alpha), cMca, 'DisplayName', 'C_{m, ca}')
plot(rad2deg(alpha), cMoy, 'DisplayName', 'C_{m, 0y}')
%plot(rad2deg(alpha), CM4, 'DisplayName', 'C_{m, 1/4}')
plot(rad2deg(alpha), cL, 'DisplayName', 'C_l')
xlabel('\alpha[º]')
ylabel('C_l, C_{m, ba}, C_{m, ca}')
legend('show', 'Location', 'northwest')
hold off


% figure
% plot(cL, cD)


%%%%%%%%%%%%%%%%% Functions

function [xw_4, yw, x_aero_ctr, y_aero_ctr, tors_ctrl, chord_ctrl, ...
    x_ctrl, y_ctrl, c_geo_mean, c_aero_mean, chord] = wing_geo(AR, sup, taper, sweep, tors, Nc, Ns, Dist)
    % Wing semiplane geometry calculations
    
    %lets assume the center wing profile is located at x,y = 0,0
    
    %Division of the y axis will be made as a distribution of points from said
    %axis, it will NOT be done dividing any of the wing chords
    
    %Position variables are stored for the right semiplane only
    
    semispan = sqrt(AR*sup)/2;
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
        yw = semispan*sin(yw);
    elseif Dist == 3
        yw = linspace(pi, 0, Ns+1);
        yw = semispan/2*(1+cos(yw));
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
    for i=1:Ns
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


    ribs = [xw_le, xw_te];
    y_ribs = [yw', yw'];


    %plots wing figure
    figure
    hold on
    plot(xw_4, yw, Color='[0.85 0.85 0.85]')
    plot(xw_le, yw, Color='black')
    plot(xw_te, yw, Color='black')
    plot(xw_ctrl, yw, Color='[0.85 0.85 0.85]')
    plot(x_ctrl,y_ctrl, LineStyle="none", Marker="o")
    %all the plot commands with a - before the y coordinate are there to
    %plot the port wing
    plot(xw_4, -yw, Color='[0.85 0.85 0.85]')
    plot(xw_le, -yw, Color='black')
    plot(xw_te, -yw, Color='black')
    plot(xw_ctrl, -yw, Color='[0.85 0.85 0.85]')
    plot(x_ctrl, -y_ctrl, LineStyle="none", Marker="o")
    for i=1:size(ribs, 1)
        plot(ribs(i, :), y_ribs(i,:), Color='black')
        plot(ribs(i, :), -y_ribs(i,:), Color='black')
    end
    xlim([0 5])
    ylim([-5 5])

end

function [gamma_real, c_L, CMoy, CMca, w_i_ctrl, cl_panel]...
    = wing_solve(x_ctrl, y_ctrl, xw_4, yw, chord, tors_ctrl, ...
    chord_ctrl, sup, AoA, u_inf, c_aero_mean, x_aero_ctr)
    % Velocity solving for the whole wing


    %Nc = size(x_ctrl, 1);
    Ns = size(x_ctrl, 1);
    
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

    fmax = 0.03;
    x_fmax = 0.3;

    dzdx = (fmax/(1-x_fmax)^2) * (2*x_fmax- 2*(3/4));
    %3/4 point at adimensional chord where the control points are
    %x_fmax fmax are hardcoded from NACA 3315

    for j=1:Ns
        b(j) = AoA+tors_ctrl(j)-dzdx;
    end
    
    b = -b;
    
    
    %circulation generated by each lifting line in each panel
    
    gamma = linsolve(v_ind_mat, b);
    
    gamma_real = gamma * u_inf;

    % Compute each panel as a lifting body

    cl_panel(:) = 2 * (gamma./chord_ctrl);
    
    %poner cuando se haga para todos los paneles el cl_perf
    
    c_L = 0;
    
    for j=1:Ns
        c_L = c_L + cl_panel(j) * (yw(j+1)-yw(j))*chord_ctrl(j);
    end
    
    c_L = 2*c_L/sup;

    % Arodynamic characteristics
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
end

function CDi = induced_drag(sup, y_ctrl, yw, gamma, u_inf, cl_panel, chord_ctrl)

    Ns = size(gamma, 1);

    % RESISTENCIA INDUCIDA DE TREFT
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

end




