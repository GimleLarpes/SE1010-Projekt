%Parametrar

clear all close all clc

v = 120/3.6;
R = 14;
m = 140;
c = 0.3;
A = 0.5;
d_f = 0.8;
d_b = 0.3;
h = 0.44;
h_1 = 0.2;
L = 1.1;
b1 = 0.16;
rho_luft = 1.21;
g = 9.81;
d_h = 0.33;
r_d = 0.11;
r_b = 0.09;
r = 0.165;
b_b = 0.28; % 0 < b_b < L/2 - b1
b_d = 0.0;
D = 0.041; %0.028
d = 0.6 * D;

step_b_d = 100;

for iter = 1:step_b_d

    b_b = b1 + (L/2-b1)/step_b_d*iter;
    b_b = 0.31 % vald b_b


b2 = L/2 - b_b - b1;



y_acc = -10;
gamma = 0.25*0; % sätt gamma = 0 om rakkörning
gamma_max = 1/v * sqrt((R*L*g)/(2*h));

if gamma == 0
    v = v;
    curve_stopper = 0;
else
    v = gamma * v; % Aktivera för kurvtagning
    curve_stopper = 1;
end


F_L = 1/2 * rho_luft * c * A * v^2;
F_D = F_L + m*y_acc;

if y_acc >= 0
    f1 = 0;
    f2 = 1;
else
    f1 = 1;
    f2 = 0;
end

F_b = f1*(-F_D * r) / (2*r_b); %breaking force, demands F_K = 0
F_K = f2*F_D * r / r_d; %Chain force, demands F_b = 0

V_b = (F_L * (h + h_1) + m*y_acc*h + m*g*d_f) / (d_f + d_b);

V_by = 1/2 * V_b * (1 + curve_stopper * (2*(v)^2*h) / (L*g*R));
V_bi = 1/2 * V_b * (1 - curve_stopper * (2*(v)^2*h) / (L*g*R));

H_by = V_by/V_b * (m*(v)^2*d_f) / (R*(d_f + d_b)) * curve_stopper;
H_bi = V_bi/V_b * (m*(v)^2*d_f) / (R*(d_f + d_b)) * curve_stopper;

R_ix = H_bi + H_by
R_yy = (((F_D/2)*b1 + (F_b *b2)/sqrt(2) - F_K * (L/2 - b1) + (F_b/sqrt(2)) * (L/2 - b1 + b_b) - F_D/2 * (L - b1)) / (L - 2*b1))
R_iy = sqrt(2) * F_b - (F_K + F_D + R_yy) 
R_yz = (H_bi * r + V_bi * b1 - (F_b *b2)/sqrt(2) - (F_b / sqrt(2)) * (b2 + 2* b_b) - V_by * (L-b1) + H_by * r ) /(L-2*b1)
R_iz = -sqrt(2) * F_b - V_bi - V_by - R_yz

Rvek(iter) = R_yy


% Kod för framtaging av moment och krafter för godtyckligt snitt
x = linspace(0, L, 100); % diskretisering av balk

for i = 1:length(x)

    xi = x(i);

    % Låt enbart krafter aktiveras då de skulle tas hänsyn till då ett
    % moment och tvärkraftsdiagram (snittning) konstrueras från en upplagd balk.
    % Detta gör att "aktiveringen" beror på xi (var på balken i x-led)


    if xi <= 0 %balkens vänsterände - ändrat!
        F_D_1 = 0;
        V_bi_ = 0;
        H_bi_ = 0;
    else
        F_D_1 = F_D;
        V_bi_ = V_bi;
        H_bi_ = H_bi;
    end

    if xi < b1
        R_ix_ = 0;
        R_iy_ = 0;
        R_iz_ = 0;
    else
        R_ix_ = R_ix;
        R_iy_ = R_iy;
        R_iz_ = R_iz;
    end

    if xi < L/2 - b_b
        F_b_1 = 0;
    else
        F_b_1 = F_b;
    end

    if xi < L/2
        F_K_ = 0;
    else
        F_K_ = F_K;
    end

    if xi < L/2 + b_b
        F_b_2 = 0;
    else
        F_b_2 = F_b;
    end

    if xi < L - b1
        R_yy_ = 0;
        R_yz_ = 0;
    else
        R_yy_ = R_yy;
        R_yz_ = R_yz;
    end

    if xi >= L % balkens högerände - ändrat!
        F_D_2 = F_D;
        V_by_ = V_by;
        H_by_ = H_by;
    else
        F_D_2 = 0;
        V_by_ = 0;
        H_by_ = 0;  
    end

    if (b1 < xi) && (xi < L-b1)
        d_ = D;
    else
        d_ = d;
    end

    % Dessa ekvationer härldes från position x=L dvs motsvarande balkens
    % högerände men gäller för godtyckligt x eftersom kratferna nu enbart
    % aktiveras om det gällande snittet har passerat det x värde för vilket
    % kraften befinner sig till vänster om snittet. Alla krafter till höger
    % om snittet blir således 0 medan de till vänster om snittet har sitt
    % "faktiska värde"

    % Tvärkraft z
    T_z(i) = -(V_bi_ + R_iz_ + 1/sqrt(2) * F_b_1 + 1/sqrt(2) * F_b_2 + R_yz_ + V_by_);

    % Moment y
    M_y(i) = - H_bi_ * r - V_bi * xi - R_iz_ * (xi - b1) - 1/sqrt(2) * F_b_1 * (xi-(L/2 - b_b)) - 1/sqrt(2) * F_b_2 * (xi-(L/2 + b_b)) - R_yz_ * (xi-(L - b1)) - H_by_*r;

    % Tvärkraft y
    T_y(i) = - F_D_1 /2 - R_iy_ + + 1/sqrt(2) * F_b_1 - F_K_ + 1/sqrt(2) * F_b_2 - R_yy_ - F_D_2 /2;

    % Moment z
    M_z(i) = - F_D_1 /2 * xi - R_iy_ * (xi - b1) + 1/sqrt(2) * F_b_1 * (xi-(L/2 - b_b)) - F_K_ * (xi-L/2) + 1/sqrt(2) * F_b_2 * (xi-(L/2 + b_b)) - R_yy_ * (xi-(L - b1)) - F_D_2 /2 * (xi - L);

    % Moment x (reaktionsmoment på torsion)
    M_x(i) = - F_D_1 /2 * r - F_b_1 *r_b + F_K_ * r_d - F_b_2 *r_b - F_D_2 /2 * r;

    % Normalkaft x
    N_x(i) = H_bi_ - R_ix_ + H_by_;

    dvek(i) = d_;

    A = pi*d_^2 / 4; % tvärsnittsarea
    I = pi * (d_/2)^4 /4; % areatröghetsmoment
    zeta1 = d_/2;
    zeta2 = -d_/2;

    sigma1(i) = N_x(i)/A + sqrt(M_y(i)^2+M_z(i)^2)/I*zeta1; % normalspänning
    sigma2(i) = N_x(i)/A + sqrt(M_y(i)^2+M_z(i)^2)/I*zeta2; % normlspänning motsatt sida

    ww = pi/16 * d_^3; % vridmotstånd?

    tau_x(i) = M_x(i)/ww; % vridskjuvpänning

    tau_y(i) = T_y(i)/A; % skjuvning y

    tau_z(i) = T_z(i)/A; % skjuvning z

    tau_yz(i) = sqrt(T_y(i)^2+T_z(i)^2)/A; % maximal skjuving i yz-planet

    von_mises(i) = sqrt ( max(sigma1(i),sigma2(i))^2 + 3*tau_x(i)^2 + 3*tau_y(i)^2 + 3*tau_z(i)^2 );

    i = i+1;

end

max_von_mises(iter) = max(von_mises)

end


volume_axle = pi*d^2 / 4 * 2 * b1 + pi * D^2 / 4 * (L-2*b1);
rho_axle = 7850;
mass_axle = rho_axle * volume_axle;

% Plot för moment (Mx, My, Mz)
figure(1)
plot(x, M_x, 'b', 'LineWidth', 1.5)
hold on
plot(x, M_y, 'r', 'LineWidth', 1.5)
plot(x, M_z, 'g', 'LineWidth', 1.5)
hold off
title('Moment längs x, y och z axlarna')
xlabel('x [m]')
ylabel('Moment [Nm]')
legend('M_x', 'M_y', 'M_z')
grid on

% Plot för krafter (Nx, Ty, Tz)
figure(2)
plot(x, N_x, 'b', 'LineWidth', 1.5)
hold on
plot(x, T_y, 'r', 'LineWidth', 1.5)
plot(x, T_z, 'g', 'LineWidth', 1.5)
hold off
title('Krafter längs x, y och z axlarna')
xlabel('x [m]')
ylabel('Krafter [N]')
legend('N_x', 'T_y', 'T_z')
grid on

% Plot för normalspänningar (sigma1, sigma2)
figure(3)
plot(x, sigma1, 'b', 'LineWidth', 1.5)
hold on
plot(x, sigma2, 'r', 'LineWidth', 1.5)
hold off
title('Normalspänningar (\sigma_1 och \sigma_2)')
xlabel('x [m]')
ylabel('Spänning [Pa]')
legend('\sigma_1', '\sigma_2')
grid on

% Plot för skjuvspänningar (tau_yz, tau_y, tau_z)
figure(4)
plot(x, tau_yz, 'b', 'LineWidth', 1.5)
hold on
plot(x, tau_y, 'r', 'LineWidth', 1.5)
plot(x, tau_z, 'g', 'LineWidth', 1.5)
hold off
title('Skjuvspänningar (\tau_{yz}, \tau_{y}, \tau_{z})')
xlabel('x [m]')
ylabel('Skjuvspänning [Pa]')
legend('\tau_{yz}', '\tau_{y}', '\tau_{z}')
grid on

% Plot för skjuvspänning längs x (tau_x)
figure(5)
plot(x, tau_x, 'b', 'LineWidth', 1.5)
title('Skjuvspänning längs x (\tau_{x})')
xlabel('x [m]')
ylabel('Skjuvspänning [Pa]')
legend('\tau_{x}')
grid on

% Plot för von mises 
figure(6)
plot(x, von_mises, 'b', 'LineWidth', 1.5)
title('Effektivspänning von mises')
xlabel('x [m]')
ylabel('Spänning [Pa]')
legend('VM')
grid on

figure(7)
plot(linspace(b1, L/2-b1, step_b_d), max_von_mises, 'LineWidth', 1.5)
title('Maximal von Mises-spänning som funktion av b_b')
xlabel('b_b [m]')
ylabel('Maximal von Mises-spänning [Pa]')
grid on

totalmax_von_mises=max(max_von_mises(2:end-1))

ns = 2.8;
sigma_max_till = 310*10^6 / ns
