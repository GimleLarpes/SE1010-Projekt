%Parametrar

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
b_b = 0.28;
b_d = 0.0;
D = 0.028;
d = 0.6 * D;

b2 = L/2 - b_b - b1;



y_acc = 0;
gamma = 0.25; % sätt gamma = 0 om rakkörning
gamma_max = 1/v * sqrt((R*L*g)/(2*h));

if gamma == 0
    v = v;
else
    v = gamma * v; % Aktivera för kurvtagning
end


F_L = 1/2 * rho_luft * c * A * v^2;
F_D = F_L + m*y_acc;
F_b = 0*(-F_D * r) / (2*r_b); %breaking force, demands F_K = 0
F_K = F_D * r / r_d; %Chain force, demands F_b = 0

V_b = (F_L * (h + h_1) + m*y_acc*h + m*g*d_f) / (d_f + d_b);

V_by = 1/2 * V_b * (1 + (2*(v)^2*h) / (L*g*R));
V_bi = 1/2 * V_b * (1 - (2*(v)^2*h) / (L*g*R));

H_by = V_by/V_b * (m*(v)^2*d_f) / (R*(d_f + d_b));
H_bi = V_bi/V_b * (m*(v)^2*d_f) / (R*(d_f + d_b));

R_ix = H_bi + H_by
R_yy = (((F_D/2)*b1 + (F_b *b2)/sqrt(2) - F_K * (L/2 - b1) + (F_b/sqrt(2)) * (L/2 - b1 + b_b) - F_D/2 * (L - b1)) / (L - 2*b1))
R_iy = sqrt(2) * F_b - (F_K + F_D + R_yy) 
R_yz = (H_bi * r + V_bi * b1 - (F_b *b2)/sqrt(2) - (F_b / sqrt(2)) * (b2 + 2* b_b) - V_by * (L-b1) + H_by * r ) /(L-2*b1)
R_iz = -sqrt(2) * F_b - V_bi - V_by - R_yz 
