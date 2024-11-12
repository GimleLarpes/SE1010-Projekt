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

% Kod för framtaging av moment och krafter för godtyckligt snitt (?)
x = linspace(0, L, 100); % diskretisering av balk

for i = 1:length(x)

    xi = x(i);

    % Låt enbart krafter aktiveras då de skulle tas hänsyn till då ett
    % moment och tvärkraftsdiagram (snittning) konstrueras från en upplagd balk.
    % Detta gör att "aktiveringen" beror på xi (var på balken i x-led).
    % Alla krafter som har "aktiveringsegenskapen" betecknas enligt
    % "F_" ifall den var "F" innan
    % om det står ett index 1 eller 2 efter en krafts underscore innebär det att
    % den är uppdelad i 2. 1 är den till vänster och 2 den till höger
    % t.ex. F_D uppdelas i F_D_1 och F_D_2


    % if x < L tekniskt sett men blir sant för alla iterationer här
    F_D_1 = F_D;
    V_bi_ = V_bi;
    H_bi_ = H_bi;

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

    if xi == L % balkens högerände
        F_D_2 = F_D;
        V_by_ = V_by;
        H_by_ = H_by;
    else
        F_D_2 = 0;
        V_by_ = 0;
        H_by_ = 0;  
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
    M_y(i) = - H_bi_ * r - V_bi * xi - R_iz_ * (xi - b1) - 1/sqrt(2) * F_b_1 * (xi-(L/2 - b_b)) - F_b_2 * (xi-(L/2 + b_b)) - R_yz_ * (xi-(L - b1)) - H_by_*r;

    % Tvärkraft y
    T_y(i) = - F_D_1 /2 - R_iy_ + + 1/sqrt(2) * F_b_1 - F_K_ + 1/sqrt(2) * F_b_2 - R_yy_ - F_D_2 /2;

    % Moment z
    M_z(i) = - F_D_1 /2 * xi - R_iy_ * (xi - b1) + 1/sqrt(2) * F_b_1 * (xi-(L/2 - b_b)) - F_K_ * (xi-L/2) + F_b_2 * (xi-(L/2 + b_b)) - R_yy_ * (xi-(L - b1)) - F_D_2 /2 * (xi - L);

    % Moment x (reaktionsmoment på torsion)
    M_x(i) = - F_D_1 /2 * r - F_b_1 *r_b + F_K_ * r_d - F_b_2 *r_b - F_D_2 /2 * r;

    % Normalkaft x
    N_x(i) = H_bi_ - R_ix_ + H_by_;

    i = i+1;
end

figure (1)
plot(x, M_x)
hold on
plot(x, M_y)
hold on
plot(x, M_z)
hold on

figure(2)
plot(x, N_x)
hold on
plot(x, T_y)
hold on
plot(x, T_z)
hold on
