clc;
clear;
load('Niculai.mat');
t = Niculai(:, 1);
u = Niculai(:, 2);
y = Niculai(:, 3);
y1 = Niculai(:, 4);

%% Identificare pentru primul semnal 90 grade
clc
close all
figure
plot(t, u);
hold on;
plot(t, y)
hold off
title('Datele primite');
grid;
legend('u', 'y');

ik1 = 22;
ik2 = 81;

k = mean(y(ik1 : ik2)) / mean(u(ik1 : ik2));

%indici pentru w calculat pe semnalu de iesire
in1 = 296;
in2 = 301;
in3 = 321;


%Calculam perioada unde se afla Wn
T = t(in3) - t(in2);
% Calculam defazajul in secunde
deltaT = t(in2) - t(in1);
% Transformam in grade
phi = (360 / T) * deltaT;
% Daca phi = 90 grade, atunci calculam Wn
Wn = (2*pi) / T;

ir1 = 211;
ir2 = 236;

% Calculam perioada in zona de rezonanta
Tr = t(ir2) - t(ir1);
% calculam omega la rezonanta
Wr = 2*pi / Tr;
zeta = sqrt(((Wr / Wn)^2 - 1) / -2);

H1Experimental = tf(k*Wn^2, [1 2*zeta*Wn Wn^2]);

% Crearea matricilor corespunzatoare sistemelor de tip Spatiul Starilor
% Ne folosim de FCOb
A = [0, 1; -Wn^2, -2*zeta*Wn];
B = [0; k*Wn^2];
C = [1, 0];
D = [0];

% Crearea sistemului propriu-zis
sys = ss(A, B, C, D);

% Simularea sistemului pe intrarea u si conditi initiale
ySim = lsim(sys, u, t, [y(1), (y(2) - y(1)) / (t(2) - t(1))]);

% Plotam rezultatul
figure
plot(t, u);
hold on;
plot(t, y);
plot(t, ySim);
hold off
title('Sistemul Simulat');
grid;
legend('u', 'y', 'ySimulat');

J = 1 / sqrt(length(t)) * norm(y-ySim);
empn = norm(y-ySim) / norm(y-mean(y)) * 100;
fprintf("\nJ = %f \n", J)
fprintf("Eroare = %f%% \n", empn)
figure
bode(H1Experimental)

%% Identificare Parametrica a semnalului I
clc
close all
% Timp esantionare
Te = t(2) - t(1);
% Date Identificare
dateIdentificare = iddata(y, u, Te);

% Pentru validarea gradului pentru polinomul B
n4sid(dateIdentificare, 1:5);

%ARMAX - validare prin autocorelatie
Harmax = armax(dateIdentificare, [2 1 2 1]);

figure
compare(dateIdentificare, Harmax);

figure
resid(dateIdentificare, Harmax);

% Trecerea din discret in continuu
HarmaxZ = tf(Harmax);
HarmaxC = d2c(HarmaxZ, "zoh");

%OE - validare prin intercorelatie
Heo = oe(dateIdentificare, [2 2 1]);
figure
compare(dateIdentificare, Heo);

figure
resid(dateIdentificare, Heo);

% Trecerea din discret in continuu
HoeZ = tf(Heo);
HoeC = d2c(HoeZ, "zoh");
rlocus(HoeC)

%% Identificare neparametrica pentru al 2-lea semnal
clc
close all
figure
plot(t, u);
hold on;
plot(t, y1)
hold off
title('Datele primite');
grid;
legend('u', 'y');

i1 = 207; % u
i2 = 209; % y

i3 = 220; % u
i4 = 222; % y

% Calculam varful maxim de magnitudine
Mr = (y1(i2) - y1(i4)) / (u(i1) - u(i3));
zeta = sqrt((Mr - sqrt(Mr^2 - 1)) / 2 / Mr);
Tr = (t(i4) - t(i2)) * 2;
Wr = 2*pi / Tr;
Wn = Wr / sqrt(1 - 2*zeta^2);

ik1 = 22;
ik2 = 81;
k = mean(y1(ik1 : ik2)) / mean(u(ik1 : ik2));

%pentru a afla zeroul, implementam o functie de transfer fara zero, pentru
%a putea calcula polii functiei
%ne folosim pe polii functiei pentru a putea calcula faza sistemului in CFI
HFaraZero = tf(Wn^2*k, [1 2*zeta*Wn Wn^2]);

icfi1 = 988;
icfi2 = 990;
icfi3 = 995;

%calculam perioada la Capatul de Frecventa Inalta
TCFI = t(icfi3) - t(icfi1);
%calculam defazajul la Capatul de Frecventa Inalta
deltaT = t(icfi2) - t(icfi1);
%transformam in grade
phi = (360 / TCFI) * deltaT;
%calculam omega la Capatul de Frecventa Inalta
WCFI = (2*pi) / TCFI;

fprintf("Faza din H(jw), cand w = %f [rad/s], este egala cu -%f grade. \n", WCFI, phi)

% dupa calcularea celor necesare, observam ca al 2 -lea semnal la pulsatia
% de 8.9760e+03 [rad/s] calculata mai sus, avem un defazaj de 102.8571
% [grade]
% Ne folosim de faza sistemului in domeniu frecvential H(jw) = -102.8571
% cu W = 8.9760e+03 [rad/s], pentru a afla Tz
Tz = tand(-phi + atand((WCFI - Wn*sqrt(1 - zeta^2)) / (Wn * zeta)) + atand((WCFI + Wn*sqrt(1 - zeta^2)) / (Wn * zeta))) / WCFI;

HCuZero = tf(k*Wn^2*[Tz 1], [1 2*zeta*Wn Wn^2]);

% Ne folosim de FCO
A = [-2*zeta*Wn 1; -Wn^2 0];
B = [k*Wn^2*Tz; k*Wn^2];
C = [1 0];
D = 0;

sys = ss(A, B, C, D);
ySim = lsim(sys, u, t, [y1(1), (y1(2) - y1(1)) / (t(2) - t(1))]);

figure
plot(t, u);
hold on;
plot(t, y1);
hold on;
plot(t, ySim);
hold off
title('Sistemul Simulat');
grid;
legend('u', 'y1', 'y1Simulat');

J = 1 / sqrt(length(t)) * norm(y1-ySim);
empn = norm(y1-ySim) / norm(y1-mean(y1)) * 100;
fprintf("\nJ = %f \n", J)
fprintf("Eroare = %f%% \n", empn)

%% Identifcare parametrica pentru semnalul al 2-lea

% Timp esantionare
Te = t(2) - t(1);
% Date Identificare
dateIdentificare = iddata(y1, u, Te);

% Pentru validarea gradului pentru polinomul B
n4sid(dateIdentificare, 1:5);

%ARMAX - validare prin autocorelatie
Harmax = armax(dateIdentificare, [2 2 10 1]);
Harmax = pem(dateIdentificare, Harmax);

figure
compare(dateIdentificare, Harmax);

figure
resid(dateIdentificare, Harmax);

% Trecerea din discret in continuu
HarmaxZ = tf(Harmax);
HarmaxC = d2c(HarmaxZ, "zoh")

%OE - validare prin intercorelatie
Hoe = oe(dateIdentificare, [2 2 1]);
figure
compare(dateIdentificare, Hoe);

figure
resid(dateIdentificare, Hoe);

HoeZ = tf(Hoe);
HoeC = d2c(HoeZ, "zoh")

%% Identificare neparametrica pentru primul semnal
clc
figure
plot(t, u);
hold on;
plot(t, y)
hold off
title('Datele primite');
grid;
legend('u', 'y');

i1 = 207; % u
i2 = 211; % y

i3 = 220; % u
i4 = 224; % y

% Calculam varful maxim de magnitudine
Mr = (y(i2) - y(i4)) / (u(i1) - u(i3));
Tr = (t(i4) - t(i2)) * 2;
zeta = sqrt((Mr - sqrt(Mr^2 - 1)) / 2 / Mr);
Wr = 2*pi / Tr;
Wn = Wr / sqrt(1 - 2*zeta^2);


ik1 = 22;
ik2 = 81;
k = mean(y(ik1 : ik2)) / mean(u(ik1 : ik2));

% Crearea matricilor corespunzatoare sistemelor de tip Spatiul Starilor
% Ne folosim de FCOb
A=[0, 1; -Wn^2, -2*zeta*Wn];
B=[0; k*Wn^2];
C=[1, 0];
D=[0];

% Crearea sistemului propriu-zis
sys = ss(A, B, C, D);

% Simularea sistemului pe intrarea u si conditi initiale
ySim = lsim(sys, u, t, [y(1), (y(2) - y(1)) / (t(2) - t(1))]);

figure
plot(t, u);
hold on;
plot(t, y);
plot(t, ySim);
hold off
title('Sistemul Simulat');
grid;
legend('u', 'y', 'ySimulat');

J = 1 / sqrt(length(t)) * norm(y-ySim);
empn = norm(y-ySim) / norm(y-mean(y)) * 100;
fprintf("\nJ = %f \n", J)
fprintf("Eroare = %f%% \n", empn)