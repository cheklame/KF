clear all;
figure(2);clf


fe = 30; % fréquence d'échantillonnage
temps = 10; % temps total de mesure en seconde
t=(0:1/fe:temps-1/fe);

 
Tx = 100/fe;
Ty = 200/fe;    
a = 4;
axy = 0.1;
ayx = 0.01;
ax = 0;
ay = 0;
f = 5;   
w = 2*pi*f;
teta = atan(Ty/Tx);


vec_etat = zeros(temps*fe, 5); % [x ; y; vx, vy,D]
vec_etat_x = zeros(temps*fe, 1); % [x ; y; vx, vy,D]

vec_etat(1,3) = 1;
vec_etat(1,4) = 1;


i = 2;

%%--real trajectory of the drone
for t = 0:1/fe:temps-2/fe
    k = i - 1;
    vec_etat(i,1) = vec_etat(k,1) + ayx* vec_etat(k,2) + vec_etat(k,3) * Tx  + vec_etat(k,5)*sin(teta);
    vec_etat(i,2) = axy* vec_etat(k,1) + vec_etat(k,2) + vec_etat(k,4) * Ty  + vec_etat(k,5)*cos(teta);
    vec_etat(i,3) = ax*vec_etat(k,1) + vec_etat(k,3);        
    vec_etat(i,4) = ay*vec_etat(k,1) + vec_etat(k,4);
    vec_etat(i,5) = a*sin(w*t) * vec_etat(k,3);

    i = i+1;
end
s = size(vec_etat);
taille = s(1);
t1=(0:1/fe:(temps/fe)*(taille-1)/10); 
plot(vec_etat(:,1),vec_etat(:,2),'-g.');
%plot(t1,vec_etat(:,1),'-g.');
hold on;

%plot(t,M2(:,2),'-g.');
%plot(vec_etat(:,1),vec_etat(:,2),'-r.');
axis([0 100 -100 100])
hold on;
%xl = a*sin(w*t1(:));
%plot(t1,xl,'-g');


%%%%-----Filtre de Kalman-------------


%%%---------Initialisation de matrices---------
%R =[];
sig1 = 0.1;
sig2 = 0.2;
A = [1 ayx Tx 0; axy 1 0 Ty; ax 0 0 1; 0 ay 0 1];
C = [ 1 0 0 0 ; 0 1 0 0 ;0 0 1 0; 0 0 0 1];           %H = [1 0 1; 0 1 0];
X = zeros(4,1);
X(3) = 2; 
X(4) = 1;
X_est = zeros(temps*fe, 4);
X_est(1,1) = X(1);
X_est(1,2) = X(2);
X_est(1,3) = X(3);
X_est(1,4) = X(4);

Q = zeros(4,4);
Q(1,1) = sig1;
Q(2,2) = sig1;
Q(3,3) = sig2;
Q(4,4) = sig2;

R = zeros(4,4);
R(1,1) = sig1;
R(2,2) = sig1;
R(3,3) = sig2;
R(4,4) = sig2;

P = zeros(4,4);
i = 1;

for t = 0:1/fe:temps-2/fe
    
%%%-----phase de prediction

Xp = A*X;
P_pred = A*P*A' + Q;

%%-----phase de mise à jour

K = P_pred * C' * inv(C * P_pred * C' + R);
P = P_pred - K * C * P_pred;
X = Xp + K * (vec_etat(i,:)' - C * Xp);


i = i+1;

end


