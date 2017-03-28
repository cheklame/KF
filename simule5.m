clear all;
figure(2);clf


fe = 30; % fr�quence d'�chantillonnage
temps = 10; % temps total de mesure en seconde
t=(0:1/fe:temps-1/fe);

 
Tx = 100/fe;
Ty = 200/fe;    
a = 4;
axy = 0.1;
ayx = 0.01;
ax = 0;
ay = 0;
f = 0;   
w = 2*pi*f;
teta = atan(Ty/Tx);


vec_etat = zeros(temps*fe, 5); % [x ; y; vx, vy,D]
vec_etat_x = zeros(temps*fe, 1); % [x ; y; vx, vy,D]

vec_etat(1,3) = 1;
vec_etat(1,4) = 1;

   
i = 2;
 
sig1 = 0.1;  %noise 1 : x,y
sig2 = 0.2;  %noise 2 : speed

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

%%---mesurement with noise-----------


   Y(:,1) = vec_etat(:,1) + sig1;
   Y(:,2) = vec_etat(:,2) + sig1;
   Y(:,3) = vec_etat(:,3) + sig2;
   Y(:,4) = vec_etat(:,4) + sig2;


%%-----Filtre de Kalman-------------


%%%---------Initialisation de matrices---------
%R =[];

A = [1 ayx Tx 0; axy 1 0 Ty; ax 0 0 1; 0 ay 0 1];
C = [ 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0; 0 0 0 1];           %H = [1 0 1; 0 1 0];
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
R(1,1) = sig1^2;
R(2,2) = sig1^2;
R(3,3) = sig2^2;
R(4,4) = sig2^2;

P = zeros(4,4);
j = 1;

for t1 = 0:1/fe:temps-2/fe
    
%%%-----phase de prediction

Xp = A*X;
P_pred = A*P*A' + Q;

%%-----phase de mise � jour

K = P_pred * C' * inv(C * P_pred * C' + R);
P = P_pred - K * C * P_pred;
X_est(j,:) = Xp + K * (Y(j,:)' - C * Xp);

j = j+1;

end

hold on;
plot(X_est(:,1),X_est(:,2),'-r.');
