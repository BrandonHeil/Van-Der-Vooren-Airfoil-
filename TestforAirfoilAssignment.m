%% Part 1 - Airfoil Geometry

clear; clc; close

% epsilon = [0.0222,.04725, 0.0721,];                                       % Defining eps. vector so for loop can change it for each airfoil
% radius = [0.551, 0.563, 0.574];                                           % Defining radius vector so for loop can change it for each airfoil
epsilon = [0.0222,0.0721,0.04725,]; % 15% last
radius = [0.551,0.574,0.563];

theta2 = zeros(1,360);                                                      % Preallocating theta2's vector size saves computing power/time
j = 1:3;
V = 1;
for i = 1:3                                                                  % index variable for e, a, and airfoil figures 1-3
    ji=j(i)+3;                                                                  % index variable for Cp vs. x figures 4-6
    k = 1.889;                                                              % T.E. angle parameter
    c = 2;                                                                  % Chord length
    tau = 20*pi/180;                                                        % T.E angle deg to rad
    a = radius(i);                                                          % Radius of circle in f-plane
    e = epsilon(i);                                                         % Thickness parameter
    l = 0.5*((a*2^k)/((1+e)^(k-1)));                                        % Half Chord length?
    theta= (1:1:360)*pi/180;                                                % Angle around circle in f-plane?

    r1 = sqrt(((a.*cos(theta) - a).^2)+((a^2).*(sin(theta)).^2));
    r2 = sqrt(((a.*cos(theta) - e*a).^2)+((a^2).*(sin(theta)).^2));
    r1(1,360) = 0;

    theta1 = atan((a.*sin(theta))./(a.*cos(theta) - a)) + pi;
    theta1(1,360)= 0.5*pi;

    % Theta2: n = 0 0-90deg, n = 1 90-270deg, n = 2 270-360deg
    theta2(1,1:89) = atan((a.*sin(theta(1,1:89)))./(a.*cos(theta(1,1:89)) - e*a)) + 0*pi;
    theta2(1,90:269) = atan((a.*sin(theta(1,90:269)))./(a.*cos(theta(1,90:269)) - e*a)) + 1*pi;
    theta2(1,270:360) = atan((a.*sin(theta(1,270:360)))./(a.*cos(theta(1,270:360)) - e*a)) + 2*pi;

    if any(theta2 < 0)                                                      % When denom. of theta2 atan term goes negative n should be 1
        index = find(theta2 < 0);
        theta2(index) = theta2(index) + pi;
    end

    if any(theta2 > 2*pi)                                                   % Accounting for when denom. in atan term goes negative and changing n back to 1
        index2 = find(theta2 > 2*pi);
        theta2(index2) = theta2(index2) - pi;
    end
                                                                            % note the atan2 function corrects this problem automatically and can save
                                                                            % computing time 
    num1 = (r1.^k)./(r2.^(k-1));
    x = num1.*(cos(k.*theta1).*cos((k-1).*theta2) + sin(k.*theta1).*sin((k-1).*theta2));
    x = x + 1;
    z = num1.*(sin(k.*theta1).*cos((k-1).*theta2) - cos(k.*theta1).*sin((k-1).*theta2));

figure(i)
plot(x,z,'k','linewidth',1.2); grid on; 
set(gca,'FontSize',13)
axis equal;
title '18% Airfoil with T.E. angle 20^{\circ} Chord Length of 2'
xlabel 'Length x'
ylabel 'Height z'
xticks (-1:.2:1);
yticks (-1:.2:1); 


%% Velocity Distribution & Pressure Distribution

    AoA = 0; % degrees
    Qinf = 1; % ft/s
    num2 = (r2.^k)./((r1.^(k-1)) + 5E-5); 

    if any(num2 > 1E12)
        index3 = find(num2 > 1E12);
        num2(index3) = (r2(index3).^k)./(r1(index3).^(k-1) + 0.000005);
    end

    A = cos((k-1).*theta1).*cos((k.*theta2)) + sin((k-1).*theta1).*sin(k.*theta2);
    B = sin((k-1).*theta1).*cos((k.*theta2)) - cos((k-1).*theta1).*sin(k.*theta2);

    D0 = a*(1-k+(k*e));
    D1 = A.*(a.*cos(theta) - D0) - B.*a.*sin(theta);
    D2 = A.*(a.*sin(theta)) + B.*(a.*cos(theta) - D0);

    u = (2.*Qinf.*num2).*((sin(AoA) - sin(AoA-theta))./(D1.^2 + D2.^2)).*(D1.*sin(theta) + D2.*cos(theta));
    w = (-2.*Qinf.*num2).*((sin(AoA) - sin(AoA-theta))./(D1.^2 + D2.^2)).*(D1.*cos(theta) - D2.*sin(theta));

    Cp = 1 - ((u.^2 + w.^2)./(Qinf^2));

figure(ji)
plot(x,Cp,'linewidth', 1.2); grid on;
axis equal;
set(gca,'FontSize',13)
title 'Pressure Coefficient (C_p) vs. Chordwise Span (x)'
xlabel 'Chordwise Span x'
ylabel 'Pressure Coefficient C_p'
xticks (-1:.2:1);
yticks (-1:.2:1);
end
%% Part 1b - C_p vs. x for 18% thickness ratio airfoil with chaning AoA

AngleofAttack = [3, 9, 15]*pi/180;
for q = 1:3
    AoA = AngleofAttack(q);

    u = (2.*Qinf.*num2).*((sin(AoA) - sin(AoA-theta))./(D1.^2 + D2.^2)).*(D1.*sin(theta) + D2.*cos(theta));
    w = (-2.*Qinf.*num2).*((sin(AoA) - sin(AoA-theta))./(D1.^2 + D2.^2)).*(D1.*cos(theta) - D2.*sin(theta));

    Cp = 1 - ((u.^2 + w.^2)./(Qinf^2));

figure(q+6);
plot(x,Cp,'linewidth',1.2); grid on;
set(gca,'FontSize',13)
title 'C_p for 18% Thickness Ratio Airfoil @ X^{\circ} AoA vs. x'
xlabel 'Chordwise Span x'
ylabel 'Pressure Coefficient C_p'
end

%% Part 2 - Panel Method 

N = 8;                                                                      % Number of panels for panel method, USE EVEN only
DELX = 4/N;                                                                 % Change in x between each panel
X = zeros(1,N/2);                                                           % Preallocate size for computing power/time
X(1) = 1;                                                                   % First X in panel method will always be 1
X(2) = x(1) - DELX;                                                         % 2nd X val calc.

for I = 3:1:(N/2 + 1)
    X(I) = X(I-1) - DELX; 
end
X = [X, flip(X(1:N/2),2)];                                                  % Airfoil is symmetrical flip first half values for next half

Z = zeros(1,N+1 - N/2);
Z(1) = 0;                                                                   % First Z will always be zero
for J = 2:(N/2 + 1)
    n = X(J);
    [val,idx] = min(abs(x-n));
    Z(J) = z(idx);
    if any(Z < 0)                                                           % Bandaid on a gunshot wound
        idx2 = find(Z < 0);
        Z(idx2) = Z(idx2)*-1;
    end
end
Z = [Z,-flip(Z(1:N/2),2)];                                                  % Airfoil is symmetrical, flip first half and mult. by neg 1 for bottom half of airfoil
plot(X,Z,'g','linewidth',1.2); grid on
xlabel 'x location'
ylabel 'z location'
title '18% t/c Airfoil with 8 Panels'
axis equal
hold on
% colocation points in the center of each panel
colocx = zeros(1,N);
colocz = zeros(1,N);
slope = zeros(1,N);
thetaN = zeros(1,N);
dels = zeros(1,N);
for Q = 1:N                                                                 % Calcs for mid point for N panels 
    dels(Q) = sqrt((X(Q+1)-X(Q))^2 + (Z(Q+1)-Z(Q))^2);                         % Panel lengths 
    colocx(Q) = ((X(Q+1)+X(Q))/2);
    colocz(Q) = ((Z(Q+1)+Z(Q))/2);
    slope(Q) = ((Z(Q+1)-Z(Q))/(X(Q+1)-X(Q)));                               % slope of each panel
    thetaN(Q) = atan2d(Z(Q+1)-Z(Q),X(Q+1)-X(Q));
    if thetaN(Q) < 0
        thetaN(Q) = 360  + thetaN(Q);
    end
end
scatter(colocx,colocz,'r','d','linewidth',1.2);
scatter(X,Z,'b','o','linewidth',1.2);
hold off

eta = zeros(2,8);
zeta = zeros(2,8);
C = zeros(N,N);
Cbar = zeros(N,N);
slope2 = zeros(N,N);
phi = zeros(N,N);
r = zeros(N,N);

%THETA = [171 173 178 196 343 1.54 6.32 8.94]*pi/180;
%PHIC = [0 352.3681, 354.0871, 358.6707, 4.31436 10.1986 16.2720 89.9542; 172.4598 0 355.8059 1.8678 0.1986 26.2986 89.9542 163.8655; 174.1787 175 0 7.66 23.26 89.9 153.55 169.59; 178 181 187 0 89.9 156 169 175; 184 190 203 269 0 172 178 181; 190 206 269 336 352 0 183 185; 196 269 333 349 358 3.94 0 187; 269 343 349 355 1.31 5.78 7.62 0];
thetaN = thetaN.*(pi/180);
for I = 1:N
    eta(1:2,I) = [colocx(I); colocz(I)];
    zeta(1:2,I) = [colocx(J); colocz(J)];
    for J = 1:N    
        zeta(1:2,J) = [colocx(J); colocz(J)];
        slope2(I,J) = ((colocz(I) - colocz(J))/(colocx(I) - colocx(J)));
        r(I,J) = sqrt((colocx(I) - colocx(J))^2 + (colocz(I) - colocz(J))^2);
        hold on; plot([eta(1,I) zeta(1,J)], [eta(2,I) zeta(2,J)],'k','linewidth',1); axis equal; grid on
        phi(I,J) = pi + atan2((colocz(J) - colocz(I)),(colocx(J) - colocx(I))); 
        if I == J 
            phi(I,J) = 0;
        end
        C(I,J) = ((sin(thetaN(I)-phi(I,J)))/(2*pi*r(I,J)))*dels(J);
        Cbar(I,J) = ((cos(thetaN(I)-phi(I,J)))/(2*pi*r(I,J)))*dels(J);
        if I == J
            C(I,J) = 0;
            Cbar(I,J) = 0;
        end
    end
end

AngleofA = 0;                                                                    % Free stream velocity in +x direction
Am = C;
Bm = zeros(1,N);
Cpnum = zeros(1,N);

for i = 1:N
    for j = 1:N
        Bm(i) = -(sin(thetaN(i)));
        if i ~= j
            Am(i,j) = C(i,j);
        else
            Am(i,j) = 0.5;
        end
    end
end
Bm = [.155 -.11 -.27 .283 .283 -.27 -.11 -.155];
Bm = Bm';
for i = 1:N
    for j = 1:N
    q = inv(Am)*Bm;
    %wtf = - sum(q(i).*C(i,j));
    vt(i) = cos(thetaN(i)) - sum(q(i).*Cbar(i,j));
    Cpnum(i) = 1 - (vt(i)).^2;
    end
end
figure(10)
plot(colocx,Cpnum)

