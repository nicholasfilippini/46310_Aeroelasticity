close all
clear all

%% Interpolation Reading - from Christian's script
thick_prof=[100,60,48,36,30.1,24.1];
[aoa,cl,cd,cm,fs,clinv,clfs]=...
    readairfoildata_stat('cylinder.txt','FFA-W3-600.txt',...
    'FFA-W3-480.txt','FFA-W3-360.txt',...
    'FFA-W3-301.txt','FFA-W3-241.txt');

%% Variable Initialisation
H = 119;            % Tower height
tr=0;               % Tower radius
ls = 7.1;           % Shaft length
r = 89.17;          % Blade radius
w = 0.757;          % Omega
bla = 3;            % Number of blades
wake =1;            % Wake filter condition
stall =0;           % Stall firlter condition
rho = 1.225;        % Density
n = 1;              % Time step
T = 10.3242*2;            % Time
v = 9;              % Velocity(m/s)
W = [0;0;0];        % Induced Wind 
Wint = [0;0;0];     % Inter. Induced Wind
mu = 0.2;           % Wind shear
dt = 0.15;          % Time increment
a = 0*pi/180;       % Yaw angle(Rad)
b = 0*pi/180;       % Tilt angle(Rad)
c = 0*pi/180;       % Cone angle(Rad)

%% Matrix Initialisation
A1 = [1,0,0; 0,cos(a),sin(a);0,-sin(a),cos(a)]; 
A2 = [cos(b),0,-sin(b);0,1,0;sin(b),0,cos(b)];
A3 = [1,0,0;0,1,0;0,0,1];
A12 = A3* A2* A1;
A21 = A12.';
A34 = [cos(c),0,-sin(c);0,1,0;sin(c),0,cos(c)];

% Position matrix
r1 = [H;0;0];
r2 = A21 * [0;0;-ls];

%% Blade Characteristics
bet = [14.50;14.43;12.55;8.89;6.38;4.67;2.89;1.21;-0.13;-1.11;-1.86;-2.08;-2.28;-2.64;-2.95;-3.18;-3.36;-3.43];
rad = [2.80;11.00;16.87;22.96;32.31;41.57;50.41;58.53;65.75;71.97;77.19;78.71;80.14;82.71;84.93;86.83;88.45;89.17];
cho = [5.38;5.45;5.87;6.18;6.02;5.42;4.70;4.00;3.40;2.91;2.54;2.43;2.33;2.13;1.90;1.63;1.18;0.60];
TC  = [100;86.05;61.10;43.04;32.42;27.81;25.32;24.26;24.10;24.10;24.10;24.10;24.10;24.10;24.10;24.10;24.10;24.10];

%% Loop Start
for i = 1:bla               % Number of blades loop
    for j = 1:length(rad)   % Section loop
        for t = 0.15:dt:T   % Time loop
            %% pitch Control

             pitch = 0;

%             pitch = 15 + (5*sin(25*t));

%             if t <100
%                 pitch = 0;
%             elseif (100<=t)<=150
%                 pitch = 2;
%             elseif t>150
%                 pitch = 0;
%             end

            if j < 18 % Tip corection condition
                
                x = rad(j); % Section Radius
                D = [(w*t);(w*t)+(2*pi/3);(w*t)+(4*pi/3)]; % Matrix for Blade angles for each balde (Rad)
                
                % Coordinate matrix calculation
                A23 = [cos(D(i)),sin(D(i)),0;-sin(D(i)),cos(D(i)),0;0,0,1];
                A14 = A34*A23*A12;
                A41 = A14.';

                % Position Matrix Calculation
                r3 = A41 * [x;0;0];
                R = r1 +r2 +r3;
                
                % Postion vector assignment
                Rx(n,i)= R(1);
                Ry(n,i)= R(2);
                Rz(n,i)= R(3);
                
                rp=(((Ry(n,i))^2)+((Rz(n,i))^2))^(0.5); % Magnitude of position vector
                
                Az(n) = w*t*180/pi; % Azimuth Angle calculation
                
                Vs = [0;0;(v*(Rx(n,i)/H)^mu)]; % Wind Shear Calculation
                
                % Tower height correction
                if Rx(n,i)>H
                    rt = tr;
                else
                    rt = 0;
                end

                % Velocity Potential Calculation
                Vr= (Rz(n,i)/rp)*Vs(3)*(1-(rt/rp)^2); % Vr calculation
                Vo= (Ry(n,i)/rp)*Vs(3)*(1+(rt/rp)^2); % V_theta calculation
                l = ((Rz(n,i)/rp)* Vr)+((Ry(n,i)/rp)*Vo); % Potential velocity in z
                k = ((Ry(n,i)/rp)* Vr)-((Rz(n,i)/rp)* Vo); % Potential velocity in y
                Vm= [0;k;l];

                V = A14*Vm; % Velocity in System 4

                Vrel =[0, V(2)+W(2,n)-(w*x*cos(c)), V(3)+W(3,n)]; % Relative Velocity

                phi= atan(Vrel(3)/-Vrel(2)); % Blade angle calculation (Rad)
                PHI(j,n,i)=phi*(180/pi);
                
                alp = (phi*180/pi) - bet(j) - pitch; % Angle of attack calculation (Deg)
                Alpha(j,n,i)=alp;
                
                thick = TC(j); % Thickness of blade assignment
                
                %% Interpolation function
                for k=1:6
                    clthick(k)=interp1(aoa(:,k),cl(:,k),alp);
                    cdthick(k)=interp1(aoa(:,k),cd(:,k),alp);
                    cmthick(k)=interp1(aoa(:,k),cm(:,k),alp);
                    fsthick(k)=interp1(aoa(:,k),fs(:,k),alp);
                    clinvthick(k)=interp1(aoa(:,k),clinv(:,k),alp);
                    clfsthick(k)=interp1(aoa(:,k),clfs(:,k),alp);
                end
                clift=interp1(thick_prof,clthick(:),thick(:)); % Coeff. Lift    
                cdrag=interp1(thick_prof,cdthick(:),thick(:)); % Coeff. Drag
                cmoment=interp1(thick_prof,cmthick(:),thick(:)); % Coeff. Moment
                fst=interp1(thick_prof,fsthick(:),thick(:)); % Seperation function
                CliftInv=interp1(thick_prof,clinvthick(:),thick(:)); % Coeff. Lift for inviscid flow
                Cliftfs=interp1(thick_prof,clfsthick(:),thick(:)); % Coeff. Lift at fully seperated flow region

                %% Stall filter
                if stall==1 % Turning on stall model
                    if n==1
                        FS = fst; % Seperation function initial value
                    else
                        tau = 4*cho(j)/magVrel;
                        FS = fst + (FS-fst)*exp(-dt/tau); % Updating sepertaion function
                        clift = (FS*CliftInv)+(Cliftfs*(1-FS)); % Coeff. lift updation
                    end
                end
                %% 

                
                cofLift(j,n,i) =clift; % Storage    
                cofDrag(j,n,i) =cdrag; % Storage
                
                magVrel= norm(Vrel); % Magnitude of Vrel
                
                % Lift force
                lift =0.5*rho*(magVrel)^2*cho(j)*clift; % Value
                LiftForce(j,n,i) =lift; % Storage
                
                % Drag force
                drag =0.5*rho*(magVrel)^2*cho(j)*cdrag; % Value
                DragForce(j,n,i) =drag; % Storage
                
                % Forces acting on the blades
                pz(j,n,i) = (lift*cos(phi))+(drag*sin(phi)); % Normal Force
                py(j,n,i) = (lift*sin(phi))-(drag*cos(phi)); % Tangential Force
                
                A = (-W(3,n)/V(3)); % Axial Induction Factor

                % Axial induction factor condition 
                if A<= 1/3
                    fg = 1; % Glauert Correction
                else
                    fg = (1/4)*(5-(3*A)); % Glauert Correction
                end

                F = (2/pi)*acos(exp((-bla/2)*((r-x)/(x*sin(abs(phi)))))); % Prandtl Tip Correction Factor
                
                Vprime = sqrt((V(2))^2+(V(3)+(fg*W(3,n)))^2); % V prime
                
                %% Wake Filter
                if wake ==1
                    Wqs(3,n+1)= (-bla*lift*cos(phi))/(4*pi*rho*x*F*Vprime); % Quasi_Steady Induced wind in Z direction
                    Wqs(2,n+1)= (-bla*lift*sin(phi))/(4*pi*rho*x*F*Vprime); % Quasi_Steady Induced wind in Y direction
                    
                    % Wake axial induction factor condition
                    if A<=0.5 
                        tau1= (1.1/(1-(1.3*A)))*(r/v);
                        tau2= (0.39-(0.26*(x/r)^2))*tau1;
                    else
                        A = 0.5;
                        tau1= (1.1/(1-(1.3*A)))*(r/v);
                        tau2= (0.39-(0.26*(x/r)^2))*tau1;
                    end
                    
                    % Backward difference 
                    Hz = Wqs(3,n+1) + (0.6*tau1*((Wqs(3,n+1)-Wqs(3,n))/dt));
                    Hy = Wqs(2,n+1) + (0.6*tau1*((Wqs(2,n+1)-Wqs(2,n))/dt));
                    
                    % Intermediate Induced Velocity
                    Wint(3,n+1) = Hz +(Wint(3,n)-Hz)*exp(-dt/tau1);
                    Wint(2,n+1) = Hy +(Wint(2,n)-Hy)*exp(-dt/tau1);
                    
                    % Updated Induced Velocity
                    W(3,n+1)= Wint(3,n+1)+ ((W(3,n)-Wint(3,n+1))*exp(-dt/tau2));
                    W(2,n+1)= Wint(2,n+1)+ ((W(2,n)-Wint(2,n+1))*exp(-dt/tau2));

                else % No Wake filter
                    % Updated Induced Velocity (Quasi Steady Value) 
                    W(3,n+1)= (-bla*lift*cos(phi))/(4*pi*rho*x*F*Vprime);
                    W(2,n+1)= (-bla*lift*sin(phi))/(4*pi*rho*x*F*Vprime);
                end
            else % Tip Correction case
                pz(j,:,:) = 0; 
                py(j,:,:) = 0;
            end
            n=n+1; % Time step increment
        end
        n= 1; % Time step intialization
    end
end

for t = 0.15:dt:T
    % Moment, Thrust, Power calculation
    Moment(n) = trapz(rad,py(:,n,1).*rad); 
    Thrust(n) = trapz(rad,pz(:,n,1));
    Power(n) =  w*Moment(n);
    time(n) = t;
    n=n+1;
end

% figure
% plot(rad,py(:,:,1))
% hold on;
% plot(rad,pz(:,:,1))

%% Power vs time graph
figure
plot(time,Power)
xlabel("Time (s)")
ylabel("Power (W)")
title("Power vs Time")
grid on
grid minor
saveas(gcf,"power-time.png")

%% Thrust vs time graph
figure
plot(time,Thrust)
xlabel("Time (s)")
ylabel("Thrust (N)")
title("Thrust vs Time")
grid on
grid minor
saveas(gcf,"thrust-time.png")

%% Cl-AOA graph
figure 
plot(Alpha(10,:),cofLift(10,:))
xlabel("Angle of attack (deg.)")
ylabel("Coefficient of lift")
title("Cl - Aoa")
grid on
grid minor
saveas(gcf,"cl-alfa.png")