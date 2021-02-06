%{
% Author: Maharshi Gurjar
% ELEC 4700 - Modeling of Integrated Devices
% Assignment 1
%}
clc; close all; clear;
set(0, 'DefaultFigureWindowStyle', 'docked')
%Define simulation envrionment and constants
M0 = 9.10938356e-31; %Rest mass of electron
Mass_n = 0.26*M0; %Effective mass of electron
T = 300; % Simulation envrionment temperature (K)
k = 1.38064852e-23; % Boltzmans constant
V_thermal = sqrt(2*k*T/Mass_n); %Thermal Veleocity
Height = 100e-9; % The height of the simulation environment
Length = 200e-9; % The lengthof the simulation environment
nElectrons = 2e3; % Total number of electrons to simulate
nPlotted_Electrons = 10; %Total number of electrons displayed 
Time_Step = Height/V_thermal/100; % Time step of simulation
Iterations = 1000; % Number of iternations to simulate
Show_Movie = 1; %Display steps control
% The mean free path is determined by multipling the thermal velocity 
... by the mean time between collisions: 
MFP = V_thermal * 0.2e-12; %Mean free path 
%The state of the electron  (postion and velocity) is stored in a array
... where each index refers to [x-position y-position v-in-x v-in-y]
Electron_State = zeros(nElectrons,4);
%Trajectories will be recorded in the array below, where there are 
... double the colums as we need both the x and y positions for each
    ... of the "to-be-plotted" electrons
Trajectories = zeros(Iterations,nPlotted_Electrons*2);
%Temperature will be recorded in the array below
Temperature = zeros(Iterations,1);
%Create a scattering probability
P_Scatterieng = 1 - exp(-Time_Step/0.2e-12);
%Create a distribution using the matlab makedist function
Velocity_PDF = makedist('Normal', 'mu', 0, 'sigma', sqrt(k*T/Mass_n));



%Setting the top/bottom of the boxes specularity
Box_Top_Specular = 1;
Box_Bottom_Specular = 1;
%Create Box-positions [x1, x2, y1,y2]
Box_pos = 1e-9*[80 120 0 40; 80 120 60 100];
%Create the state of the box (specular or 
Box_state = [0 1];


%Generate a random inital population postion and velocity
for i = 1:nElectrons
   Electron_State(i,:) = [Length*rand() Height*rand() random(Velocity_PDF) random(Velocity_PDF)];
   
   
%    if ( ((80e-9<=Electron_State(i,1))&&(Electron_State(i,1)<=120e-9)) && ...
%        ((0<=Electron_State(i,2))&&(Electron_State(i,2)<=40e-9)) )
%        Electron_State(i:1) = Length*rand();
%        Electron_State(i:2) = Height*rand();
%    elseif ( ((80e-9<=Electron_State(i,1))&&(Electron_State(i,1)<=120e-9)) && ...
%        ((50e-9<=Electron_State(i,2))&&(Electron_State(i,2)<=Height)) )
%        Electron_State(i:1) = Length*rand();
%        Electron_State(i:2) = Height*rand();
%    end
end

%Average velocity calculation
Average_Velocity = sqrt ((Electron_State(:,3).^2)/nElectrons + (Electron_State(:,4).^2)/nElectrons);

%We will now move (iterate) over time, updating the positions and direction
...while plotting the state
for i = 1:Iterations
    %The line below updates the x,y position by moving it to a new position
    ... using its current position + the velocity*(time step)
    Electron_State(:,1:2) = Electron_State(:,1:2) + Time_Step.*Electron_State(:,3:4);
    
    %Checking boundary conditions using Matlab matrix equations 
    
    %Check if and move all electrons at X=200nm Bound:
    Electron_State((Electron_State(:,1)>Length),1) = Electron_State((Electron_State(:,1)>Length),1) - Length;
 
    %Check if and move all electrons at X=0nm Bound:
    Electron_State((Electron_State(:,1)<0),1) =Electron_State((Electron_State(:,1)<0),1) + Length;
    
    %Check if and move all electrons at Y=100nm Bound:
    if (Box_Top_Specular == 1)
       Electron_State((Electron_State(:,2)>Height),4) = -1*Electron_State((Electron_State(:,2)>Height),4) ;
    elseif (Box_Top_Specular == 0)
       Electron_State((Electron_State(:,2)>Height),2) = Electron_State((Electron_State(:,2)>Height),2) - Height ;
    end
    if (Box_Bottom_Specular == 1)
        Electron_State((Electron_State(:,2)<0),4) =  -1*Electron_State((Electron_State(:,2)<0),4) ;
    elseif (Box_Bottom_Specular == 0)
        Electron_State((Electron_State(:,2)<0),2) =  Electron_State((Electron_State(:,2)<0),2) + Height;
    end
    

    
    %Add scattering
     j = rand(nElectrons,1) < P_Scatterieng;
     Electron_State(j,3:4) = random(Velocity_PDF,[sum(j),2]);
 
    % Stores the Electron [x y] posistions in the Trajectories vector
    ... for each different electron in a new coloum
    for j = 1: nPlotted_Electrons
       Trajectories(i, (j):(j+1)) = Electron_State(j,1:2); 
    end
    %To calcuatle the themal energy, Maxwell's principle of equipartion 
    ... is used,  where the final equation then becomes;
    Temperature(i) = ( sum (Electron_State(:,3).^2) + sum(Electron_State(:,4).^2)) * Mass_n / k / 2 / nElectrons;
    
    %Shows the pathing of the electron, as well as the updating trajectory
    if Show_Movie && mod(i,50)
       figure(1)
       hold off;
       plot(Electron_State(1:nPlotted_Electrons,1)./1e-9,Electron_State(1:nPlotted_Electrons,2)./1e-9,'o');
       grid on;
       axis([0 Length/1e-9 0 Height/1e-9]);
       xlabel('x (nm)');
       ylabel('y (nm)');
       title(sprintf("Plotting (%d/%d) electron at constant velocity",nPlotted_Electrons,nElectrons));
       hold on;
    end
end
figure("name","Trajectory, temperature and speed results results")
subplot(3,1,1)
hold on;
for i = 1:nPlotted_Electrons
    plot(Trajectories(:,i)./1e-9, Trajectories(:,i+1)./1e-9,'-');
end
axis([0 Length/1e-9 0 Height/1e-9]);
xlabel('x (nm)');
ylabel('y (nm)');
grid on;
title(sprintf("Trajectories of (%d/%d) electron(s) at constant velocity",nPlotted_Electrons,nElectrons));

subplot(3,1,2)
plot(Time_Step*(0:Iterations-1), Temperature);
grid on;
title('Temperature');
xlabel('Time (s)');
ylabel('Temperature (K)');

subplot(3,1,3)
Velocity = sqrt(Electron_State(:,3).^2 + Electron_State(:,4).^2);
title("Electron Speed");
histogram(Velocity);
xlabel("Speed (m/s)");
ylabel("Number of particles");
grid on;