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

%Generate a random inital population postion and velocity
for i = 1:nElectrons
   angle = rand()*2*pi; %angle of trajectory
   Electron_State(i,:) = [Length*rand() Height*rand() V_thermal*cos(angle) V_thermal*sin(angle)]; 
end

%We will now move (iterate) over time, updating the positions and direction
...while plotting the state
for i = 1:Iterations
    %The line below updates the x,y position by moving it to a new position
    ... using its current position + the velocity*(time step)
    Electron_State(:,1:2) = Electron_State(:,1:2) + Time_Step.*Electron_State(:,3:4);

    %Checks the boundary Conditions (if electron is at the bounds)
    for j = 1 : nElectrons
       if Electron_State(j,1) > Length
           Electron_State(j,1) = Electron_State(j,1) - Length;
       elseif Electron_State(j,1) <0
           Electron_State(j,1) = Electron_State(j,1) + Length;
       end
       if Electron_State(j,2) > Height
           Electron_State(j,4) = -Electron_State(j,4);
       elseif Electron_State(j,2) < 0
           Electron_State(j,4) = -Electron_State(j,4);
       end
    end
    % Stores the Electron [x y] posistions in the Trajectories vector
    ... for each different electron in a new coloum
    for j = 1: nPlotted_Electrons
       Trajectories(i, (j*2):(j*2+1)) = Electron_State(j,1:2); 
    end
    %To calcuatle the themal energy, Maxwell's principle of equipartion 
    ... is used,  where the final equation then becomes;
    Temperature(i) = ( sum (Electron_State(:,3).^2) + sum(Electron_State(:,4).^2)) * Mass_n / k / 2 / nElectrons;
    
    %Shows the pathing of the electron, as well as the updating trajectory
    ... Where the mod(i,n) controls the speed of updating the plot, for 
        ...higher levelse of n the plot updates faster 
    if Show_Movie && mod(i,5) == 0
       figure(1)
       subplot(1,1,1);
       hold off;
       plot(Electron_State(1:nPlotted_Electrons,1)./1e-9,Electron_State(1:nPlotted_Electrons,2)./1e-9,'o');
       grid on;
       axis([0 Length/1e-9 0 Height/1e-9]);
       xlabel('x (nm)');
       ylabel('y (nm)');
       title(sprintf("Plotting (%d/%d) electron at constant velocity",nPlotted_Electrons,nElectrons));
         
    end
end
figure(2)
subplot(1,1,1)
for i = 1:nPlotted_Electrons
    plot(Trajectories(:,i*2)./1e-9, Trajectories(:,i*2+1)./1e-9, '-');
end
axis([0 Length/1e-9 0 Height/1e-9]);
xlabel('x (nm)');
ylabel('y (nm)');
grid on;
title(sprintf("Trajectories of (%d/%d) electron at constant velocity",nPlotted_Electrons,nElectrons));