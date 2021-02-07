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
V_thermal = sqrt(2*k*T/Mass_n) %Thermal Veleocity
Height = 100e-9; % The height of the simulation environment
Length = 200e-9; % The lengthof the simulation environment
nElectrons = 3e3; % Total number of electrons to simulate
nPlotted_Electrons = 100; %Total number of electrons displayed 
Time_Step = Height/V_thermal/100; % Time step of simulation
Iterations = 1000; % Number of iternations to simulate
Show_Movie =0 ; %Display steps control
% The mean free path is determined by multipling the thermal velocity 
... by the mean time between collisions: 
MFP = V_thermal * 0.2e-12 %Mean free path 
%The state of the electron  (postion and velocity) is stored in a array
... where each index refers to [x-position y-position v-in-x v-in-y]
Electron_State = zeros(nElectrons,4);
%Temperature will be recorded in the array below
Temperature = zeros(Iterations,1);
%Create a scattering probability
P_Scatterieng =1 - exp(-Time_Step/0.2e-12);
%Create a distribution using the matlab makedist function
Velocity_PDF = makedist('Normal', 'mu', 0, 'sigma', sqrt(k*T/Mass_n));
%Setting the top/bottom of the boudnary specularity
Top_Specular = 1;
Bottom_Specular = 1;
%Set the box characteristics (0 = diffusive)
Top_Box_Specular = 1;
Bottom_Box_Specular = 1;
%Create Box-positions [x1, x2, y1,y2]
Box_pos = 1e-9.*[80 120 0 40; 80 120 60 100];
%Create the state of the box (specular or diffusive)
%Generate a random inital population postion and velocity
for i = 1:nElectrons
   Electron_State(i,:) = [Length*rand() Height*rand() random(Velocity_PDF) random(Velocity_PDF)];  
   Electron_State( (Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ...
       (Electron_State(:,2)<40e-9)) ,1) = Length*rand();
   Electron_State( (Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ...
       (Electron_State(:,2)>40e-9)) ,1) = Length*rand();
end
%Figure below used to test inital postion of electrons, and if movie is off
%shows the the same.
figure("Name","Electron positions")
plot(Electron_State(1:nPlotted_Electrons,1)./1e-9,Electron_State(1:nPlotted_Electrons,2)./1e-9,'o');
grid on;
axis([0 Length/1e-9 0 Height/1e-9]);
xlabel('x (nm)');
ylabel('y (nm)');
title(sprintf("Plotting (%d/%d) electron ",nPlotted_Electrons,nElectrons));
saveas(gcf,'Part_Three_Boxes.png')
hold on;
for j=1:size(Box_pos,1)
   plot([Box_pos(j, 1) Box_pos(j, 1) Box_pos(j, 2) Box_pos(j, 2) Box_pos(j, 1)]./1e-9,...
       [Box_pos(j, 3) Box_pos(j, 4) Box_pos(j, 4) Box_pos(j, 3) Box_pos(j, 3)]./1e-9, 'k-');
end
hold off;
%We will now move (iterate) over time, updating the positions and direction
...while plotting the state
for i = 1:Iterations
    %The line below updates the x,y position by moving it to a new position
    ... using its current position + the velocity*(time step)
    Electron_Prev_State = Electron_State(:,1:2);
    Electron_State(:,1:2) = Electron_State(:,1:2) + Time_Step.*Electron_State(:,3:4);
    
    %Checking boundary conditions using Matlab matrix equations/indenxing
    %Check (if) and move all electrons at X=200nm Bound:
    Electron_State((Electron_State(:,1)>Length),1) = Electron_State((Electron_State(:,1)>Length),1) - Length;
 
    %Check (if) and move all electrons at X=0nm Bound:
    Electron_State((Electron_State(:,1)<0),1) =Electron_State((Electron_State(:,1)<0),1) + Length;
    
    %Check (if) and move all electrons at Y Bounds and if specular or diffusive:
    if (Top_Specular == 1)
       Electron_State((Electron_State(:,2)>Height),4) = -1*Electron_State((Electron_State(:,2)>Height),4) ;
       Electron_State((Electron_State(:,2)>Height),2) = 2*Height - Electron_State((Electron_State(:,2)>Height),2);
    else
       Electron_State((Electron_State(:,2)>Height),4) = -random(Velocity_PDF);
       Electron_State((Electron_State(:,2)>Height),3) = random(Velocity_PDF);
    end
    if (Bottom_Specular == 1)
        Electron_State((Electron_State(:,2)<0),4) = -1*Electron_State((Electron_State(:,2)<0),4) ;
        Electron_State((Electron_State(:,2)<0),2) = -Electron_State((Electron_State(:,2)<0),2);
    else
       Electron_State((Electron_State(:,2)>Height),4) = random(Velocity_PDF);
       Electron_State((Electron_State(:,2)>Height),3) = random(Velocity_PDF);
    end

%%
    %BOX SPECULAR CASE :
    %Check if bottom Box is specular and if so then bounce 
    Electron_State( (Bottom_Box_Specular==1 & Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ((Electron_State(:,2)<40e-9 & Electron_Prev_State(:,2)<40e-9))),3)...
        =-Electron_State( (Bottom_Box_Specular==1 & Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ((Electron_State(:,2)<40e-9 & ...
        Electron_Prev_State(:,2)<40e-9))),3);
    %Check if top Box is specular and if so then bounce 
    Electron_State( (Top_Box_Specular==1 & Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ((Electron_State(:,2)>60e-9 & Electron_Prev_State(:,2)>60e-9))),3)...
        =-Electron_State( (Top_Box_Specular==1 & Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ((Electron_State(:,2)>60e-9 & Electron_Prev_State(:,2)>60e-9))),3);

    %Top Horizontal
    Electron_State( (Top_Box_Specular==1 & Electron_State(:,1)<120e-9 & Electron_State(:,1)>80e-9) & (Electron_State(:,2)>60e-9) &...
        (Electron_Prev_State(:,1)>80e-9 & Electron_Prev_State(:,1)<120e-9),4) = -Electron_State( (Top_Box_Specular==1 & Electron_State(:,1)<120e-9...
        & Electron_State(:,1)>80e-9) & (Electron_State(:,2)>60e-9) &(Electron_Prev_State(:,1)>80e-9 & Electron_Prev_State(:,1)<120e-9),4);
    %Bottom Hoizontal
    Electron_State( (Bottom_Box_Specular==1 & Electron_State(:,1)<120e-9 & Electron_State(:,1)>80e-9) & (Electron_State(:,2)<40e-9) &...
        (Electron_Prev_State(:,1)>80e-9 & Electron_Prev_State(:,1)<120e-9),4) = -Electron_State( (Bottom_Box_Specular==1 & Electron_State(:,1)<120e-9...
            & Electron_State(:,1)>80e-9) & (Electron_State(:,2)<40e-9) &(Electron_Prev_State(:,1)>80e-9 & Electron_Prev_State(:,1)<120e-9),4);
    %BOX SPECULAR CASE END :
%%
    %BOX DIFFUSIVE CASE:
    %Check if bottom Box is specular and if so then bounce 
    Electron_State( (Bottom_Box_Specular==0 & Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ((Electron_State(:,2)<40e-9 & Electron_Prev_State(:,2)<40e-9))),3)...
        =random(Velocity_PDF);
    Electron_State( (Bottom_Box_Specular==0 & Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ((Electron_State(:,2)<40e-9 & Electron_Prev_State(:,2)<40e-9))),4)...
        =random(Velocity_PDF);
    
    %Check if Top Box is specular and if so then bounce 
    Electron_State( (Top_Box_Specular==0 & Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ((Electron_State(:,2)>60e-9 & Electron_Prev_State(:,2)>60e-9))),3)...
        =random(Velocity_PDF);
    Electron_State( (Top_Box_Specular==0 & Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ((Electron_State(:,2)>60e-9 & Electron_Prev_State(:,2)>60e-9))),4)...
        =random(Velocity_PDF);

    %Bottom Horizontal
    Electron_State( (Bottom_Box_Specular==0 & Electron_State(:,1)<120e-9 & Electron_State(:,1)>80e-9) & (Electron_State(:,2)<40e-9) &...
        (Electron_Prev_State(:,1)>80e-9 & Electron_Prev_State(:,1)<120e-9),3) = random(Velocity_PDF);
    Electron_State( (Bottom_Box_Specular==0 & Electron_State(:,1)<120e-9 & Electron_State(:,4)<0 & Electron_State(:,1)>80e-9) & (Electron_State(:,2)<40e-9) &...
        (Electron_Prev_State(:,1)>80e-9 & Electron_Prev_State(:,1)<120e-9),4) = random(Velocity_PDF);
    %Top Hoizontal
    Electron_State( (Top_Box_Specular==0  & Electron_State(:,1)<120e-9 & Electron_State(:,1)>80e-9) & (Electron_State(:,2)>60e-9) &...
        (Electron_Prev_State(:,1)>80e-9 & Electron_Prev_State(:,1)<120e-9),3) = random(Velocity_PDF);
    Electron_State( (Top_Box_Specular==0 & Electron_State(:,4)>0 & Electron_State(:,1)<120e-9 & Electron_State(:,1)>80e-9) & (Electron_State(:,2)>60e-9) &...
        (Electron_Prev_State(:,1)>80e-9 & Electron_Prev_State(:,1)<120e-9),4) = -random(Velocity_PDF);
%%
    %Add scattering
    j = rand(nElectrons,1) < P_Scatterieng;
    Electron_State(j,3:4) = random(Velocity_PDF,[sum(j),2]);
%%
    %Electron shifting
    %Region 1
    Electron_State( (Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ((Electron_State(:,2)<40e-9 & Electron_Prev_State(:,2)<=40e-9)) &...
        Electron_Prev_State(:,1)<=80e-9),1) = 2*80e-9 -Electron_State( (Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 ...
        & ((Electron_State(:,2)<=40e-9 & Electron_Prev_State(:,2)<40e-9)) & Electron_Prev_State(:,1)<80e-9),1);
    %Region 2
    Electron_State( (Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 & ((Electron_State(:,2)<40e-9 & Electron_Prev_State(:,2)<=40e-9)) &...
        Electron_Prev_State(:,1)>120e-9),1) = 2*120e-9 - Electron_State( (Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 ...
        & ((Electron_State(:,2)<=40e-9 & Electron_Prev_State(:,2)<40e-9)) & Electron_Prev_State(:,1)>120e-9),1);
    %Region 3
    Electron_State( ( (Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9) & ((Electron_State(:,2)>60e-9 & Electron_Prev_State(:,2)>=60e-9)) &...
        Electron_Prev_State(:,1)<80e-9),1) = 2*80e-9 - Electron_State( (Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9...
        & ((Electron_State(:,2)>=60e-9 & Electron_Prev_State(:,2)>60e-9)) & Electron_Prev_State(:,1)<80e-9),1);
    %Region 4
    Electron_State( ( (Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 )& ((Electron_State(:,2)>60e-9 & Electron_Prev_State(:,2)>=60e-9)) &...
        Electron_Prev_State(:,1)>120e-9),1) = 2*120e-9 - Electron_State( (Electron_State(:,1)>80e-9 & Electron_State(:,1)<120e-9 &...
        ((Electron_State(:,2)>60e-9 & Electron_Prev_State(:,2)>=60e-9)) & Electron_Prev_State(:,1)>120e-9),1);
    
    %Top Horizontal
    Electron_State( (Electron_State(:,1)<120e-9 & Electron_State(:,1)>80e-9 & (Electron_State(:,2)<40e-9 & Electron_Prev_State(:,2)>=40e-9) &...
        (Electron_State(:,4)<0)),2) =  2*40e-9 + Electron_State( (Electron_State(:,1)<120e-9 & Electron_State(:,1)>80e-9 &...
        (Electron_State(:,2)<=40e-9 & Electron_Prev_State(:,2)>=40e-9) & (Electron_State(:,4)<0)),2);
    %Bottom Horizontal
    Electron_State( (Electron_State(:,1)<120e-9 & Electron_State(:,1)>80e-9 & (Electron_State(:,2)>60e-9 & Electron_Prev_State(:,2)<=60) &...
        (Electron_State(:,4)>0)),2) =  2*60e-9 - Electron_State( (Electron_State(:,1)<120e-9 & Electron_State(:,1)>80e-9 &...
        (Electron_State(:,2)>=60e-9 & Electron_Prev_State(:,2)<=60) &(Electron_State(:,4)>0)),2);
    %Electron shifting end
%%
    % Stores the Electron [x y] posistions in the Trajectories vector
    ... for each different electron in a new coloum
    for j = 1: nPlotted_Electrons
       Trajectories_x(i,j) = Electron_State(j,1); 
       Trajectories_y(i,j) = Electron_State(j,2);
    end
    %To calcuatle the themal energy, Maxwell's principle of equipartion 
    ... is used,  where the final equation then becomes;
    Temperature(i) = ( sum (Electron_State(:,3).^2) + sum(Electron_State(:,4).^2)) * Mass_n / k / 2 / nElectrons;
    
    %Shows the pathing of the electron, as well as the updating trajectory
    if Show_Movie && mod(i,20)
       figure(1)
       plot(Electron_State(1:nPlotted_Electrons,1)./1e-9,Electron_State(1:nPlotted_Electrons,2)./1e-9,'o');
       grid on;
       axis([0 Length/1e-9 0 Height/1e-9]);
       xlabel('x (nm)');
       ylabel('y (nm)');
       title(sprintf("Plotting (%d/%d) electron at constant velocity",nPlotted_Electrons,nElectrons));
       hold on;
       for j=1:size(Box_pos,1)
           plot([Box_pos(j, 1) Box_pos(j, 1) Box_pos(j, 2) Box_pos(j, 2) Box_pos(j, 1)]./1e-9,...
               [Box_pos(j, 3) Box_pos(j, 4) Box_pos(j, 4) Box_pos(j, 3) Box_pos(j, 3)]./1e-9, 'k-');
       end
       hold off;
    end
end
%%
figure("name","Trajectory, temperature and speed results results")
subplot(3,1,1)
plot(Trajectories_x(:,1:nPlotted_Electrons)./1e-9, Trajectories_y(:,1:nPlotted_Electrons)./1e-9,'.');
hold on;
for j=1:size(Box_pos,1)
   plot([Box_pos(j, 1) Box_pos(j, 1) Box_pos(j, 2) Box_pos(j, 2) Box_pos(j, 1)]./1e-9,...
       [Box_pos(j, 3) Box_pos(j, 4) Box_pos(j, 4) Box_pos(j, 3) Box_pos(j, 3)]./1e-9, 'k-');
end
hold off
axis([0 Length/1e-9 0 Height/1e-9]);
xlabel('x (nm)');
ylabel('y (nm)');
grid on;
title(sprintf("Trajectories of (%.2f/%.2f) electron(s) at constant velocity",nPlotted_Electrons,nElectrons));
subplot(3,1,2)
plot(Time_Step*(0:Iterations-1), Temperature);
grid on;
xlim([0 Time_Step*Iterations])
title(sprintf("Temperature of the region, Average Temperature: %.2f (K)",mean(Temperature)))
xlabel('Time (s)');
ylabel('Temperature (K)');
subplot(3,1,3)
Velocity = sqrt(Electron_State(:,3).^2 + Electron_State(:,4).^2);
histogram(Velocity);
title(sprintf("Electron Velocity, Average Velocity %.2d (m/s)",mean(Velocity)));
xlabel("Speed (m/s)");
ylabel("Number of particles");
grid on;


%% 
%Create density maps for the electrons and temperature
%Density without filtering :
figure("name","Electron Density")
X = [Electron_State(:,1) Electron_State(:,2)];
hist3(X,'Nbins',[200 100],'CdataMode','auto')
axis([0 Length 0 Height])
xlabel('x (nm)')
ylabel('y (nm)')
colormap;
c = colorbar;
c.Label.String = 'Electron density';
view(2)
title(sprintf("Electron density of %d electrons",nElectrons));
saveas(gcf,'Part_Three_Electron_Density.png')
%With Guassian filtering :
Density = hist3(Electron_State(:,1:2),[200 100])';
N = 20;
sigma = 3;
%Creating a Gaussian filtering matrix
[x,y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
F=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
F=F./sum(F(:));

figure("name","Electron Density (with Gaussian filtering)")
imagesc(conv2(Density,F,'same'))
xlabel('x (nm)')
ylabel('y (nm)')
c = colorbar;
c.Label.String = 'Density';
view(2)
title(sprintf("Electron density of %d electrons",nElectrons));
saveas(gcf,'Part_Three_Electron_Density_Filtered.png')
%For the temperature density plot the electrons positions were binned
Temperature_Sum_X = zeros(ceil(Length/1e-9),ceil(Height/1e-9));
Temperature_Sum_Y = zeros(ceil(Length/1e-9),ceil(Height/1e-9));
Temperature_Bin_Count = zeros(ceil(Length/1e-9),ceil(Height/1e-9));

% Look at velocities of all the particles
for i=1:nElectrons
    % bin the electrons
    x = floor(Electron_State(i,1)./1e-9);
    y = floor(Electron_State(i,2)./1e-9);
    if(x==0)
        x = 1;
    end
    if(y==0)
        y= 1;
    end
    % Add velocity components
    Temperature_Sum_X(x,y) = Temperature_Sum_X(x,y) + Electron_State(i,4)^2;
    Temperature_Sum_Y(x,y) = Temperature_Sum_Y(x,y) + Electron_State(i,3)^2;
    Temperature_Bin_Count(x,y) = Temperature_Bin_Count(x,y) + 1;
end

% Now, with the velocities added up, calculate the temperatures:
TemperatureDensity = (Temperature_Sum_X + Temperature_Sum_Y).*Mass_n./k./2./Temperature_Bin_Count;
%If somewhere in the calculation the density beame nan then make it 0
TemperatureDensity(isnan(TemperatureDensity)) = 0;
%Transpose the matrix
TemperatureDensity = TemperatureDensity';

figure("Name","Electron heat denstiy")
imagesc(TemperatureDensity)
set(gca,'YDir','normal');
title('Temperature Map');
c = colorbar;
c.Label.String = 'Temperature';
xlabel('x (nm)');
ylabel('y (nm)');

figure("Name","electron heat density (with Gaussian filtering)");
imagesc(conv2(TemperatureDensity,F,'same'));
set(gca,'YDir','normal');
title('Temperature Map');
c = colorbar;
c.Label.String = 'Temperature (K)';
xlabel('x (nm)');
ylabel('y (nm)');

saveas(gcf,'Part_Three_Heat_Map_Filtered.png')
