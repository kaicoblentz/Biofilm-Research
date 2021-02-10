%% Heun's Method, Biomass and Substrate Concetration in Bulk Liquid 

clear; clc
% Inputs

%Array Initial Conditions
xo=10; %initial biomass concentration in bulk liquid
So=25; %initial substrate concentration in bulk liquid

%Tank Parameters + Geometry
L=0.5; %[m]
W=0.5; %[m]
H=0.4; %[m]
V=L*W*H; %tank volume [m^3]
SA=(V/H)+2*((V/L)+(V/W)); %tank surface area [m^2] 
Qdot=1; %[flow rate in/out, m^3]
Sin=25; %Inflow of substrates to tank, [g/m^3]

%Biofilm Parameters
mmax=20; %max specific growth rate
Km=3; %Monod half-saturation coefficient(growth transitions from sat. to linear)
Yxs=0.5; %ratio of substrate consumed to biomass produced
Daq=2e-5; %diffusion coefficient of water(assumed at boundary) [m/s^2]
Lf=4.00E-4; %biofilm thickness [m]
LL=Lf/100; %thickness of boundary layer [m]
Co=So; %substrate concentration
Xb=20000; %g m-3	biomass density in biofilm
De=5.00E-05; %m2 d-1	effective diffusion coefficient of substrate in biofilm

%Time Constraints
tFin=2; %[s]
dt=1e-3; %Interval
N=tFin/dt; %Number of steps

%S = @(x) Sin-(V/Qdot)*x+(V^2/Qdot); %Substrate concentration wrt biomass
%x = @(S) (Qdot*Sin/V)-(S*Qdot/V)+V; %Biomass wrt substrate concetration
%m = @(x,S) ((mmax*S)/(Km+S)); %Monod Growth Kinetics



dsdt = @(x,t,S,Cs) -((((mmax*S)/(Km+S))*x)/Yxs)+((Qdot*Sin)/V)-((Qdot*S)/V)-(SA*((Daq/LL)*(Co-Cs))); 
% ^^^Substrate Concentration Change wrt time, now also considers flux
% through boundary layer of biofilm
dxdt = @(x,t,S,Cs) (((mmax*S)/(Km+S))-(Qdot/V))*x; %Biomass Concentration Change wrt time


%Preallocation
t = zeros(1,N); %Time
x = zeros(1,N); %Biomass Concentration in bulk liquid
S = zeros(1,N); %Substrate in bulk liquid

%Initial Conditions
t(1)=0;
x(1)=xo;
S(1)=So;

for i = 1:N-1
   
    % Insert loop over spacial coordinate zeta for substrate diffusion,
    % particulates, biomass growth within biofilm, etc

    t(i+1) = t(i) + dt;
    
    %Calling Biofilm Surface Substrate Concentration from 'Diffusion'
    [Sb,bflux,flux]=Diffusion(Lf,S(i),mmax,Xb,Yxs,De);
    Cs=Sb(end);
    
    xstar = x(i) + dt*dxdt(x(i),t(i),S(i),Cs);
    Sstar = S(i) + dt*dsdt(x(i),t(i),S(i),Cs);
    
    x(i+1) = x(i) + dt/2*(dxdt(x(i),t(i),S(i),Cs)+dxdt(xstar,t(i+1),Sstar,Cs));
    S(i+1) = S(i) + dt/2*(dsdt(x(i),t(i),S(i),Cs)+dsdt(xstar,t(i+1),Sstar,Cs)); 
    

end


%% plot
figure(1); clf(1)
plot(t,x)
hold on
plot(t,S)
title('Biomass and Substrate Concentrations For Filling/Draining Tank')
xlabel('Time')
ylabel('Amount of Biomass/Substrate in Tank')
legend('Biomass','Substrate')




