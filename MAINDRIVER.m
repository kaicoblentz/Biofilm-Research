 %MAINDRIVER
clear; clc
tic
%% Rewriting inputs section to run through test cases
%Call on specific test case parameters
num=1 ; %number of case, A corresponds to 1, B corresponds to 2....
param=cases(num); %structure variables are stored in

%Tank Parameters + Geometry
L=0.5; %[m]
W=0.5; %[m]
H=0.4; %[m]
SA=(param.V/H)+2*((param.V/L)+(param.V/W)); %tank surface area [m^2] 

Co=param.So; %substrate concentration

%Create initial biofilm grid
Nz=12; %Linear GridPoints in Biofilm
z=linspace(0,param.Lf,Nz); %[m] Grid of Biofilm Depth
dz=z(2)-z(1); %[m]

%Initial Boundary Conditions (in Biofilm)
Sb=zeros(1,Nz);
Sb(end)=param.So; %initially assume boundary concentration = So

%Time Constraints
tFin=20; %[s]
dt=1e-2; %Interval
N=tFin/dt; %Number of steps
outFreq=20; %Number of steps between plot updates.

%Preallocation
t = zeros(1,N); %Time
x = zeros(1,N); %Biomass Concentration in bulk liquid
S = zeros(1,N); %Substrate in bulk liquid
bflux=zeros(1,N); %Boundary Layer Flux of Biofilm Preallocate
flux=zeros(1,N); %Right hand side of power point equation to ensure matching flux
thickness=zeros(1,N); %Right hand side of power point equation to ensure matching flux

%Initial Conditions
t(1)=0;
x(1)=param.xo;
S(1)=param.So;

%Initialize plots 
outIter=outFreq-1;
plots=0;

%% Time Loop
for i = 1:N-1
    
    %Update biofilm grid as biofilm grows
    z=linspace(0,param.Lf,Nz); %[m] Grid of Biofilm Depth
    dz=z(2)-z(1); %[m]
    
    %Call on "biofilmdiffusion"
    Sbold=Sb; %store old substrate concentration array
    [Cs,Sb,bflux(i+1),flux(i+1)]=biofilmdiffusion(Sbold,S(i),Nz,dz,param);
    
    %Call on "lf"
    Lf_old=param.Lf;
    [param.Lf,Vdet]=lf(Sb,Lf_old,dt,dz,param);
    
    %Call on "tankenvironment"
    [t(i+1),x(i+1),S(i+1)]=tankenvironment(t(i),x(i),S(i),SA,Vdet,dt,Cs,Co,param);
    
    thickness(i+1)=param.Lf;
   
    %Call on desired plots from 'outputs'
    outIter=outIter+1;
    if (outIter==outFreq)
        [plots] = outputs(t(1:i+1),x(1:i+1),S(1:i+1),z,bflux(1:i+1),thickness(1:i+1),Sb,param,plots);
        outIter=0;
    end
    
end
%end
toc