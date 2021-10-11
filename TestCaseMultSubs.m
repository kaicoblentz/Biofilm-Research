% Unit tests for biofilm model
function tests = TestCaseMultSubs
clear; clc
tests = functiontests(localfunctions);
end

%% Test Diffusion
function test_mult_subs_diffusion(testCase)
% Run Test
param.mumax=2000;
param.Km=2500;
param.Yxs=[0.5 -0.278];
S=[25; 2];
param.Daq=[4e-5 6e-5];
param.De =[1e-5 1.5e-5];
Lf=5e-6;
param.LL=0;
param.Xb=20000;
param.dtol=1e-12;
param.model=1;
Yps     =1.8;

figure(1); clf(1); hold on
Nzs=[10,50,100,1000,2000]; %Grid sizes to test
error=zeros(2,length(Nzs)); %Preallocate

for i=1:length(Nzs)
    Nz=Nzs(i);
    z=linspace(0,Lf,Nz); %[m] Grid of Biofilm Depth
    dz=z(2)-z(1); %[m]
    %Sbold=linspace(0,S(1),Nz);
    Sbold=zeros(2,Nz);
    t=0;
    [Sb,~]=biofilmdiffusion_fd(Sbold,S,Nz,dz,t,param);
    figure(1)
    plot(z,Sb(1,:))
    
    Stest=S(2)+(param.Daq(1)*Yps/param.Daq(2))*(S(1)-Sb(1,:));
    figure(3)
    plot(z,Stest)
    hold on
    plot(z,Sb(2,:),'--')
    legend('Analytic','Numerical')
    
    
    % Analyze Result
    phi = sqrt(param.mumax*param.Xb*Lf*Lf/...
        (param.De(1)*param.Km*param.Yxs(1)));
    Sb_ana = S(1)*cosh(phi*z/Lf)/cosh(phi);
    
    % Error
    error(1,i)=mean(abs(Sb(1,:)-Sb_ana));
    error(2,i)=mean(abs(Sb(2,:)-Stest));
    

end
figure(1)
plot(z,Sb_ana,'--')
xlabel('z')
ylabel('Sb(z)')
title('Substrate Profiles within Biofilm')
legend(sprintf('Gridsize:%5.0f',Nzs(1)),...
       sprintf('Gridsize:%5.0f',Nzs(2)),...
       sprintf('Gridsize:%5.0f',Nzs(3)),...
       sprintf('Gridsize:%5.0f',Nzs(4)),...
       sprintf('Gridsize:%5.0f',Nzs(5)),...
       'Analytic','Location','Northwest')
set(gca,'Fontsize',20)

figure(2); clf(2)
loglog(Nzs,error(1,:),'-o')
hold on
loglog(Nzs,Nzs.^-1,'--')

xlabel('Number of grid points')
ylabel('Error')
set(gca,'Fontsize',20)

figure(3)
xlabel('z')
ylabel('Sb(z)')
title('Product Profiles within Biofilm')

% Pass/fail
tol=1e-2;
verifyLessThan(testCase,min(error),tol)
end

% %function TestCaseMultSubs(param)
% 
% Nz      =param.Nz;
% So      =param.So;
% tFin    =param.tFin;
% outFreq =param.outFreq;
% Daq     =param.Daq;
% Yps     =1.8;
% 
% N=round(tFin/param.dtmax); 
% 
% %Corresponding arrays
% t       =zeros(1,N); %Time
% x       =zeros(1,N); %Biomass Concentration in bulk liquid
% S       =zeros(2,N); %Substrate in bulk liquid
% bflux   =zeros(2,N); %Boundary Layer Flux of Biofilm Preallocate
% Lf      =zeros(1,N); %Right hand side of power point equation to ensure matching flux
% dt      =zeros(1,N); %size of each time step
% 
% %% Initial Conditions
% %Biofilm
% Sb=zeros(2,Nz);
% Sb(:,end)=So; %initially assume boundary concentration = So 
% 
% Lf(1)=param.Lfo;
% 
% %Tank
% t(1)=0;
% x(1)=param.xo;
% S(:,1)=So;
% 
% %Time
% dt(1)=param.dtmax;
% 
% %% Initialize plots 
% outIter=outFreq-1;
% plots=0; titles=0;
% 
% %% Time Loop
% i=1;
% while t(i)<tFin-dt(i)
%     
%     % Check if arrays are filling up
%     if length(bflux)==i 
%         % Compute an estimate for reamaining time steps
%         Nrem    =round((tFin-t(i))/dt(i));
%         
%         % Append time dependant arrays with estimate
%         t       =[t     zeros(1,Nrem)]; 
%         x       =[x     zeros(1,Nrem)]; 
%         S       =[S     zeros(2,Nrem)];
%         bflux   =[bflux zeros(2,Nrem)]; 
%         Lf      =[Lf    zeros(1,Nrem)]; 
%         dt      =[dt    zeros(1,Nrem)];        
%     end
%     
%     %Update biofilm grid as biofilm grows
%     z=linspace(0,Lf(i),Nz); %[m] Grid of Biofilm Depth
%     dz=z(2)-z(1); %[m]
%     
%     %Call on "biofilmdiffusion"
%     [Sb,bflux(:,i+1)]=biofilmdiffusion_fd(Sb,S(:,i),Nz,dz,t(i),bflux(:,i),param);
%     
%     %Call on "lf"
%     [Lf(i+1),Vdet]=lf(Sb,Lf(i),dt(i),dz,param);
%     
%     %Call on "tankenvironment"
%     [~,t(i+1),x(i+1),~,dt(i+1)]=tankenvironment(t(i),x(i),S(:,i),Vdet,dt(i),bflux(:,i+1),param);
%     
%     S(:,i+1) = So;
%     
%     %Call on desired plots from 'outputs'
%     outIter=outIter+1;
%     if (outIter>=outFreq)
%         Stest=So(2)+(Daq(1)*Yps/Daq(2))*(So(1)-Sb(1,:));
%         
%         figure(1); clf(1)
%         subplot(2,2,1);
%         plot(z,Stest)
%         hold on
%         plot(z,Sb(2,:))
%         legend('Analytic','Numerical')
%         xlabel('Position within Biofilm (z)')
%         ylabel('Concentration')
%         title('Comparison of Product Values')
%         
%         subplot(2,2,2)
%         plot(z,Sb(1,:))
%         ylabel('Concentration')
%         xlabel('Position within Biofilm (z)')
%         title('Substrate Concentration')
%         
%         subplot(2,2,3);
%         plot(z,Stest+28)
%         hold on
%         plot(z,Sb(2,:))
%         legend('Offset Analytic','Numerical')
%         xlabel('Position within Biofilm (z)')
%         ylabel('Concentration')
%         title('Comparison of Product Values with Offset')
%         outIter=0;
%     end
%     
%     % Update iterator
%     i=i+1;
% end
% 
% % Make final figures
% 
% 
% % Remove extra zeros if they exisit
% t=t(1:i);
% x=x(1:i);
% S=S(1:i);
% Lf=Lf(1:i);
% 
% end
