function [Lf,Vdet]=lf(Sb,Lf_old,dt,dz,param)
%This function takes the substrate concentration at a given instant in
%time and the old biofilm thickness to computes the growth velocity as well
%as the detachement velocity at a given instant in time. These results are
%used to compute a new biofilm thickness Lf

%Compute mean mu - growthrate
%muBar=mean(mu(Sb(:),param));

Kdet=param.Kdet;

%% Growth

% %Trapzoidal Integration Method
% Vg=0; %initial growth velocity
% for i=1:length(Sb)-1
%     Vg=Vg+dz*((mu(Sb(i),param)+mu(Sb(i+1),param))/2);
% end

%Midpoint Integration Method
Vg=0; %initial growth velocity
i=1;
while i<=length(Sb)-2
    Vg=Vg+dz*(mu(Sb(i),param)+mu(Sb(i+2),param));
    i=i+2;
end

%% Detachment
Vdet=Kdet*Lf_old^2; %New Velocity of mass leaving biofilm into bulk liquid

%% Biofilm Thickness
Lf=Lf_old+dt*(Vg-Vdet); %New Biofilm thickness at instant

end
