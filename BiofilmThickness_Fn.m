function [Lf,Vdet,Vg]= BiofilmThickness_Fn(S,Lf_old,mu,Kdet,mumax,dt,dz)
%This Function takes the substrate concentration at a given instant in time
%(and the corresponding equation for mu) as well as the old Biofilm
%thickness and computes a new Biofilm Thickness (Lf)

%Corresponding outputs for the growth and detachment rates are computed for
%a given instant of time.

%Biofilm Thickness
Lf=Lf_old+dt*(mu(S)*Lf_old-Kdet*Lf_old^2); %New Biofilm thickness at instant

%Detachment
Vdet=Kdet*Lf^2; %New %Velocity of mass leaving biofilm into bulk liquid

%Growth
Vg=0; %initial condition for loop
for i=1:length(S)-1
    Vg=Vg+dz*((mu(S(i))+mu(S(i+1)))/2); %trapezoidal integration method for
                                        %new growth velocity of biofilm
end
end