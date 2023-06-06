clear; 

runs = 160; % number of runs
stage =12;
ratio = 1; % ratio = GluA1 / (GluA1+GluA4)
             
              
%  group = 'Ctrl';
 group = 'KO';
 
 
% 
 
factor = 1;

D_glu = 0.4 * 1E6;% glu diffusion constant,
cleftHeight = 15;%% 28 for calyx of Held
transpDensity= 2000;%% 5000 for calyx of Held

duration = 0.2; % ms; duration of vesicle pore open 
gluPerVesicle = 3000;%% 8000 for calyx of Held
release_zone = 36; %nm; 

SpillOver = 0; % SpillOver = 0: there is no spillover
               % SpillOver = 1: include spillover with one, two, three, or
               % four neighboring synapses
                             
R_cleft = 400; % nm              
if strcmp(group, 'Ctrl')
    
    amparNo_a1 = 340;%1340; % 100 or 134, number of AMPARs，times ten for operation issue
    amparNo_a4 = 0;                           
    
elseif strcmp(group, 'KO')
    
    amparNo_a1 = 890; %88 or 89, number of AMPARs
    amparNo_a4 = 000; % 
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%      Run Simulations      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j = 1:runs
    fileNo = j;
    period = ['P',num2str(stage)]; 
    synapse_sim(gluPerVesicle, amparNo_a1, amparNo_a4, ratio, factor,...
                             fileNo,period,group,...
                             D_glu,transpDensity,cleftHeight,...                               
                             duration, release_zone, R_cleft, SpillOver);
end

beep;