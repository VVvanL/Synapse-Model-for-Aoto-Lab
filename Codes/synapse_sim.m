%
%   simulation
%

function synapse_sim(gluPerVesicle, amparNo_a1, amparNo_a4, ratio, factor,...
                                  fileNo,period, group,...
                                  D_glu,transpDensity,cleftHeight,...               
                                  duration, release_zone, R_cleft, SpillOver)


%---- globals -------------------------------------------------------------
global  matlabFileInfo  %#ok

global  amparPeakDistr  %#ok
global  amparTags       %#ok
global  amparTagsVar    %#ok
global  amparStates     %#ok
global  amparMean       %#ok
global  amparVar        %#ok
global  amparVarS       %#ok
global  gluStates       %#ok
global  gluDistrES      %#ok

global  amparState_F        %#ok
global  amparState_B        %#ok
global  amparState_D        %#ok
global  amparState_O        %#ok

%---- random seeds --------------------------------------------------------
randn('state', 0);
rand('seed', 0);

%---- parameters ----------------------------------------------------------
clear('S'); % DONT REMOVE THIS LINE

%   S.savePath              = '';
    S.workspaceFilePrefix   = 'sim';    
    S.simulation            = 'new';
        
    
%---- statistics ----------------------------------------------------------
    S.runs                  = 1;
%---- pulses --------------------------------------------------------------
    S.waitBeforeFirstPulse  = 2; % in ms
    S.pulseNo               = 1;
    S.pulseToPulseTime      = 20;% ms, 1500;
    S.gluPerVesicle         = gluPerVesicle;
    S.vesiclesPerPulse      = 1;
    S.vesicleJitter         = 0;
    S.vesFillingFraction    = [1 1];
    S.vesFusionDuration     = duration;
    
    % Kinetic schemes for AMPARs
    % Form: [x1 x2 x3 x4], where this vector contains the fractions of AMPARs
    %       with
    %  x1 = slow-GluAs scheme
    %  x2 = fast-GluAs scheme 
    
%---- time steps ----------------------------------------------------------
    S.timeStepsPerMS        = round(1000*2); %=0.5 microseconds
%---- cleft geometry control constants ------------------------------------
    S.R_Cleft               = R_cleft;  % radius, nm
    
    S.R_ES                  = S.R_Cleft + 40;   % distance between pre/post-synaptic cylinders and glial sheath
    
    
    S.R_ReleaseZone         = release_zone; % radius of release zone, 0=center 

    S.cleftHeight           = cleftHeight; % cleft height at the calyx of Held
    S.esHeight              = 10000;
    S.esBinSize             = 50; % in nm, for the distribution of glu in the ES
%---- AMPAR control constants ---------------------------------------------
    S.D_ampar_PSD           = 00;       % AMPAR diffusion constant in nm^2/ms, on the PSD % Choquet  0.2 (mu m)^2/s = 200 nm^2/ms
    S.D_ampar_outside       = 00;       % inside the reservoir
    
    a = 0.1;                            % outside/inside density ratio
    p = 1.0;    
    S.P_Reflect_Inside      = 1-p*a; % 1-p*a for a <= 0.1, otherwise 1-p   % prob to be reflected hitting the PSD boundary from inside the PSD
    S.P_Reflect_Outside     = 1-p;   % 1-p/a for a >= 0.1, otherwise 1-p    

  
    S.R_BindingToAMPAR      = 10;   %in nm, binding radius, 2.5nm yields good agreement with analytically manageable situations
%---- glu control constants -----------------------------------------------
    S.D_glu                 = D_glu;   % glu diffusion constant, Franks, Sejnowski 0.2 (mu m)^2/ms = 200000 nm^2/ms, OTHER: 0.6  (mu m)^2/ms, 0.76 (mu m)^2/ms
%---- transporter control constants ---------------------------------------
    S.D_transp              = 0;        % NOT implemented, set to 0
    S.transpDensity         = transpDensity *1/1000^2; % in number per nm^2
    S.transpBindingTime     = 5.5;      % in ms, IGNORED
    S.transpPumpingTime     = 38;       % in ms, IGNORED

    S.absorbAtGlia          = 0;        % set to 1 if glu should be absorbed upon hitting the glial sheath
    S.absorbAtCleftBD       = 0;        % set to 1 if glu should be absorbed upon leaving the cleft

    S.R_BindingToTransp     = 20;       %in nm, binding radius

    % rates for the kinetic model of a transporter, see Franks, Stevens, Sejnowsky
    S.transpRate_assoc      = 1.8  *10^7 *1/1000;
    S.transpRate_dissoc     = 1.8  *10^2 *1/1000;
    S.transpRate_pump       = 1.8  *10^2 *1/1000;
    S.transpRate_getReady   = 2.57 *10^1 *1/1000;
%--------------------------------------------------------------------------  
    S.amparNo = amparNo_a1 + amparNo_a4;
    S.ratio = ratio;
    S.SpillOver = SpillOver;
    for i = 1:2
        if i == 1 % GluA1, stands for slow-GluAs
            
            % AMPAR number on the PSD at time 0           
            S.amparNo_PSD           = amparNo_a1;
            S.amparNo_outside       = 0;
            if S.amparNo_PSD == 0
                continue
            end            
            
            S.R_PSD                 = 60;  
            S.R_Reservoir           = 140;
                       
            S.kineticsFrac_PSD      = [1 0 0 0];
            S.kineticsFrac_Outside  = [1 0 0 0];
            S.group = group;
            S.subunit = 'GluA1';
            S.factor =factor;
            S.fileName = fileNo;               
            S.period = period;
            synapse_fun(S);     % run the simulation with the parameters specified in S
        else % GluA4, stands for fast-GluAs
            
            % AMPAR number on the PSD at time 0
            S.amparNo_PSD           = amparNo_a4;
            S.amparNo_outside       = 0;
            if S.amparNo_PSD == 0
                continue
            end  
            
            S.R_PSD                 = 60; 
            S.R_Reservoir           = 140; 
                              
            S.kineticsFrac_PSD      = [0 1 0 0]; % for AMPARs initially on the PSD
            S.kineticsFrac_Outside  = [0 1 0 0];
            S.group = group;
            S.subunit = 'GluA4';
            S.factor =factor;
            S.fileName = fileNo;  
            S.period = period;
            synapse_fun(S);     % run the simulation with the parameters specified in S
        end
    end