function synapse_fun(SIM)
% LAST: SUN sep 6 10:00pm
%
% SIM structure with fields
%
%  simulation          ['new' 'skip' 'resume']   
%                                                           start new simulation (regardless of any existing results of prvious simulations)
%                                                           skip simulation if MATLAB file containing the simulation results is found
%                                                           resume old one (by loading temporary files from an interrupted run)
%                                                           Default is 'resume'.
%  getFileInfo          0 or 1                              get matlabFileInfo according to the parameters specified
%
%  savePath
%  <implement better> visualize   <<<DOES NOT WORK CURRENTLY>>>        0 or 1              0=dont visualize, 1=save .eps with amparStates and freeGlus
%


global  sim_signature

global  matlabFileInfo
global  timestep_internal   % this is for communication with printLoopStateInfo.m
global  timeSteps_internal  % this is for communication with printLoopStateInfo.m

global  amparPeakDistr      % for convenient access by caller
global  amparTags           % for convenient access by caller
global  amparTagsVar        % for convenient access by caller
global  amparStates         % for convenient access by caller
global  amparMean           % for convenient access by caller
global  amparVar            % for convenient access by caller
global  amparVarS           % for convenient access by caller
global  gluStates           % for convenient access by caller
global  gluDistrES          % for convenient access by caller

global  amparState_F        %#ok
global  amparState_B        %#ok
global  amparState_D        %#ok
global  amparState_O        %#ok

%---- additional globals for visualization routines -------------------
global  gluState_FreeES
global  gluState_FreeCleft
global  gluState_BoundToTransp
global  gluState_PumpedOut

%---- internal constants --------------------------------------------------
    kinetics_GluA1            = 1;
    kinetics_GluA4         = 2;
    
    KINETICS_NO = 2; % number of kinetic schemes

    % AMPAR states
    amparState_F = 1; %#ok % free
    amparState_B = 2; %#ok % bound
    amparState_D = 3; %#ok % desensitized
    amparState_O = 4; %#ok % open
    
    AMPARSTATES_NO = 4;% free, bound, desensitized, open

    
    MAX_AMPARSTATES_NO = 9; % number of AMPAR states (MAXIMUM used for all possible schemes)
    
    amparState_InCleft_OnPSD                =  1;   %#ok
    amparState_InCleft_OutsidePSD           =  2;   %#ok
    amparState_InES                         =  3;   %#ok    
    amparState_First_InsidePSD_Now_Outside  =  4;   %#ok
    amparState_First_OutsidePSD_Now_Inside  =  5;   %#ok
    amparState_Free_OnPSD                   =  6;   %#ok
    amparState_Free_Outside                 =  7;   %#ok
    amparState_Bound_OnPSD                  =  8;   %#ok
    amparState_Bound_Outside                =  9;   %#ok
    amparState_Desensitized_OnPSD           = 10;   %#ok
    amparState_Desensitized_Outside         = 11;   %#ok
    amparState_Open_OnPSD                   = 12;   %#ok
    amparState_Open_Outside                 = 13;   %#ok

    AMPARTAGS_NO                            = 13;   % number of AMPAR states;

    % glu states
    % NOTE that the variable "gluStates" must have no-of-possible-glu-states 
    %      entries, start with 1, because the states index an array
    gluState_FreeCleft           = 1; %#ok  % glu unbound and in happily diffusing in cleft
    gluState_FreeES              = 2; %#ok  % glu unbound and in happily diffusing in cleft
    gluState_BoundToAMPAR        = 3; %#ok  % bound to AMPAR
    gluState_BoundToTransp       = 4; %#ok  % bound to transporter
    gluState_Unreleased          = 5; %#ok  % not yet released
    gluState_Absorbed            = 6; %#ok  % absorbed---diffused out of cleft
    gluState_PumpedOut           = 7; %#ok  % absorbed---diffused out of cleft    
    gluState_UnbindingFromTransp = 8; %#ok  % use if glu has just been unbound from transporter,
    gluState_UnbindingFromAMPAR  = 9; %#ok  % use if glu has just been unbound from AMPAR
                                            % to indicate that the glu has to be re-inserted
                                            % into diffusion space correctly
    GLUSTATES_NO                 = 9; %#ok  % number of AMPAR states

    % transporter states
    transpState_unbound = 1; %#ok
    transpState_bound   = 2; %#ok
    transpState_pumping = 3; %#ok

    % further internal constants
    INVALID_INDEX = -1; % used to point to an out-of-index value


    % file prefixes
    tempFilePrefix   = 'DELETE_ME___THIS_RUN_';
    %singleRunPrefix = 'SINGLE_RUN_';
%---- set default values for fields not defined ---------------------------
    sim_signature = 'hitting (counting hits) _stochastic_ transporters, several kinetic schemes';

    %---- Runge-Kutta control constants --------------------
    if isfield(SIM, 'RK_ErrorTol') == 0
        SIM.RK_ErrorTol = 0.0001;           % to create variable
    end
    if isfield(SIM, 'RK_MaxRef') == 0
        SIM.RK_MaxRef = 1/2^25;             % to create variable
    end
    
    %---- general ------------------------------------------
    if isfield(SIM, 'timestep') == 0
        SIM.timestep = 1;                           % to create variable
    end
    if isfield(SIM, 'visualize') == 0
        SIM.visualize = 0;
    end
    if isfield(SIM, 'savePath') == 0
        SIM.savePath = cd;
    else
        if isempty(SIM.savePath) ~= 0
            SIM.savePath = cd;
        end
        if isdir(SIM.savePath) == 0
            error('____SIM.savepath is not a directory____');
        end
    end

    if isfield(SIM, 'releaseZoneCenter') == 0
        SIM.releaseZoneCenter     = 0;
    end

    if isfield(SIM, 'simulation') == 0
        SIM.simulation = 'resume';
    end
    
    if isfield(SIM, 'amparPeakRes') == 0
        SIM.amparPeakRes = 2;
    end

    if isfield(SIM, 'amparPeakRes') == 0
        SIM.amparPeakRes = 2;
    end

    if isfield(SIM, 'absorbAtGlia') == 0
        SIM.absorbAtGlia = 0;
    end
    if isfield(SIM, 'absorbAtCleftBD') == 0
        SIM.absorbAtCleftBD = 0;
    end

    if isfield(SIM, 'vesFusionDuration') == 0
        SIM.vesFusionDuration = 0;
    end

    %---- transporter kinetic model rates ---------------
        % all constants taken from franks, stevens, sejnowski "independent sources"
    if isfield(SIM, 'transpRate_assoc') == 0
        SIM.transpRate_assoc      = 1.8  *10^7 *1/1000; % in M^-1 s^-1
    end
    if isfield(SIM, 'transpRate_dissoc') == 0
        SIM.transpRate_dissoc     = 1.8  *10^2 *1/1000;
    end
    if isfield(SIM, 'transpRate_pump') == 0
        SIM.transpRate_pump       = 1.8  *10^2 *1/1000;
    end
    if isfield(SIM, 'transpRate_getReady') == 0
        SIM.transpRate_getReady   = 2.57 *10^1 *1/1000;
    end
    
    
    %---- kinetic models --------------------------------
    if isfield(SIM, 'kineticsFrac_Outside') == 0
        SIM.kineticsFrac_Outside = [1 0 0 0];
    end

    % for backward compatibility
    if isfield(SIM, 'kineticModel') ~= 0
        kineticsFracAtPulse = [0 0 0 0];
        if     strcmpi(SIM.kineticModel, 'GluA1')
            kineticsFracAtPulse(kinetics_GluA1) = 1;
        elseif strcmpi(SIM.kineticModel, 'GluA4')
            kineticsFracAtPulse(kinetics_GluA4) = 1;
        end
        SIM.kineticsFrac_PSD     = kineticsFracAtPulse;
        SIM.kineticsFrac_Outside = kineticsFracAtPulse;
    end

    
    % vesicle filling range in percent
    % if a single value is given, make it a [lower upper] array
    if isfield(SIM, 'vesFillingFraction') == 0
        SIM.vesFillingFraction = [1 1];
    elseif length(SIM.vesFillingFraction) == 1
        a = SIM.vesFillingFraction;
        if a == 0
            error('SIM.vesFillingFraction: vesicle filling fraction cannot be 0.');
        end
        SIM.vesFillingFraction = [a a];
    end
    
    SIM.vesFillingFraction = SIM.vesFillingFraction / max(SIM.vesFillingFraction);
    
%---- set internal variables ----------------------------------------------
    SIM.finalTime = (SIM.waitBeforeFirstPulse + SIM.pulseNo*SIM.pulseToPulseTime)*1;      %in ms
    SIM.timeSteps = round(SIM.finalTime*SIM.timeStepsPerMS); % +1 == first time step

    SIM.amparNo = SIM.amparNo_PSD + SIM.amparNo_outside; % all ampars together

    SIM.timeStepSize = SIM.finalTime/SIM.timeSteps;

    SIM.upperLidES =  SIM.esHeight/2 + SIM.cleftHeight/2;  % upper lid of the ES
    SIM.lowerLidES = -SIM.esHeight/2 + SIM.cleftHeight/2;  % lower lid of the ES
    SIM.upperLid   = SIM.cleftHeight;
    SIM.lowerLid   = 0;
    
    SIM.esBinNo = ceil((SIM.upperLidES-SIM.lowerLidES) / SIM.esBinSize);
    
    if SIM.transpDensity == 0
        SIM.transp2transpDist = 0;
        SIM.transpHor = 0;
        SIM.transpVer = 0;
        SIM.transpNo = 0;    
    else
        
% % %         if SIM.absorbAtGlia ~= 0  ||  SIM.absorbAtCleftBD ~= 0
% % %             error('transpDensity > 0 while absorption flags specified at the same time.');
% % %         end
        
        SIM.transp2transpDist = 1/sqrt(SIM.transpDensity); % average transp-to-transp distance in a 2D grid

        SIM.transpHor = round(2*pi*SIM.R_ES/SIM.transp2transpDist);     % horizontal num on grid
        SIM.transpVer = round((SIM.upperLidES-SIM.lowerLidES)/SIM.transp2transpDist);  % vertical num on grid

        SIM.transpNo = SIM.transpHor*SIM.transpVer;
    end

    
    % probability to to make an transp hit when hitting in the transporter
    % disk with radius SIM.R_BindingToTransp
    
    SIM.P_bindingToTransp = (10/6.023*SIM.transpRate_assoc*sqrt(SIM.timeStepSize*pi/SIM.D_glu)) / (pi*SIM.R_BindingToTransp^2); %#ok

    if SIM.P_bindingToTransp > 1
        error('transporter association rate too high (increase transp binding radius)');
    end
    
    %----prepare insertion distribution look-up table----------------------
    % NOTE: We need D_glu and timeStepSize defined at this point
    insertDistrTableLen = 2^12; % length of table. NEEDS TO BE POWER OF TWO
                                % because we want to use AND rather than MOD
                                % in C (DONT subtract 1) !!!
    
    insertDistrTable = zeros(insertDistrTableLen, 1);
    
    maxRefineStepSize = 1/insertDistrTableLen /100;
    ly = 0; % for y=0, fInv(y)=0
    uy = 3; % for y=3, fInv(y)=0.9999 has (almost) reached its full range [0, 1]
    for i = 1 : insertDistrTableLen
        RefineStepSize = uy - ly;
        while 1
            RefineStepSize = RefineStepSize/2;
            y = ly + RefineStepSize;
            fInv = sqrt(pi)* y *(1 - erf(y)) + 1 - exp(-y^2);

            if (i-1)/insertDistrTableLen <= fInv  &&  fInv < i/insertDistrTableLen
                break;
            end

            if RefineStepSize < maxRefineStepSize
                break;
            end
        end
        ly = y;
        insertDistrTable(i) = y * 2*sqrt(SIM.D_glu*SIM.timeStepSize); % adapt to D_glu, timeStepSize
    end
    

%---- book keeping variables ----------------------------------------------
%                            ----------------------------------------------
%                            ----------------------------------------------
%                            ----------------------------------------------
%                            ----------------------------------------------
    
    amparPeakDistr = zeros(SIM.amparNo*SIM.amparPeakRes, AMPARSTATES_NO, 'double');   %#ok %
    
    amparTags    = zeros(SIM.timeSteps, AMPARTAGS_NO,   'double');   %#ok % array of AMPAR states
    amparTagsVar = zeros(SIM.timeSteps, AMPARTAGS_NO,   'double');   %#ok % array of AMPAR states
   
    amparStates  = zeros(SIM.timeSteps, AMPARSTATES_NO, 'double');   %#ok % array of AMPAR states
    amparMean    = zeros(SIM.timeSteps, AMPARSTATES_NO, 'double');   %#ok % array of AMPAR states
    amparVar     = zeros(SIM.timeSteps, AMPARSTATES_NO, 'double');   %#ok % AMPAR variance (whole response)
    amparVarS    = zeros(SIM.timeSteps, AMPARSTATES_NO, 'double');   %#ok % AMPAR variance (single channel)
    
    gluStates   = zeros(SIM.timeSteps, GLUSTATES_NO,    'uint32');   %#ok % no of glus with same states
    gluDistrES  = zeros(SIM.timeSteps, SIM.esBinNo,     'uint32');   % % no of glus in vertical rings-shaped bins in ES in z-direction



%---- kinetic schemes -----------------------------------------------------

    % structure of `kinetics'
    %    name
    %    amparState_F       states in ampar.stateDistr that contribute to `free' states
    %    amparState_B
    %    amparState_D
    %    amparState_O       states in ampar.stateDistr that contribute to `open' states
    %    Q_transition       transition matrix
    %    P_transition       timestep-1 semigroup (if there is no glu bound)
    %    P_bindingToAMPAR   prob to bind upon touching by a glu (for the glu-dependent rates)
    kinetics = struct(...
                    'name',             zeros(KINETICS_NO, 60), ...
                    'amparState_F',     zeros(KINETICS_NO, MAX_AMPARSTATES_NO, 'int32'), ...
                    'amparState_B',     zeros(KINETICS_NO, MAX_AMPARSTATES_NO, 'int32'), ...
                    'amparState_D',     zeros(KINETICS_NO, MAX_AMPARSTATES_NO, 'int32'), ...
                    'amparState_O',     zeros(KINETICS_NO, MAX_AMPARSTATES_NO, 'int32'), ...
                    'Q_transition',     zeros(KINETICS_NO, MAX_AMPARSTATES_NO, MAX_AMPARSTATES_NO, 'double'), ...
                    'P_transition',     zeros(KINETICS_NO, MAX_AMPARSTATES_NO, MAX_AMPARSTATES_NO, 'double'), ...
                    'P_bindingToAMPAR', zeros(KINETICS_NO, MAX_AMPARSTATES_NO, 'double'));

    %---- GluA1 scheme --------------------------------------
    name = 'GluA1';
    kinetics.name(kinetics_GluA1, 1:length(name)) = name;
        
        % NOTE that the following numbers correspond to the indices in the Q and P matrix
        amparState_C0   = 1;
        amparState_C1   = 2;
        amparState_C2   = 3;
        amparState_D3   = 4;
        amparState_D4   = 5;
        amparState_D5   = 6;
        amparState_D6   = 7;
        amparState_D7   = 8;
        amparState_O   = 9;

    kinetics.amparState_F(kinetics_GluA1, amparState_C0 ) = 1;%#ok  F
    kinetics.amparState_B(kinetics_GluA1, amparState_C1) = 1; %#ok  B
    kinetics.amparState_B(kinetics_GluA1, amparState_C2) = 1; %#ok  B
    kinetics.amparState_D(kinetics_GluA1, amparState_D3) = 1; %#ok  D
    kinetics.amparState_D(kinetics_GluA1, amparState_D4) = 1; %#ok  D
    kinetics.amparState_D(kinetics_GluA1, amparState_D5) = 1; %#ok  D
    kinetics.amparState_D(kinetics_GluA1, amparState_D6) = 1; %#ok  D
    kinetics.amparState_D(kinetics_GluA1, amparState_D7) = 1; %#ok  D
    kinetics.amparState_O(kinetics_GluA1, amparState_O)  = 1; %#ok  O
        

%%%% 2015-Calyx of Held
    
        k1_a1  = SIM.factor*1.5E7*1         *1/1000;%#ok   % ms^-1 * M^-1 (ie, only, if singly liganded)
        k_1_a1 = SIM.factor*4.323E3*1          *1/1000;% ms^-1
        k2_a1  = SIM.factor*4E6*1              *1/1000;%#ok   % ms^-1 * M^-1 (ie, only, if singly liganded)
        k_2_a1 = SIM.factor*1.7201E4*1         *1/1000;% ms^-1
        k3_a1  = SIM.factor* 1E4*1          *1/1000;% ms^-1
        k_3_a1 = SIM.factor* 2.7E3*1         *1/1000;% ms^-1
        k4_a1  = SIM.factor* 1          *1/1000;% ms^-1
        k_4_a1 = SIM.factor*0.42183*1          *1/1000;% ms^-1
        k5_a1  = SIM.factor*1.9863E7*1        *1/1000;%#ok   % ms^-1 * M^-1 (ie, only, if singly liganded)
        k_5_a1 = SIM.factor*18.392E3*1         *1/1000;% ms^-1
        k6_a1  = SIM.factor*8.48141*1          *1/1000;% ms^-1
        k_6_a1 = SIM.factor*12.13*1          *1/1000;% ms^-1
        k7_a1  = SIM.factor* 25.85*1            *1/1000;% ms^-1
        k_7_a1 = SIM.factor* 93.287*1           *1/1000;% ms^-1
        k8_a1  = SIM.factor*1.1813E3*1           *1/1000;% ms^-1
        k_8_a1 = SIM.factor*280.35*1           *1/1000;% ms^-1
        k9_a1  = SIM.factor* 380.434*1          *1/1000;% ms^-1
        k_9_a1 = SIM.factor*1.944E1*1           *1/1000;% ms^-1
        k10_a1 = SIM.factor*0.59666*1           *1/1000;% ms^-1
        k_10_a1= SIM.factor*0.17436 *1          *1/1000;% ms^-1
        k11_a1 = SIM.factor* 2 *1             *1/1000;% ms^-1
        k_11_a1= SIM.factor*5 *1          *1/1000;% ms^-1
        
        Q = [
        %       C0         C1         C2         D3         D4         D5         D6         D7         O
            00000,       0.0,         0,         0,         0,         0,         0,         0,         0 ;    % C0
            
           k_1_a1,     00000,       0.0,     k8_a1,         0,         0,         0,         0,         0 ;    % C1
             
                0,    k_2_a1,     00000,       0.0,     k9_a1,         0,         0,         0,     k3_a1 ;    % C2
                
                0,    k_8_a1,         0,     00000,       0.0,         0,         0,         0,         0 ;    % D3
                 
                0,         0,    k_9_a1,    k_5_a1,     00000,     k6_a1,         0,         0,         0 ;    % D4
                
                0,         0,         0,         0,    k_6_a1,     00000,     k7_a1,         0,   k_10_a1 ;    % D5
                
                0,         0,         0,         0,         0,    k_7_a1,     00000,   k_11_a1,         0 ;    % D6
                
                0,         0,         0,         0,         0,         0,    k11_a1,     00000,    k_4_a1 ;    % D7
                
                0,         0,    k_3_a1,         0,         0,    k10_a1,         0,     k4_a1,     00000 ];    % O
        

        % prepare the transition matrix, where no trans        
        Q = Q - diag(diag(Q));          % set diagonal to zero
        Q = Q - diag(sum(Q, 2));        % prepare correct diagonal
    kinetics.Q_transition(kinetics_GluA1, :, :) = Q;               %#ok
    
    kinetics.P_transition(kinetics_GluA1, :, :) = expm(Q* SIM.timeStepSize);  %#ok % transition prob matrix

    kinetics.P_bindingToAMPAR(kinetics_GluA1, amparState_C0) = (10/6.023*k1_a1*sqrt(SIM.timeStepSize*pi/SIM.D_glu)) / (pi*SIM.R_BindingToAMPAR^2); %#ok
    kinetics.P_bindingToAMPAR(kinetics_GluA1, amparState_C1) = (10/6.023*k2_a1*sqrt(SIM.timeStepSize*pi/SIM.D_glu)) / (pi*SIM.R_BindingToAMPAR^2); %#ok
    kinetics.P_bindingToAMPAR(kinetics_GluA1, amparState_D3) = (10/6.023*k5_a1*sqrt(SIM.timeStepSize*pi/SIM.D_glu)) / (pi*SIM.R_BindingToAMPAR^2); %#ok

    %---- GluA4 (no TARP ligation) scheme --------------------------------
    name = 'GluA4';
    kinetics.name(kinetics_GluA4, 1:length(name)) = name;
        
        % NOTE that the following numbers correspond to the indices in the Q and P matrix

        amparState_C0   = 1;
        amparState_C1   = 2;
        amparState_C2   = 3;
        amparState_D3   = 4;
        amparState_D4   = 5;
        amparState_D5   = 6;
        amparState_D6   = 7;
        amparState_D7   = 8;
        amparState_O   = 9;

    kinetics.amparState_F(kinetics_GluA4, amparState_C0 ) = 1;%#ok  F
    kinetics.amparState_B(kinetics_GluA4, amparState_C1) = 1; %#ok  B
    kinetics.amparState_B(kinetics_GluA4, amparState_C2) = 1; %#ok  B
    kinetics.amparState_D(kinetics_GluA4, amparState_D3) = 1; %#ok  D
    kinetics.amparState_D(kinetics_GluA4, amparState_D4) = 1; %#ok  D
    kinetics.amparState_D(kinetics_GluA4, amparState_D5) = 1; %#ok  D
    kinetics.amparState_D(kinetics_GluA4, amparState_D6) = 1; %#ok  D
    kinetics.amparState_D(kinetics_GluA4, amparState_D7) = 1; %#ok  D
    kinetics.amparState_O(kinetics_GluA4, amparState_O)  = 1; %#ok  O   
        
%%%% 2015-Calyx of Held
        k1_a4  = SIM.factor*2.0E7*1.1         *1/1000;%#ok   % ms^-1 * M^-1 (ie, only, if singly liganded)
        k_1_a4 = SIM.factor*4.323E3 *1.1         *1/1000;% ms^-1
        k2_a4  = SIM.factor*4E6*1.1              *1/1000;%#ok   % ms^-1 * M^-1 (ie, only, if singly liganded)
        k_2_a4 = SIM.factor*1.7201E4*1.1         *1/1000;% ms^-1
        k3_a4  = SIM.factor* 2.5E4*1.1          *1/1000;% ms^-1
        k_3_a4 = SIM.factor* 8E3*1.1         *1/1000;% ms^-1
        k4_a4  = SIM.factor*  1E2*1.1           *1/1000;% ms^-1
        k_4_a4 = SIM.factor*105.45*1.1          *1/1000;% ms^-1
        k5_a4  = SIM.factor*1.9863E7*1.1         *1/1000;%#ok   % ms^-1 * M^-1 (ie, only, if singly liganded)
        k_5_a4 = SIM.factor*13.7937E3*1.1         *1/1000;% ms^-1
        k6_a4  = SIM.factor*282*1.1          *1/1000;% ms^-1
        k_6_a4 = SIM.factor*485.05*1.1          *1/1000;% ms^-1
        k7_a4  = SIM.factor* 51.7*1.1             *1/1000;% ms^-1
        k_7_a4 = SIM.factor* 46.644*1.1           *1/1000;% ms^-1
        k8_a4  = SIM.factor*885.99*1.1          *1/1000;% ms^-1
        k_8_a4 = SIM.factor*280.35*1.1           *1/1000;% ms^-1
        k9_a4  = SIM.factor* 380.434 *1.1         *1/1000;% ms^-1
        k_9_a4 = SIM.factor*1.944E1*1.1           *1/1000;% ms^-1
        k10_a4 = SIM.factor*5.9666*1.1           *1/1000;% ms^-1
        k_10_a4= SIM.factor*1.7436*1.1           *1/1000;% ms^-1
        k11_a4 = SIM.factor* 2000*1.1              *1/1000;% ms^-1
        k_11_a4= SIM.factor*500*1.1           *1/1000;% ms^-1


         Q = [
        %       C0         C1         C2         D3         D4         D5         D6         D7         O
            00000,       0.0,         0,         0,         0,         0,         0,         0,         0 ;    % C0
            
           k_1_a4,     00000,       0.0,     k8_a4,         0,         0,         0,         0,         0 ;    % C1
             
                0,    k_2_a4,     00000,       0.0,     k9_a4,         0,         0,         0,     k3_a4 ;    % C2
                
                0,    k_8_a4,         0,     00000,       0.0,         0,         0,         0,         0 ;    % D3
                 
                0,         0,    k_9_a4,    k_5_a4,     00000,     k6_a4,         0,         0,         0 ;    % D4
                
                0,         0,         0,         0,    k_6_a4,     00000,     k7_a4,         0,   k_10_a4 ;    % D5
                
                0,         0,         0,         0,         0,    k_7_a4,     00000,   k_11_a4,         0 ;    % D6
                
                0,         0,         0,         0,         0,         0,    k11_a4,     00000,    k_4_a4 ;    % D7
                
                0,         0,    k_3_a4,         0,         0,    k10_a4,         0,     k4_a4,     00000 ];    % O

        % prepare the transition matrix, where no trans        
        Q = Q - diag(diag(Q));          % set diagonal to zero
        Q = Q - diag(sum(Q, 2));        % prepare correct diagonal
    kinetics.Q_transition(kinetics_GluA4, :, :) = Q;               %#ok
    
    kinetics.P_transition(kinetics_GluA4, :, :) = expm(Q* SIM.timeStepSize);  %#ok % transition prob matrix
    
    % R = sqrt(10/6.023*k_plus_1*sqrt(timeStepSize/(D_glu*PI)));
    kinetics.P_bindingToAMPAR(kinetics_GluA4, amparState_C0) = (10/6.023*k1_a4*sqrt(SIM.timeStepSize*pi/SIM.D_glu)) / (pi*SIM.R_BindingToAMPAR^2); %#ok
    kinetics.P_bindingToAMPAR(kinetics_GluA4, amparState_C1) = (10/6.023*k2_a4*sqrt(SIM.timeStepSize*pi/SIM.D_glu)) / (pi*SIM.R_BindingToAMPAR^2); %#ok
    kinetics.P_bindingToAMPAR(kinetics_GluA4, amparState_D3) = (10/6.023*k5_a4*sqrt(SIM.timeStepSize*pi/SIM.D_glu)) / (pi*SIM.R_BindingToAMPAR^2); %#ok

%---- define vesicle releases (pulses) -----------------------------------------
    % structure of vesicleRelease
    %     start                   % in ms, start of vesicle release
    %     duration                % in ms, =0 immediate release
    %     gluNo                   % number released glu's per vesicle
    %     position                % release position in Cleft, =[0, 0, cleftHeight] (fixed)
    %     timeStep_Start          % (internal)--- time step at start
    %     timeStep_ReleaseRate    % (internal)--- number of glus released per time step % NOT IMPLEMENTED
    %     timeStep_End            % (internal)--- time step of end                      % NOT IMPLEMENTED

    SIM.vesicleNo = SIM.pulseNo*SIM.vesiclesPerPulse;

    vesicle = struct(...
                    'start',                zeros(1, SIM.vesicleNo), ...
                    'duration',             zeros(1, SIM.vesicleNo), ...
                    'gluNo',                zeros(1, SIM.vesicleNo), ...
                    'gluNo_Released',       zeros(1, SIM.vesicleNo), ...
                    'position',             zeros(3, SIM.vesicleNo), ...                    
                    'timeStep_Start',       zeros(1, SIM.vesicleNo), ...
                    'timeStep_ReleaseRate', zeros(1, SIM.vesicleNo), ...
                    'timeStep_End',         zeros(1, SIM.vesicleNo) );

    % define vesicle releases
    for i_pulse = 1 : SIM.pulseNo
        for i_vesiclesPerPulse = 1 : SIM.vesiclesPerPulse
            i_vesicle = (i_pulse-1)*SIM.vesiclesPerPulse + i_vesiclesPerPulse;
    
            vesicle.gluNo_Released(i_vesicle) = 0;
            
            vesicle.start(i_vesicle) = SIM.waitBeforeFirstPulse + ...
                (i_pulse-1)*SIM.pulseToPulseTime + (i_vesiclesPerPulse-1)*SIM.vesicleJitter; % in ms
            vesicle.duration(i_vesicle) = SIM.vesFusionDuration;    % in ms, =0 immediate release
            
            % NOTE: vesicle load is set to 0.
            % below, in the statistics loop, the vesicle load is chosen randomly
            vesicle.gluNo(i_vesicle) = 0;
    
            % set release position to top-center. It is randomized below
            rad = 0;
            ang = 0;
            vesicle.position(:, i_vesicle) = [rad*cos(2*pi*ang), rad*sin(2*pi*ang), SIM.upperLid]';
        end
    end

    % initialize vesicle releases
    for i_vesicle = 1 : SIM.vesicleNo
        vesicle.timeStep_Start(i_vesicle) = floor(vesicle.start(i_vesicle) / SIM.timeStepSize); % start timestep at 0
        vesicle.timeStep_ReleaseRate(i_vesicle) = 0; % NOTE: this is initialized in the statistics loop for every run

        vesicle.timeStep_End(i_vesicle) = floor((vesicle.start(i_vesicle) + vesicle.duration(i_vesicle)) / SIM.timeStepSize);
    end

%---- MATLAB file storage administration ----------------------------------

% % % if SIM.R_Reservoir == SIM.R_PSD
% % %     densOutside = -1;
% % % else
% % %     densOutside = SIM.amparNo_outside / (SIM.amparNo_PSD/(pi*SIM.R_PSD^2)* (pi*SIM.R_Reservoir^2-pi*SIM.R_PSD^2));
% % % end
% % % 
% % % fileInfo = strcat(...
% % %                 'amparNo_PSD=',             num2str(SIM.amparNo_PSD), ...
% % %                 ',ampar_dens_outside=',     num2str(densOutside, '%.3f'), ...
% % %                 ',D_amparPSD=',             num2str(SIM.D_ampar_PSD), ...
% % %                 ',D_amparoutside=',         num2str(SIM.D_ampar_outside), ...
% % %                 ',P_Refl_Ins=',             num2str(SIM.P_Reflect_Inside), ...
% % %                 ',P_Refl_Outs=',            num2str(SIM.P_Reflect_Outside), ...
% % %                 ',vesiclesPerPulse=',       num2str(SIM.vesiclesPerPulse), ...
% % %                 ',R_PSD=',                  num2str(SIM.R_PSD), ...
% % %                 ',R_Cleft=',                num2str(SIM.R_Cleft), ...
% % %                 ',h_cleft=',                num2str(SIM.cleftHeight), ...
% % %                 ',T=',                      num2str(SIM.finalTime, '%.4f'), ...
% % %                 ',runs=',                   num2str(SIM.runs));

transpDensityStr = num2str(SIM.transpDensity*1000^2, '%05.f');

if SIM.absorbAtGlia == 1
    transpDensityStr = 'absGlia';
end
if SIM.absorbAtCleftBD == 1
    transpDensityStr = 'absCleftBD';
end


fileInfo = strcat(...
                '',                  num2str(SIM.period),...
                '_Subunit=',              SIM.subunit,...
                ',amparNo=',             num2str(SIM.amparNo), ...
                ',gluNo=',              num2str(SIM.gluPerVesicle),...
                ',amparNoNC=',        num2str(SIM.amparNo_PSD),...
                ',amparNoOutside=',        num2str(SIM.amparNo_outside),...
                ',VesicleDuration=',   num2str(SIM.vesFusionDuration),...
                ',transpDensity=',      transpDensityStr, ...
                ',R_NC=',              num2str(SIM.R_PSD), ...
                ',R_SC=',              num2str(SIM.R_Reservoir), ...
                ',R_release=',          num2str(SIM.R_ReleaseZone),...
                ',R_ES=',               num2str(SIM.R_ES),...
                ',R_Cleft=',            num2str(SIM.R_Cleft),...
                ',D_glu=',               num2str(SIM.D_glu));

            
matlabFileInfo = fullfile(SIM.savePath, strcat('matlab_', SIM.workspaceFilePrefix, '_', fileInfo, '.mat'));

if isfield(SIM, 'getFileInfo') ~= 0
    if SIM.getFileInfo ~= 0
        return;
    end
end

fid = fopen(matlabFileInfo);
if fid ~= -1  &&  strcmpi(SIM.simulation, 'new') == 0
    fclose(fid);
    fprintf('>>>> skipping sim (set SIM.simulation = ''new'' to rerun)...\n');
    
    % @@@ clear memory, ie, book keeping variables.
    return;
end

runsStart = 1;
fileList = dir(strcat(tempFilePrefix, '*.mat'));
if ~isempty(fileList)
    for i = 1 : length(fileList)
%         [pathstr, name, ext, versn] = fileparts(fileList(i).name); %#ok
        [pathstr, name, ext] = fileparts(fileList(i).name); %#ok
        fileNum = sscanf(name, strcat('%*', num2str(length(tempFilePrefix)), 'c%d'));
        runsStart = max([fileNum runsStart]);
    end
    if runsStart < 1  ||  runsStart > SIM.runs
        runsStart = 1;
    else
        fprintf('>>>> resuming simulation from temp file');
        tempName = fileList(runsStart).name;
        load(fileList(runsStart).name);
        % load twice
        load(tempName);
        clear('tempName'); % make sure that `tempName' does not appear in the mat file
        runsStart = runsStart + 1;
        fprintf(', starting run %d\n', runsStart);
    end
end

runsBegin = runsStart;
clear('runsStart'); % make sure that `runsStart' DOES NOT appear in a temp file

%    XXX                         XXX
%    XXX                         XXX
%    XXX  AXXXXM  MM  MAA   AAM  XXX
%    XXX  MM      MM  MMAA AAMM  XXX
%    XXX  VXXXXA  MM  MM VXV MM  XXX
%  VVVVV      MM  MM  MM  V  MM  VVVVV
%   VVVV      MM  MM  MM     MM  VVVV
%    VVV  MXXXXV  MM  MM     MM  VVV
%     VV                         VV
%      V                         V

%---- statistics loop ----------------------------------------------------------

fprintf('-------- start -----                            %s\n', datestr(now)); % INFO OUTPUT

for i_runs = runsBegin : SIM.runs
    fprintf('run: %3.0f (of %3.0f)                %% -----', i_runs, SIM.runs); % INFO OUTPUT


    %---- AMPAR initialization -------------------------------------------------
    fprintf('\b\b\b\b\bAMPAR'); % INFO OUTPUT
    
    dt_ampar_PSD     = sqrt(2*SIM.D_ampar_PSD    *SIM.timeStepSize); %#ok %Brownian motion variance per timestep
    dt_ampar_outside = sqrt(2*SIM.D_ampar_outside*SIM.timeStepSize); %#ok %Brownian motion variance per timestep
    % array of all AMPARs,
    % AMPAR is a structure with fields
    %   traj          discretized trajectory
    %   stateDistr    current state distribution (of the kinetics Markov chain)
    %   kinetics      kinetic model used for every AMPAR
    %   ligand1       index of glu liganding first
    %   ligand2       index of glu liganding second
    %   tag           position of the AMPAR (PSD, outside PSD but in cleft, outside cleft)
    %   tag_initial   position of the AMPAR at time 0
    %   T             waiting time (exponentially distributed by state)
    ampar = struct(...
        'traj',         zeros(3, 1, SIM.amparNo), ... % 3=dimension, 2d diffusion but keep z-coord constant
        'stateDistr',   zeros(MAX_AMPARSTATES_NO, SIM.amparNo, 'double'), ...
        'kinetics',     zeros(1, SIM.amparNo, 'int8'), ...
    ...%'ligand1',      ones (1, SIM.amparNo, 'int32')*INVALID_INDEX, ...
    ...%'ligand2',      ones (1, SIM.amparNo, 'int32')*INVALID_INDEX, ...
        'hits',         zeros(1, SIM.amparNo, 'int32'), ...
        'currentState', zeros(1, SIM.amparNo, 'uint32'), ...
        'tag',          zeros(1, SIM.amparNo, 'uint8'), ...
        'tag_initial',  zeros(1, SIM.amparNo, 'uint8'), ...
        'T',            zeros(1, SIM.amparNo));

    % distribute AMPARs uniformly on PSD.
    if SIM.SpillOver == 1
        NND = 200;  % nearst-neighbor distance between central AMPARs and synapses nearby
                    % 93 +93 is around 200 nm, 93 is the radius of the second circle (SSD). 
    else % SIM.SpillOver == 0
        NND = 0;
    end
    
%     for i_ampar = 1 : SIM.amparNo_PSD
%         rad = sqrt(unifrnd(0, SIM.R_PSD^2)); % (Take square root of unif. distributed random value)
%         ang = unifrnd(0, 1);
%         ampar.traj(:, i_ampar) = [NND + rad*cos(2*pi*ang), rad*sin(2*pi*ang), SIM.lowerLid]';
%         ampar.tag(i_ampar)     = amparState_InCleft_OnPSD;
%     end
% 
%     % distribute AMPARs uniformly outside PSD
%     for i_ampar = 1 : SIM.amparNo_outside
%         rad = sqrt(unifrnd(SIM.R_PSD^2, SIM.R_Reservoir^2));
%         ang = unifrnd(0, 1);
%         ampar.traj(:, i_ampar +SIM.amparNo_PSD) = [NND + rad*cos(2*pi*ang), rad*sin(2*pi*ang), SIM.lowerLid]';
%         
%         if rad < SIM.R_Cleft
%             ampar.tag(i_ampar +SIM.amparNo_PSD) = amparState_InCleft_OutsidePSD;
%         else
%             ampar.tag(i_ampar +SIM.amparNo_PSD) = amparState_InES;
%         end
%     end


if strcmp(SIM.group,'Ctrl')
    if strcmp(SIM.subunit,'GluA1')
        %%%%%%% A1 ctrl
            for i_ampar = 1 : 416
                rad = sqrt(unifrnd(0, 34^2)); % (Take square root of unif. distributed random value)
                ang = unifrnd(0, 1);
                ampar.traj(:, i_ampar) = [NND + rad*cos(2*pi*ang), rad*sin(2*pi*ang), SIM.lowerLid]';
                ampar.tag(i_ampar)     = amparState_InCleft_OnPSD;
            end
            for i_ampar = 1 : 584
                rad = sqrt(unifrnd(34^2, 93^2)); % (Take square root of unif. distributed random value)
                ang = unifrnd(0, 1);
                ampar.traj(:, i_ampar + 0) = [NND + rad*cos(2*pi*ang), rad*sin(2*pi*ang), SIM.lowerLid]';
                ampar.tag(i_ampar + 0)     = amparState_InCleft_OnPSD;
            end
            for i_ampar = 1:340
                rad = sqrt(unifrnd(93^2, 273^2)); % (Take square root of unif. distributed random value)
                ang = unifrnd(0, 1);
                ampar.traj(:, i_ampar + 0) = [NND + rad*cos(2*pi*ang), rad*sin(2*pi*ang), SIM.lowerLid]';
                ampar.tag(i_ampar + 0)     = amparState_InCleft_OnPSD;  
            end

    end

elseif strcmp(SIM.group,'KO')
    if strcmp(SIM.subunit,'GluA1')
        %%%%%%% A1 KO
            for i_ampar = 1 : 312 %332
                rad = sqrt(unifrnd(0, 36^2)); % (Take square root of unif. distributed random value)
                ang = unifrnd(0, 1);
                ampar.traj(:, i_ampar) = [NND + rad*cos(2*pi*ang), rad*sin(2*pi*ang), SIM.lowerLid]';
                ampar.tag(i_ampar)     = amparState_InCleft_OnPSD;
            end
            for i_ampar = 1 : 548
                rad = sqrt(unifrnd(36^2, 83.5^2)); % (Take square root of unif. distributed random value)
                ang = unifrnd(0, 1);
                ampar.traj(:, i_ampar + 312) = [NND + rad*cos(2*pi*ang), rad*sin(2*pi*ang), SIM.lowerLid]';
                ampar.tag(i_ampar + 312)     = amparState_InCleft_OnPSD;
            end
            for i_ampar = 1:30 %10
                rad = sqrt(unifrnd(83.5^2, 196^2)); % (Take square root of unif. distributed random value)
                ang = unifrnd(0, 1);
                ampar.traj(:, i_ampar + 860) = [NND + rad*cos(2*pi*ang), rad*sin(2*pi*ang), SIM.lowerLid]';
                ampar.tag(i_ampar + 860)     = amparState_InCleft_OnPSD;  
            end   

    end       
end

    % save initial tags
    for i_ampar = 1 : SIM.amparNo
        ampar.tag_initial(i_ampar) = ampar.tag(i_ampar);
    end

    % initialize AMPARs: choose kinetic model and set initial state accordingly
    SIM.kineticsFrac_PSD     = SIM.kineticsFrac_PSD     / sum(SIM.kineticsFrac_PSD);     % normalize
    SIM.kineticsFrac_Outside = SIM.kineticsFrac_Outside / sum(SIM.kineticsFrac_Outside); % normalize
    for i_ampar = 1 : SIM.amparNo
        
        if ampar.tag(i_ampar) == amparState_InCleft_OnPSD
            kFrac = SIM.kineticsFrac_PSD;
        else
            kFrac = SIM.kineticsFrac_Outside;
        end
        
        
        % set kinetic model for every AMPAR
        u = unifrnd(0, 1);
        ll = 0;
        ampar.kinetics(i_ampar) = kinetics_GluA1; % default (numerical reasons, DONT CHANGE this line)
        for i_kinetics = 1 : KINETICS_NO
            lr = ll + kFrac(i_kinetics);
            if ll <= u  &&  u <= lr
                ampar.kinetics(i_ampar) = i_kinetics;
                break;
            end
            ll = lr;
        end
        
        states_F = find(kinetics.amparState_F(ampar.kinetics(i_ampar), :));
        % set initial states to `free' state(s)
        ampar.stateDistr(states_F, i_ampar) = 1/length(states_F);
        state = floor(unifrnd(0, length(states_F)))+1;
        if state > length(states_F)
            state = length(states_F);
        end
        ampar.currentState(i_ampar) = states_F(state); % choose a free state with uniform distribution
    end

    %---- transporter initialization -------------------------------------------
    fprintf('\b\b\b\b\bTRANS'); % INFO OUTPUT
    dt_transp = sqrt(2*SIM.D_transp*SIM.timeStepSize); %#ok %Brownian motion variance per timestep
    % array of all transporters,
    % transp is a structure with fields
    %   pos           position of transporter on glia cell
    %   state         state (bound, unbound, pumping)
    %   ligand        index of glu liganding
    %   T             if state=pumping: time of pumping until transpPumpingTime is reached
    %                 if state=bound:   time of pumping until transpPumpingTime is reached
    transp = struct(...
        'pos',         zeros(3, SIM.transpHor, SIM.transpVer), ...
        'state',       ones (   SIM.transpHor, SIM.transpVer, 'uint8') .* transpState_unbound, ...
        'ligand',      zeros(   SIM.transpHor, SIM.transpVer, 'uint32'), ...
        'pumpingTime', zeros(   SIM.transpHor, SIM.transpVer), ...
        'T',           zeros(   SIM.transpHor, SIM.transpVer));

    % uniform distribution (on a grid)
    for i_hor = 1 : SIM.transpHor
        for i_ver = 1 : SIM.transpVer
            %hei = transpOccupiedHeight*unifrnd(-1/2, 1/2) + cleftHeight/2;
            %ang = 2*Pi*unifrnd(0, 1);
            hei = (i_ver-1)*SIM.transp2transpDist + SIM.lowerLidES + SIM.transp2transpDist/2;
            ang =  i_hor   *2*pi/SIM.transpHor;
            transp.pos(:, i_hor, i_ver) = [SIM.R_ES*cos(ang), SIM.R_ES*sin(ang), hei];
            transp.state( i_hor, i_ver) = transpState_unbound;
        end
    end


    %---- vesicle initialization (continued) ------------------------------
        %---- vesicle release site initialization -------------------------
        for i_vesicle = 1 : SIM.vesicleNo
            rad = sqrt(unifrnd(0, SIM.R_ReleaseZone^2));
            ang = unifrnd(0, 1);
            vesicle.position(:, i_vesicle) = [SIM.releaseZoneCenter + rad*cos(2*pi*ang), rad*sin(2*pi*ang), SIM.upperLid]';
        end

        %---- vesicle load initialization ---------------------------------
        for i_vesicle = 1 : SIM.vesicleNo
            vesicle.gluNo_Released(i_vesicle) = 0;

            vesicle.gluNo(i_vesicle) = round( unifrnd(SIM.vesFillingFraction(1), SIM.vesFillingFraction(2)) * SIM.gluPerVesicle );

            % NOTE THAT this is a decimal number, not necessarily an integer
            if vesicle.duration(i_vesicle) == 0
                vesicle.timeStep_ReleaseRate(i_vesicle) = vesicle.gluNo(i_vesicle);
            else
                vesicle.timeStep_ReleaseRate(i_vesicle) = vesicle.gluNo(i_vesicle) / (vesicle.duration(i_vesicle) / SIM.timeStepSize);
            end
        end
    %---- glu initialization ---------------------------------------------------
    fprintf('\b\b\b\b\bGLU  '); % INFO OUTPUT

    SIM.gluNo = 0; %#ok     % released glu's

    gluNo_AllGlus = 0;  % number of glu's that is going to be released
    for i_vesicle = 1 : SIM.vesicleNo
        gluNo_AllGlus = gluNo_AllGlus + vesicle.gluNo(i_vesicle);
    end

    dt_glu = sqrt(2*SIM.D_glu*SIM.timeStepSize); %#ok %glu Brownian motion variance per timestep
    % array of all glu's,
    % glu is a structure with fields
    %   traj          trajectory
    %     MODIFICATION: save just current time step
    %   boundToAMPAR  index of AMPAR in ampar array (in case of binding)
    %   exitTime      time of exiting (in ms) (by absorption, pumping, ...)
    %!!!!4
    glu = struct(...
        'traj',            zeros(3, gluNo_AllGlus), ... % 3= dimension, 3d diffusioon
        'state',           ones (1, gluNo_AllGlus, 'uint8') *gluState_Unreleased, ...
        'boundToAMPAR',    ones (1, gluNo_AllGlus, 'uint32')*INVALID_INDEX, ...
        'bindingPosition', zeros(3, gluNo_AllGlus), ...
        'birthTime',       zeros(1, gluNo_AllGlus), ...
        'exitTime',        ones (1, gluNo_AllGlus) .* SIM.timeSteps);  %#ok
        
    %---- simulate -------------------------------------------------------------
    %   VVVVVVV    VVVVVVV    VVVVVVV    VVVVVVV    VVVVVVV    VVVVVVV   %
    %    VVVVV      VVVVV      VVVVV      VVVVV      VVVVV      VVVVV    %
    %     VVV        VVV        VVV        VVV        VVV        VVV     %
    %      V          V          V          V          V          V      %
    fprintf('\b\b\b\b\bRUN :'); % INFO OUTPUT
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SIM.timestep = 1;
    timestep_internal = SIM.timestep;   % this is for printLoopStateInfo.m
    timeSteps_internal = SIM.timeSteps; % this is for printLoopStateInfo.

    amparStatesOld = amparStates;       % save amparStates, diff after sim gives new sample

    % call external C function
    if     SIM.absorbAtCleftBD == 1
        synapse_absorb_at_cleftbd_C;
    elseif SIM.absorbAtGlia    == 1
        synapse_absorb_at_glia_C;
    else
        synapse_C;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    amparVar = amparVar + (amparStates-amparStatesOld) .^2;
    
    % save single run file
    %%%%% singleRunPrefix = 'SINGLE_RUN_';
    %%%%% amparStates_ThisRun = amparStates - amparStatesOld; %#ok
    %%%%% save(strcat(singleRunPrefix, num2str(i_runs), '.mat'), 'amparStates_ThisRun');
    %%%%% clear('amparStates_ThisRun');
    
    % save amparPeakDistr (if SIM.amparPeakRes > 0)
    if SIM.amparNo > 0
        if SIM.amparPeakRes > 0
            m = round(max(amparStates - amparStatesOld)*SIM.amparPeakRes);
            for i = 1 : AMPARSTATES_NO
                if m(i) < 1
                    m(i) = 1;
                end
                if m(i) > SIM.amparNo*SIM.amparPeakRes
                    m(i) = SIM.amparNo*SIM.amparPeakRes;
                end
                amparPeakDistr(m(i), i) = amparPeakDistr(m(i), i) + 1;
            end
        end
    end
    
    %---- save run just finished ------------------------------------------
    matlabFile = strcat(tempFilePrefix, num2str(i_runs), '.mat');
    save(matlabFile);

    fprintf('        %s\n', datestr(now)); % INFO OUTPUT
    
end % for i_runs

fprintf('-------- end\n'); % INFO OUTPUT



%----calculate AMPAR states --- mean and variance -------------------------
% So far, amparVar/amparVarS contains sum_1^SIM.runs x_n ^2
amparMean = amparMean/SIM.runs;
amparVar  = amparVar /SIM.runs - amparMean .^2;
amparVarS = amparVarS/SIM.runs - amparMean .^2;

%----calculate AMPAR tags --- variance
% So far, amparTagsVar contains sum_1^SIM.runs x_n ^2
for timestep = 1 : SIM.timeSteps
    amparTagsVar(timestep, :) = (amparTagsVar(timestep, :) - ...
                                 2* amparTags(timestep, :) .* amparTags(timestep, :)/SIM.runs + ...
                                 (amparTags(timestep, :).^2/SIM.runs ))/SIM.runs;
end

%---- save MATLAB file ----------------------------------------------------
% 
% fprintf('Saving %s ... ', matlabFileInfo);
% save(matlabFileInfo);

if SIM.SpillOver ==1
    char = [num2str(SIM.subunit),'_',num2str(SIM.period),'_',num2str(SIM.fileName),'_SpillOver.mat'];
else % SIM.SpillOver == 0
    char = [num2str(SIM.subunit),'_',num2str(SIM.period),'_',num2str(SIM.fileName),'.mat'];
end

fprintf('Saving %s ... \n ReleaseZone = %f', char, SIM.R_ReleaseZone);
fprintf('# on NC =  %f\n', SIM.amparNo_PSD);
save(char,'amparStates');


%---- save visualization files --------------------------------------------
if SIM.visualize ~= 0
    V_ampar_states(SIM);
    V_free_glus(SIM);
end
%---- delete temporary files (single-run files)----------------------------
for i_runs = 1 : SIM.runs
    matlabFile = strcat(tempFilePrefix, num2str(i_runs), '.mat');
    delete(matlabFile);
end

fprintf('Done\n');
