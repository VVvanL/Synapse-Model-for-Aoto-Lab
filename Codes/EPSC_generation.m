% clear;
hold off;
% %%% Parameters can be changed

runs = 160;
% SSD_Num = 10; % eEPSC for WT, KO_2 and KO_4
SSD_Num = 8.2; % eEPSC for KO with reduced SSDs and Nrxn3 KO

period = 12 % P4, P8, P12, P18, P30
ratio = 1; % GluA1 / (GluA1 + GluA4)
              % P4 : 0.85;  P8 : 0.7; P12: 0.45
              % P16: 0.42;  P30: 0.25
SpillOver = 0;
NND_No = 0; % number of synapse @NND

%%%% Calculation of mEPSC %%%%
timestepSize = 10 * 2000;% 10 ms * 2000 timeStepSize per microsecond
baseline = 2 * 2000; % 2 ms
dt = 0.5; % 0.5us, time step
range = 1000; % the range of fitting, 
              % it should be changed sometimes, because the fitting
              % equation can't be got when the curve larger than zero

%%%% AMPAR attribute
Pr = 1; % release probability
Vampa = 0; % reversal potential

%%%%% The membrane is holding in -70 mV.
V_a1 = -70; %mV, initial volt
V_a4 = -70; %mV, initial volt

I_total_a1 = zeros(20000,runs);% Store the traces of GluA1 in each loop
I_total_a4 = zeros(20000,runs);% Store the traces of GluA4 in each loop
decayTime_var = zeros(1, runs);
riseTime_var  = zeros(1, runs);
amp_var = zeros(1, runs);
for i = 1: runs
    run = i
    %%%%  slow GluAs  %%%%
    if ratio == 0
        I_total_a1 = zeros(20000,runs);
        continue
    else
        name = ['GluA1_P',num2str(period),'_',num2str(i) '.mat'];
        load(name);
        
%         open_state_total_a1 = amparStates(:,4);
        open_state_total_a1 = amparStates(:,4)/10 *SSD_Num;
        open_state_a1 = open_state_total_a1;
        
       
        g_a1 = 0.031;% nanoSieman, Robert A. and Howe J., 2003        

        if SpillOver == 1
            name_SpillOver = ['GluA1_P',num2str(period),'_',num2str(i) '_SpillOver.mat'];
            load(name_SpillOver);
            
            open_state_a1 = open_state_a1 + NND_No * amparStates(:,4)/10; % AMPAR open number     
        end
                  
        I_2_a1 = zeros(timestepSize ,1);
        for k = 1:timestepSize
            I_2_a1(k) = (g_a1*open_state_a1(k)) * (V_a1-Vampa);%Ampa excitatory current in picoAmpere.
        end            
    
        I_total_a1(:,i) = I_2_a1;
    end
 

    %%%%  fast GluAs  %%%%
    if ratio ==1
        I_total_a4 = zeros(20000,runs);
        continue
    else
        name = ['GluA4_P',num2str(period),'_',num2str(i), '.mat'];
        load(name);        
        
%         open_state_total_a4 = amparStates(:,4);  
        open_state_total_a4 = amparStates(:,4)/10;       
        open_state_a4 = open_state_total_a4; % AMPAR open number
        
        g_a4 = 0.045;% nanoSieman, Robert A. and Howe J., 2003  
        
        if SpillOver == 1
            name_SpillOver = ['GluA4_P',num2str(period),'_',num2str(i) '_SpillOver.mat'];
            load(name_SpillOver);          
            
            open_state_a4 = open_state_a4 + NND_No * amparStates(:,4); % open probability on PSD        
        end
        
        I_2_a4 = zeros(timestepSize ,1);
        for k = 1:timestepSize
            I_2_a4(k) = (g_a4 * open_state_a4(k)) * (V_a4-Vampa);%Ampa excitatory current in picoAmpere.
        end
        I_total_a4(:,i) = I_2_a4;
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Calculate the 10-90% rise time and decay time %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = timestepSize - baseline; 
        X = (1:tmp)'; 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GluA1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_2_a1 = mean(I_total_a1,2);
%         I_2_a1 = sum(I_total_a1,2); % eEPSCs
%         I_2_a1 = I_total_a1*1; % eEPSCs
        I_fit1 = I_2_a1;
        I_fit1(1:baseline) = []; % Remove the '0' part

        I_max_a1 = min(I_fit1); % Actually it's the peak value of current              
        index_max_a1 = find(I_fit1 == I_max_a1); 

        % extract rise phase
        X_rise_a1 = X; X_rise_a1((index_max_a1+1):tmp) = [];
        I_rise_a1 = I_fit1; I_rise_a1((index_max_a1+1):tmp) = [];
        % extract decay phase
        X_decay_a1 = X; X_decay_a1(1:(index_max_a1-1)) = [];
        I_decay_a1 = I_fit1; I_decay_a1(1:(index_max_a1-1)) = [];


        f_decay_a1 = fit(X_decay_a1(1:range),I_decay_a1(1:range),'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
        I_37_a1 = 0.37 * I_max_a1; % 37% of the current peak
        syms x ;
        fun1_a1 = f_decay_a1.a * exp(f_decay_a1.b * x) + f_decay_a1.c * exp(f_decay_a1.d * x) == I_37_a1;  % V = V0 * exp(-t/tau) = a * exp(b * x)
        decayTime_a1 = double(vpasolve(fun1_a1,x,1000) - index_max_a1)/2000;

        f_rise_a1 = fit(X_rise_a1(:),I_rise_a1,'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
        I_10_a1 = 0.1 * I_max_a1; % 10% of the current peak
        I_90_a1 = 0.9 * I_max_a1; % 90% of the current peak
        syms x ;
        fun10_a1 = f_rise_a1.a * exp(f_rise_a1.b * x) + f_rise_a1.c * exp(f_rise_a1.d * x) == I_10_a1;  % V = V0 * exp(-t/tau) = a * exp(b * x)
        fun90_a1 = f_rise_a1.a * exp(f_rise_a1.b * x) + f_rise_a1.c * exp(f_rise_a1.d * x) == I_90_a1;
        riseTime10_a1 = double(vpasolve(fun10_a1,x));
        riseTime90_a1 = double(vpasolve(fun90_a1,x));

        riseTime_a1 = double(riseTime90_a1-riseTime10_a1)/2000;


        
%         %%%%%%%%%%%%%%%%%%%% GluA4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         I_2_a4 = mean(I_total_a4,2);
% %         I_2_a4 = sum(I_total_a4,2);
% %         I_2_a4 = I_total_a4 * 1;
%         I_fit_a4 = I_2_a4;
%         I_fit_a4(1:baseline) = []; % Remove the '0' part    
% 
%         I_max_a4 = min(I_fit_a4);% Actually it's the peak value of current
%         index_max_a4 = find(I_fit_a4 == I_max_a4);
% 
% %        extract rise phase
%         X_rise_a4 = X; X_rise_a4((index_max_a4 + 1) : tmp) = [];
%         I_rise_a4 = I_fit_a4; I_rise_a4 ((index_max_a4 + 1) : tmp) = [];
% %         extract decay phase
%         X_decay_a4 = X; X_decay_a4(1:(index_max_a4-1)) = [];
%         I_decay_a4 = I_fit_a4; I_decay_a4(1:(index_max_a4-1)) = [];                
% 
% 
%         f_decay_a4 = fit(X_decay_a4(1:range),I_decay_a4(1:range),'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
%         I_37_a4 = 0.37 * I_max_a4; % 37% of the current peak
%         syms x ;
%         fun1_a4 = f_decay_a4.a * exp(f_decay_a4.b * x) + f_decay_a4.c * exp(f_decay_a4.d * x) == I_37_a4;  % V = V0 * exp(-t/tau) = a * exp(b * x)
%         decayTime_a4 = double(vpasolve(fun1_a4,x,1000) - index_max_a4)/2000;
% 
%         f_rise_a4 = fit(X_rise_a4(:),I_rise_a4,'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
%         I_10_a4 = 0.1 * I_max_a4; % 10% of the current peak
%         I_90_a4 = 0.9 * I_max_a4; % 90% of the current peak
%         syms x ;
%         fun10_a4 = f_rise_a4.a * exp(f_rise_a4.b * x) + f_rise_a4.c * exp(f_rise_a4.d * x) == I_10_a4;  % V = V0 * exp(-t/tau) = a * exp(b * x)
%         fun90_a4 = f_rise_a4.a * exp(f_rise_a4.b * x) + f_rise_a4.c * exp(f_rise_a4.d * x) == I_90_a4;
%         riseTime10_a4 = double(vpasolve(fun10_a4,x));
%         riseTime90_a4 = double(vpasolve(fun90_a4,x));
% 
%         riseTime_a4 = double(riseTime90_a4 - riseTime10_a4)/2000; 


        %%%%%%%%%%%%%%%%%%%%%%%%%%% Combine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = timestepSize - baseline; 
        X = (1:tmp)'; 
        
        if ratio == 0
            I_2_a1 = zeros(20000,runs);
        elseif ratio ==1
            I_2_a4 = zeros(20000,runs);
        end
        
        I_total = I_2_a1 + I_2_a4;
        I_fit = I_total;
        I_fit(1:baseline) = []; % Remove the '0' part    

        I_max = min(I_fit);% Actually it's the peak value of current
        index_max = find(I_fit == I_max);

        % extract rise phase
        X_rise = X; X_rise((index_max + 1) : tmp) = [];
        I_rise = I_fit; I_rise((index_max + 1) : tmp) = [];
       % extract decay phase
        X_decay = X; X_decay(1:(index_max-1)) = [];
        I_decay = I_fit; I_decay(1:(index_max-1)) = [];                

        
%         f_decay = fit(X_decay(1:range),I_decay(1:range),'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
%         I_37 = 0.37 * I_max; % 37% of the current peak
%         syms x ;
%         fun = f_decay.a * exp(f_decay.b * x) + f_decay.c * exp(f_decay.d * x) == I_37;  % V = V0 * exp(-t/tau) = a * exp(b * x)
%         decayTime = double(vpasolve(fun,x,2000) - index_max)/2000;
% 
%         f_rise = fit(X_rise(:),I_rise,'exp2','Robust','LAR','Algorithm','Levenberg-Marquardt'); % double exponential fit, for decay time                                                   
%         I_10 = 0.1 * I_max; % 10% of the current peak
%         I_90 = 0.9 * I_max; % 90% of the current peak
%         syms x ;
%         fun10 = f_rise.a * exp(f_rise.b * x) + f_rise.c * exp(f_rise.d * x) == I_10;  % V = V0 * exp(-t/tau) = a * exp(b * x)
%         fun90 = f_rise.a * exp(f_rise.b * x) + f_rise.c * exp(f_rise.d * x) == I_90;
%         riseTime10 = double(vpasolve(fun10,x));
%         riseTime90 = double(vpasolve(fun90,x));
% 
%         riseTime = double(riseTime90 - riseTime10)/2000; 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        plot(I_total_a1,'Color',[0.65,0.65,0.65]); hold on;
        plot(I_2_a1,'b');hold on;  
        plot(I_total_a4,'Color',[0.65,0.65,0.65]); hold on;
        plot(I_2_a4,'r');hold on; 
        plot(I_total,'k'); hold off; 
        
%         decay_Text_a1 = ['GluA1 decay time = ', num2str(decayTime_a1)];
%         text(12000,-1,decay_Text_a1);        
%         rise_Text_a1 = ['GluA1 rise time = ', num2str(riseTime_a1)];
%         text(12000,-2,rise_Text_a1);
%         amplitude_a1 = ['GluA1 amplitude = ', num2str(I_max_a1)];
%         text(12000,-3,amplitude_a1); 
% 
%         decay_Text_a4 = ['GluA4 decay time = ', num2str(decayTime_a4)];
%         text(12000,-5,decay_Text_a4);        
%         rise_Text_a4 = ['GluA4 rise time = ', num2str(riseTime_a4)];
%         text(12000,-6,rise_Text_a4);
%         amplitude_a4 = ['GluA4 amplitude = ', num2str(I_max_a4)];
%         text(12000,-7,amplitude_a4);         
                
%         decay_Text = ['GluA1&GluA4 decay time = ', num2str(decayTime)];
%         text(12000,-1,decay_Text);
%         rise_Text = ['GluA1&GluA4 rise time = ', num2str(riseTime)];
%         text(12000,-3,rise_Text);
%         amplitude = ['Combined amplitude = ', num2str(I_max)];
%         text(12000,-5,amplitude); 

        decay_Text = ['GluA1 decay time = ', num2str(decayTime_a1)];
        text(12000,-1,decay_Text);
        rise_Text = ['GluA1 rise time = ', num2str(riseTime_a1)];
        text(12000,-3,rise_Text);
        amplitude = ['GluA1 amplitude = ', num2str(I_max_a1)];
        text(12000,-5,amplitude);
        xlim([0 timestepSize]);
        set(gca,'XTick',[0 4000 8000 12000 16000 20000],'XTickLabel',[0 2 4 6 8 10]);

%         title('Combined ')
        xlabel('time (ms)');ylabel('EPSC (pA)');
        box off;
        
beep;