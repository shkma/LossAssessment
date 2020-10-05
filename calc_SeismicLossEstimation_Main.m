% Created by: Shoma Kitayama (shomakit@buffalo.edu)
% Last modified: 05.Oct.2020

% Seismic loss assessment of seismically isolated building
% The following exanple code implements seismic loss assessment procedure
% outlined in the manuscript for seismically isolated building with SCBF
% and with TFP-1 (notations are explained in the manuscript).
% The code below was prepared for analysis based on Conditional Spectra procedure but 
% can be modified to use for Incremental Dynamic Analysis procedure.

% NOTE: Before running the code below, please download the data of results of nonlinear response history
% analysis from the following website (GoogleDrive) and locate the folder in the same place as the matlab files are located:
% https://drive.google.com/file/d/1pNAoIicndCVjONx4P9fMqTAlJawJbsfy/view?usp=sharing

clear, clc

NumIM       =  10;    % Number of intensity measure (return period)
NumGM       =  40;    % Number of ground motions per each RP (Return Period)
NumStory    =  6;     % Number of stories of building
NumFloor    =  7;     % Number of floors of building including ground floor and roof
g           = 386.0;  % Gravity acceleration (in/sec^2)

lamda       = 0.03;   % Discount rate (=3%) per Hwang and Lignos (2017a, b)

n           = 50.;    % Lifetime duration of buildings.. (yrs.)

y_drift_Col = 0.05;   % Story drift ratio that causes collapse
Disp_limit  = 26.9;   % 32.6; % 38.3; % inch; for isolator displacement, D_Ultimate (collapse)

RP          = [43, 144, 289, 475, 949, 1485, 2475, 3899, 7462, 10000]; % Return period

% Specified for different systems --------------------------------------------------------------------------------
system        = 4;    % '1'=Nonisolated SCBF;                   '2'=Nonisolated SMRF;
                      % '3'=Isolated SCBF (RI=1) (Lower Bound); '4'=Isolated SCBF (RI=2) (Lower Bound);
                      % '5'=Isolated SMRF (RI=1) (Lower Bound); '6'=Isolated SMRF (RI=2) (Lower Bound);
                      % '7'=Isolated Upper Bound
                       
iso_size      = 1;    % '1'=TFP-1 or DC-1;
                      % '3'=TFP-3;
                       
wall_or_none  = 0;    % '0'=no moat wall;
                      % '1'=moat wall;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if     (system == 1) || (system == 2)
        ImpRows_Drift       = 2:7;  % Rows that contain the seismic response (max response is taken from these rows; change this per structural system & EDP)
        ImpRows_AccelAll    = 2:8;  % Rows that contain the seismic response (max response is taken from these rows; change this per structural system & EDP)
elseif (system == 3) || (system == 4) || (system == 5) || (system == 6) || (system == 7)
        ImpRows_Drift       = 3:8;  % Rows that contain the seismic response (max response is taken from these rows; change this per structural system & EDP)
        ImpRows_AccelAll    = 3:9;  % IMPORTANT: PGA is not considered for isolated structures,, elevator is sensitive to PFA at first floor (ie., isolated)
        ImpRows_TFPdisp     = 2:3;  % Rows that contain the seismic response (max response is taken from these rows; change this per structural system & EDP)
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Values used to calculate collapse probability (Baker JW,2015) ---------
num_gms                   = [NumGM,NumGM,NumGM,NumGM,NumGM,NumGM,NumGM,NumGM,NumGM,NumGM];
returnPeriod              = [43,144,289,475,949,1485,2475,3899,7462,10000];
num_collapse              = zeros( 1, NumIM );      % Vector for having number of collapse 
num_failure_StoryDrift    = zeros( 1, NumIM );      % Vector for having number of failure of story drift
num_failure_NumericalProb = zeros( 1, NumIM );      % Vector for having number of numerical failure
MaxDriftVector            = zeros( NumGM, NumIM );  % Vector that install maximum drift values to be checked how many collapse happens out of 40 GMs
MaxResDriftVector         = zeros( NumGM, NumIM );  % Vector that install maximum residual drift values to be checked how many collapse happens out of 40 GMs
MaxAccelFloorVector       = zeros( NumGM, NumIM );  % Vector that install maximum floor accel values to be checked how many collapse happens out of 40 GMs
MaxAccelRoofVector        = zeros( NumGM, NumIM );  % Vector that install maximum roof accel values to be checked how many collapse happens out of 40 GMs
MaxAccelAllVector         = zeros( NumGM, NumIM );  % Vector that install maximum roof accel values to be checked how many collapse happens out of 40 GMs

n_NumericalProb_matrix    = zeros( NumGM, NumIM );  % Matrix that install identification of if there is a numerical problem or not
n_Collapse_matrix         = zeros( NumGM, NumIM );  % Matrix that install identification of if there is a collapse or not
n_StoryDrift_matrix       = zeros( NumGM, NumIM );  % Matrix that install identification of if there is a excessive story drift or not

% Store Sa(T1) for different systems ------------------------------------
if system == 1      % Nonisolated SCBF
    Sa_T1     = [0.31,0.63,0.88,1.09,1.39,1.60,1.85,2.19,2.63,2.85];
    Period_T  = 0.524;
    Sa_MCE_T1 = 1.500;        % Spectral acceleration at T1 or TM
    beta_MDL  = 0.2;          % Modeling uncertainyy
    beta_TD   = 0.2;          % Test data uncertainyy
    beta_DR   = 0.2;          % Design requirement
    HazardCurve_Name = 'SeismicHazardData_0.524sec.txt';  % Specify name of file containing seismic hazard data
elseif system == 2  % Nonisolated SMRF
    Sa_T1     = [0.16,0.36,0.54,0.68,0.90,1.06,1.25,1.44,1.74,1.88];
    Period_T  = 1.186;
    Sa_MCE_T1 = 0.756302521;  % Spectral acceleration at T1 or TM
    beta_MDL  = 0.2;          % Modeling uncertainyy
    beta_TD   = 0.2;          % Test data uncertainyy
    beta_DR   = 0.1;          % Design requirement
    HazardCurve_Name = 'SeismicHazardData_1.186sec.txt';  % Specify name of file containing seismic hazard data
elseif (system == 3) || (system == 4) || (system == 5) || (system == 6)  % Isolated (lower bound)
    Sa_T1     = [0.02,0.06,0.10,0.16,0.24,0.31,0.39,0.44,0.53,0.58];
    Period_T  = 3.660;
    Sa_MCE_T1 = 0.246;  % Spectral acceleration at T1 or TM
    beta_MDL  = 0.2;    % Modeling uncertainyy
    beta_TD   = 0.2;    % Test data uncertainyy
    beta_DR   = 0.1;    % Design requirement
    HazardCurve_Name = 'SeismicHazardData_3.660sec.txt';  % Specify name of file containing seismic hazard data
elseif system == 7  % Isolated
    Sa_T1     = [0.04,0.10,0.17,0.24,0.35,0.42,0.50,0.57,0.69,0.75];
    Period_T  = 2.990;
    Sa_MCE_T1 = 0.301;  % Spectral acceleration at T1 or TM
    beta_MDL  = 0.2;    % Modeling uncertainyy
    beta_TD   = 0.2;    % Test data uncertainyy
    beta_DR   = 0.1;    % Design requirement
    HazardCurve_Name = 'SeismicHazardData_2.990sec.txt';  % Specify name of file containing seismic hazard data
end

% Construct collapse fragility curve --------------------------------------
% & store maximum residual drift along the height -------------------------
% & store maximum story drift at each story -------------------------------
% & store maximum PGA, PFA & PRF at each floor ----------------------------
fprintf(strcat('Compute collapse fragility curve - - - - - - - - - - - - - -  \n'));
for i_RP = 1 : NumIM
    dirname  = strcat('Response_analysis_data_CS_SCBF_RI2.0_Iso1_ASCE16/ReturnPeriodID_', num2str(i_RP), '_dynamicData');
    
    % Display the progress of processing data to generate collapse fragility curve
    fprintf('~ Progress (Seismic loss assessment): %s th RP out of %s RP is currently being processed \n', num2str(i_RP), num2str(NumIM));

    % Location of folder that contains the data for analysis (change this depending on location of analysis files)
    if     i_RP == 1  % (RP=   43 yrs)
       DataFolder = strcat('PEER_GroundMotionData/PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=43');
    elseif i_RP == 2  % (RP=  144 yrs)
       DataFolder = strcat('PEER_GroundMotionData/PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=144');
    elseif i_RP == 3  % (RP=  289 yrs)
       DataFolder = strcat('PEER_GroundMotionData/PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=289');
    elseif i_RP == 4  % (RP=  475 yrs)
       DataFolder = strcat('PEER_GroundMotionData/PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=475');
    elseif i_RP == 5  % (RP=  949 yrs)
       DataFolder = strcat('PEER_GroundMotionData/PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=949');
    elseif i_RP == 6  % (RP=  1485 yrs)
       DataFolder = strcat('PEER_GroundMotionData/PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=1485');
    elseif i_RP == 7  % (RP=  2475 yrs)
       DataFolder = strcat('PEER_GroundMotionData/PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=2475');
    elseif i_RP == 8  % (RP=  3899 yrs)
       DataFolder = strcat('PEER_GroundMotionData/PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=3899');
    elseif i_RP == 9  % (RP=  7462 yrs)
       DataFolder = strcat('PEER_GroundMotionData/PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=7462');
    elseif i_RP == 10 % (RP= 10000 yrs)
       DataFolder = strcat('PEER_GroundMotionData/PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=10000');
    end
        
        for i_gm = 1 : NumGM
            data_Drift      = load( strcat(char(dirname), '/GM_', num2str(i_gm), '_DriftStories.out') );  % Load story drift data (also for check if there's a collapse)
            data_Accel      = load( strcat(char(dirname), '/GM_', num2str(i_gm), '_AccelStories.out') );  % Load floor accel data
            data_TFPdisp    = load( strcat(char(dirname), '/GM_', num2str(i_gm), '_TFPdisp.out') );       % Load triple FP disp data.. Seismic frame
            endTime_Anal    = data_Drift(end,1);  % End (last) time of the time history response data (to be used to check if analysis completes or not; if not=collapse)
            
            % Collapse if defined as exceeding limit of story drift ratio (5% for SCBF) or
            % if there is numerical stability problem (analysis terminates)
                
                % % Load 'NumberGMpoint.txt' to be used in parametric study
                load( strcat(strcat(char(DataFolder), '/', 'NumberGMpoint.txt') ) ); % Load data of number of steps in acceleration history
                % % Load 'TimeStepGM.txt' to be used in parametric study
                load( strcat(strcat(char(DataFolder), '/', 'TimeStepGM.txt') ) );    % Load data of time steps in acceleration history
                % % Calculate max duration of ground motion
                endTime_GM = NumberGMpoint(i_gm) * TimeStepGM(i_gm); % Max duration of ground motion
                
                % Load horizontal displacmeent of isolator (to be checked if it exceeds D_Ultimate)
                MaxTFPdisp = max( max( abs(data_TFPdisp(:, ImpRows_TFPdisp)) ) );   % Obtain max horizontal TFP response (TFPdisp at SF)
                
                % Collapse is defined either by exceeding y_drift_Col OR D_Ultimate in isolator OR Uplift_limit of uplift
                temp_MaxDrift = max( max( abs(data_Drift(:, ImpRows_Drift)) ) );    % Obtain max response (PSDR)
                
                % % For Collapse - - - - -
                    if endTime_Anal < (endTime_GM - 0.0) || MaxTFPdisp > Disp_limit || temp_MaxDrift > y_drift_Col % || MaxTFPdisp > Disp_limit || MaxUpliftdisp > Uplift_limit
                       n_Collapse_matrix(i_gm,i_RP) = 1;                            % identification that there was a collapse
                    else
                       n_Collapse_matrix(i_gm,i_RP) = 0;                            % identification that there wasn't a collapse
                    end
                
                
                % % For max. residual story drift along the height of the
                % building given IM=im (=return period of i_RP)  - - - - -
                MaxResDriftVector(i_gm,i_RP)= max( abs(data_Drift(end, ImpRows_Drift)) ); % Obtain max response (RSDR)
                
                % % For peak seismic response at each story or floor - - - - -
                
                for i_story = 1 : length(ImpRows_Drift)
                    % % For story based peak response - - - - - (story drift)
                    if     (system == 1) || (system == 2) % i.e., non-isolated buildings
                        eval(['MaxDriftVector_Story_',num2str(i_story),'(i_gm,i_RP)', '= max( max( abs(data_Drift(:, i_story+1)) ) )']) % Obtain max PSDR at i_story
                    elseif (system == 3) || (system == 4) || (system == 5) || (system == 6) || (system == 7) % i.e., isolated building
                        eval(['MaxDriftVector_Story_',num2str(i_story),'(i_gm,i_RP)', '= max( max( abs(data_Drift(:, i_story+2)) ) )']) % Obtain max PSDR at i_story
                    end
                end
                
                for i_floor = 1 : length(ImpRows_AccelAll)
                    % % For floor (story) based peak response - - - - - (floor & roof acceleration)
                    if     (system == 1) || (system == 2) % i.e., non-isolated buildings
                        eval(['MaxAccelVector_Floor_',num2str(i_floor),'(i_gm,i_RP)', '= max( max( abs(data_Accel(:, i_floor+1)) ) )/g']) % Obtain max PGA, PFA, PRA at i_floor
                    elseif (system == 3) || (system == 4) || (system == 5) || (system == 6) || (system == 7) % i.e., isolated building
                        eval(['MaxAccelVector_Floor_',num2str(i_floor),'(i_gm,i_RP)', '= max( max( abs(data_Accel(:, i_floor+2)) ) )/g']) % Obtain max PGA, PFA, PRA at i_floor
                    end
                end
                
        end
        
    % Count number of GMs that causes collapse (out of 40) in current return period
    num_collapse(1, i_RP) = sum(n_Collapse_matrix(:,i_RP) > 0);  % Count and save in "num_collapse" (used later for maximum likelifood method)
            
end


% Estimate params of collapse fragility function using MLE method (equation 11 in Baker, 2015, Earthquake Spectra) --------------------------
[theta_hat_mle, beta_hat_mle] = fn_mle_pc(Sa_T1, num_gms, num_collapse);
% Compute collapse fragility curve (function) PC|IM using estimated parameters --------------------------
Sa_Range = 0.001:0.001:50;   % IM levels to plot fragility function
PC_Sa    = normcdf((log(Sa_Range/theta_hat_mle))/beta_hat_mle); % compute fragility function using Eq. 1 and estimated parameters
PC_Sa_write = [Sa_Range', PC_Sa'];  % Write out the output in text file
dlmwrite(strcat('PC_Sa.txt'), PC_Sa_write, 'delimiter', '\t', 'precision', 8); % Output the data of fragility curves

PC_IM = zeros(NumIM,2);

for i = 1:NumIM
    PC_IM(i,1)=i; 
    PC_IM(i,2)=normcdf((log(Sa_T1(i)/theta_hat_mle))/beta_hat_mle); 
    
        if num_collapse(1, i) == NumGM % Added 27June2019, to consider Probability of collapse of 100% (all GM caused collapse)
           PC_IM(i,2) = 1.;
        end
    
end

% Compute probability that the building is being considered to be demolished --------------------------
% First, compute probability density function of maximum residual story drift ratio, fRSDR|IM, along the height --------------------------
x_RSDR_pdf = [0.0002:0.0002:0.06]';
    % Assumed fragility curve for decision of demolition of building based
    % on Ramirez and Miranda (2012)...
    theta_hat_Demolish = 0.015; beta_hat_Demolish = 0.3; % From Ramirez and Miranda (2012)
    PD_RSDR = normcdf((log(x_RSDR_pdf/theta_hat_Demolish))/beta_hat_Demolish); % compute fragility function using Eq. 1 and estimated parameters
    vectorPD_RSDR = [x_RSDR_pdf, PD_RSDR];
    
    PD_IM_NC = zeros(NumIM,2);
    vector_mu_sigma_IM = zeros(NumIM,3);
    
for i_RP = 1 : NumIM
    
    MaxResDriftVector_NoCollapse      = MaxResDriftVector( n_Collapse_matrix(:,i_RP) < 1, i_RP ); % Remove any cases that collapse occurs
    MaxResDriftVector_NoCollapse_Log  = log(MaxResDriftVector_NoCollapse);
    mu_RSDR_IM                        = mean(MaxResDriftVector_NoCollapse);     % Mean of RSDR
    sigma_lnRSDR_IM                   = std(MaxResDriftVector_NoCollapse_Log);  % Standard deviation of lnPSDR
            
            if isscalar(MaxResDriftVector_NoCollapse) == 1
               sigma_lnRSDR_IM = 0.009999;
            end

            if isempty(MaxResDriftVector_NoCollapse) == 1
               mu_RSDR_IM      = y_drift_Col; % 9.999; Changed per correction in downtime study (19.Apr.2020)
               sigma_lnRSDR_IM = 0.009999; % 0.09999; Changed per correction in downtime study (19.Apr.2020)
            end
            
    vector_mu_sigma_IM(i_RP,1)        = i_RP;
    vector_mu_sigma_IM(i_RP,2)        = mu_RSDR_IM;
    vector_mu_sigma_IM(i_RP,3)        = sigma_lnRSDR_IM;
    
    vectorPDF_RSDR_IM_col1  = x_RSDR_pdf;                                        % column=1 of vector of PDF of RSDR at IM=im=i_RP
    vectorPDF_RSDR_IM_col2  = lognpdf(x_RSDR_pdf,log(mu_RSDR_IM),sigma_lnRSDR_IM);    % column=2 of vector of PDF of RSDR at IM=im=i_RP
    vectorPDF_RSDR_IM       = [vectorPDF_RSDR_IM_col1, vectorPDF_RSDR_IM_col2];
    eval(['vectorPDF_RSDR_IM_',num2str(i_RP), '= [vectorPDF_RSDR_IM_col1, vectorPDF_RSDR_IM_col2]']);
    
    PD_RSDR_x_fRSDR_IM_col1 = x_RSDR_pdf;
    PD_RSDR_x_fRSDR_IM_col2 = vectorPD_RSDR(:,2) .* vectorPDF_RSDR_IM(:,2);
    
    PD_IM_NC(i_RP,1) = i_RP;
    PD_IM_NC(i_RP,2) = trapz(PD_RSDR_x_fRSDR_IM_col1, PD_RSDR_x_fRSDR_IM_col2); % Probability that building is being considered to be demolished
    
end


% Compute mean values of loss for each source of loss ------------------------------------------------------------------------------------------
% Source of losses = "Structural repair loss";
%                    "Non-structural repair loss (drift)";
%                    "Non-structural repair loss (acc)"
%                    "Demolish loss";
%                    "Collapse loss"

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Structural repair loss"  - - - - - - - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %

x_PSDR_pdf = [0.00002:0.00002:3.000]'; % [0.0002:0.0002:3.000]'; Corrected on 19.Apr.2020

sigma_story_sigma_m_sigma_n_integral_Eq4=zeros(NumIM,2);

for i_RP = 1 : NumIM
    
    sigma_m_sigma_n_integral_Eq4 = zeros(NumStory,2); % Reset the vector for each i_RP
    
    for i_story = 1 : NumStory   % consider each story
        
        [m] = info_num_Components_Structural(i_story, system);  % from outside of main MATLAB file.. from i_story, obtain total # of damageable component IDs
        
        MaxPSDR_Vector_NoCollapse = eval(['MaxDriftVector_Story_',num2str(i_story),'( n_Collapse_matrix(:,i_RP) < 1, i_RP )']); % Added 25June2019
        mean_PSDR     = eval(['mean(MaxPSDR_Vector_NoCollapse','(:,1))']); % Added 25June2019
        std_log_PSDR  = eval(['std(log(MaxPSDR_Vector_NoCollapse','(:,1)))']);
            
            if isscalar(MaxPSDR_Vector_NoCollapse) == 1
               std_log_PSDR = 0.009999;
            end
            
            if isempty(MaxPSDR_Vector_NoCollapse) == 1
               mean_PSDR    = y_drift_Col;
               std_log_PSDR = 0.009999;
            end

        PDF_PSDR_IM   = lognpdf(x_PSDR_pdf,log(mean_PSDR),std_log_PSDR);
        
        sigma_n_integral_Eq4 = zeros(length(m),2);  % Reset the vector for each i_story
        
      for i_m = m  % consider number of types of damageable components at a story.. note: 1 to m corresponds to component ID
          
          [n] = info_num_DamageStates_Structural(i_m);  % then, from damageable each component ID, obtain number of damage states
          
          integral_Eq4=zeros(n+1,2);  % Reset the vector for each i_m
          
         for i_n = 0 : n  % consider number of damage states a component may experience
             
             [PDS_ij_EDP, xm_Cost, numCompPerStory] = info_Comp_Fragility_Structural(i_n, i_m, x_PSDR_pdf, system, i_story); % then, from damageable each component ID, obtain number of damage states
             
             integral_Eq4(i_n+1,1) = i_n+1;
             integral_Eq4(i_n+1,2) = trapz(x_PSDR_pdf, numCompPerStory*xm_Cost*(PDS_ij_EDP.*PDF_PSDR_IM));
            
         end
         
         sigma_n_integral_Eq4(i_m,1) = i_m;
         sigma_n_integral_Eq4(i_m,2) = sum(integral_Eq4(:,2));
         
      end
      
      sigma_m_sigma_n_integral_Eq4(i_story,2) = i_story;
      sigma_m_sigma_n_integral_Eq4(i_story,2) = sum(sigma_n_integral_Eq4(:,2));
      
    end
    
    sigma_story_sigma_m_sigma_n_integral_Eq4(i_RP,2) = i_RP;
    sigma_story_sigma_m_sigma_n_integral_Eq4(i_RP,2) = sum(sigma_m_sigma_n_integral_Eq4(:,2));
    
end

mu_LR_R_IM_NC_Structural = sigma_story_sigma_m_sigma_n_integral_Eq4; % mean of losses because of repairs for structural components for IM=im


% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Non-structural repair loss (Drift)"  - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %

x_PSDR_pdf = [0.00002:0.00002:3.000]'; % [0.0002:0.0002:3.000]'; Corrected on 19.Apr.2020

sigma_story_sigma_m_sigma_n_integral_Eq4=zeros(NumIM,2);

for i_RP = 1 : NumIM
    
    sigma_m_sigma_n_integral_Eq4 = zeros(NumStory,2); % Reset the vector for each i_RP
    
    for i_story = 1 : NumStory   % consider each story
        
        [m] = info_num_Components_NonStructural_Drift(i_story);  % from outside of main MATLAB file.. from i_story, obtain total # of damageable component IDs
        
        MaxPSDR_Vector_NoCollapse = eval(['MaxDriftVector_Story_',num2str(i_story),'( n_Collapse_matrix(:,i_RP) < 1, i_RP )']); % Added 25June2019
        mean_PSDR     = eval(['mean(MaxPSDR_Vector_NoCollapse','(:,1))']); % Added 25June2019
        std_log_PSDR  = eval(['std(log(MaxPSDR_Vector_NoCollapse','(:,1)))']);
            
            if isscalar(MaxPSDR_Vector_NoCollapse) == 1
               std_log_PSDR = 0.009999;
            end

            if isempty(MaxPSDR_Vector_NoCollapse) == 1
               mean_PSDR    = y_drift_Col; % 9.999; Changed per correction in downtime study (19.Apr.2020)
               std_log_PSDR = 0.009999;    % 0.09999; Changed per correction in downtime study (19.Apr.2020)
            end

        PDF_PSDR_IM   = lognpdf(x_PSDR_pdf,log(mean_PSDR),std_log_PSDR);
        
        sigma_n_integral_Eq4 = zeros(length(m),2); % Reset the vector for each i_story
        
      for i_m = m            % consider number of types(?) of damageable components at a story.. note: 1 to m corresponds to component ID
        
          [n] = info_num_DamageStates_NonStructural_Drift(i_m);  % then, from damageable each component ID, obtain number of damage states
              
          integral_Eq4=zeros(n+1,2); % Reset the vector for each i_m
          
         for i_n = 0 : n         % consider number of damage states a component may experience
             
             [PDS_ij_EDP, xm_Cost, numCompPerStory] = info_Comp_Fragility_NonStructural_Drift(i_n, i_m, x_PSDR_pdf); % then, from damageable each component ID, obtain number of damage states
             
             integral_Eq4(i_n+1,1) = i_n;
             integral_Eq4(i_n+1,2) = trapz(x_PSDR_pdf, numCompPerStory*xm_Cost*(PDS_ij_EDP.*PDF_PSDR_IM));
            
         end
         
         sigma_n_integral_Eq4(i_m,1) = i_m;
         sigma_n_integral_Eq4(i_m,2) = sum(integral_Eq4(:,2));
         
      end
      
      sigma_m_sigma_n_integral_Eq4(i_story,2) = i_story;
      sigma_m_sigma_n_integral_Eq4(i_story,2) = sum(sigma_n_integral_Eq4(:,2));
      
    end
    
    sigma_story_sigma_m_sigma_n_integral_Eq4(i_RP,2) = i_RP;
    sigma_story_sigma_m_sigma_n_integral_Eq4(i_RP,2) = sum(sigma_m_sigma_n_integral_Eq4(:,2));
    
end

mu_LR_R_IM_NC_NonStructural_Drift = sigma_story_sigma_m_sigma_n_integral_Eq4; % mean of losses because of repairs for structural components for IM=im


% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Non-structural repair loss (Accel)"  - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %

x_Accel_pdf = [0.0002:0.0002:5.000]';

sigma_story_sigma_m_sigma_n_integral_Eq4=zeros(NumIM,2);

for i_RP = 1 : NumIM
    
    sigma_m_sigma_n_integral_Eq4 = zeros(length(ImpRows_AccelAll),2); % Reset the vector for each i_RP
    
    for i_floor = 1 : NumFloor   % consider each floor
        
        [m] = info_num_Components_NonStructural_Accel(i_floor);  % from outside of main MATLAB file.. from i_floor, obtain total # of damageable component IDs
        
        MaxAccel_Vector_NoCollapse = eval(['MaxAccelVector_Floor_',num2str(i_story),'( n_Collapse_matrix(:,i_RP) < 1, i_RP )']); % Added 25June2019
        mean_Accel     = eval(['mean(MaxAccel_Vector_NoCollapse','(:,1))']); % Added 25June2019
        std_log_Accel  = eval(['std(log(MaxAccel_Vector_NoCollapse','(:,1)))']);

            if isscalar(MaxAccel_Vector_NoCollapse) == 1
               std_log_Accel = 0.009999;
            end

            if isempty(MaxAccel_Vector_NoCollapse) == 1
               mean_Accel    = x_Accel_pdf(end);   % 9.999; % Corrected: 22.Apr.2019
               std_log_Accel = 0.009999;           % 0.09999; Corrected: 22.Apr.2019
            end
            
        PDF_Accel_IM   = lognpdf(x_Accel_pdf,log(mean_Accel),std_log_Accel);
        
        sigma_n_integral_Eq4 = zeros(length(m),2); % Reset the vector for each i_floor
        
      for i_m = m   % consider number of types(?) of damageable components at a floor.. note: 1 to m corresponds to component ID
        
          [n] = info_num_DamageStates_NonStructural_Accel(i_m);   % then, from damageable each component ID, obtain number of damage states
              
          integral_Eq4=zeros(n+1,2);   % Reset the vector for each i_m
          
         for i_n = 0 : n   % consider number of damage states a component may experience
             
             [PDS_ij_EDP, xm_Cost, numCompPerStory] = info_Comp_Fragility_NonStructural_Accel(i_n, i_m, x_Accel_pdf); % then, from damageable each component ID, obtain number of damage states
             
             integral_Eq4(i_n+1,1) = i_n+1;
             integral_Eq4(i_n+1,2) = trapz(x_Accel_pdf, numCompPerStory*xm_Cost*(PDS_ij_EDP.*PDF_Accel_IM));
             
         end
         
         sigma_n_integral_Eq4(i_m,1) = i_m;
         sigma_n_integral_Eq4(i_m,2) = sum(integral_Eq4(:,2));
         
      end
      
      sigma_m_sigma_n_integral_Eq4(i_floor,2) = i_floor;
      sigma_m_sigma_n_integral_Eq4(i_floor,2) = sum(sigma_n_integral_Eq4(:,2));
      
    end
    
    sigma_story_sigma_m_sigma_n_integral_Eq4(i_RP,2) = i_RP;
    sigma_story_sigma_m_sigma_n_integral_Eq4(i_RP,2) = sum(sigma_m_sigma_n_integral_Eq4(:,2));
    
end

mu_LR_R_IM_NC_NonStructural_Accel = sigma_story_sigma_m_sigma_n_integral_Eq4; % mean of losses because of repairs for structural components for IM=im



% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Demolish loss"   - - - - - - - - - - - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% "Demolish loss" is assumed to be equal to the total replacement cost of
% the building (Hwang and Lignos, 2017a,b).

if     system == 1   % non-isolated SCBF
          mu_LD_D_IM_NC = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) ); % "$1880.00/m^2" is from (Hwang and Lignos, 2017 in ASCE - SCBF)
elseif system == 2   % non-isolated SMRF
          mu_LD_D_IM_NC = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) ); % "$2690.97/m^2" is from (Hwang and Lignos, 2017 in EESD - SMRF)
          
elseif system == 3   % isolated SCBF, RI=1
          mu_LD_D_IM_NC = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) ) + 1378766.00;
          
elseif (system == 4) && (iso_size == 1) && (wall_or_none == 0)   % isolated SCBF, RI=2
          mu_LD_D_IM_NC = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) ) + 1378766.00 - 521664.00;
elseif (system == 4) && (iso_size == 1) && (wall_or_none == 1)   % isolated SCBF, RI=2
          mu_LD_D_IM_NC = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) ) + 1550000.00 - 521664.00;
elseif (system == 4) && (iso_size == 3) && (wall_or_none == 0)   % isolated SCBF, RI=2
          mu_LD_D_IM_NC = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) ) + 1447704.00 - 521664.00;

elseif system == 5   % isolated SMRF, RI=1
          mu_LD_D_IM_NC = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) ) + 1378766.00 + 243728.00;
          
elseif (system == 6) && (iso_size == 1) && (wall_or_none == 0)   % isolated SMRF, RI=2
          mu_LD_D_IM_NC = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) ) + 1378766.00 - 000000.00;
elseif (system == 6) && (iso_size == 1) && (wall_or_none == 1)   % isolated SMRF, RI=2
          mu_LD_D_IM_NC = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) ) + 1550000.00 - 000000.00;
elseif (system == 6) && (iso_size == 3) && (wall_or_none == 0)   % isolated SMRF, RI=2
          mu_LD_D_IM_NC = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) ) + 1447704.00 - 000000.00;

end

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Collapse loss"   - - - - - - - - - - - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% "Collapse loss" is assumed to be equal to the total replacement cost of
% the building (Hwang and Lignos, 2017a,b).

if     system == 1   % non-isolated SCBF
          mu_Lc_C        = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) ); % "$1880.00/m^2" is from (Hwang and Lignos, 2017 in ASCE - SCBF)
          mu_Lc_C_nonIso = mu_Lc_C;
          
elseif system == 2   % non-isolated SMRF
          mu_Lc_C        = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) ); % "$2690.97/m^2" is from (Hwang and Lignos, 2017 in EESD - SMRF)
          mu_Lc_C_nonIso = mu_Lc_C;
          
elseif system == 3   % isolated SCBF, RI=1
          mu_Lc_C        = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) ) + 1378766.00;
          mu_Lc_C_nonIso = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) );
          
elseif (system == 4) && (iso_size == 1) && (wall_or_none == 0)   % isolated SCBF, RI=2
          mu_Lc_C        = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) ) + 1378766.00 - 521664.00;
          mu_Lc_C_nonIso = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) );
          
elseif (system == 4) && (iso_size == 1) && (wall_or_none == 1)   % isolated SCBF, RI=2
          mu_Lc_C        = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) ) + 1550000.00 - 521664.00;
          mu_Lc_C_nonIso = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) );
          
elseif (system == 4) && (iso_size == 3) && (wall_or_none == 0)   % isolated SCBF, RI=2
          mu_Lc_C        = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) ) + 1447704.00 - 521664.00;
          mu_Lc_C_nonIso = (NumFloor-1) * 1880.00 * ( 21*(9.144*9.144) );

elseif system == 5   % isolated SMRF, RI=1
          mu_Lc_C        = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) ) + 1378766.00 + 243728.00;
          mu_Lc_C_nonIso = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) );
          
elseif (system == 6) && (iso_size == 1) && (wall_or_none == 0)   % isolated SMRF, RI=2
          mu_Lc_C        = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) ) + 1378766.00 - 000000.00;
          mu_Lc_C_nonIso = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) );
          
elseif (system == 6) && (iso_size == 1) && (wall_or_none == 1)   % isolated SMRF, RI=2
          mu_Lc_C        = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) ) + 1550000.00 - 000000.00;
          mu_Lc_C_nonIso = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) );
          
elseif (system == 6) && (iso_size == 3) && (wall_or_none == 0)   % isolated SMRF, RI=2
          mu_Lc_C        = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) ) + 1447704.00 - 000000.00;
          mu_Lc_C_nonIso = (NumFloor-1) * 2690.97 * ( 21*(9.144*9.144) );

end

% ------------------------------------------------------------------------------------------
% MULTIPLY PROBABILITY OF DEMOLISH & COLLAPSE ----------------------------------------------
% ------------------------------------------------------------------------------------------

% Source of losses = "Structural repair loss";
%                    "Non-structural repair loss (drift)";
%                    "Non-structural repair loss (acc)"
%                    "Demolish loss";
%                    "Collapse loss"

Loss_Structural_RP           =  zeros(NumIM+1,3);
Loss_NonStructural_Drift_RP  =  zeros(NumIM+1,3);
Loss_NonStructural_Accel_RP  =  zeros(NumIM+1,3);
Loss_Demolish_RP             =  zeros(NumIM+1,3);
Loss_Collapse_RP             =  zeros(NumIM+1,3);
Loss_Total_RP                =  zeros(NumIM+1,3);

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Structural repair loss"  - - - - - - - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
Loss_Structural_RP(2:end,2)          = mu_LR_R_IM_NC_Structural(:,2)          .* (1.-PD_IM_NC(:,2)) .* (1.-PC_IM(:,2));

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Non-structural repair loss (drift)"  - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
Loss_NonStructural_Drift_RP(2:end,2) = mu_LR_R_IM_NC_NonStructural_Drift(:,2) .* (1.-PD_IM_NC(:,2)) .* (1.-PC_IM(:,2));

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Non-structural repair loss (acc)"  - - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
Loss_NonStructural_Accel_RP(2:end,2) = mu_LR_R_IM_NC_NonStructural_Accel(:,2) .* (1.-PD_IM_NC(:,2)) .* (1.-PC_IM(:,2));

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Demolish loss"   - - - - - - - - - - - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
Loss_Demolish_RP(2:end,2) = mu_LD_D_IM_NC * PD_IM_NC(:,2) .* (1.-PC_IM(:,2));

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Collapse loss"   - - - - - - - - - - - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
Loss_Collapse_RP(2:end,2) = mu_LD_D_IM_NC * PC_IM(:,2);

for i = 1:NumIM
   
    Loss_Structural_RP(i+1,1)           =  Sa_T1(i);  Loss_Structural_RP(i+1,3)          = RP(i);
    Loss_NonStructural_Drift_RP(i+1,1)  =  Sa_T1(i);  Loss_NonStructural_Drift_RP(i+1,3) = RP(i);
    Loss_NonStructural_Accel_RP(i+1,1)  =  Sa_T1(i);  Loss_NonStructural_Accel_RP(i+1,3) = RP(i);
    Loss_Demolish_RP(i+1,1)             =  Sa_T1(i);  Loss_Demolish_RP(i+1,3)            = RP(i);
    Loss_Collapse_RP(i+1,1)             =  Sa_T1(i);  Loss_Collapse_RP(i+1,3)            = RP(i);
    Loss_Total_RP(i+1,1)                =  Sa_T1(i);  Loss_Total_RP(i+1,3)               = RP(i);
    
end

% ----------------------------------------------------------------------------
% COMPUTE TOTAL LOSS FOR IM=im  ----------------------------------------------
% ----------------------------------------------------------------------------
% Source of losses = "Structural repair loss";
Loss_Total_RP(2:end,2) = Loss_Structural_RP(2:end,2) + Loss_NonStructural_Drift_RP(2:end,2) + Loss_NonStructural_Accel_RP(2:end,2) + Loss_Demolish_RP(2:end,2) + Loss_Collapse_RP(2:end,2);

% ----------------------------------------------------------------------------
% OUTPUT VULNERABILITY FUNCTIONS  --------------------------------------------
% ----------------------------------------------------------------------------

% As-is - - - - -
dlmwrite(strcat('Loss_Structural_RP.txt'),               Loss_Structural_RP,                  'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Loss_NonStructural_Drift_RP.txt'),      Loss_NonStructural_Drift_RP,         'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Loss_NonStructural_Accel_RP.txt'),      Loss_NonStructural_Accel_RP,         'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Loss_Demolish_RP.txt'),                 Loss_Demolish_RP,                    'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Loss_Collapse_RP.txt'),                 Loss_Collapse_RP,                    'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Loss_Total_RP.txt'),                    Loss_Total_RP,                       'delimiter', '\t', 'precision', 8);

% Normalized by total repair cost - - - - -
dlmwrite(strcat('Norm_Loss_Structural_RP.txt'),          [Loss_Structural_RP(:,1),          Loss_Structural_RP(:,2)/mu_Lc_C,           Loss_Structural_RP(:,3)],          'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Norm_Loss_NonStructural_Drift_RP.txt'), [Loss_NonStructural_Drift_RP(:,1), Loss_NonStructural_Drift_RP(:,2)/mu_Lc_C,  Loss_NonStructural_Drift_RP(:,3)], 'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Norm_Loss_NonStructural_Accel_RP.txt'), [Loss_NonStructural_Accel_RP(:,1), Loss_NonStructural_Accel_RP(:,2)/mu_Lc_C,  Loss_NonStructural_Accel_RP(:,3)], 'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Norm_Loss_Demolish_RP.txt'),            [Loss_Demolish_RP(:,1),            Loss_Demolish_RP(:,2)/mu_Lc_C,             Loss_Demolish_RP(:,3)],            'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Norm_Loss_Collapse_RP.txt'),            [Loss_Collapse_RP(:,1),            Loss_Collapse_RP(:,2)/mu_Lc_C,             Loss_Collapse_RP(:,3)],            'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Norm_Loss_Total_RP.txt'),               [Loss_Total_RP(:,1),               Loss_Total_RP(:,2)/mu_Lc_C,                Loss_Total_RP(:,3)],               'delimiter', '\t', 'precision', 8);


% ----------------------------------------------------------------------------
% COMPUTE EXPECTED ANNUAL LOSS (EAL) -----------------------------------------
% ----------------------------------------------------------------------------
fprintf(strcat('Compute Expected Annual Loss (EAL) - - - - -  \n'));

% Following code includes the integration of seismic hazard curve and loss
% vulnerability curve. The integration scheme may be similar to the computation
% of mean annual frequency of exceedance of seismic demand parameters in
% Table 6-5 & 6-6 in the following document:
% Kumar M, Whittaker AS, Constantinou MC. (2015). "Seismic Isolation of Nuclear Power Plants
% using Sliding Bearings" Technical Report MCEER-15-0006. December 27,
% 2015.

HazardCurveData  =  load(strcat('Seismic_Hazard_Data/',char(HazardCurve_Name))); % Load selected seismic hazard data (differs per each period - system)

Sa_1=zeros(length(Sa_T1)+1,1); Sa_2=zeros(length(Sa_T1)+1,1);

Sa_1(2:end,1)    =  Sa_T1(:);
Sa_2(1:end-1,1)  =  Sa_T1(:);
Sa_2(end,1)      =  99.9;

Sa_Ave           = (Sa_1+Sa_2)/2.0;

for choice = 1:6

switch choice
    case 1  % Loss_Structural_RP
        Loss_Vulnerability  =  Loss_Structural_RP(:,2);
        type_of_EAL         =  'Loss_Structural';
    case 2  % Loss_NonStructural_Drift_RP
        Loss_Vulnerability  =  Loss_NonStructural_Drift_RP(:,2);
        type_of_EAL         =  'Loss_NonStructural_Drift';
    case 3  % Loss_NonStructural_Accel_RP
        Loss_Vulnerability  =  Loss_NonStructural_Accel_RP(:,2);
        type_of_EAL         =  'Loss_NonStructural_Accel';
    case 4  % Loss_Demolish_RP
        Loss_Vulnerability  =  Loss_Demolish_RP(:,2);
        type_of_EAL         =  'Loss_Demolish';
    case 5  % Loss_Collapse_RP
        Loss_Vulnerability  =  Loss_Collapse_RP(:,2);
        type_of_EAL         =  'Loss_Collapse';
    case 6  % Loss_Total_RP
        Loss_Vulnerability  =  Loss_Total_RP(:,2);
        type_of_EAL         =  'Loss_Total';
end
    
    MeanLoss_Ave = mean([Loss_Vulnerability(1:end-1)';Loss_Vulnerability(2:end)'])'; % cdf('Lognormal', Sa_Ave(:,1), log(medianSa), beta);  
    
    lamb_Sa_1=zeros(length(Sa_T1),1); lamb_Sa_2=zeros(length(Sa_T1),1);
    
    for s = 1 : NumIM
        [c, index]        = min(abs(HazardCurveData(:,1)-Sa_T1(s)));
        lamb_Sa_1(s, 1)   = HazardCurveData(index,2);
    end
        lamb_Sa_2(1:end-1,1) = lamb_Sa_1(2:end,1);
        lamb_Sa_2(  end,  1) = 1.e-9;   % arbitrary very small frequency of occurance

    Delta_lamb_Sa         =  abs(lamb_Sa_2 - lamb_Sa_1);
    lamb_F_i              =  Delta_lamb_Sa .* MeanLoss_Ave;  
    lamb_F                =  sum(lamb_F_i);           % Expected Annual Loss (EAL)
    EAL                   =  lamb_F;                  % Expected Annual Loss (EAL)
    Normalized_EAL        =  lamb_F/mu_Lc_C;          % Expected Annual Loss (EAL)
    
    dlmwrite(strcat('EAL_',            char(type_of_EAL), '.txt'),            EAL, 'delimiter', '\t', 'precision', 8); % for PSDR
    eval(['EAL_',                      char(type_of_EAL), '= EAL']);                  % Obtain EAL
    dlmwrite(strcat('Normalized_EAL_', char(type_of_EAL), '.txt'), Normalized_EAL, 'delimiter', '\t', 'precision', 8); % for PSDR
    eval(['Normalized_EAL_',           char(type_of_EAL), '= Normalized_EAL']);       % Obtain Normalized_EAL

end


% ----------------------------------------------------------------------------
% COMPUTE Expected Loss (=EL) Over Time (=t) ---------------------------------
% ----------------------------------------------------------------------------

t  = 0.2:0.2:100;
EL = zeros(length(t), 2);

Cp = mu_Lc_C - mu_Lc_C_nonIso; % Cost premium for isolated buildings relative to corresponding non-isolated buildings.

EL(:,1) = t';
EL(:,2) = (1.0 - exp(-lamda*t'))/lamda * EAL + Cp;

dlmwrite(strcat('EL.txt'), EL, 'delimiter', '\t', 'precision', 8);


