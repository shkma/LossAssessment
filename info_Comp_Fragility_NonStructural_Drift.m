function [PDS_ij_EDP, xm_Cost, numCompPerStory] = info_Comp_Fragility_NonStructural_Drift(i_n, i_m, x_PSDR_pdf)

% This function file returns from story number IDs (for each story):
% n = number of damage states a component may experience (this is per m)

% i_m = ID of component
       
    % % % % % % % % % % % % % % % %
    if i_m == 1
    % % % % % % % % % % % % % % % %
        
        numCompPerStory = 21*(9.144*9.144)  /  6.0; % total area of each floor divided by 6m^2 (Hwang and Lignos, 2017)
          
       if i_n == 0  % Drywall partition
          
          xm_EDP=0.0039; beta_EDP=0.17; xm_Cost=0.%   90.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.0039; beta_EDP=0.17; xm_Cost=90.;
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.0085; beta_EDP=0.23;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 2
          
          xm_EDP=0.0085; beta_EDP=0.23; xm_Cost=530.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end
       
       
    % % % % % % % % % % % % % % % %
    elseif i_m == 2
    % % % % % % % % % % % % % % % %
        
        numCompPerStory = 21*(9.144*9.144)  /  6.0; % total area of each floor divided by 6m^2 (Hwang and Lignos, 2017)
        
       if i_n == 0  % Drywall finish

          xm_EDP=0.0039; beta_EDP=0.17; xm_Cost=0.%   90.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.0039; beta_EDP=0.17; xm_Cost=90.;
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.0085; beta_EDP=0.23;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 2
          
          xm_EDP=0.0085; beta_EDP=0.23; xm_Cost=250.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end
       
       
    % % % % % % % % % % % % % % % %
    elseif i_m == 3
    % % % % % % % % % % % % % % % %
        
        numCompPerStory = (9.144*3.0*4.0 + 9.144*8.0) * 3.6576 * 6.0   /   2.81278   /6.0; % For each story.. See SMRF paper by Hwang
        
       if i_n == 0  % Exterior glazing

          xm_EDP=0.04; beta_EDP=0.36; xm_Cost=0.%   440.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.04; beta_EDP=0.36; xm_Cost=440.;
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.046; beta_EDP=0.33;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 2
          
          xm_EDP=0.046; beta_EDP=0.33; xm_Cost=440.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end

    end
    
end


