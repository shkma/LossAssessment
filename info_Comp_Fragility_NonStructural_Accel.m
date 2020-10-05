function [PDS_ij_EDP, xm_Cost, numCompPerStory] = info_Comp_Fragility_NonStructural_Accel(i_n, i_m, x_Accel_pdf)

% This function file returns from story number IDs (for each story):
% n = number of damage states a component may experience (this is per m)

% i_m = ID of component

    % % % % % % % % % % % % % % % %
    if i_m == 1
    % % % % % % % % % % % % % % % %
    
        numCompPerStory = 21*(9.144*9.144) /232.0; % total area divided by 232 m^2
        
       if i_n == 0  % Corrugated slab

          xm_EDP=0.35; beta_EDP=0.40; xm_Cost=0.%   3542.;
          F_DS_ij    = normcdf((log(x_Accel_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.35; beta_EDP=0.40; xm_Cost=3542.;
          F_DS_i1    = normcdf((log(x_Accel_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.55; beta_EDP=0.40;
          F_DS_i2    = normcdf((log(x_Accel_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 2
               
          xm_EDP=0.55; beta_EDP=0.40; xm_Cost=29337.;
          F_DS_i1    = normcdf((log(x_Accel_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.80; beta_EDP=0.40;
          F_DS_i2    = normcdf((log(x_Accel_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 3
          
          xm_EDP=0.80; beta_EDP=0.40; xm_Cost=55200.;
          F_DS_ij    = normcdf((log(x_Accel_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end
       
       
    % % % % % % % % % % % % % % % %
    elseif i_m == 2
    % % % % % % % % % % % % % % % %
        
        numCompPerStory = 21*(9.144*9.144) / (3.66*3.66); % total area of each floor divided by 3.66^2m^2 (each sprinkler protects 3.66*3.66m^2)
          
       if i_n == 0  % Automatic sprinklers
          
          xm_EDP=0.32; beta_EDP=1.40; xm_Cost=0.%   900.;
          F_DS_ij    = normcdf((log(x_Accel_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.32; beta_EDP=1.40; xm_Cost=900.;
          F_DS_ij    = normcdf((log(x_Accel_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end
       
       
    % % % % % % % % % % % % % % % %
    elseif i_m == 3
    % % % % % % % % % % % % % % % %
        
        numCompPerStory = 2.0; % total number of elevators (same as: Hwang and Lignos, 2017)
        
       if i_n == 0  % Drywall finish

          xm_EDP=0.50; beta_EDP=0.28; xm_Cost=0.%   868.;
          F_DS_ij    = normcdf((log(x_Accel_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.50; beta_EDP=0.28; xm_Cost=868.;
          F_DS_ij    = normcdf((log(x_Accel_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end
       
       
    end
    
end


