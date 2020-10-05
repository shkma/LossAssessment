function [PDS_ij_EDP, xm_Cost, numCompPerStory] = info_Comp_Fragility_Structural(i_n, i_m, x_PSDR_pdf, system, i_story)

% This function file returns from story number IDs (for each story):
% n = number of damage states a component may experience (this is per m)

% i_m = ID of component
    
    % % % % % % % % % % % % % % % %
    if i_m == 1
    % % % % % % % % % % % % % % % %
        
        numCompPerStory = 16; % below each SCBF
        
       if i_n == 0
          
          xm_EDP=0.04; beta_EDP=0.40; xm_Cost=0.; %   19224.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.04; beta_EDP=0.40;
          
              if     system == 1 || system == 3 || system == 4 % Non-isolated SCBF or Isolated SCBF (RI=1 or 2)
                     xm_Cost=19224.;
              elseif system == 2 || system == 6                % Non-isolated SMRF or Isolated SMRF (RI=2)
                     xm_Cost=20082.;
              elseif system == 5                               % Isolated SMRF (RI=1)
                     xm_Cost=21363.;
              end
                 
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.07; beta_EDP=0.40;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 2
               
          xm_EDP=0.07; beta_EDP=0.40;
          
              if     system == 1 || system == 3 || system == 4 % Non-isolated SCBF or Isolated SCBF (RI=1 or 2)
                     xm_Cost=27263.;
              elseif system == 2 || system == 6                % Non-isolated SMRF or Isolated SMRF (RI=2)
                     xm_Cost=29395.;
              elseif system == 5                               % Isolated SMRF (RI=1)
                     xm_Cost=32567.;
              end
          
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.10; beta_EDP=0.40; 
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 3
          
          xm_EDP=0.10; beta_EDP=0.40;
          
              if     system == 1 || system == 3 || system == 4 % Non-isolated SCBF or Isolated SCBF (RI=1 or 2)
                     xm_Cost=32423.;
              elseif system == 2 || system == 6                % Non-isolated SMRF or Isolated SMRF (RI=2)
                     xm_Cost=36657.;
              elseif system == 5                               % Isolated SMRF (RI=1)
                     xm_Cost=41890.;
              end
          
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end
       
    % % % % % % % % % % % % % % % %
    elseif i_m == 2
    % % % % % % % % % % % % % % % %
        
        numCompPerStory = 32;
        
       if i_n == 0

          xm_EDP=0.04; beta_EDP=0.40; xm_Cost=0.; %   9446.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.04; beta_EDP=0.40;
          
          if     i_story == 3
              
              if     system == 1 || system == 3 || system == 4 % Non-isolated SCBF or Isolated SCBF (RI=1 or 2)
                     xm_Cost=9446.;
              elseif system == 2 || system == 5 || system == 6 % Non-isolated SMRF or Isolated SMRF (RI=1 or 2)
                     xm_Cost=(10246.     + 9446.)/2.0;   % Corrected to account for size of elements in gravity frames
              end
              
          elseif i_story == 5
              
              if     system == 1 || system == 3 || system == 4 % Non-isolated SCBF or Isolated SCBF (RI=1 or 2)
                     xm_Cost=9446.;
              elseif system == 2 || system == 5 || system == 6 % Non-isolated SMRF or Isolated SMRF (RI=1 or 2)
                     xm_Cost=9446.;
              end
              
          end
              
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.07; beta_EDP=0.40;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 2
          
          xm_EDP=0.07; beta_EDP=0.40;
                    
          if     i_story == 3
          
              if     system == 1 || system == 3 || system == 4 % Non-isolated SCBF or Isolated SCBF (RI=1 or 2)
                     xm_Cost=11246.;
              elseif system == 2 || system == 5 || system == 6 % Non-isolated SMRF or Isolated SMRF (RI=1 or 2)
                     xm_Cost=(13012.     + 11246.) / 2.0;   % Corrected to account for size of elements in gravity frames
              end
              
          elseif i_story == 5
              
              if     system == 1 || system == 3 || system == 4 % Non-isolated SCBF or Isolated SCBF (RI=1 or 2)
                     xm_Cost=11246.;
              elseif system == 2 || system == 5 || system == 6 % Non-isolated SMRF or Isolated SMRF (RI=1 or 2)
                     xm_Cost=11246.;
              end
              
          end
          
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.10; beta_EDP=0.40;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 3
          
          xm_EDP=0.10; beta_EDP=0.40;
                    
          if     i_story == 3
          
              if     system == 1 || system == 3 || system == 4 % Non-isolated SCBF or Isolated SCBF (RI=1 or 2)
                     xm_Cost=38473.;
              elseif system == 2 || system == 5 || system == 6 % Non-isolated SMRF or Isolated SMRF (RI=1 or 2)
                     xm_Cost=42533.;
              end
              
          elseif i_story == 5
              
              if     system == 1 || system == 3 || system == 4 % Non-isolated SCBF or Isolated SCBF (RI=1 or 2)
                     xm_Cost=38473.;
              elseif system == 2 || system == 5 || system == 6 % Non-isolated SMRF or Isolated SMRF (RI=1 or 2)
                     xm_Cost=38473.;
              end
              
          end
          
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end
       
       
    % % % % % % % % % % % % % % % %
    elseif i_m == 3
    % % % % % % % % % % % % % % % %
        
        numCompPerStory = 32;
        
       if i_n == 0  % Column (W<=W27)
          
          xm_EDP=0.03; beta_EDP=0.30; xm_Cost=0.%   16033.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.03; beta_EDP=0.30; xm_Cost=16033.;
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.04; beta_EDP=0.30;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 2
          
          xm_EDP=0.04; beta_EDP=0.30; xm_Cost=25933.;
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.05; beta_EDP=0.30;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 3
          
          xm_EDP=0.05; beta_EDP=0.30; xm_Cost=25933.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end
       
       
    % % % % % % % % % % % % % % % %
    elseif i_m == 4
    % % % % % % % % % % % % % % % %
        
        numCompPerStory = 24;
        
             % Determining different global slenderness ratio x_GP (=KL/r; Lignos and Karamanci, 2013) 
             % for different systems and different stories (K=1; L=all same; r=varies) - - - - - - - - -
             K=1.0; L = 0.70*sqrt((15*12)^2+(12*12)^2);% Unit=in
             
               if     (i_story == 1 && system == 1) || (i_story == 2 && system == 1)
                  r = 2.89;
               elseif (i_story == 3 && system == 1) || (i_story == 4 && system == 1)
                  r = 2.32;
               elseif (i_story == 5 && system == 1)
                  r = 1.96;
               elseif (i_story == 6 && system == 1)
                  r = 1.69;
                  
               elseif (i_story == 1 && system == 3) || (i_story == 2 && system == 3)
                  r = 2.89;
               elseif (i_story == 3 && system == 3) || (i_story == 4 && system == 3)
                  r = 2.49;
               elseif (i_story == 5 && system == 3)
                  r = 1.96;
               elseif (i_story == 6 && system == 3)
                  r = 1.48;
                  
               elseif (i_story == 1 && system == 4) || (i_story == 2 && system == 4)
                  r = 2.18;
               elseif (i_story == 3 && system == 4) || (i_story == 4 && system == 4)
                  r = 1.96;
               elseif (i_story == 5 && system == 4)
                  r = 1.61;
               elseif (i_story == 6 && system == 4)
                  r = 1.34;
               end
               
               x_GP = K*L/r;
             %   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                   
       if (i_story == 1 && system == 1) || (i_story == 2 && system == 1) || (i_story == 1 && system == 3) || (i_story == 2 && system == 3) %  for first story non-isolated SCBF & isolated SCBF w. Ri=1, thicker HSS is more expensive for DS2 & DS3.
            
           if i_n == 0  % Round HSS (60kg/m < brace weight < 147kg/m); Round HSS

              xm_EDP=0.41  /100.; beta_EDP=0.51; xm_Cost=0.; % added "/100." on 19.Apr.2020
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (KL/r) - - - - -
                    xm_GP=63.6; beta_GP=0.46;
                    F_DS_ij_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
                    
              PDS_ij_EDP = 1.0 - F_DS_ij*F_DS_ij_GP;   % j=0, i.e. no damage

           elseif i_n == 1

              xm_EDP=0.41  /100.; beta_EDP=0.51; xm_Cost=29983.;  % added "/100." on 19.Apr.2020
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.96  /100.; beta_EDP=0.45;                  % added "/100." on 19.Apr.2020
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (KL/r) - - - - -
                    xm_GP=63.6; beta_GP=0.46;
                    F_DS_i1_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    xm_GP=66.1; beta_GP=0.45;
                    F_DS_i2_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
              
              PDS_ij_EDP = F_DS_i1*F_DS_i1_GP - F_DS_i2*F_DS_i2_GP;    % 1<=j<=n, i.e. some damage

           elseif i_n == 2

              xm_EDP=0.96  /100.; beta_EDP=0.45; xm_Cost=47115.;     % added "/100." on 19.Apr.2020
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=2.75  /100.; beta_EDP=0.51;                     % added "/100." on 19.Apr.2020
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=66.1; beta_GP=0.45;
                    F_DS_i1_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    xm_GP=68.9; beta_GP=0.40;
                    F_DS_i2_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
              
              PDS_ij_EDP = F_DS_i1*F_DS_i1_GP - F_DS_i2*F_DS_i2_GP;   % 1<=j<=n, i.e. some damage

           elseif i_n == 3

              xm_EDP=2.75  /100.; beta_EDP=0.51; xm_Cost=47882.;      % added "/100." on 19.Apr.2020
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=68.9; beta_GP=0.40;
                    F_DS_ij_GP    = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
                    
              PDS_ij_EDP = F_DS_ij*F_DS_ij_GP;             % j=n, i.e. biggest damage

           end
       
       else  % for others, the price will be lower.
           
           if i_n == 0  % Round HSS (brace weight < 60kg/m); Round HSS

              xm_EDP=0.41  /100.; beta_EDP=0.51; xm_Cost=0.%   29983.;   % added "/100." on 19.Apr.2020
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=63.6; beta_GP=0.46;
                    F_DS_ij_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
                    
              PDS_ij_EDP = 1.0 - F_DS_ij*F_DS_ij_GP;   % j=0, i.e. no damage

           elseif i_n == 1

              xm_EDP=0.41  /100.; beta_EDP=0.51; xm_Cost=29983.;   % added "/100." on 19.Apr.2020
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.96  /100.; beta_EDP=0.45;                   % added "/100." on 19.Apr.2020
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=63.6; beta_GP=0.46;
                    F_DS_i1_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    xm_GP=66.1; beta_GP=0.45;
                    F_DS_i2_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
              
              PDS_ij_EDP = F_DS_i1*F_DS_i1_GP - F_DS_i2*F_DS_i2_GP;    % 1<=j<=n, i.e. some damage

           elseif i_n == 2

              xm_EDP=0.96  /100.; beta_EDP=0.45; xm_Cost=37014.;   % added "/100." on 19.Apr.2020
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=2.75  /100.; beta_EDP=0.51;                   % added "/100." on 19.Apr.2020
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=66.1; beta_GP=0.45;
                    F_DS_i1_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    xm_GP=68.9; beta_GP=0.40;
                    F_DS_i2_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
              
              PDS_ij_EDP = F_DS_i1*F_DS_i1_GP - F_DS_i2*F_DS_i2_GP;  % 1<=j<=n, i.e. some damage

           elseif i_n == 3

              xm_EDP=2.75  /100.; beta_EDP=0.51; xm_Cost=36480.;     % added "/100." on 19.Apr.2020
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=68.9; beta_GP=0.40;
                    F_DS_ij_GP    = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
                    
              PDS_ij_EDP = F_DS_ij*F_DS_ij_GP;             % j=n, i.e. biggest damage

           end
       
       end
           
           
    % % % % % % % % % % % % % % % %
    elseif i_m == 5
    % % % % % % % % % % % % % % % %
        
        numCompPerStory = 8;
        
        if system == 1  || system == 3  || system == 4 % if SCBF (isolated OR non-isolated) 
            
        
           if i_n == 0  % Moment connection; one-sided; <= W27

              xm_EDP=0.03; beta_EDP=0.30; xm_Cost=0.%   16033.;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage

           elseif i_n == 1

              xm_EDP=0.03; beta_EDP=0.30; xm_Cost=16033.;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.04; beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage

           elseif i_n == 2

              xm_EDP=0.04; beta_EDP=0.30; xm_Cost=25933.;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.05; beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage

           elseif i_n == 3

              xm_EDP=0.05; beta_EDP=0.30; xm_Cost=25933.;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage

           end
       
           
        elseif system == 2  || system == 5  || system == 6 % if SMRF (isolated OR non-isolated) 
            
        
           if i_n == 0  % RBS connection; one-sided; <= W27

              xm_EDP=0.01;   beta_EDP=0.17; xm_Cost=0.%   16033.;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage

           elseif i_n == 1

              xm_EDP=0.01;   beta_EDP=0.17; xm_Cost=0.;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.0216; beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage

           elseif i_n == 2

              xm_EDP=0.0216; beta_EDP=0.30; xm_Cost=16033.;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.05; beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage

           elseif i_n == 3

              xm_EDP=0.05;   beta_EDP=0.30; xm_Cost=25933.;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage

           end
       
        end
           
           
    % % % % % % % % % % % % % % % %
    elseif i_m == 6
    % % % % % % % % % % % % % % % %
        
       numCompPerStory = 8;  % only located at odd-numbered stories (= 1, 3, 5)
        
       if system == 1  || system == 3  || system == 4 % if SCBF (isolated OR non-isolated) 

            
           if i_n == 0  % Moment connection; two-sided; <= W27

              xm_EDP=0.03; beta_EDP=0.30; xm_Cost=0.%   30400.;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage

           elseif i_n == 1

              xm_EDP=0.03; beta_EDP=0.30; xm_Cost=30400.;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.04; beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage

           elseif i_n == 2

              xm_EDP=0.04; beta_EDP=0.30; xm_Cost=47000.;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.05; beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage

           elseif i_n == 3

              xm_EDP=0.05; beta_EDP=0.30; xm_Cost=47000.;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage

           end
       
           
        elseif system == 2  || system == 5  || system == 6 % if SMRF (isolated OR non-isolated) 

            
           if i_n == 0  % Moment connection; two-sided; <= W27

              xm_EDP=0.01; beta_EDP=0.17; xm_Cost=0.%   30400.;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage

           elseif i_n == 1

              xm_EDP=0.01;   beta_EDP=0.17; xm_Cost=0.;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.0216; beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage

           elseif i_n == 2

              xm_EDP=0.0216; beta_EDP=0.30; xm_Cost=26567.;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.05;   beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage

           elseif i_n == 3

              xm_EDP=0.05;   beta_EDP=0.30; xm_Cost=46999.;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
              
           end
          
       end
       
       
    % % % % % % % % % % % % % % % %
    elseif i_m == 7
    % % % % % % % % % % % % % % % %
        
        numCompPerStory = 80;
        
       if i_n == 0  % Shear tab connections
          
          xm_EDP=0.04; beta_EDP=0.40; xm_Cost=0.%   12107.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.04; beta_EDP=0.40; xm_Cost=12107.;
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.08; beta_EDP=0.40;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 2
          
          xm_EDP=0.08; beta_EDP=0.40; xm_Cost=12357.;
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.11; beta_EDP=0.40;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 3
          
          xm_EDP=0.11; beta_EDP=0.40; xm_Cost=12357.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end

       
    % % % % % % % % % % % % % % % %
    elseif i_m == 8
    % % % % % % % % % % % % % % % %
    
        numCompPerStory = 21*(9.144*9.144); % total area of each floor (including area of elevator shafts for simplicity)
        
       if i_n == 0  % Corrugated slab

          xm_EDP=0.00375; beta_EDP=0.13; xm_Cost=0.;%   180.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.00375; beta_EDP=0.13; xm_Cost=180.;
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.01; beta_EDP=0.22; 
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 2
               
          xm_EDP=0.01; beta_EDP=0.22; xm_Cost=330.;
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.05; beta_EDP=0.35;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 3
          
          xm_EDP=0.05; beta_EDP=0.35; xm_Cost=570.;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end
   
    end
    
end


