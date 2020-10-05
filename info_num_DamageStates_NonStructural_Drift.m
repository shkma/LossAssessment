function [n] = info_num_DamageStates_NonStructural_Drift(i_m)

% This function file returns from story number IDs (for each story):
% n = number of damage states a component may experience (this is per m)

% i_m = ID of component

    switch i_m

        case 1  % i_m = 1

            n = 2;

        case 2  % i_m = 2

            n = 2;

        case 3  % i_m = 3

            n = 2;
            
    end

end

