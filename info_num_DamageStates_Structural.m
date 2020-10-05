function [n] = info_num_DamageStates_Structural(i_m)

% This function file returns from story number IDs (for each story):
% n = number of damage states a component may experience (this is per m)

% i_m = ID of component

    switch i_m

        case 1  % i_m = 1

            n = 3;

        case 2  % i_m = 2

            n = 3;

        case 3  % i_m = 3

            n = 3;

        case 4  % i_m = 4

            n = 3;

        case 5  % i_m = 5

            n = 3;

        case 6  % i_m = 6

            n = 3;

        case 7  % i_m = 7

            n = 3;

        case 8  % i_m = 8;  25June2019 "corrugated slab" is "structural" component. (previously used as "non-structural - drift")

            n = 3;
    end

end


