function [m] = info_num_Components_Structural(i_story, system)

% This function file returns from story number IDs (for each story):
% m = number of damageable components at a story.. note: 1 to m corresponds to component ID


if system == 1 % '1'=Nonisolated SCBF - - - - - - - - - -
    
    
    switch i_story

        case 1  % i_story = 1

%           m = [1,  3,4,5,6,7,8];
            m = [1,  3,4, 7,8];
            
        case 2  % i_story = 2

%           m = [    3,4,5,6,7,8];
            m = [    3,4,    7,8];
            
        case 3  % i_story = 3

%           m = [  2,3,4,5,6,7,8];
            m = [  2,3,4,    7,8];

        case 4  % i_story = 4

%           m = [    3,4,5,6,7,8];
            m = [    3,4,    7,8];
            
        case 5  % i_story = 5

%           m = [  2,3,4,5,6,7,8];
            m = [  2,3,4,    7,8];

        case 6  % i_story = 6

%           m = [    3,4,5,6,7,8];
            m = [    3,4,    7,8];
            
    end
    
    
elseif system == 2 % '2'=Nonisolated SMRF - - - - - - - - - -
    
    
    switch i_story
        
        case 1  % i_story = 1

            m = [1,  3,  5,6,7,8];

        case 2  % i_story = 2

            m = [    3,  5,6,7,8];

        case 3  % i_story = 3

            m = [  2,3,  5,6,7,8];

        case 4  % i_story = 4

            m = [    3,  5,6,7,8];

        case 5  % i_story = 5

            m = [  2,3,  5,6,7,8];


        case 6  % i_story = 6

            m = [    3,  5,6,7,8];
    end
    
    
elseif system == 3 || system == 4 % '3'=Isolated SCBF (RI=1) (Lower Bound); '4'=Isolated SCBF (RI=2) (Lower Bound); - - - - - - - - - -

    
    switch i_story

        case 1  % i_story = 1

%           m = [1,  3,4,5,6,7,8];
            m = [1,  3,4,    7,8];
            
        case 2  % i_story = 2

%           m = [    3,4,5,6,7,8];
            m = [    3,4,    7,8];
            
        case 3  % i_story = 3

            m = [  2,3,4,5,6,7,8];

        case 4  % i_story = 4

%           m = [    3,4,5,6,7,8];
            m = [    3,4,    7,8];
            
        case 5  % i_story = 5

%           m = [  2,3,4,5,6,7,8];
            m = [  2,3,4,    7,8];

        case 6  % i_story = 6

%           m = [    3,4,5,6,7,8];
            m = [    3,4,    7,8];
            
    end
    
    
elseif system == 5 || system == 6 % '5'=Isolated SMRF (RI=1) (Lower Bound); '6'=Isolated SMRF (RI=2) (Lower Bound); - - - - - - - - - -

    
    switch i_story

        case 1  % i_story = 1

            m = [1,  3,  5,6,7,8];

        case 2  % i_story = 2

            m = [    3,  5,6,7,8];

        case 3  % i_story = 3

            m = [  2,3,  5,6,7,8];

        case 4  % i_story = 4

            m = [    3,  5,6,7,8];

        case 5  % i_story = 5

            m = [  2,3,  5,6,7,8];


        case 6  % i_story = 6

            m = [    3,  5,6,7,8];
    end
  
end
