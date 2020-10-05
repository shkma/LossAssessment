function [m] = info_num_Components_NonStructural_Accel(i_floor)

% This function file returns from floor number IDs (for each floor):
% m = number of damageable components at a floor.. note: 1 to m corresponds to component ID

switch i_floor
    
    case 1  % i_floor = 1 = Ground (floor)
        
        m = [    3];

    case 2  % i_floor = 2
        
        m = [1,2, ];
        
    case 3  % i_floor = 3
        
        m = [1,2  ];
        
    case 4  % i_floor = 4
        
        m = [1,2  ];
        
    case 5  % i_floor = 5
        
        m = [1,2  ];
        
    case 6  % i_floor = 6
        
        m = [1,2  ];
        
    case 7  % i_floor = 7 = Roof
        
        m = [1,2  ];
        
end


