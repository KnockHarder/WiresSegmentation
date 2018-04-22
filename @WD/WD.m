classdef WD
    %WD 
    %   图像增强使用形态学运算
    
    properties
    end
    
    methods
    end
    
    methods( Static )
        labels = vesselCluster( inImg )
        
        function enImg = contrastEnhance( inImg )
            enImg = WD.enhance( inImg, 5000 );
        end
        
    end
    
    methods( Static, Access = private )
        enImg = enhance( inImg, P);
    end
end

