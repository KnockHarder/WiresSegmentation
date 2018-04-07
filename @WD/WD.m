classdef WD
    %WD 
    %   图像增强使用形态学运算
    
    properties
    end
    
    methods
    end
    
    methods( Static )
        enImg = areaCut( inImg, minValue, maxValue, phases );
        
        function enImg = contrastEnhance( inImg, alpha )
            enImg = WD.enhance( inImg, 5000, alpha );
        end
    end
    
    methods( Static, Access = private )
        enImg = enhance( inImg, P, alpha );
    end
end

