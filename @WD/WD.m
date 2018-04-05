classdef WD
    %WD 
    %   图像增强使用形态学运算
    
    properties
    end
    
    methods
    end
    
    methods( Static )
        enImg = areaCut( inImg, minValue, maxValue, phases );
        enImg = contrastEnhance( inImg, P, alpha );
    end
        
end

