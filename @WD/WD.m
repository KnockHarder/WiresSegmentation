classdef WD
    %WD 
    %   图像增强使用形态学运算
    
    properties
    end
    
    methods
    end
    
    methods( Static )
        enhI = vesselEnhanceLDE(param);
        labels = vesselCluster( inImg )
        
        function enImg = contrastEnhance( inImg )
            param_LDE.origin = inImg;
            param_LDE.Img = WD.enhance( inImg, 5000 );
            param_LDE.detector_orientation  = -90:15:90;
            param_LDE.offset_factor         = 0.6;
            enImg = WD.vesselEnhanceLDE( param_LDE );
        end
        
    end
    
    methods( Static, Access = private )
        enImg = enhance( inImg, P);
    end
end

