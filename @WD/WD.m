classdef WD
    %WD 
    %   ͼ����ǿʹ����̬ѧ����
    
    properties
    end
    
    methods
    end
    
    methods( Static )
        [LDE, pointsImg, labelImg] = distanceMethod( grayImg, enImg, spacing, bar );        
        enImg = contrastEnhance( inImg, pixNum );
    end
    
    methods( Static, Access = private )
        enImg = enhance( inImg, P);
    end
end

