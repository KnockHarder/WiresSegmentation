classdef WD
    %WD 
    %   ͼ����ǿʹ����̬ѧ����
    
    properties
    end
    
    methods
    end
    
    methods( Static )
        labelImg = localGrowing( grayImg, enImg );
        
        function enImg = contrastEnhance( inImg, pixNum )
            enImg = WD.enhance( inImg, pixNum );
        end
        
    end
    
    methods( Static, Access = private )
        enImg = enhance( inImg, P);
    end
end

