classdef WD
    %WD 
    %   ͼ����ǿʹ����̬ѧ����
    
    properties
    end
    
    methods
    end
    
    methods( Static )
        enImg = areaCut( inImg, minValue, maxValue, phases );
        
        function enImg = contrastEnhance( inImg )
            enImg = WD.enhance( inImg, 5000);
%             figure, imhist( enImg );
        end
    end
    
    methods( Static, Access = private )
        enImg = enhance( inImg, P);
    end
end

