classdef WD
    %WD 
    %   ͼ����ǿʹ����̬ѧ����
    
    properties
    end
    
    methods
    end
    
    methods( Static )
        [labelImg, labels] = vesselCluster( orImg, enImg )
        
        function enImg = contrastEnhance( inImg )
            enImg = WD.enhance( inImg, 5000 );
        end
        
    end
    
    methods( Static, Access = private )
        enImg = enhance( inImg, P);
    end
end

