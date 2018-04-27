function [labelImg, labels] = vesselGrowing( oriImg, enImg )
%VESSELCLUSTER 此处显示有关此函数的摘要
%   此处显示详细说明
    assert( ndims(enImg) == 2 );

    I = enImg;
    [m, n] = size( I );

    all_theta = -90:15:90;
    l_theta = length( all_theta );
    sz = 6; assert( m > sz * 3  &&  n > sz * 3 );
    filters = getLDE( sz, 0.6, all_theta );
    evidences = zeros( m, n, l_theta ); 
    varTresh = 0.2;
    for i = 1 : l_theta
        evidences( :, :, i ) = imfilter( I, filters(:,:,i), 'same', 'replicate' );
    end

    [maxEvi, directMax] = max(evidences, [], 3);
    weightImg = oriImg .* maxEvi;

    variances = var(evidences, 1, 3);
    variances = variances * l_theta;
    labelImg = zeros(m,n);
    labelImg( variances > varTresh ) = -1;
    labelImg( variances <= varTresh & weightImg > 0.1 ) = -2;
    
%     figure, imshow( imresize( label2rgb( labelImg * -1, @jet, [.5,.5,.5] ), 0.5 ) );

    hitsThresh = 10;
    smooth = 0.25;
    sigma = 2;
    width = 2;
    filters_band = filters( sz+1-width:sz+1+width, sz+1-width:sz+1+width, : );
    masks_band = filters_band > 0;
    [X,Y] = meshgrid(-width:width);
    G = exp(-(X.^2+Y.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
    G = G / G( width+1, width+1 );

    labelImg( 1:width, : ) = 0;
    labelImg( :, 1:width ) = 0;
    labelImg( m-width+1:m, : ) = 0;
    labelImg( :, n-width+1:n ) = 0;
    labelHits = zeros(m*n, 1);
%     maxp = sum( sum( labelImg < 0 ) );
% 
%     [x,y] = getNextPoint( labelImg, 1, 1, -1 );
%     count = 0;
    [x,y] = find( labelImg == -1, 1 );
    x_pre = -1; y_pre = -1;
    labelImg(x,y) = 1;
    maxLabel = 1;
    labelHits(maxLabel) = 0;
    while ~isempty(x)
%         count = count + 1;
%         if rem( count, 1000 ) == 0
%             count, 
%         end
%         assert( count <= maxp );
        
        labelImg(x,y) = maxLabel;
        labelHits(maxLabel) = labelHits(maxLabel) + 1;
        directIJ = directMax( x, y );
        labels_nei = labelImg( x-width:x+width, y-width:y+width );
        directs_nei = directMax( x-width:x+width, y-width:y+width );

        mask = labels_nei < 0;
        mask = mask .* masks_band(:,:,directIJ);
        D = directs_nei + mask * 0.5;
        D = abs(D - directIJ);
        
        while ~isempty( find(mask > 0, 1) )
            errors = mask .* D .* G ;
            errors( errors == 0 ) = inf;
            [minCol, X] = min( errors );
            [~,y_end] = min( minCol );
            x_end = X(y_end);
            x_next = x + ( x_end - width - 1);
            y_next = y + ( y_end - width - 1 );
            assert( labels_nei(x_end,y_end) == labelImg(x_next,y_next) );
            
            if isequal( [x_pre,y_pre], [-1, -1] )
                theta = 0;
            else
                A = [x,y] - [x_pre,y_pre];
                A = A / sqrt( A(1)^2 + A(2)^2 );
                B = [x_next,y_next] - [x,y];
                B = B / sqrt( B(1)^2 + B(2)^2 );
                theta = acos( A * B' );
            end

            if theta < smooth * pi
                if labelImg(x_next, y_next) < -1
                    directMax(x_next, y_next) = directMax(x,y);
                end

                dx = x_next - x;
                dy = y_next - y;

                if abs(dx) > abs(dy)
                    k = dy / dx;
                    for i_x = x+1 : x_next-1
                        i_y = k * (i_x - x ) + y;
                        i_y = floor( i_y + 0.5 );
                        if labelImg(i_x, i_y) <= 0
                            labelImg(i_x, i_y) = maxLabel;
                            labelHits(maxLabel) = labelHits(maxLabel) + 1;
                        end
                    end
                else
                    k = dx / dy;
                    for i_y = y+1 : y_next-1
                        i_x = k * (i_y - y ) + x;
                        i_x = floor( i_x + 0.5 );
                        if( labelImg(i_x, i_y) <= 0 )
                            labelImg(i_x, i_y) = maxLabel;
                            labelHits(maxLabel) = labelHits(maxLabel) + 1;
                        end
                    end
                end

                x_pre = x ;  y_pre = y ;
                x = x_next;  y = y_next;
                break;
            else
                mask( x_end, y_end ) = 0;
            end
        end
        
        if isempty( find(mask > 0, 1) )
%             [x,y] = getNextPoint( labelImg, x, y, -1 );
            [x,y] = find( labelImg == -1, 1 );
            x_pre = -1;  y_pre = -1;
            maxLabel = maxLabel + 1;
            labelHits(maxLabel) = 0;
            continue;
        end
    end

    labelImg( labelImg < 0 ) = 0;

    labelHits = labelHits(1:maxLabel);
    idx =  1:length(labelHits);
    for badLabel = idx( labelHits < hitsThresh )
        labelImg( labelImg == badLabel ) = 0;
    end
    labels = idx( labelHits >= hitsThresh );
end


function filters = getLDE( sz,d, orientations )
    all_theta = orientations;
    l_theta = length( all_theta );
    
    sigma = 1;
    filters = zeros( sz*2 + 1, sz*2 + 1, l_theta );
    [X,Y] = meshgrid( -sz:sz );
    
    G       = exp( -(X.^2+Y.^2) / (2*sigma^2) ) / ( sqrt(2*pi) * sigma );
    Gxx     = (X.^2-sigma^2).*G/(sigma^4);
    Gyy     = (Y.^2-sigma^2).*G/(sigma^4);
    Gxy     = (X.*Y).*G/(sigma^4);
    
    for i = 1 : l_theta
        rot = all_theta(i);
        
        X_d_f = X - d*sind(rot);
        Y_d_f = Y + d*cosd(rot);
        G_d_f = exp(-(X_d_f.^2+Y_d_f.^2)/(2*sigma^2))/( sqrt(2*pi) * sigma );
        Gxx_d_f = (X_d_f.^2-sigma^2).*G_d_f/(sigma^4);
        Gyy_d_f = (Y_d_f.^2-sigma^2).*G_d_f/(sigma^4);
        Gxy_d_f = (X_d_f.*Y_d_f).*G_d_f/(sigma^4);

        X_d_b = X - d*sind(rot);
        Y_d_b = Y - d*cosd(rot);
        G_d_b = exp(-(X_d_b.^2+Y_d_b.^2)/(2*sigma^2))/( sqrt(2*pi) * sigma );
        Gxx_d_b = (X_d_b.^2-sigma^2).*G_d_b/(sigma^4);
        Gyy_d_b = (Y_d_b.^2-sigma^2).*G_d_b/(sigma^4);
        Gxy_d_b = (X_d_b.*Y_d_b).*G_d_b/(sigma^4);
        
        R_d = Gxx*(cosd(rot))^2 + Gyy*(sind(rot))^2 + Gxy*(sind(2*rot));
        R_b = Gxx_d_b*(cosd(rot))^2+(Gyy_d_b)*(sind(rot))^2 + (Gxy_d_b)*(sind(2*rot));
        R_f = Gxx_d_f*(cosd(rot))^2+(Gyy_d_f)*(sind(rot))^2 + (Gxy_d_f)*(sind(2*rot));
        
        filters(:,:,i) = -sigma ^ 1.5 * ( R_d + R_b + R_f );
    end
end

% function [x_next,y_next] = getNextPoint( labelImg, x, y, label )
%     assert( ndims(labelImg) == 2 );
%     [m,n] = size( labelImg );
%     assert( x <=m  &  y <= n );
% 
%     x_next = x + 1;
%     y_next = y;
%     while 1
%         if x_next > m
%             x_next = 1;
%             y_next = y_next + 1;
%         end
%         if y_next > n
%             x_next = [];
%             y_next = [];
%             break;
%         end
% 
%         if labelImg(x_next,y_next) == label
%             break;
%         else
%             x_next = x_next + 1;
%         end
%     end
% end