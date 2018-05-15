function labelImg = vesselCluster( oriImg, bwImg )
%VESSELCLUSTER 此处显示有关此函数的摘要
%   此处显示详细说明
enImg = oriImg .* bwImg;
thresh = min( enImg(enImg > 0) );
thresh = thresh * 0.9;
display( thresh );
enImg = (enImg-thresh)/(1-thresh);
% figure, imshow( enImg );

I = enImg;
[m, n] = size( I );
inc = 15;
step = 90/inc;
assert( step == 6 );
all_theta = inc:inc:180;
l_theta = length( all_theta );
sz = 6; assert( m > sz * 3  &&  n > sz * 3 );
filters = getLDE( sz, 0.6, all_theta );

evidences = zeros( m, n, l_theta ); 
varTresh = 0.2;
for i = 1 : l_theta
    evidences( :, :, i ) = imfilter( I, filters(:,:,i), 'same', 'replicate' );
end
variances = var(evidences, 1, 3);
variances = variances * l_theta;
[maxEvi, ~] = max(evidences, [], 3);
wid = sz;
pointLabel = zeros(m,n);
pointLabel( variances > varTresh  & bwImg) = 1;
pointLabel( variances <= varTresh & bwImg ) = 2;
pointLabel( 1:wid, : ) = 0;
pointLabel( :, 1:wid ) = 0;
pointLabel( m-wid+1:m, : ) = 0;
pointLabel( :, n-wid+1:n ) = 0;
% labelImg = pointLabel;
% return;
I = oriImg;
for i = 1 : l_theta
    evidences( :, :, i ) = imfilter( I, filters(:,:,i), 'same', 'replicate' );
end
variances = var(evidences, 1, 3);
variances = variances * l_theta;
[maxEvi, directV] = max(evidences, [], 3);
directH = rem( directV + step, step*2 );
pointLabel( variances > varTresh & pointLabel == 2 ) = 1;

sigma = 1;
[X,Y] = meshgrid(-wid:wid);
G = exp(-(X.^2+Y.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
G = G / max(max(G));

hitsThresh = 50;
labelImg = zeros( m,n );
maxLabel = 0;
labelHits = zeros(m*n, 1);
clusterCenter = zeros( m*n, 2, 2 );
for     ind_p = find( pointLabel == 1 )'
    [x,y] = ind2sub( [m,n], ind_p );
    maxLabel = maxLabel + 1;
    labelImg(x,y) = maxLabel;
    clusterCenter(maxLabel,1,:) = [x,y];
    clusterCenter(maxLabel,2,:) = [x,y];
    labelHits(maxLabel) = 1;
end
for i = 1 :1 :maxLabel
    assert( sum(sum(labelImg == i)) == labelHits(i) );
end 
masks_band = getLDE( wid, 0.5, all_theta ) > 0;

preI = labelImg + 1;
whilecount = 0;
while   ~isequal( preI, labelImg )
    if  whilecount > 50
       break;
    end
    whilecount = whilecount + 1,
    
    preI = labelImg;
    for     label = 1: 1: maxLabel
        if  labelHits(label) == 0
            continue;
        end
        
        for     i = 1: 2
            x = clusterCenter(label,i,1);
            y = clusterCenter(label,i,2);
            labels_nei = labelImg( x-wid:x+wid, y-wid:y+wid );
            dV_nei = directV( x-wid:x+wid, y-wid:y+wid );
            dH_nei = directH( x-wid:x+wid, y-wid:y+wid );
            evi_nei = maxEvi( x-wid:x+wid, y-wid:y+wid );
            dV = directV( x, y );
            dH = directH( x, y );
            evi = maxEvi( x, y );
            
            diffV = dV_nei + 0.5 - dV;
            diffV = abs( floor(diffV + sign(diffV)) );
            diffH = dH_nei + 0.5 - dH;
            diffH = abs( floor(diffH + sign(diffH)) );
            diff = min( diffV, diffH );
            errG = (1-G) .* (diff.^2) .* abs(evi_nei - evi);
            
            mask = pointLabel( x-wid:x+wid, y-wid:y+wid ) > 0;
            mask( labels_nei == label ) = 0;
            mask = mask .* masks_band(:,:,dV);
            while   ~isempty( find( mask > 0, 1 ) )
                errors = mask .* errG;
                errors( errors == 0 ) = inf;
                [minCol, X] = min( errors );
                [~,y_end] = min( minCol );
                x_end = X(y_end);
                x_find = x + ( x_end - wid - 1);
                y_find = y + ( y_end - wid - 1 );
                
                if  diff(x_end,y_end) > step / 2
                    mask( x_end, y_end ) = 0;
                    continue;
                end
                
                label_find = labelImg( x_find, y_find );
                if  label_find > 0
                    source_find = squeeze( clusterCenter(label_find,:,:) );
                    idx_source = 0;
                    for     ii = 1 : 2
                        if  isequal( [x_find,y_find], source_find(ii,:) )
                            idx_source = ii;
                        end
                    end
                    
                    if  idx_source == 0
                        mask( x_end, y_end ) = 0;
                        continue;
                    end
                    
                    labelImg( labelImg == label_find ) = label;
                    labelHits( label ) = labelHits(label) + labelHits(label_find);
                    labelHits( label_find ) = 0;
                    clusterCenter( label, i, : ) = source_find( 3-idx_source, : );
                    break;
                else
                    labelImg(x_find, y_find) = label;
                    labelHits( label ) = labelHits( label ) + 1;
                    break;
                end
            end
        end
    end
end

for i = 1 :1 :maxLabel
    assert( sum(sum(labelImg == i)) == labelHits(i) );
end
goodLabels = find( labelHits > hitsThresh );
temp = labelImg;
labelImg = zeros(m,n);
for i = 1: 1: length(goodLabels)
    labelImg( temp == goodLabels(i) ) = i;
end
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