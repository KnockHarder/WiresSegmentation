function [labelImg, labels] = vesselCluster( inImg )
%VESSELCLUSTER 此处显示有关此函数的摘要
%   此处显示详细说明
assert( ndims(inImg) == 2 );

I = inImg;
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

variances = var(evidences, 1, 3);
variances = variances * l_theta;
[maxEvi, directMax] = max(evidences, [], 3);

labelImg = zeros(m,n);
labelImg( variances > varTresh ) = -1;
labelImg( variances <= varTresh & maxEvi > 0.1 ) = -2;

maxLabel = 0;
labelHits = zeros(m*n, 1);
hitsThresh = 20;
bendError = 25;

labelImg( 1:sz, : ) = 0;
labelImg( :, 1:sz ) = 0;
labelImg( m-sz+1:m, : ) = 0;
labelImg( :, n-sz+1:n ) = 0;

sigma = 1;
X = 1 : (2 * sz + 1) ^ 2;
GX = exp( -(X.^2) / (2*sigma^2) ) / ( sqrt(2*pi) * sigma );

for i = sz + 1 : m - sz
    for j = sz + 1 : n -sz
        if labelImg(i,j) >= 0
            continue;
        end
        
        labels_neighbor = labelImg( i-sz:i+sz, j-sz:j+sz );
        directs_neighbor = directMax( i-sz:i+sz, j-sz:j+sz );
        directIJ = directMax(i,j);
        
        idx_band =  filters(:,:,directIJ) > 0 & labels_neighbor > 0 ;
        labels_band = unique( labels_neighbor( idx_band ) );
        
        if isempty(labels_band)
            if labelImg(i,j) == -1
                maxLabel = maxLabel + 1;
                labelImg(i,j) = maxLabel;
                labelHits(maxLabel) = 1;
            else
                labelImg(i,j) = 0;
            end
        else
            l_LB = length(labels_band);
            errors = zeros(l_LB, 1);
            for jj = 1 : l_LB
                label = labels_band(jj);
                [X,Y] = find( labels_neighbor == label );
                D = (X - sz - 1) .^ 2  +  (Y - sz - 1) .^ 2;
                
                Dlen = length(D);
                [~,idx_a] = sort( D, 'ascend' );
                for jjj = 1 : Dlen
                    idx = idx_a(jjj);
                    DD = D * 0;
                    DD(jjj) = directMax( X(idx), Y(idx) ) - directIJ;
                end
                errors(jj) = sum( sum( (DD .^ 2) .* GX(1:Dlen) ) );
                errors(jj) = errors(jj) / sum( GX(1:Dlen) ) * D( idx_a(1) );
            end
            
            [minError,idx] = min(errors);
            if minError < bendError
                label = labels_band(idx);
                labelImg(i,j) = label;
                labelHits(label) = labelHits(label) + 1;
            else
                if labelImg(i,j) == -1
                    maxLabel = maxLabel + 1;
                    labelImg(i,j) = maxLabel;
                    labelHits(maxLabel) = 1;
                else
                    labelImg(i,j) = 0;
                end
            end
        end
        
    end
end

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