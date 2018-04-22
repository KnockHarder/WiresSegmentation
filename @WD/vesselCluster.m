function labels = vesselCluster( inImg )
%VESSELCLUSTER 此处显示有关此函数的摘要
%   此处显示详细说明
assert( ndims(inImg) == 2 );

I = inImg;
[m, n] = size( I );

thresh = 0.2;
sz = 6;
assert( m > sz * 3  &&  n > sz * 3 );
I( 1:sz, : ) = 0;
I( m-sz+1:m, :) = 0;
I( :, 1:sz ) = 0;
I( :, n-sz+1:n) = 0;

labels = ones( m, n );
labels = labels * -1;
labels(  I == 0  ) = 0;
nextLabel = 1;

all_theta = -90:15:90;
l_theta = length( all_theta );
filters = getLDE( sz, 0.6, all_theta );

evidences = zeros( l_theta, 1 );
for i = sz +1 : m - sz
    for j = sz +1 : n - sz
        if labels(i,j) == -1
            sam = I( i-sz:i+sz, j-sz:j+sz );
            samLables = labels( i-sz:i+sz, j-sz:j+sz );
            for ii = 1 : l_theta
                evidences(ii) = sum( sum( sam .* filters(:,:,ii) ) );    
            end
            [maxEvi, max_idx] = max( evidences );

            if maxEvi - I(i,j) < thresh
                labels(i,j) = 0;
            else
                maxFilter = filters( :, :, max_idx );
                pred = find( maxFilter > 0 & samLables > 0 );
                
                if( isempty( pred ) )
                    labels(i,j) = nextLabel;
                    nextLabel = nextLabel + 1;
                else
                    predLables = samLables( pred );
                    labels(i,j) = mode( predLables );
                end
            end
        end
    end
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