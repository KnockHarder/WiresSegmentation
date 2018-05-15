function enImg = enhance( inImg, P)
%ENHANCE 此处显示有关此函数的摘要
%   此处显示详细说明
assert( ndims( inImg ) == 2 );

[x, y] = size( inImg );
enImg = zeros( x, y );

bw = edge( inImg, 'sobel' );
se = strel( 'disk', 2 );
bw = imclose( bw, se );
% figure, imshow( bw );
se = strel( 'disk', 1 );
er = imerode( bw, se );

wiresImg = inImg .* er;
histValue = imhist( wiresImg );
halfPost = histValue( ceil(end/2): end );
dctV = dct( halfPost );
k = 30;
idctV = idct( dctV(1:1:k) );
[~,peak] = max( abs(idctV) );
lambda = peak/30/2 + 0.3;

enImg = (inImg > lambda) .* bw;
end

