function enImg = contrastEnhance( inImg, P, alpha )
%ENHANCE 此处显示有关此函数的摘要
%   此处显示详细说明
assert( ndims( inImg ) == 2 );
assert( alpha > 1 );

[x, y] = size( inImg );
enImg = zeros( x, y );
values = zeros( 1, x * y );
len = uint64(1);

bw = edge( inImg, 'sobel' );
se = strel( 'disk', 2 );
bw = imclose( bw, se );
mask = bwareaopen( bw, P );
for i = 1 : x
    for j = 1 : y
        values( len ) = inImg(i,j) * mask(i,j);
        if values( len ) ~= 0
            len = len + 1;
        end
    end
end
values = values( 1:len );
lambda = median( values );

for i = 1 : x
    for j = 1 : y
        pixel = inImg(i,j);
        if abs( pixel - lambda ) > eps
            if pixel < lambda
                enImg(i,j) = lambda * pixel / ( alpha * lambda - ...
                    (alpha - 1 ) * pixel );
            else
                enImg(i,j) = 1 - ( 1 - pixel ) / ( ...
                    ( alpha - 1 ) * ( pixel - lambda ) / ( 1 - lambda )...
                    + 1 );
            end
        end
    end
end     

end

