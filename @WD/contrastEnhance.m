function enImg = contrastEnhance( inImg, P)
%ENHANCE 此处显示有关此函数的摘要
%   此处显示详细说明
assert( ndims( inImg ) == 2 );

[x, y] = size( inImg );
enImg = zeros( x, y );
interested = zeros( 1, x * y );
len = uint64(1);

bw = edge( inImg, 'sobel' );
se = strel( 'disk', 2 );
bw = imclose( bw, se );
% figure, imshow( bw );
mask = bwareaopen( bw, P );
bw = mask;
for i = 1 : x
    for j = 1 : y
        if mask(i,j) > 0
            interested( len ) = inImg(i,j) * mask(i,j);
            len = len + 1;
        end
    end
end
interested = interested( 1:len );
interested = sort( interested );
maxValue = interested( len );
nValues = length( unique( interested ) );
frequence = zeros( 1, nValues );
index = zeros( 1, nValues );
k = interested(1);
k_index = 1;
for i = 1 : len
    if interested(i) == k
        frequence( k_index ) = frequence( k_index ) + 1;
    else
        k = interested(i);
        k_index = k_index + 1;
    end
end
index(1) = 1;
for i = 2 : nValues
    index(i) = index(i - 1) + frequence( i - 1 );
end

[modeF, freqIndex] = max( frequence );
modeI = index( freqIndex );
modeV = interested( modeI );
midV = interested( round( len / 2 ) );
while  abs(modeV - midV) > 0.01
    prePart = interested( 1 : (modeI - 1) );
    postPart = interested( (modeI + modeF) : len );
    interested = [ prePart postPart ];
    len = len - modeF;
    
    prePart = frequence( 1 : (freqIndex - 1 ) );
    postPart = frequence( (freqIndex + 1) : nValues );
    frequence = [ prePart postPart ];
    prePart = index( 1 : (freqIndex - 1 ) );
    postPart = index( (freqIndex + 1) : nValues );
    index = [ prePart postPart ];
    nValues = nValues - 1;
    for i = freqIndex : nValues
        index(i) = index(i) - modeF;
    end
    
    [modeF, freqIndex] = max( frequence );
    modeI = index( freqIndex );
    modeV = interested( modeI );
    midV = interested( round( len / 2 ) );
end
lambda = modeV;

for i = 1 : x
    for j = 1 : y
%         pixel = inImg(i,j) * bw(i,j);
        pixel = inImg(i,j);
        if abs( pixel - lambda ) > eps
            if pixel < lambda 
                enImg(i,j) = 0;
            else
                enImg(i,j) = power( ( pixel - lambda ) / ...
                    ( maxValue - lambda ), 3 );
            end
        end
    end
end     

end

