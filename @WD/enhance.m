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
histValue( 1: floor(end/2) ) = 0;
sumHist = histValue;
for     i = 2: 1: length(sumHist)
    sumHist(i) = sumHist(i) + sumHist(i-1);
end
[~,modeV] = max( histValue );
midV = find( sumHist >= sumHist(end)/2, 1 );

while   modeV ~= midV
    histValue(modeV) = 0;
    sumHist = histValue;
    for     i = 2: 1: length(sumHist)
        sumHist(i) = sumHist(i) + sumHist(i-1);
    end
    [~,modeV] = max( histValue );
    midV = find( sumHist >= sumHist(end)/2, 1 );
end
lambda = (modeV-10)/255;
enImg = (inImg > lambda) .* bw;
end

