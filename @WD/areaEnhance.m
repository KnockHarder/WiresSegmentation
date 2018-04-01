function enImg = areaEnhance( inImg )
%AREAENHANCE 去除背景噪声，提取出铁丝网的区域
%   由于金属丝丝为极丝的状结构，找出边缘后，执行开操作，仍呈现出丝状；
%	由于金属丝相互交错，呈网状，且几乎都是网络的部分，因此可以使用
%areaopen操作，将其他边缘过滤掉

bw = edge( inImg, 'sobel' );
se = strel( 'disk', 2 );
bw = imclose( bw, se );

enImg = bwareaopen( bw, 5000 );

end

