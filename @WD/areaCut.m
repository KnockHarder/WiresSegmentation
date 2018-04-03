function enImg = areaCut( inImg, minValue, maxValue, phases )
%AREAENHANCE 去除背景噪声，提取出铁丝区域
%  找出边缘，去除背景
%  使用imclose，去除噪声
%  使用areaopen，去除小块区域
bw = edge( inImg, 'sobel' );
se = strel( 'disk', 2 );
bw = imclose( bw, se );

[x,y] = size( bw );
imgs = zeros( x, y, phases );

step = round( (maxValue - minValue) / phases );
for  i = 1 : phases
  imgs( :,:,i ) = bwareaopen( bw, minValue + step ); 
  N(i) = sum( sum( imgs(:,:,i) > 0 ) );
end

D = abs( diff(N) ); 
[~,peak] = min( D );
enImg = squeeze( imgs(:,:,peak) );

end

