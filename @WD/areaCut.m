function enImg = areaCut( inImg, minValue, maxValue, phases )
%AREAENHANCE ȥ��������������ȡ����˿����
%  �ҳ���Ե��ȥ������
%  ʹ��imclose��ȥ������
%  ʹ��areaopen��ȥ��С������
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

