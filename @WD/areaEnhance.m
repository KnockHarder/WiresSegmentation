function enImg = areaEnhance( inImg )
%AREAENHANCE ȥ��������������ȡ����˿��������
%   ���ڽ���˿˿Ϊ��˿��״�ṹ���ҳ���Ե��ִ�п��������Գ��ֳ�˿״��
%	���ڽ���˿�໥��������״���Ҽ�����������Ĳ��֣���˿���ʹ��
%areaopen��������������Ե���˵�

bw = edge( inImg, 'sobel' );
se = strel( 'disk', 2 );
bw = imclose( bw, se );

enImg = bwareaopen( bw, 5000 );

end

