function enhI = vesselEnhanceLDE(param)
%VESSELENHANCELDE Enhance the vessels via LDE

origin = param.origin;
I     = param.Img;                      %-- original signal
theta = param.detector_orientation;     %-- orientation space for detectors    
sigma = 1;                    %-- scale space
k     = param.offset_factor;            %-- d = k*sigma

d = k*sigma;
enhI = LDE(I,theta,d,sigma);

eg = edge( origin, 'sobel' );
str = strel( 'disk', 20 );
mask = imclose( eg, str );
mask = bwareaopen( mask, 5000 );
enhI = enhI .* mask;

end



function [resp] = LDE(I,all_theta,d,sigma)
%LDE Perform local directional evidence filtering
%   all_theta -- orientation for detection filter
%   psi       -- orientations of evidence filter, for each value of theta
%   d         -- d = k*sigma
l_theta             = length(all_theta);
ind                 = 0;
siz     = round(6*sigma);
[mI,nI] = size( I );
[X,Y]   = meshgrid(-siz:siz);
G       = exp( -(X.^2+Y.^2) / (2*sigma^2) ) / ( sqrt(2*pi) * sigma );
Gxx     = (X.^2-sigma^2).*G/(sigma^4);
Gyy     = (Y.^2-sigma^2).*G/(sigma^4);
Gxy     = (X.*Y).*G/(sigma^4);

for ii = 1 : l_theta
    rot = all_theta(ii);
    substack = zeros( mI, nI);

    ind = ind+1;
    
    if d == -1
        X_d = 0;
        Y_d = 0;
        G_d_f = 0;
        Gxx_d_f = 0;
        Gyy_d_f = 0;
        Gxy_d_f = 0;
        Gxx_d_b = 0;
        Gyy_d_b = 0;
        Gxy_d_b = 0;
    else
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
    end

    R_d = Gxx*(cosd(rot))^2 + Gyy*(sind(rot))^2 + Gxy*(sind(2*rot));
    R_b = Gxx_d_b*(cosd(rot))^2+(Gyy_d_b)*(sind(rot))^2 + (Gxy_d_b)*(sind(2*rot));
    R_f = Gxx_d_f*(cosd(rot))^2+(Gyy_d_f)*(sind(rot))^2 + (Gxy_d_f)*(sind(2*rot));

    I_d = imfilter(I,-sigma^1.5*R_d,'same','replicate');   % local detector
    I_b = imfilter(I,-sigma^1.5*R_b,'same','replicate');   % boosting -- backward
    I_f = imfilter(I,-sigma^1.5*R_f,'same','replicate');   % boosting -- forward
    I_theta = (I_d + 1*I_b + 1*I_f) ;                         % superposition
    
    t = find(I_theta<0);
    I_theta(t) = 0;
    stack(:,:,ii) = I_theta;
end

resp = max(stack,[],3);

end

