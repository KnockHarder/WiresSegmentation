function labelImg = localGrowing( grayImg, enImg )
    step = 6;
    inc = 90 / step;
    all_theta = inc:inc:180;
    sz = 6;
    
    I = enImg;
    [m, n] = size( I );
    l_theta = length( all_theta );
    assert( l_theta == 2 * step );
    assert( m > sz * 3  &&  n > sz * 3 );
    filters = getLDE( sz, 0.6, all_theta );
    evidences = zeros( m, n, l_theta ); 
    varTresh = 0.2;
    for i = 1 : l_theta
        evidences( :, :, i ) = imfilter( I, filters(:,:,i), 'same', 'replicate' );
    end
    [~, directV] = max(evidences, [], 3);
    
    variances = var(evidences, 1, 3);
    variances = variances * l_theta;   
    eI = edge( grayImg, 'sobel' );
    str = strel( 'disk', 2 );
    bw = imclose( eI, str );
%     figure, imshow(bw);
    allwires = (bw .* I) > 0;
%     figure, imshow( allwires );    
    skI = bwmorph( allwires, 'thin', inf ) | bwmorph( bw, 'thin', inf );    
    labelImg = zeros(m,n);
    labelImg( variances > varTresh  &  skI == 1 ) = -1;
    labelImg( variances <= varTresh & skI == 1 ) = -2;
    
    inc = 90 / step;
    all_theta = inc:inc:180;
    masks_band = filters > 0;
    minlen = 2*sz + 1;
    smooth = 0.25;
    [X,Y] = meshgrid(-sz:sz);
    sigma = 2;
    G = exp(-(X.^2+Y.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
    
    lines = [];
    preI = labelImg * 0;
    ep = zeros( 2, 2 );
    pre_ep = zeros( 2, 2 );
    whilecount = 0; %%%%%%%%%%%%%%
    while ~isequal( preI, labelImg )
        %%%%%%%%%
        if whilecount == 1
            break;
        end
        whilecount = whilecount + 1,
%         if rem( whilecount, 1 ) == 0
%             WI = preI;
%             WI( WI < 0 ) = 0;
%             CI = label2rgb( WI, @jet, [.5, .5, .5 ] );
%             figure,imshow( CI );
%             
%             WI = labelImg;
%             WI( WI < 0 ) = 0;
%             CI = label2rgb( WI, @jet, [.5, .5, .5 ] );
%             figure,imshow( CI );
%             
%             WI = labelImg - preI;
%             WI(WI ~= 0) = 1;
%             CI = label2rgb( WI, @jet, [.5, .5, .5 ] );
%             figure,imshow( CI );
%         end
        %%%%%%%%%
        preI = labelImg;
        [labelImg, lines, directV] = vesselGrowing( labelImg, lines, masks_band,...
            directV, sz, minlen);
        directH = rem( directV + l_theta/2, l_theta );
        directH( directH == 0 ) = l_theta;
        %%%%%
%         labelImg( labelImg < 0 ) = 0;
%         return;
        %%%%%
        l_labels = length( lines );
        for label = 1 :1 :l_labels
            add_seg = cell( m+n, 2 );
            count_seg = [0 0];
            theLine = lines{label};
            l_line = length( theLine );
            if( l_line < minlen )
                labelImg( labelImg == label ) = -1;
                lines{label} = [];
                continue;
            end
            %%%%%%%%%%%%%%
            assert( length(theLine) == sum(sum(labelImg == label)) );
            ep(1,:) = theLine{1};         pre_ep(1,:) = theLine{sz};
            ep(2,:) = theLine{l_line};    pre_ep(2,:) = theLine{l_line-sz+1};
            
            for ii = 1 : 2
                x = ep(ii,1);           y = ep(ii, 2);
                x_pre = pre_ep(ii,1);   y_pre = pre_ep(ii,2);
                %%%%%%%%%%
%                 wwcount = 0;
                %%%%%%%%%%
                while 1
                    %%%%%%%
%                     wwcount = wwcount + 1;
%                     if rem(wwcount, 10) == 0
%                         WI = labelImg;
%                         WI( WI < 0 ) = 0;
%                         CI = label2rgb( WI, @jet, [.5, .5, .5 ] );
%                         figure, imshow( CI );
%                     end
                    %%%%%%%%%
                    dV = directV( x, y );
                    dH = directH( x, y );
                    labels_nei = labelImg( x-sz:x+sz, y-sz:y+sz);
                    dV_nei = directV( x-sz:x+sz, y-sz:y+sz);
                    dH_nei = directH( x-sz:x+sz, y-sz:y+sz);
                    
                    mask = labels_nei ~= 0  &  labels_nei ~= label;
                    mask = mask .* masks_band(:,:,dV);
                    diffV = dV_nei + 0.5 - dV;
                    diffV = abs( floor(diffV + sign(diffV)) );
                    diffH = dH_nei + 0.5 - dH;
                    diffH = abs( floor(diffH + sign(diffH)) );
                    diff = min( diffV, diffH );
                    errG = (1-G) .* (diff.^2);
                    
                    while ~isempty( find(mask > 0, 1) )
                        errors = mask .* errG;
                        errors( errors == 0 ) = inf;
                        [minCol, X] = min( errors );
                        [~,y_end] = min( minCol );
                        x_end = X(y_end);
                        x_next = x + ( x_end - sz - 1);
                        y_next = y + ( y_end - sz - 1 );

                        if labelImg(x_next, y_next) > 0
                            label_n = labelImg(x_next, y_next); 
                            line_n = lines{label_n};
                            idx = find( cellfun( @(x) isequal(x, [x_next, y_next]), line_n ) );
                            assert( ~isempty(idx) ); % for debuging
                            
                            lln = length(line_n);
                            if lln < minlen
                                labelImg( labelImg == label_n ) = -1;
                                lines{label_n} = [];
                                continue;
                            end
                            %%%%%%%%%%%%%%%%%
                            assert( length(line_n) == sum(sum(labelImg == label_n)) );
                            if idx > sz  &&  idx < lln - sz
                                pre_n = line_n{idx - sz};
                                next_n = line_n{idx + sz};

                                A = [ pre_n(1) - x_next, pre_n(2) - y_next ];
                                B = [ next_n(1) - x_next, next_n(2) - y_next ];
                                C = [ x - x_next, y - y_next ];
                                alpha = acos( dot(-A,C)/norm(A)/norm(C) );
                                beta = acos( dot(-B,C)/norm(B)/norm(C) );
                                theta = acos( dot(-A,B)/norm(A)/norm(B) );
                                
                                angles = [ alpha, beta, theta ];
                                [minAng,idx_minAng] = min( angles );
                                if  idx_minAng == 3  ||  minAng > smooth * pi
                                    mask( x_end, y_end ) = 0;
                                else
                                    dx = x_next - x;
                                    dy = y_next - y;
                                    delta_x = sign(dx);
                                    delta_y = sign(dy);
                                    if abs(dx) > abs(dy)
                                        k = dy / dx;
                                        for i_x = x+delta_x :delta_x :x_next-delta_x
                                            i_y = k * (i_x - x ) + y;
                                            i_y = floor( i_y + 0.5 );
                                            if labelImg(i_x, i_y) <= 0
                                                if norm(i_x-x_pre, i_y-y_pre) >= sz
                                                    labels_nei = labelImg( i_x-sz:i_x+sz, i_y-sz:i_y+sz );
                                                    [directV(i_x,i_y), directH(i_x,i_y)] = ...
                                                        adjustDir( labels_nei, label, masks_band, l_theta );
                                                end
                                                labelImg(i_x, i_y) = label;
                                                count_seg(ii) = count_seg(ii) + 1;
                                                add_seg{ count_seg(ii), ii } = [i_x, i_y];
                                            end
                                        end
                                    else
                                        k = dx / dy;
                                        for i_y = y+delta_y :delta_y :y_next-delta_y
                                            i_x = k * (i_y - y ) + x;
                                            i_x = floor( i_x + 0.5 );
                                            if( labelImg(i_x, i_y) <= 0 )
                                                if norm(i_x-x_pre, i_y-y_pre) >= sz
                                                    labels_nei = labelImg( i_x-sz:i_x+sz, i_y-sz:i_y+sz );
                                                    [directV(i_x,i_y), directH(i_x,i_y)] = ...
                                                        adjustDir( labels_nei, label, masks_band, l_theta );
                                                end
                                                labelImg(i_x, i_y) = label;
                                                count_seg(ii) = count_seg(ii) + 1;
                                                add_seg{ count_seg(ii), ii } = [i_x, i_y];
                                            end
                                        end
                                    end
                                    
                                    if idx_minAng == 1
                                        idx_ep = 1;
                                        delta = -1;
                                    else
                                        idx_ep = lln;
                                        delta = 1;
                                    end
                                    count = count_seg(ii);
                                    count_seg(ii) = count + (idx_ep-idx)/delta + 1;
                                    add_seg( count+1:count_seg(ii), ii) = ...
                                        line_n(idx:delta:idx_ep);
                                    for jjj = idx : delta: idx_ep
                                        labelImg( line_n{jjj}(1), line_n{jjj}(2) ) = label;
                                    end
                                    
                                    lines{label_n} = line_n( idx-delta:-delta:lln-idx_ep+1 );

                                    x = line_n{idx_ep}(1);
                                    y = line_n{idx_ep}(2);
                                    x_pre = line_n{ idx_ep - delta * (sz-1) }(1);
                                    y_pre = line_n{ idx_ep - delta * (sz-1) }(2);
                                    break;
                                end
                            else
                                if idx <= sz
                                    next_n = line_n{idx + sz};
                                    B = [ next_n(1) - x_next, next_n(2) - y_next ];
                                    C = [ x - x_next, y - y_next ];
                                    gama = acos( dot(-B,C)/norm(B)/norm(C) );
                                    idx_ep = lln;
                                    delta = 1;
                                else
                                    pre_n = line_n{idx - sz};
                                    A = [ pre_n(1) - x_next, pre_n(2) - y_next ];
                                    C = [ x - x_next, y - y_next ];
                                    gama = acos( dot(-A,C)/norm(A)/norm(C) );
                                    idx_ep = 1;
                                    delta = -1;
                                end
                                if gama < pi * smooth
                                    count = count_seg(ii);
                                    count_seg(ii) = count + (idx_ep-idx)/delta + 1;
                                    add_seg( count+1:count_seg(ii), ii) = ...
                                        line_n(idx:delta:idx_ep);
                                    for jjj = idx:delta:idx_ep
                                        labelImg( line_n{jjj}(1), line_n{jjj}(2) ) = label;
                                    end
                                    
                                    lines{label_n} = line_n( idx-delta:-delta:lln-idx_ep+1 );

                                    x = line_n{idx_ep}(1);
                                    y = line_n{idx_ep}(2);
                                    x_pre = line_n{ idx_ep - delta * (sz-1) }(1);
                                    y_pre = line_n{ idx_ep - delta * (sz-1) }(2);
                                    break;
                                else
                                    mask( x_end, y_end ) = 0;
                                end
                            end
                        else
                            A = [x,y] - [x_pre,y_pre];
                            B = [x,y] - [x_next,y_next];
                            theta = acos( dot(-A,B)/norm(A)/norm(B) );

                            if theta < smooth * pi
                                dx = x_next - x;
                                dy = y_next - y;
                                delta_x = sign(dx);
                                delta_y = sign(dy);

                                if abs(dx) > abs(dy)
                                    k = dy / dx;
                                    for i_x = x+delta_x :delta_x :x_next-delta_x
                                        i_y = k * (i_x - x ) + y;
                                        i_y = floor( i_y + 0.5 );
                                        if labelImg(i_x, i_y) <= 0
                                            if norm(i_x-x_pre, i_y-y_pre) >= sz
                                                labels_nei = labelImg( i_x-sz:i_x+sz, i_y-sz:i_y+sz );
                                                [directV(i_x,i_y), directH(i_x,i_y)] = ...
                                                        adjustDir( labels_nei, label, masks_band, l_theta );
                                            end
                                            labelImg(i_x, i_y) = label;
                                            count_seg(ii) = count_seg(ii) + 1;
                                            add_seg{ count_seg(ii), ii } = [i_x, i_y];
                                        end
                                    end
                                else
                                    k = dx / dy;
                                    for i_y = y+delta_y :delta_y :y_next-delta_y
                                        i_x = k * (i_y - y ) + x;
                                        i_x = floor( i_x + 0.5 );
                                        if( labelImg(i_x, i_y) <= 0 )
                                            if norm(i_x-x_pre, i_y-y_pre) >= sz
                                                labels_nei = labelImg( i_x-sz:i_x+sz, i_y-sz:i_y+sz );
                                                [directV(i_x,i_y), directH(i_x,i_y)] = ...
                                                        adjustDir( labels_nei, label, masks_band, l_theta );
                                            end
                                            labelImg(i_x, i_y) = label;
                                            count_seg(ii) = count_seg(ii) + 1;
                                            add_seg{ count_seg(ii), ii } = [i_x, i_y];
                                        end
                                    end
                                end
                                if norm(B) > sz
                                    x_pre = x ;  y_pre = y ;
                                end
                                x = x_next;  y = y_next;
                                
                                if labelImg(x, y) < -1
                                    labels_nei = labelImg( x-sz:x+sz, y-sz:y+sz );
                                    [directV(x,y), directH(x,y)] = adjustDir( labels_nei, ...
                                        label, masks_band, l_theta );
                                end
                                labelImg(x,y) = label;
                                count_seg(ii) = count_seg(ii) + 1;
                                add_seg{ count_seg(ii), ii } = [x, y];
                                break;
                            else
                                mask( x_end, y_end ) = 0;
                            end
                        end
                    end
                    
                    if isempty( find( mask > 0, 1 ) )
                        break;
                    end
                end
            end
            lines{label} = [ add_seg( count_seg(1):-1:1, 1 ); theLine; ...
                        add_seg( 1:count_seg(2), 2 ) ];
        end
    end
    labelImg( labelImg < 0 ) = 0;
    labels = unique( labelImg );
    labels = labels( labels ~= 0 );
    for i = 1 : length(labels)
        lb = labels(i);
        labelImg( labelImg == lb ) = i;
    end
end

function [labelImg, lines, directV] = vesselGrowing( labelImg, lines, masks_band, ...
    directV, width, minlen)
%VESSELCLUSTER 此处显示有关此函数的摘要
%   此处显示详细说明
    hitsThresh = minlen;
    [~,~,l_theta] = size( masks_band );
    directH = rem( directV + l_theta/2, l_theta );
    directH( directH == 0 ) = l_theta;
    smooth = 0.25;
    sigma = 1;
    [X,Y] = meshgrid(-width:width);
    G = exp(-(X.^2+Y.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
    G = G / G( width+1, width+1 );
    
    [m,n] = size( labelImg );
    labelImg( 1:width, : ) = 0;
    labelImg( :, 1:width ) = 0;
    labelImg( m-width+1:m, : ) = 0;
    labelImg( :, n-width+1:n ) = 0;
    labelHits = zeros(m*n, 1);
    ori_lines = lines;
    l_lines = length( lines );
    maxLabel = l_lines + 1;
    lines = cell( m*n, 1 );
    lines(1:1:l_lines) = ori_lines(:);
%     maxp = sum( sum( labelImg < 0 ) );
% 
%     [x,y] = getNextPoint( labelImg, 1, 1, -1 );
%     count = 0;
    [x,y] = find( labelImg == -1, 1 );
    x_pre = -1; y_pre = -1;
    labelImg(x,y) = maxLabel;
    labelHits(maxLabel) = 0;
    lines{maxLabel} = cell( m+n, 1 );
    while ~isempty(x)
%         count = count + 1;
%         if rem( count, 1000 ) == 0
%             count, 
%         end
%         assert( count <= maxp );
        if norm(x-x_pre, y-y_pre) >= width
            labels_nei = labelImg( x-width:x+width, y-width:y+width );
            [directV(x,y), directH(x,y)] = ...
                adjustDir( labels_nei, maxLabel, masks_band, l_theta );
        end
        labelImg(x,y) = maxLabel;
        labelHits(maxLabel) = labelHits(maxLabel) + 1;
        lines{maxLabel}{labelHits(maxLabel)} = [x,y];
        dV = directV( x, y );
        dH = directH( x, y );
        labels_nei = labelImg( x-width:x+width, y-width:y+width );
        dV_nei = directV( x-width:x+width, y-width:y+width );
        dH_nei = directH( x-width:x+width, y-width:y+width );

        mask = labels_nei < 0;
        mask = mask .* masks_band(:,:,dV);
        diffV = dV_nei + 0.5 - dV;
        diffV = abs( floor( diffV + sign(diffV) ) );
        diffH = dH_nei + 0.5 - dH;
        diffH = abs( floor( diffH + sign(diffH) ) );
        diff = min( diffV, diffH );
        
        while ~isempty( find(mask > 0, 1) )
            errors = mask .* (1-G) .* (diff.^2);
            errors( errors == 0 ) = inf;
            [minCol, X] = min( errors );
            [~,y_end] = min( minCol );
            x_end = X(y_end);
            x_next = x + ( x_end - width - 1);
            y_next = y + ( y_end - width - 1 );
            assert( labels_nei(x_end,y_end) == labelImg(x_next,y_next) );
            
            if isequal( [x_pre,y_pre], [-1, -1] )
                theta = 0;
            else
                A = [x,y] - [x_pre,y_pre];
                B = [x_next,y_next] - [x,y];
                theta = acos( dot(A,B)/norm(A)/norm(B) );
            end

            if theta < smooth * pi
                dx = x_next - x;
                dy = y_next - y;
                delta_x = sign(dx);
                delta_y = sign(dy);

                if abs(dx) > abs(dy)
                    k = dy / dx;
                    for i_x = x+delta_x :delta_x :x_next-delta_x
                        i_y = k * (i_x - x ) + y;
                        i_y = floor( i_y + 0.5 );
                        if labelImg(i_x, i_y) <= 0
                            if norm(i_x-x_pre, i_y-y_pre) >= width
                                labels_nei = labelImg( i_x-width:i_x+width, i_y-width:i_y+width );
                                [directV(i_x,i_y), directH(i_x,i_y)] = ...
                                    adjustDir( labels_nei, maxLabel, masks_band, l_theta );
                            end
                            labelImg(i_x, i_y) = maxLabel;
                            labelHits(maxLabel) = labelHits(maxLabel) + 1;
                            lines{maxLabel}{labelHits(maxLabel)} = [i_x, i_y];
                        end
                    end
                else
                    k = dx / dy;
                    for i_y = y+delta_y :delta_y :y_next-delta_y
                        i_x = k * (i_y - y ) + x;
                        i_x = floor( i_x + 0.5 );
                        if( labelImg(i_x, i_y) <= 0 )
                            if norm(i_x-x_pre, i_y-y_pre) >= width
                                labels_nei = labelImg( i_x-width:i_x+width, i_y-width:i_y+width );
                                [directV(i_x,i_y), directH(i_x,i_y)] = ...
                                    adjustDir( labels_nei, maxLabel, masks_band, l_theta );
                            end
                            labelImg(i_x, i_y) = maxLabel;
                            labelHits(maxLabel) = labelHits(maxLabel) + 1;
                            lines{maxLabel}{labelHits(maxLabel)} = [i_x, i_y];
                        end
                    end
                end
                if  theta == 0 || norm(B) > width
                    x_pre = x ;  y_pre = y ;
                end
                x = x_next;  y = y_next;
                break;
            else
                mask( x_end, y_end ) = 0;
            end
        end
        
        if isempty( find(mask > 0, 1) )
            lines{maxLabel} = lines{maxLabel}( 1:labelHits(maxLabel) );
            [x,y] = find( labelImg == -1, 1 );
            x_pre = -1;  y_pre = -1;
            maxLabel = maxLabel + 1;
            labelHits(maxLabel) = 0;
            lines{maxLabel} = cell( m+n, 1 );
        end
    end

    labelHits = labelHits(1:maxLabel);
    lines = lines(1:maxLabel);
    labels = find( labelHits >= hitsThresh );
    badLabels = find( labelHits < hitsThresh );
    
    lines = lines( labelHits >= hitsThresh );
    
    for lb = badLabels'
        labelImg( labelImg == lb ) = -1;
    end
    for i = 1 : length(labels)
        labelImg( labelImg == labels(i) ) = i;
    end
end

function [dV, dH] = adjustDir( labels_nei, label, masks_band, l_theta)
    assert( rem(l_theta, 2) == 0 );
    preP = labels_nei == label;
    scores = zeros( l_theta, 1 );
    for i = 1 : l_theta
        scores(i) = sum(sum( preP .* masks_band(:,:,i) ));
    end
    [~,dV] = max(scores);
    dH = rem( dV + l_theta/2, l_theta );
    dH = dH + (dH == 0)*l_theta;
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