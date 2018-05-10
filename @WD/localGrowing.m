function labelImg = localGrowing( grayImg, enImg, iter )
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
    
    variances = var(evidences, 1, 3);
    variances = variances * l_theta;   
    eI = edge( grayImg, 'sobel' );
    str = strel( 'disk', 2 );
    bw = imclose( eI, str );
    allwires = (bw .* I) > 0;

    skI = bwmorph( allwires, 'thin', inf ) | bwmorph( bw, 'thin', inf );    
    labelImg = zeros(m,n);
    labelImg( variances <= varTresh & skI == 1 ) = -2;
    inds = find( variances > varTresh  &  skI == 1 );
    lines = cell( m*n, 1 );
    l_eachLine = zeros(m*n, 1);
    for i = 1 : length(inds)
        idx = inds(i);
        [x,y] = ind2sub( [m,n], idx );
        labelImg(x,y) = i;
        lines{i} = {[x,y]};
        l_eachLine(i) = 1;
    end
    l_lines = length(inds);
    
    [~, directV] = max(evidences, [], 3);
    directH = rem( directV + step, l_theta );
    directH( directH == 0 ) = l_theta;
    masks_band = filters > 0;
    minlen = sz + 1;
    smooth = 0.25;
    [X,Y] = meshgrid(-sz:sz);
    sigma = 2;
    G = exp(-(X.^2+Y.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
    
    add_seg = cell( m+n, 1 );
    preI = labelImg * 0;
    whilecount = 0; %%%%%%%%%%%%%%
    assert( minlen > sz );
    while ~isequal( preI, labelImg )
        %%%%%%%%%
        if whilecount > iter
            break;
        end
        whilecount = whilecount + 1;
        %%%%%%%%%
        preI = labelImg;
        for label = 1 :1 :l_lines
            if l_eachLine(label) == 0
                continue;
            end
            for ii = -1 : 2 :1
                theLine = lines{label};
                l_theline = l_eachLine(label);
                
                if ii == -1
                    idx = 1;
                    x = theLine{idx}(1);
                    y = theLine{idx}(2);
                else
                    idx = l_theline;
                    x = theLine{idx}(1);
                    y = theLine{idx}(2);
                end
                if l_theline < minlen
                    x_pre = -1;
                    y_pre = -1;
                else
                    idx = idx - ii * (sz -1);
                    x_pre = theLine{idx}(1);
                    y_pre = theLine{idx}(2);
                end
                
                count_seg = 0;
                while 1
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
                        x_find = x + ( x_end - sz - 1);
                        y_find = y + ( y_end - sz - 1 );

                        if labelImg(x_find, y_find) > 0
                            label_f = labelImg(x_find, y_find); 
                            line_f = lines{label_f};
                            idx_f = find( cellfun( @(x) isequal(x, [x_find, y_find]), line_f ), 1 );
                            
                            len_f = l_eachLine(label_f);
                            if len_f < 2
                                labelImg(x_find, y_find) = label;
                                count_seg = count_seg + 1;
                                add_seg{count_seg} = [x_find, y_find];
                                l_eachLine(label_f) = 0;
                                break;
                            end
                            
                            A = [x,y] - [x_pre, y_pre];
                            B = [x,y] - [x_find,y_find];
                            if x_pre == -1
                                alpha = 0;
                            else
                                alpha = acos( dot(-A,B)/norm(A)/norm(B) );
                            end
                            if alpha > pi * smooth
                                mask(x_end,y_end) = 0;
                                continue;
                            end
                            
                            idx_pre = max( 1, idx_f - sz );
                            pre_f = line_f{idx_pre};
                            idx_next = min( len_f, idx_f + sz );
                            next_f = line_f{idx_next};
                            A = pre_f - [x_find, y_find];
                            B = next_f - [x_find, y_find];
                            C = [ x - x_find, y - y_find ];
                            if idx_f == 1 || idx_f == len_f
                                if idx_f == 1
                                    alpha = inf;
                                    beta = 0;
                                    theta = inf;
                                else
                                    alpha = 0;
                                    beta = inf;
                                    theta = inf;
                                end
                            else                                
                                theta_n = all_theta( directV(x_find,y_find) );
                                if idx_pre == 1
                                    A1 = [ -cos(theta_n), sin(theta_n) ];
                                    A2 = -A1;
                                    if dot(A1, -A) >= 0
                                        A = A1;
                                    else
                                        A = A2;
                                    end 
                                end
                                if idx_next == len_f
                                    B1 = [ -cos(theta_n), sin(theta_n) ];
                                    B2 = -B1;
                                    if dot(B1,B) >= 0
                                        B = B1;
                                    else
                                        B = B2;
                                    end
                                end
                                
                                alpha = acos( dot(-A,C)/norm(A)/norm(C) );
                                beta = acos( dot(-B,C)/norm(B)/norm(C) );
                                theta = acos( dot(-A,B)/norm(A)/norm(B) );
                            end
                                
                            angles = [ alpha, beta, theta ];
                            [~,idx_min] = min( angles );
                            if  idx_min == 3
                                mask( x_end, y_end ) = 0;
                                continue;
                            end
                            
                            dx = x_find - x;
                            dy = y_find - y;
                            delta_x = sign(dx);
                            delta_y = sign(dy);
                            if abs(dx) > abs(dy)
                                k = dy / dx;
                                for i_x = x+delta_x :delta_x :x_find-delta_x
                                    i_y = k * (i_x - x ) + y;
                                    i_y = floor( i_y + 0.5 );
                                    if labelImg(i_x, i_y) <= 0
                                        if labelImg(i_x, i_y) < -1
                                            lbn_n = labelImg( i_x-sz:i_x+sz, i_y-sz:i_y+sz );
                                            [directV(i_x, i_y), directH(i_x, i_y)] = adjustDir(...
                                                lbn_n, label, masks_band, l_theta );
                                        end
                                        labelImg(i_x, i_y) = label;
                                        count_seg = count_seg + 1;
                                        add_seg{count_seg} = [i_x, i_y];
                                    end
                                end
                            else
                                k = dx / dy;
                                for i_y = y+delta_y :delta_y :y_find-delta_y
                                    i_x = k * (i_y - y ) + x;
                                    i_x = floor( i_x + 0.5 );
                                    if labelImg(i_x, i_y) <= 0 
                                        if labelImg(i_x, i_y) < -1
                                            lbn_n = labelImg( i_x-sz:i_x+sz, i_y-sz:i_y+sz );
                                            [directV(i_x, i_y), directH(i_x, i_y)] = adjustDir(...
                                                lbn_n, label, masks_band, l_theta );
                                        end
                                        labelImg(i_x, i_y) = label;
                                        count_seg = count_seg + 1;
                                        add_seg{count_seg} = [i_x, i_y];
                                    end
                                end
                            end

                            if idx_min == 1
                                count = count_seg;
                                count_seg = count + idx_f;
                                add_seg( count+1:1:count_seg ) = line_f(idx_f:-1:1);
                                for jjj = 1 : 1 : idx_f
                                    labelImg( line_f{jjj}(1), line_f{jjj}(2) ) = label;
                                end
                                x = line_f{1}(1);
                                y = line_f{1}(2);

                                lines{label_f}(1:1:len_f-idx_f) = line_f( idx_f+1:1:len_f );
                                l_eachLine(label_f) = len_f - idx_f;
                            else
                                count = count_seg;
                                count_seg = count + len_f - idx_f + 1;
                                add_seg( count+1:1:count_seg ) = line_f(idx_f:1:len_f);
                                for jjj = idx_f : 1 : len_f
                                    labelImg( line_f{jjj}(1), line_f{jjj}(2) ) = label;
                                end
                                x = line_f{len_f}(1);
                                y = line_f{len_f}(2);

                                l_eachLine(label_f) = idx_f - 1;
                            end
                            break;
                        else
                            A = [x,y] - [x_find,y_find];
                            if x_pre == -1
                                theta = 0;
                            else
                                B = [x,y] - [x_pre,y_pre];
                                theta = acos( dot(-A,B)/norm(A)/norm(B) );
                            end

                            if theta < smooth * pi
                                dx = x_find - x;
                                dy = y_find - y;
                                delta_x = sign(dx);
                                delta_y = sign(dy);

                                if abs(dx) > abs(dy)
                                    k = dy / dx;
                                    for i_x = x+delta_x :delta_x :x_find-delta_x
                                        i_y = k * (i_x - x ) + y;
                                        i_y = floor( i_y + 0.5 );
                                        if labelImg(i_x, i_y) <= 0
                                            if labelImg(i_x, i_y) < -1
                                                labels_nei = labelImg( i_x-sz:i_x+sz, i_y-sz:i_y+sz );
                                                [directV(i_x,i_y), directH(i_x,i_y)] = ...
                                                        adjustDir( labels_nei, label, masks_band, l_theta );
                                            end
                                            labelImg(i_x, i_y) = label;
                                            count_seg = count_seg + 1;
                                            add_seg{ count_seg } = [i_x, i_y];
                                        end
                                    end
                                else
                                    k = dx / dy;
                                    for i_y = y+delta_y :delta_y :y_find-delta_y
                                        i_x = k * (i_y - y ) + x;
                                        i_x = floor( i_x + 0.5 );
                                        if( labelImg(i_x, i_y) <= 0 )
                                            if labelImg(i_x, i_y) < -1
                                                labels_nei = labelImg( i_x-sz:i_x+sz, i_y-sz:i_y+sz );
                                                [directV(i_x,i_y), directH(i_x,i_y)] = ...
                                                        adjustDir( labels_nei, label, masks_band, l_theta );
                                            end
                                            labelImg(i_x, i_y) = label;
                                            count_seg = count_seg + 1;
                                            add_seg{ count_seg } = [i_x, i_y];
                                        end
                                    end
                                end
                                x = x_find;  y = y_find;
                                
                                if labelImg(x, y) < -1
                                    labels_nei = labelImg( x-sz:x+sz, y-sz:y+sz );
                                    [directV(x,y), directH(x,y)] = adjustDir( labels_nei, ...
                                        label, masks_band, l_theta );
                                end
                                labelImg(x,y) = label;
                                count_seg = count_seg + 1;
                                add_seg{ count_seg } = [x, y];
                                break;
                            else
                                mask( x_end, y_end ) = 0;
                            end
                        end
                    end
                    
                    if count_seg > sz
                        point_pre = add_seg{ count_seg - sz };
                        x_pre = point_pre(1);
                        y_pre = point_pre(2);
                    end               
                    if isempty( find( mask > 0, 1 ) )
                        break;
                    end
                end
                
                if ii == -1
                    lines{label}(1:1:count_seg) = add_seg( count_seg:-1:1 );
                    lines{label}(count_seg+1:1:count_seg+l_theline) = ...
                        theLine(1:1:l_theline);
                    l_eachLine(label) = count_seg + l_theline;
                else
                    lines{label}(l_theline+1:1:l_theline+count_seg) = ...
                        add_seg( 1:1:count_seg );
                    l_eachLine(label) = count_seg + l_theline;
                end
            end
        end
        
        rem_labels = find( l_eachLine(1:1:l_lines) > 0 );
        l_lines = length( rem_labels );
        for i = 1 : 1 :l_lines
            label = rem_labels(i);
            theLine = lines{label};
            l_theline = l_eachLine(label);
            for idx = 1 : 1 : l_theline
                x = theLine{idx}(1);
                y = theLine{idx}(2);
                labelImg(x,y) = i;
            end
            lines{i}(1:1:l_theline) = lines{label}(1:1:l_theline);
            l_eachLine(i) = l_eachLine(label);
        end
    end
    
    labelImg = labelImg * 0;
    minlen = 20;
    rem_labels = find( l_eachLine(1:1:l_lines) > minlen );
    l_lines = length( rem_labels );
    for i = 1 : 1 :l_lines
        label = rem_labels(i);
        theLine = lines{label};
        l_theline = l_eachLine(label);
        for idx = 1 : 1 : l_theline
            x = theLine{idx}(1);
            y = theLine{idx}(2);
            labelImg(x,y) = i;
        end
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