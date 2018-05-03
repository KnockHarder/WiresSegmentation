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
    for i = 1 : length(inds)
        idx = inds(i);
        [x,y] = ind2sub( [m,n], idx );
        labelImg(x,y) = i;
        lines{i} = {[x,y]};
    end
    lines = lines(1:length(inds));
    
    [~, directV] = max(evidences, [], 3);
    directH = rem( directV + step, l_theta );
    directH( directH == 0 ) = l_theta;
    masks_band = filters > 0;
    minlen = 2*sz + 1;
    smooth = 0.25;
    [X,Y] = meshgrid(-sz:sz);
    sigma = 2;
    G = exp(-(X.^2+Y.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
    
    preI = labelImg * 0;
    whilecount = 0; %%%%%%%%%%%%%%
    while ~isequal( preI, labelImg )
        %%%%%%%%%
        if whilecount == 10
            break;
        end
        whilecount = whilecount + 1,
        %%%%%%%%%
        preI = labelImg;
        l_labels = length( lines );
        for label = 1 :1 :l_labels
            if isempty( lines{label} )
                continue;
            end
            for ii = -1 : 2 :1
                theLine = lines{label};
                l_line = length( theLine );
                assert( length(theLine) == sum(sum(labelImg == label)) );
                
                if ii == -1
                    idx = 1;
                    x = theLine{idx}(1);
                    y = theLine{idx}(2);
                else
                    idx = l_line;
                    x = theLine{idx}(1);
                    y = theLine{idx}(2);
                end
                if l_line < minlen
                    x_pre = -1;
                    y_pre = -1;
                else
                    idx = idx - ii * (sz -1);
                    x_pre = theLine{idx}(1);
                    y_pre = theLine{idx}(2);
                end
                
                add_seg = cell( m+n, 1 );
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
                        x_next = x + ( x_end - sz - 1);
                        y_next = y + ( y_end - sz - 1 );

                        if labelImg(x_next, y_next) > 0
                            label_n = labelImg(x_next, y_next); 
                            line_n = lines{label_n};
                            idx_n = find( cellfun( @(x) isequal(x, [x_next, y_next]), line_n ) );
                            assert( ~isempty(idx_n) ); % for debuging
                            
                            lln = length(line_n);
                            if lln < minlen
                                labelImg( labelImg == label_n ) = -1;
                                lines{label_n} = [];
                                continue;
                            end
                            %%%%%%%%%%%%%%%%%
                            assert( length(line_n) == sum(sum(labelImg == label_n)) );
                            if idx_n > sz  &&  idx_n < lln - sz
                                A = [ x - x_next, y - y_next ];
                                if x_pre == -1
                                    gama = 0;
                                else
                                    B = [ x - x_pre, y - y_pre ];
                                    gama = acos( dot(-A,B)/norm(A)/norm(B) );
                                end
                                if gama > pi * smooth
                                    mask( x_end, y_end ) = 0;
                                    continue;
                                end
                                
                                pre_n = line_n{idx_n - sz};
                                next_n = line_n{idx_n + sz};
                                A = pre_n - [x_next, y_next];
                                B = next_n - [x_next, y_next];
                                C = [ x - x_next, y - y_next ];
                                
                                alpha = acos( dot(-A,C)/norm(A)/norm(C) );
                                beta = acos( dot(-B,C)/norm(B)/norm(C) );
                                theta = acos( dot(-A,B)/norm(A)/norm(B) );
                                
                                angles = [ alpha, beta, theta ];
                                [~,idx_min] = min( angles );
                                if  idx_min == 3
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
                                        for i_y = y+delta_y :delta_y :y_next-delta_y
                                            i_x = k * (i_y - y ) + x;
                                            i_x = floor( i_x + 0.5 );
                                            if( labelImg(i_x, i_y) <= 0 )
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
                                        idx_ep = 1;
                                        delta = -1;
                                    else
                                        idx_ep = lln;
                                        delta = 1;
                                    end
                                    count = count_seg;
                                    count_seg = count + (idx_ep-idx_n)/delta + 1;
                                    add_seg( count+1:count_seg ) = line_n(idx_n:delta:idx_ep);
                                    for jjj = idx_n : delta: idx_ep
                                        labelImg( line_n{jjj}(1), line_n{jjj}(2) ) = label;
                                    end
                                    
                                    lines{label_n} = line_n( idx_n-delta:-delta:lln-idx_ep+1 );

                                    x = line_n{idx_ep}(1);
                                    y = line_n{idx_ep}(2);
                                    x_pre = line_n{ idx_ep - delta * (sz-1) }(1);
                                    y_pre = line_n{ idx_ep - delta * (sz-1) }(2);
                                    break;
                                end
                            else
                                if idx_n <= sz
                                    next_n = line_n{idx_n + sz};
                                    B = next_n - [x_next, y_next];
                                    C = [ x - x_next, y - y_next ];
                                    gama = acos( dot(-B,C)/norm(B)/norm(C) );
                                    idx_ep = lln;
                                    delta = 1;
                                else
                                    pre_n = line_n{idx_n - sz};
                                    A = pre_n - [x_next, y_next];
                                    C = [ x - x_next, y - y_next ];
                                    gama = acos( dot(-A,C)/norm(A)/norm(C) );
                                    idx_ep = 1;
                                    delta = -1;
                                end
                                if gama < pi * smooth
                                    count = count_seg;
                                    count_seg = count + (idx_ep-idx_n)/delta + 1;
                                    add_seg( count+1:count_seg ) = ...
                                        line_n(idx_n:delta:idx_ep);
                                    for jjj = idx_n:delta:idx_ep
                                        labelImg( line_n{jjj}(1), line_n{jjj}(2) ) = label;
                                    end
                                    
                                    lines{label_n} = line_n( idx_n-delta:-delta:lln-idx_ep+1 );

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
                            A = [x,y] - [x_next,y_next];
                            if x_pre == -1
                                theta = 0;
                            else
                                B = [x,y] - [x_pre,y_pre];
                                theta = acos( dot(-A,B)/norm(A)/norm(B) );
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
                                    for i_y = y+delta_y :delta_y :y_next-delta_y
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
                                A = [x_next,y_next] - [x,y];
                                if norm(A) > sz
                                    x_pre = x ;  y_pre = y ;
                                end
                                x = x_next;  y = y_next;
                                
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
                    
                    if isempty( find( mask > 0, 1 ) )
                        break;
                    end
                end
                
                if ii == -1
                    lines{label} = [ add_seg( count_seg:-1:1 ); lines{label} ];
                else
                    lines{label} = [ lines{label}; add_seg( 1:1:count_seg ) ];
                end
            end
        end
        
        labels = unique( labelImg );
        labels = labels( labels > 0 );
        l_labels = length( labels );
        for i = 1 : 1 :l_labels
            label = labels(i);
            theLine = lines{label};
            for idx = 1 : 1 : length(theLine)
                x = theLine{idx}(1);
                y = theLine{idx}(2);
                labelImg(x,y) = i;
            end
            lines{i} = theLine;
        end
        tempLines = lines;
        lines = cell( m*n, 1 );
        lines(1:l_labels) = tempLines(1:l_labels);
        inds = find( labelImg == -1 );
        for i = inds'
            [x,y] = ind2sub( size(labelImg), i );
            l_labels = l_labels + 1;
            labelImg(x,y) = l_labels;
            lines{l_labels} = {[x,y]};
        end
        lines = lines( 1:1:l_labels );
    end
    
    labelImg( labelImg < 0 ) = 0;
    labels = unique( labelImg );
    labels = labels( labels ~= 0 );
    l_labels = length(labels);
    l_lines = zeros( l_labels, 1 );
    for i = 1 : l_labels
        l_lines(i) = length( lines{i} );
    end
    lines = lines( l_lines > minlen );
    l_labels = length(lines);
    labelImg = zeros( size(labelImg) );
    for i = 1 : 1: l_labels
        theLine = lines{i};
        for idx = 1 : 1 : length(theLine)
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