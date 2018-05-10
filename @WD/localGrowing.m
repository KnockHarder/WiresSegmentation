function labelImg = localGrowing( grayImg, enImg )
    step = 6;
    inc = 90 / step;
    all_theta = inc:inc:180;
    sz = 6;
    rad_band = 4;
    
    I = enImg;
    [m, n] = size( I );
    l_theta = length( all_theta );
    assert( l_theta == 2 * step );
    assert( m > sz * 3  &&  n > sz * 3 );
    filters = getLDE( sz, 0.6, all_theta );
    evidences = zeros( m, n, l_theta ); 
    varTresh = 0.2;
    for idx_theta = 1 : l_theta
        evidences( :, :, idx_theta ) = imfilter( I, filters(:,:,idx_theta), 'same', 'replicate' );
    end
    variances = var(evidences, 1, 3);
    variances = variances * l_theta;   
    eI = edge( grayImg, 'sobel' );
    str = strel( 'disk', 2 );
    bw = imclose( eI, str );
    allwires = (bw .* I) > 0;

    skI = bwmorph( allwires, 'thin', inf ) | bwmorph( bw, 'thin', inf );  
    inds = find( variances > varTresh  &  skI == 1 );
    [X,Y] = meshgrid(-sz:sz);
    mask_disk = (X.^2 + Y.^2) <= sz^2;
    dots = cell( length(inds), 1 );
    l_dots = 0;
    sub2idx_dot = zeros(m,n);
    for idx_ind = 1 : length(inds)
        theInd = inds(idx_ind);
        [x,y] = ind2sub( [m,n], theInd );
        if  sum(sum( sub2idx_dot( x-sz:x+sz, y-sz:y+sz ) .* mask_disk )) == 0
            l_dots = l_dots + 1;
            dots{l_dots} = [x,y];
            sub2idx_dot(x,y) = l_dots;
        end
    end
    
    ind_ap = find( allwires | imfilter(sub2idx_dot,ones(3,3)) );
    l_ap = length( ind_ap );
    sub2idx_ap = zeros(m,n);
    for idx_ap = 1: 1: l_ap
        [x,y] = ind2sub( [m,n], ind_ap(idx_ap) );
        sub2idx_ap(x,y) = idx_ap;
    end
    %%%%%%%%%%%
%     bw_dots = sub2idx_dot > 0;
%     figure, imshow( bw_dots );
%     bw_ap = sub2idx_ap > 0;
%     figure, imshow( bw_ap );
%     bw_dots = imfilter( bw_dots, ones(3,3) ) > 0;
%     figure, imshow( bw_dots );
%     figure, imshow( bw_dots & (1-bw_ap) );
    %%%%%%%%%%%
    I = grayImg;
    for idx_theta = 1 : l_theta
        evidences( :, :, idx_theta ) = imfilter( I, filters(:,:,idx_theta), 'same', 'replicate' );
    end
    [maxEvi,~] = max( evidences, [], 3 );
    maxEvi = maxEvi / max(max(maxEvi)) * m;
    [DX,DY,DE] = get8diff( maxEvi );
    
    
    I_nei9 = ones(3,3);
    dis_d2p = zeros( l_dots, l_ap );
    par_d2p = zeros( l_dots, l_ap );
%     l_dots,
    [X,Y] = meshgrid(-rad_band:rad_band, -rad_band:rad_band);
    mask_disk = X.^2 + Y.^2 <= rad_band^2;
    for idx_dot = 1: 1: l_dots
        %%%%
        if  rem( idx_dot, 200 ) == 0
            display(idx_dot);
        end
        %%%%
        x_dot = dots{idx_dot}(1);
        y_dot = dots{idx_dot}(2);
        idx_sourceP = sub2idx_ap( x_dot, y_dot );
        distence = inf - zeros(1, l_ap);
        distence(idx_sourceP) = 0;
        parents = zeros(1, l_ap );
        visited = zeros( 1, l_ap );
        
        for i = 1 : 8
            x_i = x_dot + DX(i);
            y_i = y_dot + DY(i);
            idx_i = sub2idx_ap( x_i, y_i );
            if  idx_i < 1
                continue;
            end
            distence(idx_i) = DE(x_dot, y_dot, i);
            parents(idx_i) = idx_sourceP;
        end
        visited(idx_sourceP) = 1;
        temp = distence;
        temp(idx_sourceP) = inf;
        [~,idx_u] = min( temp );
        ind_u = ind_ap(idx_u);
        [x_u,y_u] = ind2sub( [m,n], ind_u );
        for i = 1: l_ap - 1
            for j = 1 : 8
                x_i = x_u + DX(j);
                y_i = y_u + DY(j);
                idx_i = sub2idx_ap( x_i, y_i );
                if  idx_i < 1
                    continue;
                end
                if  visited(idx_i) < 1
                    dis_new = distence(idx_u) + DE(x_u,y_u,j);
                    if  dis_new < distence(idx_i)
                        distence(idx_i) = dis_new;
                        parents(idx_i) = idx_u;
                    end
                end
            end
            
            visited( idx_u ) = 1;
            while 1
                temp = distence;
                temp( visited == 1 ) = inf;
                [min_dis,idx_u] = min( temp );
                if min_dis == inf
                    break;
                end
                
                ind_u = ind_ap(idx_u);
                [x_u, y_u] = ind2sub( [m,n], ind_u );
                if  sub2idx_dot( x_u, y_u ) > 0
                    diskNei = sub2idx_ap( x_u-rad_band:x_u+rad_band, y_u-rad_band:y_u+rad_band) & mask_disk;
                    connection = zeros( 2*rad_band + 1, 2*rad_band + 1 );
                    connection( rad_band+1, rad_band+1 ) = 1;
                    for i_con = 1: 1: rad_band
                        temp = imfilter( connection, I_nei9 );
                        temp = temp & diskNei;
                        
                        if isequal( temp, connection )
                            break;
                        end
                        connection = temp;
                    end
                    [X_con,Y_con] = find( connection > 0 );
                    DX_con = X_con - rad_band - 1;
                    DY_con = Y_con - rad_band - 1;
                    for i_con = 1: 1: length(DX_con)
                        x_con = x_u + DX_con(i_con);
                        y_con = y_u + DY_con(i_con);
                        idx_con = sub2idx_ap( x_con, y_con );
                        visited(idx_con) = 1;
                    end
                else
                    break;
                end
            end
            if min_dis == inf
                break;
            end
        end
        %%%%%%%%%%%%%
%         temp = distence;
%         temp( idx_sourceP ) = inf;
%         t_findD = find( temp < inf );
%         t_findP = find( parents > 0 );
%         assert( isequal( t_findD, t_findP ) );
        %%%%%%%%%%%%%%%%%%%%%%%%
        dis_d2p(idx_dot, :) = distence(:);
        par_d2p(idx_dot, :) = parents(:);
    end
    
    dis_dd = inf - zeros( l_dots, l_dots );
    for idx_headD = 1 : l_dots
        x_start = dots{idx_headD}(1);
        y_start = dots{idx_headD}(2);
        idx_startP = sub2idx_ap( x_start, y_start );
        for idx_endD = 1 : l_dots
            x_end = dots{idx_endD}(1);
            y_end = dots{idx_endD}(2);
            idx_endP = sub2idx_ap( x_end, y_end );
            if  dis_d2p( idx_headD, idx_endP ) == 0 
                dis_dd( idx_headD,idx_endD ) = 0;
                continue;
            end
            if dis_d2p( idx_headD, idx_endP ) == inf
                continue;
            end
                
            max_theta = 0;
            idx_midP = par_d2p( idx_headD, idx_endP );
            while idx_midP ~= idx_startP
%                 idx_midP, tempdis = dis_d2p(idx_headD, idx_endP),
                ind_mid = ind_ap( idx_midP );
                [x_mid, y_mid] = ind2sub( [m,n], ind_mid );
                A = [x_start, y_start] - [x_mid, y_mid];
                B = [x_end, y_end] - [x_mid, y_mid];
                if norm(A) > 3 && norm(B) > 3
                    theta = acos( dot(-A,B)/norm(A)/norm(B) );
                    if max_theta < theta
                        max_theta = theta;
                    end
                end
                idx_midP = par_d2p( idx_headD, idx_midP );
            end

            param_k = max_theta * step / pi + 1;
            dis_dd(idx_headD,idx_endD) = dis_d2p(idx_headD,idx_endP) * param_k;
        end
    end
    
    old_dis = dis_dd;
    dis_dd = min( dis_dd, dis_dd' );
    for     i = 1: 1: l_dots
        dis_dd(i,i) = inf;
    end
    connection = zeros(l_dots+2,l_dots+2);
    while   1
        count_con = sum( connection );
        sourceDot = find( count_con < 2, 1 );
        if  isempty( sourceDot )
            break;
        end
        
        [minValue,dot_minDis] = min( dis_dd(sourceDot,:) );
        if  minValue == inf
            connection(sourceDot, l_dots+1) = 1;
            connection(l_dots+1, sourceDot) = 1;
            
            connection(sourceDot, l_dots+2) = 1;
            connection(l_dots+2, sourceDot) = 1;
            continue;
        end
        
        dot_con = find( connection(sourceDot,:) > 0 );
        if  ~isempty(dot_con)
            assert( dot_con ~= l_dots + 1 );
            A = dots{dot_con} - dots{sourceDot};
            B = dots{dot_minDis} - dots{sourceDot};
            gama = acos( dot(-A,B)/norm(A)/norm(B) );
        else
            gama = 0;
        end
        if  gama >= pi/2
            dis_dd(sourceDot,dot_minDis) = inf;
            dis_dd(dot_minDis,sourceDot) = inf;
            continue;
        end
        
        dots_comp = find( connection(dot_minDis,:) > 0 );
        dots_comp = dots_comp( dots_comp <= l_dots );
        if   ~isempty(dots_comp)
            dotA = dots_comp(1);
            dotB = 0;
            dotC = sourceDot;
            A = dots{dotA} - dots{dot_minDis};
            C = dots{dotC} - dots{dot_minDis};
            alpha = acos( dot(-A,C)/norm(A)/norm(C) );
            if  length(dots_comp) > 1
                dotB = dots_comp(2);
                B = dots{dotB} - dots{dot_minDis};
                beta = acos( dot(-B,C)/norm(B)/norm(C) );
                gama = acos( dot(-A,B)/norm(A)/norm(B) );
            else
                beta = inf;
                gama = inf;
            end
            [~, idx_min] = min( [alpha,beta,gama] );
            if  idx_min == 3
                dis_dd(sourceDot,dot_minDis) = inf;
                dis_dd(dot_minDis,sourceDot) = inf;
                continue;
            else
                connection(sourceDot,dot_minDis) = 1;
                connection(dot_minDis,sourceDot) = 1;
                dis_dd(sourceDot,dot_minDis) = inf;
                dis_dd(dot_minDis,sourceDot) = inf;
                if  idx_min == 1
                    if  dotB > 0
                        connection(dot_minDis,dotB) = 0;
                        connection(dotB,dot_minDis) = 0;
                    end
                else
                    connection(dot_minDis,dotA) = 0;
                    connection(dotA,dot_minDis) = 0;
                end
            end
        else
            connection(sourceDot,dot_minDis) = 1;
            connection(dot_minDis,sourceDot) = 1;
            dis_dd(sourceDot,dot_minDis) = inf;
            dis_dd(dot_minDis,sourceDot) = inf;
        end
    end
    
    maxLabel = 0;
    labelHist = zeros(l_dots,1);
    labelImg = zeros(m,n);
    while   1
        headDot = find( connection(l_dots+1,:) > 0, 1 );
        if  isempty( headDot )
            break;
        end
        
        maxLabel = maxLabel + 1;
        connection(l_dots+1,headDot) = 0;
        connection(headDot,l_dots+1) = 0;
        startD = headDot;
        endD = find( connection(headDot,:) > 0, 1 );
        while  endD < l_dots + 1
            if  old_dis(startD, endD) == inf
                tempD = startD;
                startD = endD;
                endD = tempD;
            end
            startP = sub2idx_ap( dots{startD}(1), dots{startD}(2) );
            x = dots{endD}(1);
            y = dots{endD}(2);
            idx = sub2idx_ap( x, y );  
            while   idx ~= startP
                ind = ind_ap( idx );
                [x, y] = ind2sub( [m,n], ind );
                labelImg( x, y ) = maxLabel;
                labelHist(maxLabel) = labelHist(maxLabel) + 1;
                idx = par_d2p( startD, idx );
                assert( idx > 0 );
            end
            connection(startD,endD) = 0;
            connection(endD,startD) = 0;
            startD = endD;
            endD = find( connection(endD,:) > 0 );
        end
    end
    
    useable = find(labelHist > 20 );
    temp = labelImg;
    labelImg = zeros(m,n);
    for     i = 1: 1: length(useable)
        labelImg( temp == useable(i) ) = i;
    end
end

function [DX,DY,DEvi] = get8diff( evidences )
    [X,Y] = meshgrid(-1:1);
    nei8 = find( X~=0 | Y~=0 );
    DX = X(nei8);
    DY = Y(nei8);
    
    [m,n] = size( evidences );
    assert( m > 10 && n > 10 );
    
    temp = zeros(m,n);
    temp(1,:) = inf;
    temp(2:m,:) = evidences(1:m-1,:);
    diff_t = abs(temp - evidences);
    
    temp(1:m-1,:) = evidences(2:m,:);
    temp(m,:) = inf;
    diff_b = abs( temp - evidences );
    
    temp(:,1) = inf;
    temp(:,2:n) = evidences(:,1:n-1);
    diff_l = abs( temp - evidences );
    
    temp(:,1:n-1) = evidences(:,2:n);
    temp(:,n) = inf;
    diff_r = abs( temp - evidences );
    
    temp(1,:) = inf;
    temp(:,1) = inf;
    temp( 2:m, 2:n ) = evidences( 1:m-1, 1:n-1 );
    diff_lt = abs( temp - evidences ) * sqrt(2);
    
    temp(m,:) = inf;
    temp(:,1) = inf;
    temp(1:m-1,2:n) = evidences(2:m,1:n-1);
    diff_lb = abs( temp - evidences ) * sqrt(2);
    
    temp(1,:) = inf;
    temp(:,n) = inf;
    temp(2:m,1:n-1) = evidences(1:m-1,2:n);
    diff_rt = abs( temp - evidences ) * sqrt(2);
    
    temp(m,:) = inf;
    temp(:,n) = inf;
    temp(1:m-1,2:n) = evidences(2:m,1:n-1);
    diff_rb = abs( temp - evidences ) * sqrt(2);
    
    DEvi = zeros(m,n,8);
    DEvi(:,:,1) = diff_lt;
    DEvi(:,:,2) = diff_l;
    DEvi(:,:,3) = diff_lb;
    DEvi(:,:,4) = diff_t;
    DEvi(:,:,5) = diff_b;
    DEvi(:,:,6) = diff_rt;
    DEvi(:,:,7) = diff_r;
    DEvi(:,:,8) = diff_rb;
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