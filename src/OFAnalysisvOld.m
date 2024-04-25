% Parameters
margin = 20; % min cut distance from image frame

alpha = 80  % optical flow parameters
ofs = 5
qscale = 4;
delta = 10 % time resolution - compute flow every delta frame
w = 12
%% Data 
% folder with input data
folder = ".\Data"
% result folder
resF = fullfile(folder,'Results')
mkdir(resF)

% pattern of files to be analysed
pattern = 'DV*'; % select files that match the pattern

% files
fn = dir(fullfile(folder, pattern));

for k = 3:length(fn)
    if (isempty(strfind(fn(k).name, 'PMT')))
        
    % read image with cut position
    cutIm = imread(fullfile(fn(k).folder,fn(k).name));
    % read data
    ch1V = Read3d(fullfile(fn(k).folder,strcat(strtok(fn(k).name, '.'), 'PMT - PMT [560-] _C1.ome.tif')));
    
    % Pre-processing
    
    % remove saturated spots
    ch1V = min(ch1V, quantile(ch1V(:),0.99));
    ch1Im = ch1V(:,:,1); 
    % destripe
    im_filtered = Destripe3d(ch1V);
    % denoise
    for z = 1:size(ch1V,3)
        A(:,:,z)=medfilt2(im_filtered(:,:,z));
    end

    % find cut
    mm = imbinarize(medfilt2(cutIm(:,:,2)-cutIm(:,:,1),[25 25]));
    stat = regionprops(mm,'Area','BoundingBox');
    [v pos] = max([stat.Area]);
    partcutIm = cutIm(floor(stat(pos).BoundingBox(2):stat(pos).BoundingBox(2)+stat(pos).BoundingBox(4)-1),floor(stat(pos).BoundingBox(1):stat(pos).BoundingBox(1)+stat(pos).BoundingBox(3)-1),2);
    cut = (partcutIm==0)+ (partcutIm==max(max(partcutIm==0)));
    partcutImR = imresize(partcutIm, size(ch1Im));
    cut(1:margin,:)=0; cut(end-margin:end,:)=0; cut(:,1:margin)=0; cut(:,end-margin:end)=0;
    cutR = imresize(cut, size(ch1Im));
    figure; imshowpair(cutR>0,ch1Im)

    % create mask from cut
    figure;
    [H,T,R] = hough(cutR);
    P  = houghpeaks(H,1,'threshold',ceil(0.3*max(H(:))));
    x = T(P(:,2)); y = R(P(:,1));

    % Find lines and plot them
    lines = houghlines(cutR,T,R,P,'FillGap',35,'MinLength',25);
    figure, 
    %imshow(cut), hold on
    imshow(sum(A,3),[]); hold on
    max_len = 0;
    for ll = 1:length(lines)
       xy = [lines(ll).point1; lines(ll).point2];
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

       % Determine the endpoints of the longest line segment
       len = norm(lines(ll).point1 - lines(ll).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end

    % highlight the longest line segment
    plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');
    mask0 =  0*cutR;
    mask0(xy(1,2),xy(1,1)) =1;
    mask0(xy(2,2),xy(2,1)) =1;
    mask0 = imdilate(mask0, strel('disk',7));
    
    mask1 = 0*cutR;
    id = LineFind(xy_long(1,1), xy_long(1,2), xy_long(2,1), xy_long(2,2)); 
    mask1(sub2ind(size(cutR),id(:,2),id(:,1))) = 1;
    
    mask1 = imdilate(bwmorph(mask1,'spur',5),strel('disk',5));
    mask2 = 0*cutR;
    mask2(sub2ind(size(cutR),id(:,2),id(:,1))) = 1;
    %mask2 = bwmorph(bwmorph(mask2,'spur',20),'dilate',20);
    mask2 = imdilate(bwmorph(mask2,'spur',30),strel('disk',30));
    
    % compute orientation of the cit
    Or = regionprops(mask1, 'Orientation');

    % Since 0,0 i top ;egft in the image it seems the code below is
    % necessary to compute normal;
    Or.Orientation = -Or.Orientation
    if Or.Orientation>0
        gamma = deg2rad(Or.Orientation-90);
    else
        gamma = deg2rad(Or.Orientation+90);
    end
    % create ROI label image with the two sides of the cut
    Lb = bwlabel(mask2 > (mask1+mask0));
    figure; imagesc(Lb)
    
    maskstat = regionprops(Lb, 'Centroid', 'PixelIdxList');
    if max(Lb(:))==2
        % LB==1 for leftmost x centroid
        if maskstat(1).Centroid(1)>maskstat(2).Centroid(1)
            auxLb = Lb;
            auxLb(Lb==1) = 2;
            auxLb(Lb==2) = 1;
            maskstat = regionprops(auxLb, 'Centroid', 'PixelIdxList');
        end
        
        [bp1y bp1x] = find(bwperim(mask1)>0);
        [bp2y bp2x] = find(bwperim(mask2)>0);

        w =12;
        delta = 3
        clear DA DetA BA noiseA waveA

        screen_dim=get(0, 'MonitorPositions');
        mainGUI_pos=[10 10 1920/2 1200/2];

        clear n
        uv = [];
        close all
   
        filename = fullfile(resF, strcat(strtok(fn(k).name, '.'), '.csv')); 
        fileID = fopen(filename,'wt');
        fprintf(fileID, '%s\n', 'File, Time, mean_all_disp_left, mean_all_disp_right, mean_bright_disp_left,mean_bright_disp_right');
        for z = delta:delta:min(30,size(A,3)-delta )%size(A,3)-delta 
            z
            close all
            tic
            %uv{z} = estimate_flow_hs(A(:,:,z),A(:,:,z+delta));

            %[uu vv] = LucasKanade(sum(A(:,:,z-delta+1:z),3),sum(A(:,:,z+1:z+delta),3), w);
            aa = cat(3, sum(A(:,:,z-delta+1:z),3), sum(A(:,:,z+1:z+delta),3));
            [D magnitude] = OFBrox(aa, 0.25, alpha,5);
            %print(gca, fullfile(resF, strcat(strtok(fn(k).name, '.'), '_Brox',sprintf('%03d', z+1), '.png')),'-dpng');
            uu = D{1}(:,:,1);
            vv = D{1}(:,:,2);
            
            % compute the flow orthogonal to the cut
            Vperp = uu(:)*cos(gamma)+vv(:)*sin(gamma);

            all_left = mean(Vperp(maskstat(1).PixelIdxList));
            all_right = mean(Vperp(maskstat(2).PixelIdxList));

            [DA(:,:,z) DetA(:,:,z) BA noiseA waveA] = FindPeakWav(sum(A(:,:,z+1:z+delta),3), 0.05, 0, 3, 2, 1); 
            DetA(:,:,z) = bwareaopen(DetA(:,:,z),5);
            DetA(:,:,z) = imclose(DetA(:,:,z), strel('disk',3));
            aux = DetA(:,:,z);  

            id = find(aux(maskstat(1).PixelIdxList)>0);
            sel_left = mean(Vperp(maskstat(1).PixelIdxList(id)));
            id = find(aux(maskstat(2).PixelIdxList)>0);
            sel_right = mean(Vperp(maskstat(2).PixelIdxList(id)));

            h= figure('WindowState','maximized');%'units','normalized','outerposition',[0 0 1 1]);
            reset(h)
            clf
            imshowpair(aa(:,:,1),aa(:,:,2)); hold on; axis equal tight

            [xx yy] = meshgrid(1:ofs:size(A,2),1:ofs:size(A,1));
            u2 = interp2(reshape(Vperp*cos(gamma), size(aux)), xx, yy);
            v2 = interp2(reshape(Vperp*sin(gamma), size(aux)), xx, yy);
            h1 = quiver(xx,yy,u2,v2, 0, 'm', 'AutoScale', 'on','AutoScaleFactor', 3);
            hold on
            hU = get(h1,'UData') ;
            hV = get(h1,'VData') ;
            set(h1,'UData',qscale*hU,'VData',qscale*hV)

            plot(xy_long(:,1),xy_long(:,2),'LineWidth',0.75, 'Color', 'c', 'LineStyle', '--');%'Color',[163,193,173]/255); % DarkTeal
            plot(bp2x, bp2y, 'c.','MarkerSize',6);%, 'Color', clr(2,:));%, 'Marker','-')
            print(h, fullfile(resF, strcat(strtok(fn(k).name, '.'), '_Proj',sprintf('%03d', z+1), '.png')),'-dpng');
            fprintf(fileID, '%s, %d, %f, %f, %f, %f \n',strtok(fn(k).name, '.'), z+1, all_left,all_right,sel_left,sel_right);    
            pause(0.3)
        end
        fclose(fileID);
    end
    end
end


%%
summary = [];
summaryMeans = [];
mkdir(resF);
%imwrite(cut, strcat(resF,filesep,'cut.tif'), 'tif');

sA = A;
B = A;


clear DA DetA BA noiseA waveA

screen_dim=get(0, 'MonitorPositions');
mainGUI_pos=[10 screen_dim(4)-200 650 160];

clear n
uv = [];
close all
h= figure('pos',mainGUI_pos); clf
im1 = sum(B(:,:,1:1+delta-1),3);
im1 = im1/max(im1(:));
imagesc(im1); hold on; axis equal tight
savedc = caxis();
set(gcf, 'Position', get(0,'Screensize'));

for z = 1+delta-1:delta:min(24,size(A,3)-delta )%size(A,3)-delta 
    clf
    z
    close all
    h= figure('units','normalized','outerposition',[0 0 1 1]);
    %set(gcf, 'Position', get(0,'Screensize'));
    
    im1 = sum(B(:,:,z+1:z+delta),3);
    im1 = im1/max(im1(:));
       
    %imagesc(sum(B(:,:,z+1:z+delta),3)); hold on; axis equal tight
    imagesc(imadjust(im1, [0.35 0.75], [0 0.9])); hold on; axis equal tight
    
    %imagesc(sum(sA(:,:,z-delta+1:z),3)); hold on; axis equal tight
    colormap gray
    caxis (savedc)
    tic
    %uv{z} = estimate_flow_hs(A(:,:,z),A(:,:,z+delta));
 
    [u v] = LucasKanade(sum(sA(:,:,z-delta+1:z),3),sum(sA(:,:,z+1:z+delta),3), w);
    aa = cat(3, sum(sA(:,:,z-delta+1:z),3), sum(sA(:,:,z+1:z+delta),3));
    [D magnitude] = OFBrox(aa, 0.25, alpha,5);
    uu = D{1}(:,:,1);
    vv = D{1}(:,:,2);
    Dperp = sqrt(uu(:).^2 + vv(:).^2)*sin(deg2rad(Or.Orientation+90));
    gamma = deg2rad(Or.Orientation+90);
    Vperp = uu(:)*cos(gamma)+vv(:)*sin(gamma);
    
    % Projection on the perpeniduclar direction to the cut 
    vcut =  [xy_long(2,1)-xy_long(1,1), xy_long(2,2)-xy_long(1,2) ];
    vpcut = [-vcut(2) vcut(1)];
    
    
    for p = 1:length(u(:))
        a =atan2(v(p),u(p)) - atan2(vcut(2), vcut(1));
        if a< -pi
            a = a+2*pi;
        end
        if a> pi
            a = a-2*pi;
        end  

        perp = norm([u(p) v(p)])* sin(a)*vpcut/norm(vpcut);
        u(p)  = perp(1);
        v(p)  = perp(2);
    end

    if (pos == z+1)
         u = u/(delta+1);
         v = v/(delta+1);
    else
        u = u/(delta);
        v = v/(delta);
    end
    uv{z} = cat(3,u,v);       

    toc
   
    % detection on original. 5% false positives
    [DA(:,:,z) DetA(:,:,z) BA noiseA waveA] = FindPeakWav(sum(B(:,:,z+1:z+delta),3), 0.05, 0, 3, 2, 1); 
 
    % sparse mesh - for visualisation
    ofs = 3
    [xx yy] = meshgrid(1:ofs:size(A,2),1:ofs:size(A,1));
    u2 = interp2(uv{z}(:, :, 1), xx, yy);
    v2 = interp2(uv{z}(:, :, 2), xx, yy);
    % dense mesh
    [x y] = meshgrid(1:size(A,2),1:size(A,1));
    uu = uv{z}(:, :, 1);
    vv = uv{z}(:, :, 2);

    n(:,:,z)=sqrt(uv{z}(:,:,1).^2+uv{z}(:,:,2).^2);
    
    % find signal
    DetA(:,:,z) = bwareaopen(DetA(:,:,z),5);
    DetA(:,:,z) = imclose(DetA(:,:,z), strel('disk',3));
    aux = DetA(:,:,z);   

    id = find(aux(sub2ind(size(aux), yy(:), xx(:)))>0);
    %h1 = quiver(xx(id),yy(id),u2(id),v2(id), 0, 'Color',[88,166,24]/255 );  % Core Orange
    %h1 = quiver(xx(id),yy(id),u2(id),v2(id), 0, 'Color',[138,226,52]/255 );% clr(5,:));%[20,195,24]/255 );  % former green - all arrows
    
    h = quiver(xx,yy,u2,v2, 0, 'Color',[138,226,52]/255 );% clr(5,:));%[20,195,24]/255 );  % former green - all arrows
     
    h1 = quiver(xx,yy,u2,v2, 0, 'Color',[138,226,52]/255 );% clr(5,:));%[20,195,24]/255 );  % former green - all arrows
    hU = get(h1,'UData') ;
    hV = get(h1,'VData') ;
    set(h1,'UData',qscale*hU,'VData',qscale*hV)

    if (pos == z+1) & withcyan
        u2b = interp2(uv{z-delta}(:, :, 1), xx, yy);
        v2b = interp2(uv{z-delta}(:, :, 2), xx, yy);        
        h2 = quiver(xx(id)+1,yy(id)+1,u2b(id),v2b(id), 0, 'c');
        hU = get(h2,'UData') ;
        hV = get(h2,'VData') ;
        set(h2,'UData',qscale*hU,'VData',qscale*hV)
    end
%     pause(0.01);
%     print(h, strcat(resF,filesep,'D',fn(1:end-9), sprintf('%03d', z+1), '.png'),'-dpng');
    if z >= pos-delta
        plot(xy_long(:,1),xy_long(:,2),'LineWidth',0.75, 'Color', 'w', 'LineStyle', '--');%'Color',[163,193,173]/255); % DarkTeal
        %plot(bp1x, bp1y, 'y.')
        %plot(bp2x, bp2y, 'w.', 'MarkerSize',0.8);%, 'Color', clr(2,:));%, 'Marker','-') % cut
        %plot(bp2x, bp2y, '.', 'Color', [0,62,114]/255, 'Marker','.')
    end
     plot(bp2x, bp2y, 'w.','MarkerSize',6);%, 'Color', clr(2,:));%, 'Marker','-')
    % DetA - signal. Intersect with mask2. Angle with xy_long(:,1),xy_long(:,2)
    newid = find((mask2 & DetA(:,:,z)) > 0); 

    
    for m = 1:length(newid)    
       newid(m) 
       if (norm([uu(newid(m)) vv(newid(m))])== 0)
            angles(m) = 0;
       else
            angles(m) =atan2(vv(newid(m)),uu(newid(m))) - atan2(vcut(2), vcut(1));
            if angles(m)< -pi
                angles(m) = angles(m)+2*pi;
            end
            if angles(m)> pi
                angles(m) = angles(m)-2*pi;
            end
            angles(m) = angles(m)*180/pi;
       end
       [y1 x1] = ind2sub(size(mask2), newid(m));
       

       % Pointing away from the cut. Line of equation: aa*x-y+cc=0
       if (xy_long(2,1)~=xy_long(1,1))
            aa = (xy_long(2,2)-xy_long(1,2))/(xy_long(2,1)-xy_long(1,1));           
       else
           aa = 0;  
       end
       bb = -1;
       cc =  - xy_long(2,1)*(xy_long(2,2)-xy_long(1,2))/(xy_long(2,1)-xy_long(1,1)) + xy_long(2,2);
       d1 = abs(dot([aa bb cc], [x1 y1 1]))/sqrt(aa*aa+bb*bb); 
       d2 = abs(dot([aa bb cc], [x1+uu(newid(m)) y1+vv(newid(m)) 1]))/sqrt(aa*aa+bb*bb);
       away = d1<d2;
       % column 6 is an indicator if above or under the line. 
       %>0 means under the line
       summary = [summary; z+1 y1 x1 angles(m) norm([uu(newid(m)) vv(newid(m))])...
           (y1 > x1*(xy_long(2,2)-xy_long(1,2))/(xy_long(2,1)-xy_long(1,1)) - xy_long(2,1)*(xy_long(2,2)-xy_long(1,2))/(xy_long(2,1)-xy_long(1,1)) + xy_long(2,2)) away]; 
    end
    if q<0
        id2 = find((summary(:,1)==z+1) &(abs(summary(:,4))>minangle) & (abs(summary(:,4))<maxangle)  &(summary(:,7)>0));
    else
        id2 = find((summary(:,1)==z+1) &(abs(summary(:,4))>minangle) & (abs(summary(:,4))<maxangle) &(summary(:,5)>q) &(summary(:,7)>0));
    end
   
   
    if (z+delta<pos )
         id4 = find((summary(:,1)==z+1)  & (summary(:,6)>0));
    else
        id4 = find((summary(:,1)==z+1)  & (summary(:,6)>0) );
    end
    mmx(1) = mean(uu(sub2ind(size(mask2), summary(id4,2),summary(id4,3))));
    mmy(1) = mean(vv(sub2ind(size(mask2), summary(id4,2),summary(id4,3))));
    mx(1) = mean(summary(id4,3));
    my(1) = mean(summary(id4,2));
    
        
    if (z+delta<pos )
         id4b = find((summary(:,1)==z+1)  & (summary(:,6)==0));
    else
        id4b = find((summary(:,1)==z+1)  & (summary(:,6)==0) );
    end
    
    mmx(2) = mean(uu(sub2ind(size(mask2), summary(id4b,2),summary(id4b,3))));
    mmy(2) = mean(vv(sub2ind(size(mask2), summary(id4b,2),summary(id4b,3))));
    mx(2) = mean(summary(id4b,3));
    my(2) = mean(summary(id4b,2));

    a = [0 0];
    if length(id4)>0
    if norm([mmx(1) mmy(1)])> 0
        a(1) =atan2(mmy(1),mmx(1)) - atan2(vcut(2), vcut(1));
        if angles(m)< -pi
            a(1) = a(1)+2*pi;
        end
        if a(1)> pi
           a(1) = a(1)-2*pi;
        end
        a(1) = a(1)*180/pi;
    end
    else
        a1 = 0; mx(1) = 0; my(1) = 0; mmx(1) = 0; mmy(1)= 0;
    end
    
    if length(id4b)>0
    if norm([mmx(2) mmy(2)])> 0
        a(2) =atan2(mmy(2),mmx(2)) - atan2(vcut(2), vcut(1));
        if angles(m)< -pi
            a(2) = a(2)+2*pi;
        end
        if a(2)> pi
           a(2) = a(2)-2*pi;
        end
        a(2) = a(2)*180/pi;
    end
    else
        a2 = 0; mx(2) = 0; my(2) = 0; mmx(2) = 0; mmy(2)= 0;
    end
    
    
    % Pointing away from the cut. Line of equation: aa*x-y+cc=0
   if (xy_long(2,1)~=xy_long(1,1))
        aa = (xy_long(2,2)-xy_long(1,2))/(xy_long(2,1)-xy_long(1,1));           
   else
       aa = 0;  
   end
   aa = (xy_long(2,2)-xy_long(1,2))/(xy_long(2,1)-xy_long(1,1));
   bb = -1;
   cc =  - xy_long(2,1)*(xy_long(2,2)-xy_long(1,2))/(xy_long(2,1)-xy_long(1,1)) + xy_long(2,2);
   d1 = abs(dot([aa bb cc], [mx(1) my(1) 1]))/sqrt(aa*aa+bb*bb); 
   d2 = abs(dot([aa bb cc], [mx(1)+mmx(1) my(1)+mmy(1) 1]))/sqrt(aa*aa+bb*bb);
   away1 = d1<d2;
   
   if (xy_long(2,1)~=xy_long(1,1))
        aa = (xy_long(2,2)-xy_long(1,2))/(xy_long(2,1)-xy_long(1,1));           
   else
       aa = 0;  
   end
   bb = -1;
   cc =  - xy_long(2,1)*(xy_long(2,2)-xy_long(1,2))/(xy_long(2,1)-xy_long(1,1)) + xy_long(2,2);
   d1 = abs(dot([aa bb cc], [mx(2) my(2) 1]))/sqrt(aa*aa+bb*bb); 
   d2 = abs(dot([aa bb cc], [mx(2)+mmx(2) my(2)+mmy(2) 1]))/sqrt(aa*aa+bb*bb);
   away2 = d1<d2;


    summaryMeans = [summaryMeans; z+1 my(1) mx(1) a(1) norm([mmx(1) mmy(1)])  mmx(1) mmy(1) length(id4) 1 ...
                    my(2) mx(2)  a(2) norm([mmx(2) mmy(2)]) mmx(2) mmy(2) length(id4b) 0 away1 away2]; 
    
    if ~isempty(find(isnan(summaryMeans(end,:))))  
        disp('stop')
    end
                
    %h3 = quiver(mx,my,mmx,mmy, 0, 'Color', [227,114,34]/255, 'LineWidth',1.5);  % Dark Purple
    pink = [199, 21, 133]/255 ; % dark (now) medium violet
    pink = [218 112 214]/255 ; % orchid
    pink = [1 0 1]; % magenta
    
    h3 = quiver(mx, my, mmx, mmy, 0, 'Color', pink, 'MaxHeadSize', 0.4, 'LineWidth',1.2);
    hU = get(h3,'UData') ;
    hV = get(h3,'VData') ;
    set(h3,'UData',qscale*hU,'VData',qscale*hV)

end

save( strcat(resF,filesep,'DataD',fn(1:end-4),'.mat'), 'uv', 'n');

% write header, name of the file
%csvwrite(strcat(resF,filesep,'Summary',fn(1:end-4),'.csv'),summary);

filename = strcat(resF,filesep,'Summary',fn(1:end-4),'.csv'); 
fileID = fopen(filename,'wt');
fprintf(fileID, '%s\n', fn(1:end-4));
fprintf(fileID, '%s\n', 'Time, y, x, Angle, D, BelowCut, Away');

for index = 1:length(summary)    
    
    fprintf(fileID, '%d, %d, %d, %f, %f,  %d, %d \n', summary(index,1),summary(index,2),summary(index,3),summary(index,4),summary(index,5),...
        summary(index,6),summary(index,7));    
end
fclose(fileID);
