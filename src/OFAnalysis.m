% Parameters

margin = 20;    % min cut distance from image frame
alpha = 80      % optical flow parameters
ofs = 5         % offset for visualisation of quivers
qscale = 4;
delta = 3       % time resolution - compute flow every delta frame
w = 12          % window size for basic Lukas Kanade
codeBrox =1     % select the code in Brox if ~=0
totT = 30       % analyse first tot timepoinmts at most
%% Data 

% folder with input data
folder = "..\Data"

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
        
        clear DA DetA BA noiseA waveA

        screen_dim=get(0, 'MonitorPositions');
        mainGUI_pos=[10 10 1920/2 1200/2];

        clear n
        uv = [];
        close all
   
        filename = fullfile(resF, strcat(strtok(fn(k).name, '.'), '.csv')); 
        fileID = fopen(filename,'wt');
        fprintf(fileID, '%s\n', 'File, Time, mean_all_disp_left, mean_all_disp_right, mean_bright_disp_left,mean_bright_disp_right');
        for z = delta:delta:min(totT,size(A,3)-delta )%size(A,3)-delta 
            z
            close all
            tic
            %uv{z} = estimate_flow_hs(A(:,:,z),A(:,:,z+delta));
            aa = cat(3, sum(A(:,:,z-delta+1:z),3), sum(A(:,:,z+1:z+delta),3));
            % select optical flow method
            if codeBrox                 
                [D magnitude] = OFBrox(aa, 0.25, alpha,5);
                %print(gca, fullfile(resF, strcat(strtok(fn(k).name, '.'), '_Brox',sprintf('%03d', z+1), '.png')),'-dpng');
                uu = D{1}(:,:,1);
                vv = D{1}(:,:,2);
            else
                [uu vv] = LucasKanade(aa, w);
            end
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

