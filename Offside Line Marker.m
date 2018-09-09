clear all
close all
clc
%% Reading the Video
obj = VideoReader('match.mp4');
v = VideoWriter('v_1','MPEG-4');
open(v);
start_frame = 39;

%% Detecting the Offside Line
for frame_index = 1 : 600 %fill in the appropriate number
    disp(frame_index)
    if(exist('img'))
        prev_img = img;
    end
    img = readFrame(obj);
    
    %Skip a few frames
    if(frame_index<start_frame)
        continue
    end
    
    %Calculate vanishing point. User input required here. Please mark the
    %lines in the GUI
    
    if(frame_index==start_frame)
       
        imshow(img)
        [x,y] = getpts;
        points = [x,y];
        close all
        m = zeros(size(points,1)/2,1);
        c = zeros(size(points,1)/2,1);
        k = 1;
        vp = zeros(2,1);
        for j = 1:2:size(points,1)
            m(k) = (points(j+1,2) - points(j,2))/ (points(j+1,1) - points(j,1));
            c(k) = - points(j,1) * m(k) + points(j,2);
            k = k+1;
        end
        count = 0;
        for p = 1 : size(points,1)/2

           for q = p+1 : size(points,1)/2
               count = count + 1;
               A = [-m(p),1;-m(q),1];
               b = [c(p);c(q)];
               vp = vp + inv(A)*b;
           end
        end
        vp = int16(vp / count);
        disp(vp)
%       vp = [935;-1115];
        continue
    end
    
    %Do tracking for 19/20 frames. Only 1 in every 20 frames we do the
    %actual detection. Rest 19 frames are tracked using KLT algorithm. The
    %below code is for the tracking.
    %Note here that for the first time detection runs and creates variable
    %S which contains all the bounding boxes and Team_Ids which contains
    %corresponding teams.
    if(frame_index>start_frame+1 && mod(frame_index,20)~=0)
    f = figure('visible','off');
        imshow(img)
        left_most = 9999;
        for i = 1:size(S,1)
            BB = S(i).BoundingBox;
            if(( BB(1)+(BB(3)/2)<115 || BB(1)+(BB(3)/2)>130) && (BB(2)+(BB(4)/2)<990 || BB(2)+(BB(4)/2)>1010))
                if(S(i).BoundingBox(1)<1)
                    S(i).BoundingBox(1) = 1;
                    BB(1) = S(i).BoundingBox(1);
                end
                if(S(i).BoundingBox(2)<1)
                    S(i).BoundingBox(2) = 1;
                    BB(2) = S(i).BoundingBox(2);
                end
                if(S(i).BoundingBox(1)+BB(3)>size(img,2))
                    S(i).BoundingBox(3) = size(img,2)-S(i).BoundingBox(1);
                    BB(3) = S(i).BoundingBox(3);
                end
                if(S(i).BoundingBox(2)+BB(4)>size(img,1))
                    S(i).BoundingBox(4) = size(img,1)-S(i).BoundingBox(2);
                    BB(4) = S(i).BoundingBox(4);
                end
                points = detectMinEigenFeatures(rgb2gray(prev_img),'ROI',S(i).BoundingBox);
                if(size(points,1) ==0)
                    disp('ERROR in points here')
                    continue
                end
                pointImage = insertMarker(prev_img,points.Location,'+','Color','white');
                tracker = vision.PointTracker('MaxBidirectionalError',1);
                initialize(tracker,points.Location,prev_img);         
                frame = img;
                [points, validity] = step(tracker,frame);
                mean_x = mean(points(:,1));
                mean_y = mean(points(:,2));
                S(i).BoundingBox(1) = floor(mean_x - BB(3)/2);
                S(i).BoundingBox(2) = floor(mean_y - BB(4)/2);
                S(i).BoundingBox(3) = BB(3);
                S(i).BoundingBox(4) = BB(4);
                img1 = insertMarker(frame,points(validity, :),'+');
                hold on;
                rectangle('Position',[S(i).BoundingBox(1),S(i).BoundingBox(2),S(i).BoundingBox(3),S(i).BoundingBox(4)],...
                        'LineWidth',2,'EdgeColor','red')
                if(Team_Ids(i)==1)
                text(BB(1)-2, BB(2)-2,'D_T');
                end
                if(Team_Ids(i)==2)
                    text(BB(1)-2, BB(2)-2,'A_T');
                end
                %Calculating the last defender on the left side using
                %vanishing point. Same can be done symmetrically to the
                %right hand side as well.
                
                x1 = floor(BB(1) + BB(3)/2);
                y1 = floor(BB(2) + BB(4));
                ly = size(img,1);
                slope = (vp(2) - y1)/(vp(1) - x1);
                y_int = - x1 * slope + y1;
                lx = (ly - y_int)/slope;
                if(lx<left_most && Team_Ids(i) == 1)
                 left_most = lx;
                end
            end
        end
        plot([left_most,vp(1)],[ly ,vp(2)],'y','LineWidth',1)
        fig = getframe(gcf);
        writeVideo(v,fig);
        close(gcf);
        continue;
    end
    
    
%% Actual Detection starts (one every 20 frames).
%Pre Processing the image to grayscale

    BW_img = rgb2gray(img);
    Edge_img = edge(BW_img,'sobel');
%% Removing the TOP Boundary
    start_angle = 89;
    end_angle = 89.99;
    theta_resolution = 0.01;
    
    [hou,theta,rho] = hough(Edge_img(1:floor(size(Edge_img,1)/2),:), 'Theta', start_angle:theta_resolution:end_angle);
    peaks = houghpeaks(hou,2,'threshold',ceil(0.3*max(hou(:))));
    lines = houghlines(Edge_img(1:floor(size(Edge_img,1)/2),:),theta,rho,peaks,'FillGap',5,'MinLength',7);    
    min_row = lines(1).point1(2);
    xy_long = [lines(1).point1; lines(1).point2];
    for k = 1:length(lines) 
        xy = [lines(k).point1; lines(k).point2];
        row_index = lines(k).point1(2);
        if (row_index < min_row)
            min_row = row_index; 
            xy_long = xy;
            index = k;
        end
    end
    img(1:xy_long(:,2),:,:)=0;
    BW_img(1:xy_long(:,2),:,:)=0;
    Edge_img(1:xy_long(:,2),:,:)=0;
   
%% Determining the actual play area
    
    %Find all dominant greens
    indg = find(fuzzycolor(im2double(img),'green')<0.1);
    n = size(img,1)*size(img,2);

    imggreen = img;
    imggreen([indg;indg+n;indg+2*n]) = 0;

    mask = imbinarize(rgb2gray(imggreen));
    mask = imfill(mask,'holes'); 
    mask_open = bwareaopen(mask,300);
    mask_open = imfill(mask_open,'holes'); 
    Conn_Comp_green = bwconncomp(mask_open,8);
    S_green = regionprops(Conn_Comp_green,'BoundingBox','Area');

    [~,max_ind_green] = max([S_green.Area]);
    bb_max_green = S_green(max_ind_green).BoundingBox;
    
    %Get a new `valid' image, which contatins only the actual play area.
    img_valid = img;
    max_h = size(img,1);
    if(bb_max_green(1)>1)
        for row = 1:max_h
            x_curr  = bb_max_green(1) + ((bb_max_green(1)-vp(1))/(max_h-vp(2))) * (row-max_h);
            x_curr = floor(x_curr);
            img_valid(row,1:x_curr,:) = 0;
        end
    else
        current_ind = 1;
        while(1)
            current_value = mask_open(current_ind,1);
            if(current_value ==1)
                break
            end
            current_ind = current_ind + 1;
        end
        for row = 1:current_ind
            x_curr  = bb_max_green(1) + ((bb_max_green(1) -vp(1))/(current_ind-vp(2))) * (row-current_ind);
            x_curr = floor(x_curr);
            img_valid(row,1:x_curr,:) = 0;
        end

    end
    
%% Determining the players and Team_Ids
    
    %Team Red, replace color with the color of the defending team.
    indg = find(fuzzycolor(im2double(img_valid),'red')<0.1);
    n = size(img,1)*size(img,2);

    img_team_read = img_valid;
    img_team_read([indg;indg+n;indg+2*n]) = 0;

    mask = imbinarize(rgb2gray(img_team_read));
    mask = imfill(mask,'holes'); 
    mask_open = bwareaopen(mask,20);
    mask_open = imfill(mask_open,'holes'); 
    S_E_D = strel('disk',15);
    mask_open = imdilate(mask_open,S_E_D);

    Conn_Comp_team_red = bwconncomp(mask_open,8);
    S_team_red = regionprops(Conn_Comp_team_red,'BoundingBox','Area');
    
    %Team Blue, replace color with the color of the defending team.
    indg = find(fuzzycolor(im2double(img_valid),'blue')<0.1);
    n = size(img,1)*size(img,2);

    img_team_read = img_valid;
    img_team_read([indg;indg+n;indg+2*n]) = 0;

    mask = imbinarize(rgb2gray(img_team_read));
    mask = imfill(mask,'holes'); 
    mask_open = bwareaopen(mask,20);
    mask_open = imfill(mask_open,'holes'); 
    S_E_D = strel('disk',15);
    mask_open = imdilate(mask_open,S_E_D);

    Conn_Comp_team_blue = bwconncomp(mask_open,8);
    S_team_blue = regionprops(Conn_Comp_team_blue,'BoundingBox','Area');

    %Getting all players/teamids in one list
    S = [S_team_red; S_team_blue];
    Team_Ids = [ones(size(S_team_red,1),1); 2*ones(size(S_team_blue,1),1)];

%% Mark the bounding boxes
    f = figure('visible','off');    
    figure();
    imshow(img)
    hold on;
    left_most = 9999;
    for i =1:size(S,1)
        BB = S(i).BoundingBox;
        %Accounting for static UI elements with team logos which disrupt
        %the detection. We know the locations of these static UI elements.
        %In real life, this need not be done as we can work with the raw
        %camera feed directly. Remove/change this portion to suit your needs.
        
         if(( BB(1)+(BB(3)/2)<115 || BB(1)+(BB(3)/2)>130) || (BB(2)+(BB(4)/2)<990 || BB(2)+(BB(4)/2)>1010)) 
            
            if(Team_Ids(i)==1)
                text(BB(1)-2, BB(2)-2,'D');
                BB(4)  = 1.5*BB(4);
                S(i).BoundingBox(4) = BB(4);
            end
            if(Team_Ids(i)==2)
                text(BB(1)-2, BB(2)-2,'A');
            end
            rectangle('Position',[BB(1),BB(2),BB(3),BB(4)],...
            'LineWidth',2,'EdgeColor','red')
        
            x1 = floor(BB(1) + BB(3)/2);
            y1 = floor(BB(2) + BB(4));
            ly = size(img,1);
            slope = (vp(2) - y1)/(vp(1) - x1);
            y_int = - x1 * slope + y1;
            lx = (ly - y_int)/slope;
            if(lx<left_most && Team_Ids(i) == 1)
             left_most = lx;
            end
            
         end
         
    end
    %Plot the offside line, this is currently done for the left most
    %player. Same can be repeated for the right most player as well.
    plot([left_most,vp(1)],[ly ,vp(2)],'c','LineWidth',1)
    fig = getframe(gcf);
    writeVideo(v,fig);
    close(gcf);
end
close(v);
