clear all
close all
clc
%% Reading the Video
obj = VideoReader('match.mp4');
v = VideoWriter('v_1','MPEG-4');
open(v);
for frame_index = 1 : 600 %fill in the appropriate number
    disp(frame_index)
    if(exist('img'))
        prev_img = img;
    end
    img = readFrame(obj);
    
    if(frame_index<39)
        continue
    end
    

    if(frame_index==39)
%         imshow(img)
%         [x,y] = getpts;
%         points = [x,y];
%         close all
%         m = zeros(size(points,1)/2,1);
%         c = zeros(size(points,1)/2,1);
%         k = 1;
%         vp = zeros(2,1);
%         for j = 1:2:size(points,1)
%             m(k) = (points(j+1,2) - points(j,2))/ (points(j+1,1) - points(j,1));
%             c(k) = - points(j,1) * m(k) + points(j,2);
%             k = k+1;
%         end
%         count = 0;
%         for p = 1 : size(points,1)/2
% 
%            for q = p+1 : size(points,1)/2
%                count = count + 1;
%                A = [-m(p),1;-m(q),1];
%                b = [c(p);c(q)];
%                vp = vp + inv(A)*b;
%            end
%         end
%         vp = int16(vp / count);
%         disp(vp)
       vp = [935;-1115];
        continue
    end
    
    if(frame_index>40 && mod(frame_index,20)~=0)
%         f = figure('visible','off');    
%         figure();
%         imshow(img)
%         hold on;
%         for i =1:size(S,1)
%             BB = S(i).BoundingBox;
%              if(A(i).Area < 4000 && A(i).Area > 650 && BB(3) < 110 && BB(4) < 110 && BB(3) > 20 && BB(4) > 20 && Team_Ids(i)~=0  && BB(1)>=temp_max)
%                 rectangle('Position',[BB(1),BB(2),BB(3),BB(4)],...
%                 'LineWidth',2,'EdgeColor','red')
%                 if(Team_Ids(i)==1)
%                     text(BB(1)-2, BB(2)-2,'A');
%                 end
%                 if(Team_Ids(i)==2)
%                     text(BB(1)-2, BB(2)-2,'D');
%                 end
% 
%                 if(Team_Ids(i)==3)
%                     text(BB(1)-2, BB(2)-2,'GK');
%                 end
% 
%              end
%         end
%         plot([left_most,vp(1)],[ly ,vp(2)],'y','LineWidth',1)
    %f = figure('visible','off');
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
    
    
    
    BW_img = rgb2gray(img);
    Edge_img = edge(BW_img,'sobel');
%% Removing the TOP Boundary
    start_angle = 89;
    end_angle = 89.99;
    theta_resolution = 0.01;
    
    [hou,theta,rho] = hough(Edge_img(1:floor(size(Edge_img,1)/2),:), 'Theta', start_angle:theta_resolution:end_angle);
    peaks = houghpeaks(hou,2,'threshold',ceil(0.3*max(hou(:))));
    lines = houghlines(Edge_img(1:floor(size(Edge_img,1)/2),:),theta,rho,peaks,'FillGap',5,'MinLength',7);
%     [hou,theta,rho] = hough(Edge_img, 'Theta', start_angle:theta_resolution:end_angle);
%     peaks = houghpeaks(hou,2,'threshold',ceil(0.3*max(hou(:))));
%     lines = houghlines(Edge_img,theta,rho,peaks,'FillGap',5,'MinLength',7);
%     
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
   
%%  Pre-Processing the image (Dilation - Opening - Erosion)
%     S_E_D = strel('disk',5);
%     Dilate_img = imdilate(Edge_img,S_E_D);
%     S_E_E = strel('line',15,90);
%     Open_img = bwareaopen(Dilate_img,1000);
%     Erode_img = imerode(Open_img,S_E_E);
%% Removing the mid-line
%     start_angle = -45;
%     end_angle = 45;
%     theta_resolution = 0.01;
%     
%     [H,T,R] = hough(Erode_img, 'Theta', start_angle:theta_resolution:end_angle);
%     P  = houghpeaks(H,10,'threshold',ceil(0.3*max(H(:))));
%     lines = houghlines(Erode_img,T,R,P,'FillGap',5,'MinLength',7);
%     pts = findline_pts(lines,Erode_img);
%     for i = 1:size(pts,1)
%         Erode_img(pts(i,2),pts(i,1))=0;
%     end
%     
%     Conn_Comp = bwconncomp(Erode_img,8);
%     S = regionprops(Conn_Comp,'BoundingBox');  
%     A = regionprops(Conn_Comp,'Area');
%     Erode_rgb(:,:,1)=Erode_img;
%     Erode_rgb(:,:,2)=Erode_img;
%     Erode_rgb(:,:,3)=Erode_img;
%     Final = uint8(Erode_img).*img;
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


    S = [S_team_red; S_team_blue];

    Team_Ids = [ones(size(S_team_red,1),1); 2*ones(size(S_team_blue,1),1)];

    
    
    


















%% Get the Attacker Seed
%     if(frame_index==1)
%         disp('Mark the attacker');
%         figure();imshow(img);
%         rect  = getrect;
%         Reqd_Img = Final(rect(2):rect(2)+rect(4)/2,rect(1):rect(1)+rect(3),:);
%         Attack_ID = mean(Reqd_Img(find(Reqd_Img>0)));
%     end
%     if(size(S,1)==0)
%         continue
%     end
%%  Mark the Team-Players
%     Team_Ids = get_team_ids(Final,S,A);
%% Mark the bounding boxes
    %f = figure('visible','off');    
    figure();
    imshow(img)
    hold on;
    left_most = 9999;
    for i =1:size(S,1)
        BB = S(i).BoundingBox;
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
    plot([left_most,vp(1)],[ly ,vp(2)],'c','LineWidth',1)
    fig = getframe(gcf);
    writeVideo(v,fig);
    close(gcf);
end
close(v);
