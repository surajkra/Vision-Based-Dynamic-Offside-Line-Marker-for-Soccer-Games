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
    
    if(frame_index<439)
        continue
    end
    

    if(frame_index==439)
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
        continue
    end
    
    if(frame_index~=601 && mod(frame_index,10)~=0)
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
    f = figure('visible','off');    
        imshow(img)
        for i = 1:size(S,1)
            BB = S(i).BoundingBox;
            if(A(i).Area < 4000 && A(i).Area > 650 && BB(3) < 110 && BB(4) < 110 && BB(3) > 20 && BB(4) > 20 && Team_Ids(i)~=0  && BB(1)>=temp_max )
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
                text(BB(1)-2, BB(2)-2,'A');
                end
                if(Team_Ids(i)==2)
                    text(BB(1)-2, BB(2)-2,'D');
                end

                if(Team_Ids(i)==3)
                    text(BB(1)-2, BB(2)-2,'GK');
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
    S_E_D = strel('rectangle',[16,9]);
    Dilate_img = imdilate(Edge_img,S_E_D);
    S_E_E = strel('line',15,90);
    Open_img = bwareaopen(Dilate_img,1000);
    Erode_img = imerode(Open_img,S_E_E);
%% Removing the mid-line
    start_angle = -45;
    end_angle = 45;
    theta_resolution = 0.01;
    
    [H,T,R] = hough(Erode_img, 'Theta', start_angle:theta_resolution:end_angle);
    P  = houghpeaks(H,10,'threshold',ceil(0.3*max(H(:))));
    lines = houghlines(Erode_img,T,R,P,'FillGap',5,'MinLength',7);
    pts = findline_pts(lines,Erode_img);
    for i = 1:size(pts,1)
        Erode_img(pts(i,2),pts(i,1))=0;
    end
    Conn_Comp = bwconncomp(Erode_img,8);
    S = regionprops(Conn_Comp,'BoundingBox');  
    A = regionprops(Conn_Comp,'Area');
    Erode_rgb(:,:,1)=Erode_img;
    Erode_rgb(:,:,2)=Erode_img;
    Erode_rgb(:,:,3)=Erode_img;
    Final = uint8(Erode_img).*img;

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
    Team_Ids = get_team_ids(Final,S,A);
%% Mark the bounding boxes
    f = figure('visible','off');    
    figure();
    imshow(img)
    hold on;
    temp_indices = find(Team_Ids == 3);
    temp_BB = S(temp_indices);
    temp_len = size(temp_BB);
    temp_max = 0;
    for t = 1:temp_len
        if(temp_BB(t).BoundingBox(1)>temp_max)
            temp_max = temp_BB(t).BoundingBox(1);
        end
    end
    left_most = 9999;
    for i =1:size(S,1)
        BB = S(i).BoundingBox;
         if(A(i).Area < 4000 && A(i).Area > 650 && BB(3) < 110 && BB(4) < 110 && BB(3) > 20 && BB(4) > 20 && Team_Ids(i)~=0  && BB(1)>=temp_max)
            rectangle('Position',[BB(1),BB(2),BB(3),BB(4)],...
            'LineWidth',2,'EdgeColor','red')
            if(Team_Ids(i)==1)
                text(BB(1)-2, BB(2)-2,'A');
            end
            if(Team_Ids(i)==2)
                text(BB(1)-2, BB(2)-2,'D');
            end
            
            if(Team_Ids(i)==3)
                text(BB(1)-2, BB(2)-2,'GK');
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
end
close(v);
