clear all
close all
clc
%% Reading the Video
obj = VideoReader('cropped_match_1.mp4');
v = VideoWriter('v_1','MPEG-4');
open(v);
for frame_index = 1 : 500 %fill in the appropriate number
    disp(frame_index)
    if(exist('img'))
        prev_img = img;
    end
    img = readFrame(obj);    
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
   
%%  Pre-Processing the image (Dilation - Opening - Erosion)
    
    if(frame_index~=1 && mod(frame_index,10)~=0)
        f = figure('visible','off');    
        imshow(img)
        for i = 1:size(S,1) 
            if(S(i).BoundingBox(3) < 110 && S(i).BoundingBox(4) < 110 && S(i).BoundingBox(3) > 20 && S(i).BoundingBox(4) > 20 )
            BB = S(i).BoundingBox;
            if(S(i).BoundingBox(1)+BB(3)>size(img,2))
                S(i).BoundingBox(3) = size(img,2)-S(i).BoundingBox(1);
                BB(3) =S(i).BoundingBox(3);
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
            end
        end
        fig = getframe(gcf);
        writeVideo(v,fig);
        close(gcf)
    else
        S_E_D = strel('disk',5);
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
        if(frame_index==1)
            disp('Mark the attacker');
            figure();imshow(img);
            rect  = getrect;
            Reqd_Img = Final(rect(2):rect(2)+rect(4)/2,rect(1):rect(1)+rect(3),:);
            Attack_ID = mean(Reqd_Img(find(Reqd_Img>0)));
        end
        if(size(S,1)==0)
            continue
        end
    %%  Mark the Team-Players
        Team_Ids = get_team_ids(Final,S,A);
    %% Mark the bounding boxes
        f = figure('visible','off');    
    %     figure();
        imshow(img)
        hold on;
        for i =1:size(S,1)
            BB = S(i).BoundingBox;
             if(A(i).Area < 4000 && A(i).Area > 650 && BB(3) < 110 && BB(4) < 110 && BB(3) > 20 && BB(4) > 20 )
                rectangle('Position',[BB(1),BB(2),BB(3),BB(4)],...
                'LineWidth',2,'EdgeColor','red')
                if(Team_Ids(i)==1)
                    text(BB(1)-2, BB(2)-2,'A');
                else
                    text(BB(1)-2, BB(2)-2,'D');
                end
             end
        end
        fig = getframe(gcf);
        writeVideo(v,fig);
        close(gcf);
   end
end
close(v);
