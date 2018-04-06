function    [pts] = findline_pts(lines,img)
len_x = size(img,2);
len_y = size(img,1);
pts= [];
count=1;
for k=1:size(lines,2)
    len = norm(lines(k).point1 - lines(k).point2);
    if(len>100)
       x1 = lines(k).point1(1);
       y1 = lines(k).point1(2);
       x2 = lines(k).point2(1);
       y2 = lines(k).point2(2); 
       if(x1<x2)
           start_x = x1;
           end_x = x2;
       else
           start_x = x2;
           end_x = x1;
       end
       for x = start_x:end_x
           y =floor(y2 + (x-x2)*(y2-y1)/(x2-x1));
           if(y<=len_y-5 && x<=len_x-5 && y>5 && x>5)
               [A,B] = meshgrid(x-5:x+5,y-5:y+5);
                c = cat(2,A',B');
                pts  = [pts;reshape(c,[],2)]; 
           end
       end 
    end
end
end

