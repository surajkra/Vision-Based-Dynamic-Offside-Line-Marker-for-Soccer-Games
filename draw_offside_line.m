list = dir('I:\CVGA videos\*.mp4');
video = VideoReader(strcat('I:\CVGA videos\',list(9).name));
for i = 1:500
a = video.readFrame();
end
imshow(a)
%select lines and press enter
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
imshow(a)
hold on 
[x1,y1] = getpts;
plot([x1,vp(1)],[y1 ,vp(2)],'y','LineWidth',3)

