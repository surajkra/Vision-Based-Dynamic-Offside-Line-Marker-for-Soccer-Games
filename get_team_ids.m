function [Team_IDs] = get_team_ids(Final,S,A)
    Team_IDs=zeros(size(S,1),1);
    for i =1:size(S,1)
        BB = S(i).BoundingBox;
        if(A(i).Area < 4000 && A(i).Area > 650 && BB(3) < 110 && BB(4) < 110 && BB(3) > 20 && BB(4) > 50 && BB(2)<950 )
            
            Reqd_Img = Final(BB(2):floor(BB(2)+BB(4)/2),BB(1):floor(BB(1)+BB(3)),:);
            indg = find(fuzzycolor(im2double(Reqd_Img),'red')<0.2);
            n = size(Reqd_Img,1)*size(Reqd_Img,2);
            Reqd_Img([indg;indg+n;indg+2*n]) = 0;
            
            Reqd_Img2 = Final(BB(2):floor(BB(2)+BB(4)/2),BB(1):floor(BB(1)+BB(3)),:);
            indg2 = find(fuzzycolor(im2double(Reqd_Img2),'blue')<0.2);
            n = size(Reqd_Img2,1)*size(Reqd_Img2,2);
            Reqd_Img2([indg2;indg2+n;indg2+2*n]) = 0;
            
            
            Reqd_Img3 = Final(BB(2):floor(BB(2)+BB(4)/2),BB(1):floor(BB(1)+BB(3)),:);
            indg3 = find(fuzzycolor(im2double(Reqd_Img3),'green')<0.5);
            n = size(Reqd_Img3,1)*size(Reqd_Img3,2);
            Reqd_Img3([indg3;indg3+n;indg3+2*n]) = 0;
            
            if(nnz(Reqd_Img)>100)
                Team_IDs(i) = 1;
            end
            
            if(nnz(Reqd_Img2)>10)
                Team_IDs(i) = 2;
            end
            
            if(nnz(Reqd_Img3)>50 && nnz(Reqd_Img)<50 && nnz(Reqd_Img2)<20 && BB(1)<700)
                Team_IDs(i) = 3;
            end
            
        else
            Team_IDs(i)=0;
        end
    end
end

