function [Team_IDs] = get_team_ids(Final,S,A)
    for i =1:size(S,1)
%         Team_IDs=zeros(size(S,1),1);
        BB = S(i).BoundingBox;
        if(A(i).Area < 4000 && A(i).Area > 650 && BB(3) < 110 && BB(4) < 110 )
            
            Reqd_Img = Final(BB(2):floor(BB(2)+BB(4)/2),BB(1):floor(BB(1)+BB(3)),:);
            indg = find(fuzzycolor(im2double(Reqd_Img),'neutral')<0.75);
            n = size(Reqd_Img,1)*size(Reqd_Img,2);
            Reqd_Img([indg;indg+n;indg+2*n]) = 0;
            if(nnz(Reqd_Img)>150)
                Team_IDs(i) = 1;
            else
                Team_IDs(i) = 2;
            end
        else
            Team_IDs(i)=0;
        end
    end
end

