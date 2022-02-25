function [tp,ind]=PickTime(data,t,dt)



for j=1:length(data(1,:))
    
    index=0;
    
    
    for i=2:(length(data(:,1))-1)
      
            k=((data(i+1,j)-2*data(i,j)+data(i-1,j))/(dt^2));
            if(k>2 || k<-2);
                index=i;
      
            break
        end
    end
    ind(1,j)=index;
    tp(1,j)=t(index);
end
    
end 


