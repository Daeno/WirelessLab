function [xn,yn]=findrange(xn,yn,range)
for i=1:length(yn),
   if xn(i)<=range(1,1), xn(i)=range(1,1); end
   if xn(i)>=range(1,2), xn(i)=range(1,2); end
   if yn(i)<=range(2,1), yn(i)=range(2,1); end
   if yn(i)>=range(2,2), yn(i)=range(2,2); end
end

end