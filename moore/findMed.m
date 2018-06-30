function [backMed] = findMed(extA,updn,forwardMed,iter)

n=0;

if(updn > 0) % High Limit
    for i=1:length(extA)
        if(extA(i) > forwardMed)
            n = n+1;
            newA(n) = extA(i);
        end
    end
else         % low limit
    for i=1:length(extA)
        if(extA(i) < forwardMed)
            n = n+1;
            newA(n) = extA(i);
        end
    end
end

newMed = median(newA);

if(iter < 1)
    backMed = newMed;
else
    backMed = findMed(newA,updn,newMed,iter-1);
end
