function [outvec] = vecize(inMat)
% "vectorizes" the matrix Mat by by placing the columns tip to tail

m = length(inMat(:,1));
n = length(inMat(1,:));

for i = 1:n
    k = (i-1)*m+1;
    outvec(k:k+m-1) = inMat(:,i);
end