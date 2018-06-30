function outA = Zstretch(inA,iter)
% Zstretch Performs the streching previously incorporated in EandR.m after
% the fact

A = inA;

%% Convert to vector
nrow = length(A(:,1));
ncol = length(A(1,:));

extA = zeros(ncol*nrow,1);
extA(1:nrow) = A(:,1);
for i=1:ncol-1
    extA(i*nrow+1:(i+1)*nrow) = A(:,i+1);
end

%% Determine the lower and upper limits
firstMed = median(extA); 

lowLimit = findMed(extA,-1,firstMed,iter);
highLimit = findMed(extA,1,firstMed,iter);

%% Exclude out of range pixels
for i=1:nrow
    for j=1:ncol
        if(A(i,j) > highLimit)
            A(i,j) = highLimit;
        end
        if(A(i,j) < lowLimit)
            A(i,j) = lowLimit;
        end
    end
end

%% Perform Stretch
A = A-min(min(A));
A = A./max(max(A)).*256;

outA = A;
