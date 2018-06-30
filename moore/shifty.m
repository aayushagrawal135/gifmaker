%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% PHY 495S ASSIGNMENT #2 %
% %
% The Elastic Thickness of the Earth?s Crust %
% %
% by: John Moores (c) February - April 2003 %
% Supervisor: R. Bailey %
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSFORMED MATRIX TREATMENT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ M ] = shifty(MAP,nx,nxq,nyq)
%% Function to reorganize MAP into M by concatenating
%% to a matrix of size nxq by nyq since we require only
%% information from one of quadrants (1,3) and one of (2,4)
%% we select quadrants 1 and 4. The elements of 4 are reversed
%% and added to the elements of quandrant 1 (for averaging purposes)

M = zeros(nyq,nxq);

M(:,1) = MAP(1:nyq,1).*2;

for j = 2:nxq
    if( nx-(j-2) > nxq )
        M(:,j) = MAP(1:nyq,j)+MAP(1:nyq,nx-(j-2));
    else
    M(:,j) = MAP(1:nyq,j).*2;
    end
end 