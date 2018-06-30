%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% PHY 495S ASSIGNMENT #2 %
% %
% The Elastic Thickness of the Earth?s Crust %
% %
% by: John Moores (c) February - April 2003 %
% Supervisor: R. Bailey %
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONCENTRIC AVERAGING %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ B,S ] = mAVG(T, bn,latlim)
% Calculates concentric circular averag es of matrices treated by shifty.m
% (i.e. assumes an overlapping 1st and 2nd quadrants). Furthermore assumes a geographic
% projection with the y-axis of the matrix running North-South and the x-axis running East West
% designed to divide the matrix along concentric circles centered at the first element
% into several bins averaged over bn matrix elements
%
% Modified to accept non-square matrices and to compensate the x-coordinates for the
% change in period at the average y lattitude

%clear all
%T = [ones(50,100); zeros(50,100) ]; %% sample input matrix (used only for test cases)
%bn = 1; %% averaging (bin) width ? elemental, but does not have to be an integral number

%a = 1/cos((latlim(2)+latlim(1))/2*pi/180); %% lattitude scale factor applied to x-domain
a = 1;

sx = length(T(1,:));
sy = length(T(:,1));
k = 0;
while( bn*k < (sqrt((sx*a)^2 + sy^2)+bn) )  %% Creates vector of widths
    R(k+1) = (bn*k)^2;
    B(k+1) = k+1;
    S(k+1) = 0;
    N(k+1) = 0;
    k=k+1;
end 

for y = 1:sy
r = ((1-1)*a)^2+(y-1)^2;
while( r < R(k) ) %% selects correct contour to be one BELOW current value
k = k-1; %% at BEGINNING of matrix row (thus all other moves are greater by one contour)
end
for x = 1:sx
r = ((x-1)*a)^2+(y-1)^2;
if( r >= R(k+1) ) %% Advances contour by one if appropriate
k = k+1;
end
S(k) = S(k) + T(y,x);
N(k) = N(k) + 1;
end
end
S = S(1:length(S)-2); %% concatenation to remove extra elements
B = B(1:length(S));
N = N(1:length(S));
S = S./N;
%figure
%plot((B-1).*bn,S) 