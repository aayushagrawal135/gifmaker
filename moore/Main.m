clc
%clearvars
close all

sclockroot = 111110000;     % change 1's to first 5 digits of sclock 
soln = 1;    % Sol
framen = 8; % # of images in sequence



offset = [42924];   % Size of header. For reading in images. Offset = image size - (height*width*bit depth)

Height = 511;   % Navcam pixel height after binning
Width = 511;    % Navcam pixel width after binning
%n = 1;

Amean = zeros(Height,Width,1,1);
sclock = sclockroot;
while( sclock < (sclockroot+20000) && i < framen)

    fid = -1;
    sclock = sclock +1;
    while( sclock < (sclockroot+20000) && fid == -1 )
        modspec = 'RAD_M0550538NCAM00536M'; % adjust for each sol
        fname = strcat('NRB_',sprintf('%i',sclock),modspec,'1.IMG');
        fid = fopen(fname);
        sclock = sclock +1;
    end
    sclock = sclock-1;
    fname2 = strcat('NRB_',sprintf('%i',sclock),modspec,'2.IMG');
    fname3 = strcat('NRB_',sprintf('%i',sclock),modspec,'3.IMG');
    if(fopen(fname2) ~= -1)
        fname = fname2;
    end
    if(fopen(fname3) ~= -1)
        fname = fname3;
    end
    i=i+1;
    disp(fname)
    A = multibandread(fname,[Height,Width,1],'uint16',offset,'bsq','ieee-be');  %read in images
    B(:,:,1,i) = double(A);
    Amean = Amean + B(:,:,1,i);  % assemble mean frame
end



if(sclock >= sclockroot+19999)
    disp('Error, Maximum Travel Exceeded')
end

Amean = Amean./framen;  %calculate mean frame
iter = 3;


for i=1:framen
    X(:,:,1,i) = B(:,:,1,i)-Amean(:,:,1,1); % subtract mean frame from each image
    %Y(:,:,1,i) = Zstretch(X(:,:,1,i),iter);
    
    A = X(:,:,1,i);

        % Convert to vector
        nrow = length(A(:,1));
        ncol = length(A(1,:));

        extA = zeros(ncol*nrow,1);
        extA(1:nrow) = A(:,1);
        for i=1:ncol-1
            extA(i*nrow+1:(i+1)*nrow) = A(:,i+1);
        end

        % Determine the lower and upper limits
        firstMed = median(extA); 

        lowLimit = findMed(extA,-1,firstMed,iter);
        highLimit = findMed(extA,1,firstMed,iter);

        % Exclude out of range pixels
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

        % Perform Stretch
        A = A-min(min(A));
        A = A./max(max(A)).*256;

        Y(:,:,1,i) = A;
    
    
    Z(:,:,1,i) = cast(Y(:,:,1,i),'uint8');
end

imwrite(Z,sprintf('sol0%i_%i_Zenith.gif',soln,sclock),'gif','DelayTime',0.10,'LoopCount',inf); % create cloud animation

rowN = 362;

overallmS = zeros(1,rowN);   %363
overallRmS = zeros(1,rowN);
overallCmS = zeros(1,rowN);
C  = zeros(1,rowN);


% FIND HIGH AND LOW RADIANCE

for currN = 1:framen

     if(currN == 4) %  use middle frame of sequence (4,5,6)
        uplimit = 25;  
        lowlimit = -30;


        highPointy = 418;   %high pixel y
        highPointx = 177;   %high pixel x
        
        lowPointy = 120;    %low pixel y
        lowPointx = 60;     %low pixel x
        hlpRadius = 5;      %pixel radius for averagign
        %initialize
        highTotal = 0;
        lowTotal = 0;
        nPointPoints = 0;
        

        for i = (-1*hlpRadius):1:hlpRadius
            for j = (-1*hlpRadius):1:hlpRadius
                if( i^2+j^2 <= hlpRadius^2 )
                    highTotal = highTotal + X(i+highPointy,j+highPointx,1,currN);
                    lowTotal = lowTotal + X(i+lowPointy,j+lowPointx,1,currN);
                    nPointPoints = nPointPoints+1;
                    hpvec(nPointPoints) = X(i+highPointy,j+highPointx,1,currN);
                    lpvec(nPointPoints) = X(i+lowPointy,j+lowPointx,1,currN);
                end
            end
        end
        
        highTotal = highTotal/nPointPoints;
        lowTotal = lowTotal/nPointPoints;
       

        ST(1) = (highTotal-lowTotal)/mean(mean(B(10:511,:,1,currN)));
        ST(2) = (highTotal-lowTotal);
        ST(3) = std(hpvec);
        ST(4) = std(lpvec);
        
        disp(ST)
        
        for i = 1:511
            for j = 1:511
        
                if( X(i,j,1,currN) > uplimit )
                    X(i,j,1,currN) = uplimit;
                end
                if( X(i,j,1,currN) < lowlimit )
                    X(i,j,1,currN) = lowlimit;
                end
            end
        end
        
        figure(1) %% sufrace plot of currN
        surf(1:511,10:511,X(10:511,:,1,currN))
        view(2)
        title(sprintf('Time-Variable Component of Zenith Movie on Sol %i, Frame %i',soln,currN'),'FontSize',18)
        xlabel('Pixels horizontally along NavCam Frame','FontSize',16)
        ylabel('Pixels vertically along NavCam Frame','FontSize',16)
        set(gca,'FontSize',16); shading interp; axis tight; colorbar  
    end


%% ENTER THE SPECTRAL DOMAIN!
map = X(:,:,1,currN);

m = length(map(:,1));
n = length(map(1,:));

for i = 1:n
    k = (i-1)*m+1;
    outvec(k:k+m-1) = map(:,i);
end

clear m n k;


meanr = mean(outvec);


stdr = (0.005./3).*max(max(B(10:511,:,1,currN))); % synthetic noise distribution based on camera SNR
maxr = meanr + 2.*stdr; 
minr = meanr - 2.*stdr;

randmap = meanr + stdr.*randn(length(map(1,:)),length(map(:,1)));

nx = length(map(1,:));
ny = length(map(:,1));
nxq = floor(nx/2)+1; %% nyquist elements in x and y
nyq = floor(ny/2)+1; 


for i = 1:nx %% welch window processing
    for j = 1:ny
        welch(j,i) = (1-(((i-1)-(nx-1)/2)/((nx-1)/2))^2 ) *(1-(((j-1)-(ny-1)/2)/((ny-1)/2))^2 );
    end
end 


map = map.*welch;   %% Applying the window
randmap = randmap.*welch;  

MAP = fft2(map);    %% 2D Fourier transform
RANDMAP = fft2(randmap);
    

M = zeros(nyq,nxq);


%%% superimposes quadrants 1 and 2
M(:,1) = MAP(1:nyq,1).*2;

for j = 2:nxq
    if( nx-(j-2) > nxq )
        M(:,j) = MAP(1:nyq,j)+MAP(1:nyq,nx-(j-2));
    else
    M(:,j) = MAP(1:nyq,j).*2;
    end
end 



%%% superimposes quadrants 1 and 2
RM = zeros(nyq,nxq);
RM(:,1) = RANDMAP(1:nyq,1).*2;

for j = 2:nxq
    if( nx-(j-2) > nxq )
        RM(:,j) = RANDMAP(1:nyq,j)+RANDMAP(1:nyq,nx-(j-2));
    else
    RM(:,j) = RANDMAP(1:nyq,j).*2;
    end
end 


PM = M.*conj(M);        %% map power
PRM = RM.*conj(RM);     %%random power
CPM = M.*conj(RM);      %% cross power


tau=511;
omega = 2*pi/tau; %% frequency domain

% AVERAGING PROTOCOL
bn = 1; %% matrix elements per average

%[mB,mS] = mAVG(PM,ea,1); %% concentric circular averaging (compensated for latitude) 

a = 1;

sx = length(PM(1,:));
sy = length(PM(:,1));
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
S(k) = S(k) + PM(y,x);
N(k) = N(k) + 1;
end
end
S = S(1:length(S)-2); %% concatenation to remove extra elements
mB = B(1:length(S));
N = N(1:length(S));
mS = S./N;

clear a sx sy k R B S N r 

a = 1;

%[RmB,RmS] = mAVG(PRM,bn,1);

sx = length(PRM(1,:));
sy = length(PRM(:,1));
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
S(k) = S(k) + PM(y,x);
N(k) = N(k) + 1;
end
end
S = S(1:length(S)-2); %% concatenation to remove extra elements
RmB = B(1:length(S));
N = N(1:length(S));
RmS = S./N;

clear a sx sy k R B S N r 

%[CmB,CmS] = mAVG(CPM,bn,1);

a = 1;

%[RmB,RmS] = mAVG(PRM,bn,1);

sx = length(CPM(1,:));
sy = length(CPM(:,1));
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
S(k) = S(k) + PM(y,x);
N(k) = N(k) + 1;
end
end
S = S(1:length(S)-2); %% concatenation to remove extra elements
CmB = B(1:length(S));
N = N(1:length(S));
CmS = S./N;


overallmS = overallmS+log10(mS);
overallRmS = overallRmS+log10(RmS);
overallCmS = overallCmS+log10(CmS); 

disp('inloop')

end



overallmS = overallmS./framen;
overallRmS = overallRmS./framen;


figure(7)
semilogx(((mB-1).*ea+1).*omega,overallmS.^2,'r')
hold on
semilogx(((RmB-1).*ea+1).*omega,overallRmS.^2,'b')
hold off

figure(8)
semilogx(((RmB-1).*ea+1).*omega,overallmS.^2./overallRmS.^2,'g')
grid on

