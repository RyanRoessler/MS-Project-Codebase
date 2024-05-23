% TAMSD function
function [XX,YY,ZZ]=TAMSD(x,y,DeltaMax)
% This function computes the TA-MSD of trajectory x,y
% if DeltaMax is not specified, it computes MSD up to delta=n-1

n = length(x);

if (nargin < 3)
    DeltaMax=n-1;
end % Default value for DeltaMax

for i = 1:DeltaMax
   XX(i) = mean((x((i+1):n)-x(1:(n-i))).^2);
   YY(i) = mean((y((i+1):n)-y(1:(n-i))).^2);
end
ZZ=XX+YY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes TA-EA-MSD and EA-MSD
% The data should be as XY columns
% TA stands for time-averaged
% EA stands for ensemble averaged
%
file_name='C:\Users\ryanr\OneDrive\Desktop\Krapf Research\Tracking\Capture_18.txt';
tracks_headers = 0; % Number of header rows in file
M=15; % number of lag times to use (do not exceed 1/5 of minlength for good results)
minlength=64; %minimum trajectory length to use in calculation
t = 0.2; % frame time in s
%N=101;% number of frames to use (not used)
%%
A0 = importdata(file_name, '\t', tracks_headers);
A0 = A0.data;
%A0=A0(1:N,:);
%%
% filter trajectories to have more than minlength frames
%
lengths=lengthsTrajFunction(A0);
len=NaN(1,length(lengths)*2);
len(1:2:length(len)-1)=lengths;
len(2:2:length(len))=lengths;
A0(:,len<minlength)=[];
%%
[p1,p2]=size(A0);
XX=NaN(M,p2/2);  %2 because it is XY
YY=XX;
ZZ=XX;
xav=XX;
yav=XX;

for i=1:p2/2	%2 because it is XY
%i=8; %trajectory number

x=A0(:,(i-1)*2+1);	%2 because it is XY
y=A0(:,(i-1)*2+2);	%2 because it is XY

x=InterpNans(x); %Replaces NaNs with values according to linear interpolation
y=InterpNans(y);
x(isnan(x))=[]; %Remove NaNs
y(isnan(y))=[];

[XX(:,i),YY(:,i),ZZ(:,i)]=TAMSD(x,y,M); %TA MSD

%EAMSD
xa=(x-x(1)).^2; %square displacements
ya=(y-y(1)).^2; %square displacements
xa=xa(2:M+1);
ya=ya(2:M+1);
xav(:,i)=xa;
yav(:,i)=ya;

end

ETMSDx=mean(XX,2);
ETMSDy=mean(YY,2);
EAMSDx=mean(xav,2);
EAMSDy=mean(yav,2);

% ETMSD=mean(ZZ,2);
% EMSD=mean(xt+yt,2);
% EMSD=EMSD(2:end);
% S=std(xt+yt,0,2)./sqrt(p2/3);
% S=S(2:end);

% The rest plots the results
lagt=t:t:M*t;
figure ('Name', 'Log-log plot of MSD')
loglog(lagt,ETMSDx*1e-6,'ob')
hold on
loglog(lagt,ETMSDy*1e-6,'o')
loglog(lagt,EAMSDx*1e-6)
loglog(lagt,EAMSDy*1e-6)
hold off
xlabel('t_{lag} (s)')
ylabel('MSD')
legend('EA-TA-MSDx', 'EA-TA-MSDy','EA-MSDx', 'EA-MSDy')
%saveas(gcf,'ETAMSD in log-log plot.png')

figure ('Name', 'Linear plot of ETAMSD')
plot(lagt,ETMSDx*1e-6,'ob')
hold on
plot(lagt,ETMSDy*1e-6,'o')
plot(lagt,EAMSDx*1e-6)
plot(lagt,EAMSDy*1e-6)
hold off
xlabel('t_{lag} (s)')
ylabel('MSD')
legend('EA-TA-MSDx', 'EA-TA-MSDy','EA-MSDx', 'EA-MSDy')
%saveas(gcf,'ETAMSD plot.png')

figure ('Name', 'Log-log plot of TAMSDx for individual particles')
loglog(lagt,XX*1e-6)
xlabel('t_{lag} (s)')
ylabel('MSD')
%saveas(gcf,'Individual TAMSDx.png')

figure ('Name', 'Linear plot of TAMSDx for individual particles')
plot(lagt,XX*1e-6)
xlabel('t_{lag} (s)')
ylabel('MSD')
%saveas(gcf,'Individual TAMSDy.png')


% time=0.05:0.05:(N-1)*0.05;
% 
% figure
% loglog(time,ETMSD*1e-6,'r');
% hold on
% loglog(time,(EMSD)*1e-6);
% loglog(time,(EMSD-1.96*S)*1e-6);
% loglog(time,(EMSD+1.96*S)*1e-6);
% hold off
% xlim([0 3.5])
% 
% figure
% plot(time,ETMSD*1e-6,'r');
% hold on
% plot(time,(EMSD)*1e-6);
% plot(time,(EMSD-1.96*S)*1e-6);
% plot(time,(EMSD+1.96*S)*1e-6);
% xlim([0 3.5])

%xlim([0 0.17])
%dlmwrite('stdevEA_MSD.txt',S,'\t');
%dlmwrite('T:\projects\nanobio\Carsten\170727\Tracks\EA_TA_MSDX5.txt',[lagt',ETMSDx*1e-6,ETMSDy*1e-6],'\t');
%dlmwrite('U:\MATLAB\trajectory analysis\test data\TA_MSDX5.txt',[lagt',XX*1e-6],'\t'); %saves in um^2
%dlmwrite('U:\MATLAB\trajectory analysis\test data\TA_MSDY5.txt',[lagt',YY*1e-6],'\t');
%dlmwrite('EA_MSD.txt',EMSD,'\t');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Computes TA-EA-MSD and EA-MSD
%
% ns: number of frames to use
% filter: 1 to use min length, 0 not to
% M: number of lag times. i.e. lags = 1,2,3, ..., M in frames
%
% ETMSDx: Ensemble-time averaged MSD in x or y
% ETMSD: Total ensemble-time averaged MSD (each column is one ns) 
% XX: TA-MSD in x for each trajectory - resets for each ns
% sem: standard error of the mean of MSD
%
tracks_headers = 0; % Number of header rows in file
directory='C:\Users\ryanr\OneDrive\Desktop\Krapf Research\RNA data analysis\Trajectories\All puro data concatenated\+20-30min\';
file_name='+20-30min_min64.txt'; % XY file (wiht large min length)
ns=[8 16 32 64]; % realization times
%
M=max(ns);		% number of lag times to use
filter=0;  % 1 to filter for min length, 0 for no filter
minlength=64; % minimum trajectory length to use in calculation if filter=1.
               % if filter=0, minlength is not used.

correct=1; %1/0.193; %correct A0 to um
t = 0.1; % frame time
sigma2 = 0.0035; % 2*sigma^2 errors for MSD correction
               
               
%% Do not change anything after this line, except the legend in the plot (line 98)
%
save_file_name=strcat(directory,'MSD',file_name); 
file_name=strcat(directory,file_name);

if tracks_headers>0
    A0 = importdata(file_name, '\t', tracks_headers);
    A0 = A0.data;
else
    A0 = importdata(file_name, '\t'); %when data does not have header
end

A0=A0*correct;

%% Filter trajectories to have more than minlength frames
% 
if filter==1
	lengths=lengthsTrajFunction(A0); %lengths of all trajectories
	len=NaN(1,length(lengths)*2);
	len(1:2:length(len)-1)=lengths; %X
	len(2:2:length(len))=lengths;   %Y
	A0(:,len<minlength)=[];
end

%% Initializations

[p1,p2]=size(A0);
p2=p2/2; %2 because it is XY. p2 is number of trajectories

XX=NaN(M,p2);  %TA-MSD x
YY=XX;	% TA-MSD y
%RR=XX;	% TA-MSD R

ETMSDx=NaN(M,length(ns));
ETMSDy=NaN(M,length(ns));

count=0;

%% Compute TA-MSD
for n=ns
    count=count+1;
    for i=1:p2	%2 because it is XY
        %i=8; %trajectory number
        
        x=A0(:,(i-1)*2+1);	%2 because it is XY
        y=A0(:,(i-1)*2+2);	%2 because it is XY
        
        x=InterpNans(x); %Replaces NaNs with values according to linear interpolation
        y=InterpNans(y);
        
        x(isnan(x))=[];
        y(isnan(y))=[];
        
        [XX(1:n,i),YY(1:n,i),~]=TAMSDv2_toN(x,y,n,n);
    end
    
    ETMSDx(:,count)=mean(XX,2);
    ETMSDy(:,count)=mean(YY,2);
    semETx(:,count)=(std(XX,0,2)./sqrt(p2))';
    semETy(:,count)=(std(YY,0,2)./sqrt(p2))';
end

%% Compute EA-MSD

X0=A0(:,1:2:p2*2-1);
Y0=A0(:,2:2:p2*2);
xt=(X0(:,:)-X0(1,:)).^2;
yt=(Y0(:,:)-Y0(1,:)).^2;
EAMSDx=mean(xt(1:M+1,:),2);
EAMSDy=mean(yt(1:M+1,:),2);
semEAx=(std(xt(1:M+1,:),0,2)./sqrt(p2));
semEAy=(std(yt(1:M+1,:),0,2)./sqrt(p2));
EAMSDx(1,:)=[];
EAMSDy(1,:)=[];
semEAx(1,:)=[];
semEAy(1,:)=[];

%% Compute MSD for r
EAMSDx=EAMSDx-sigma2;
EAMSDy=EAMSDy-sigma2;

EAMSD=EAMSDx+EAMSDy;
semEA=sqrt(semEAx.^2+semEAy.^2);

ETMSDx=ETMSDx-sigma2;
ETMSDy=ETMSDy-sigma2;

ETMSD=ETMSDx+ETMSDy;
semET=sqrt(semETx.^2+semETy.^2);

%% Plots the results
% 
lagt=t:t:M*t;
figure ('Name', 'Log-log plot of ETAMSD')
loglog(lagt,ETMSD)
hold on
errorbar(lagt,EAMSD,semEA)
hold off
%legend('8','16','32','64','EAMSD')
legend('16','32','64','128','EAMSD')
ylabel('MSD ({\mu}m^2)')
xlabel('t_{lag}(s)')

lagt=lagt';
%origin=[lagt EAMSDx semEAx EAMSDy semEAy EAMSD semEA ETMSDx ETMSDy ETMSD];
origin=[lagt EAMSD semEA ETMSD];
%label=["lagt" "EAMSDx" "semEAx" "EAMSDy" "semEAy" "EAMSD" "semEA" "ETMSDx8" "ETMSDx16" "ETMSDx32" "ETMSDx64" "ETMSDy8" "ETMSDy16" "ETMSDy32" "ETMSDy64" "ETMSD8" "ETMSD16" "ETMSD32" "ETMSD64"];
%label=["lagt" "EAMSDx" "semEAx" "EAMSDy" "semEAy" "EAMSD" "semEA" "ETMSDx8" "ETMSDx16" "ETMSDx32" "ETMSDx64" "ETMSDx128" "ETMSDy8" "ETMSDy16" "ETMSDy32" "ETMSDy64"  "ETMSDy128" "ETMSD8" "ETMSD16" "ETMSD32" "ETMSD64" "ETMSD128"];
label=["lagt" "EAMSD" "semEA" "ETMSD8" "ETMSD16" "ETMSD32" "ETMSD64" "ETMSD128"];


%saveas(gcf,'ETAMSD in log-log plot.png')


% time=0.05:0.05:(N-1)*0.05;
% 
% figure
% loglog(time,ETMSD*1e-6,'r');
% hold on
% loglog(time,(EMSD)*1e-6);
% loglog(time,(EMSD-1.96*S)*1e-6);
% loglog(time,(EMSD+1.96*S)*1e-6);
% hold off
% xlim([0 3.5])
% 
% figure
% plot(time,ETMSD*1e-6,'r');
% hold on
% plot(time,(EMSD)*1e-6);
% plot(time,(EMSD-1.96*S)*1e-6);
% plot(time,(EMSD+1.96*S)*1e-6);
% xlim([0 3.5])

%% Saves the results
% save('U:\SBME courses\ECE 681 Random Walks\Zach\experimental data\Nav16 data neurons\agingMSD260.mat')
%sem=sem';
%ETMSD=ETMSD';
% dlmwrite('U:\MATLAB\trajectory analysis\test data\TA_MSDY5.txt',[lagt',YY*1e-6],'\t');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes Power Spectral Density (PSD)
% The data should be as XY columns

directory='C:\Users\ryanr\OneDrive\Desktop\Krapf Research\RNA data analysis\2022-11-09\Removed - dn6, th0.55\';
file_name='nTracks_T0_Conc128.txt';
tracks_headers = 0; % Number of header rows in file
time1=0.1; %t_step (frames in s)
minlength=128; %minimum trajectory length to use in calculation
Ns = [2^4 2^5 2^6 2^7];  % set of times to assess aging
leg = {'t=16','32','64','128'}; % same times as Ns for figure legend

%%
% No need to change anything after this line

%% Read data
%
file_name=strcat(directory,file_name);
if tracks_headers>0
    X0 = importdata(file_name, '\t', tracks_headers);
    X0 = X0.data;
else
    X0 = importdata(file_name, '\t'); %when data does not have header
end

m=size(X0,2);
%% Select X coordinate and filter trajectories for minlength frames
% 
lengths=lengthsTrajFunction(X0); %lengths of all trajectories
X0=X0(:,1:2:m-1); %select X coordinate
X0(:,lengths<minlength)=[];
[frames,m]=size(X0);
% Interpolate and remove initial NaNs
for j=1:m
    X0(:,j)=InterpNans(X0(:,j));
    x1=X0(:,j);
    x1(isnan(x1))=[];
    X0(1:length(x1),j)=x1;
    if length(x1)<frames
        X0(length(x1)+1:frames,j)=NaN(frames-length(x1),1);
    end
end

%% 
%initialize whole arrays
p3=Ns(length(Ns))/2+1;  %number of frequencies
pmeanx=NaN(p3,length(Ns));
freqs=NaN(p3,length(Ns));
%
count=0;
for N = Ns
    display(N);
    j=(0:N/2)';
    w=j./(N*time1); %frequencies
    x=X0(1:N,:)-X0(1,:); %subtract first frame
    [pxx,f] = periodogram(x,[],w,1/time1);
    pxx=pxx/2;   % need to change to pxx/2
    f=f*2*pi;
        
    p3=length(f);
    count=count+1;
    pmeanx(1:p3,count)=mean(pxx')';  %mean PSD
    freqs(1:p3,count)=f;
end

pmeanx=pmeanx*2;

loglog(freqs,pmeanx)
legend(leg)
xlabel('\omega')
ylabel('\langle S(\omega,t)\rangle')
%title(['m=' num2str(m) ', \alpha=' num2str(alpha) ', H=' num2str(H)])
%load gong
%sound(y,Fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Segment trajectory state only function
function state=SegmentTrajectoryStateOnly(Sd,t,thresh)

% Segment trajectories according to threshold
%
% th   threshold
%
% Sd   time series (data), e.g., local convex hull
% t    time array (start from time>1 because of window size)
% Sd and t must have same length
% 
% Confined: state=1
% Free: state=0
% 

%%
winn=t(1)*2-1; % window size
f1=(find(Sd<thresh)); % immobile points
state=zeros(1,length(Sd));
state(f1)=1;
if winn>1
    state=[ones(1,t(1)-1)*state(1) state ones(1,t(1)-1)*state(end)]'; % fill beginning and end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convex Hull function for single trajectory
function [t,Sd]=ConvexHullSingleTraj(x,y)

%------------------------------------------------------
% max. diameter of Local Convex Hull (LCH), single traj
%------------------------------------------------------
% dt   time increment
% dn   (# vertices in one side for convex hull)
% x    x coordinates (1:N)
% y    y coordinates (1:N)
%------------------------------------------------------

dn=6;
dt=1;

%%
x(isnan(x))=[];
y(isnan(y))=[];
N    = numel(x);
t    = [];
Sd   = [];
%ab=0;
for i=dn+1:N-dn
     sx = x(i-dn:i+dn);
     sy = y(i-dn:i+dn);
     Dt = convhull(sx,sy);
     nn = numel(Dt);
     ab = 0;
     for jj=1:nn
         pf = Dt(jj+1:nn);
         aa = max(sqrt((sx(Dt(jj))-sx(pf)).^2+(sy(Dt(jj))-sy(pf)).^2));
         if (aa > ab) ab=aa; end
     end
     Sd = [Sd,ab];
     t  = [t,dt*i];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script segments trajectories according to 
% the local convex hull and it removes the first part of the trajectory up to the 
% first data point that belongs to a free start (state UP or 1). 

% Trajectories are first segmented into 2 states called FREE (UP,1) and CONFINED (DOWN,0).
% Trajectories without free parts are removed.
% Trajectories that start within the free state are left unchanged.

tracks_headers = 0; % Number of header rows in file, 0 for no header

directory='Z:\krapflab\Ryan\Analysis and concatenated files\Puromycin and Cycloheximide\2023-05-16 (cycloheximide)\All (both incubators, 23-05-16)\Ethanol +19-31min\';
file_name='Kipper_EtOH_+19-31min_min128.txt'; % XY file (high min length)
minlength=13; %minimum length of trajectory to analyze

thresh=0.37;

  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not change anything after this line

savefilename=strcat(directory,'T0_',file_name); %modified file name to save
file_name=strcat(directory,file_name);
if tracks_headers>0
    A0 = importdata(file_name, '\t', tracks_headers);
    A0 = A0.data;
else
    A0 = importdata(file_name, '\t'); %when data does not have header
    %A0 = A0.data % I'VE ADDED THIS LINE
end


[p1,N]=size(A0);% N is twice the number of trajectories

x=A0(:,1:2:N-1); %trajectory x coordinates
y=A0(:,2:2:N);   %trajectory y coordinates

N=N/2; % N now is number of trajectories
removals=[]; %array of trajectories to remove
% first=[]; %first free point
% firsttot=[]; %first index

zeroindex=find(y==0);
zeroindex=ceil(zeroindex/p1);
x(:,zeroindex)=[];
y(:,zeroindex)=[];
[p1,N]=size(x);

for i=1:N
    x1=x(:,i); %specific trajectory
    y1=y(:,i);
    
    if length(x1(~isnan(x1)))>=minlength % do not analyze if <minlength
        x1=InterpNans(x1)'; % Replace inner NaNs with interpolated values
        y1=InterpNans(y1)';
        
        [t,Sd]=ConvexHullSingleTraj(x1,y1);
        state=SegmentTrajectoryStateOnly(Sd,t,thresh);
		firstpt=find(state==0,1); %first free point
		firstindex=find(~isnan(x1),1); %index of trajectory first point
		if isempty(firstpt) % no free part, remove trajectory
			removals=[removals;2*i-1 2*i];
		elseif firstpt>1 %remove first confined part
            x1(1:firstpt-1)=[];
            x1=[x1;nan(firstpt-1,1)];
            y1(1:firstpt-1)=[];
            y1=[y1;nan(firstpt-1,1)];
            
			A0(:,2*i-1:2*i)=[x1 y1];
            %first=[first firstpt];
            %firsttot=[firsttot firstindex];
        end
%         if isempty(firstindex)
%             firstindex=0;
%         end
    else
        removals=[removals;2*i-1 2*i];
    end
    
    
end

A0(:,removals)=[];

%first=[first' firsttot'];

dlmwrite(savefilename,A0,'delimiter','\t','precision','%.3f');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script segments trajectories according to 
% the local convex hull. 
% Trajectories are segmented into 2 states called UP (1) and DOWN (0).
% It provides four matrices of only UP and DOWN states for x and y
% and the state matrix (0 or 1) of the original data
%
% Trajectories are segmented into 2 states called CONFINED (UP,1) and FREE (DOWN,0).
%
% The code only works on files where all the trajectories are the same length
% It is useful for the output of ConcatenateTrajFiles.m
%
% state:    State matrix
% L0: durations of free state
% L1: durations of immobilizations
% freeX, freeY: trajectories in free state (confined are NaNs)
% confinedX, confinedY: trajectories in confined state

close all % closes all figure windows
tracks_headers = 0; % Number of header rows in file

directory='C:\Users\ryanr\OneDrive\Desktop\Krapf Research\RNA data analysis\Trajectories\05-16-23 (cyclo)\All\Cyclo +19-31min\';
file_name='nnCyclo_+19-31min_min64.txt'; % XY file (high min length)
lag = 0.1;		% lag time(in frames) for displacements
pixelsize = 0.130; % Pixel size in um  
                % use 1 to keep pixels 
minlength=64; %minimum length of trajectory to analyze
thresh=0.37; %threshold for convex hull

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not change anything after this line

file_name=strcat(directory,file_name);
if tracks_headers>0
    A0 = importdata(file_name, '\t', tracks_headers);
    A0 = A0.data;
else
    A0 = importdata(file_name, '\t'); %when data does not have header
end


[p1,N]=size(A0);% N is number of trajectories*2

x=A0(:,1:2:N-1); %trajectory x coordinates
y=A0(:,2:2:N);

N=N/2; % N now is number of trajectories
state=NaN(p1,N); % state matrix (all trajectories), 1 or 0
L0=[]; %array of dwell times of free state
L0traj=[]; %%% I'M ADDING FOR L > 40 ANALYSIS %%%
L1=[]; %array of dwell times of immobilizations
ConvHull=[];

for i=1:N
    x1=x(:,i);
    y1=y(:,i);
    
    if length(x1(~isnan(x1)))>=minlength
        x1=InterpNans(x1)';
        y1=InterpNans(y1)';
        
        [t,Sd]=ConvexHullSingleTraj(x1,y1);
        [~,changeup,lengths0,lengths1,st]=SegmentTrajectory(Sd,t,thresh);
        L0  = [L0;i;i;i;i;i;lengths0]; % I'M CHANGING FOR L > 40 ANALYSIS
        L1  = [L1;lengths1];

        ConvHull = [ConvHull;Sd'];
        state(~isnan(x1),i)=st; %State matrix (1 or 0)
        %
        x1(isnan(x1))=[];
        y1(isnan(y1))=[];
        ximm=x1(changeup); %immobilization coordinates
        yimm=y1(changeup);
        
    end
end

confinedstate=state;
freestate=1-state;
confinedstate(state==0)=NaN;
freestate(state==1)=NaN;
freeX=freestate.*x;
freeY=freestate.*y;
confinedX=confinedstate.*x;
confinedY=confinedstate.*y;

%I'm adding
% b=[]
% for i=1:N
% a=find(L0>40)
%     if ~isempty(a)
%         b=[a;i]
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Segment trajectories according to threshold, then plot
%
% Sd   time series (data), e.g., local convex hull
% x,y  trajectory
%
% th   threshold

th=0.37; % threshold
%x=0.193*x
%y=0.193*y

%%
%
winn=t(1)*2-1;
f1=(find(Sd<th));
state=zeros(1,length(Sd));
state(f1)=1;
if winn>1
    state=[ones(1,t(1)-1)*state(1) state ones(1,t(1)-1)*state(end)]';
end

upx=x;
upy=y;
upx(state==0)=NaN;
upy(state==0)=NaN;
figure(4)
%plot(x,y)
plot(x,y,'b','Linewidth',1.5) %%% divided by 0.193 to get original pixels in imageJ %%%
hold on
plot(upx,upy,'r','Linewidth',1.5)
hold off
xlabel('X (\mum)','FontSize',20);
ylabel('Y (\mum)','FontSize',20);
set(gca,'FontSize',20)
daspect([1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this script should create an array of stepsizes from a file of
%trajectories so that I can find zeros and do a mass deletion of bad
%trajectories from collinear error
close all

directory='C:\Users\ryanr\OneDrive\Desktop\Krapf Research\RNA data analysis\2022-11-09\Removed - dn6, th0.55\';
file_name='nTracks_T0_Conc64.txt'; % XY file (high min length)
tracks_headers = 1; % Number of header rows in file, 0 for no header

savefilename=strcat(directory,'diff_',file_name); %modified file name to save
file_name=strcat(directory,file_name);
if tracks_headers>0
    A0 = importdata(file_name, '\t', tracks_headers);
    A0 = A0.data;
else
    A0 = importdata(file_name, '\t'); %when data does not have header
end

stepsize = diff(A0,2); %takes difference twice so that trajectories with random, single zeros don't get deleted

[p1,N]=size(stepsize);% N is twice the number of trajectories

x=stepsize(:,1:2:N-1); %trajectory x coordinates
y=stepsize(:,2:2:N);   %trajectory y coordinates
N=N/2;

zeroindex=find(y==0);
zeroindex=ceil(zeroindex/p1);
x(:,zeroindex)=[];
y(:,zeroindex)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Alpha, K function
function [alpha,K]=AlphaK(MSD,deltas)

% Computes alpha and K from MSD=K t^alpha
%
% deltas is number of times to use in fit

MSD=MSD(1:deltas);
times=1:deltas;
Plog = polyfit(log(times),log(MSD),1);  %linear fit (log-log scales)
alpha = Plog(1);
K = exp(Plog(2))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes alpha and K from MSD=2K t^alpha
% Also computes D and sigma2 where MSD=2Dt+sigma2
%
% N is number of times to use
%
clear all
close all
%tracks_headers = 1; % Number of header rows in file it can also be 2
file_name='Z:\krapflab\Ryan\Analysis and concatenated files\Snehal miRNA\13 fps, single plane\DKalpha_Cell4_min64.txt'; %File with MSD organized as 
								  %Time in first column
                                  %MSDs in all other columns
N=11;		% number of times to use in fit (8 points)

A0 = importdata(file_name);
%A0 = importdata(file_name, '\t', tracks_headers);
%A0 = A0.data;
times=A0(1:N,1);
A0=A0(1:N,2:end); % First column is time
[p1,p2]=size(A0);

alpha=NaN(p2,1);  
K=alpha;
D=alpha;
sigma2=alpha;

for i=1:p2
%i=5; %trajectory number for tests
MSD=A0(:,i);
Plog = polyfit(log(times),log(MSD),1);  %linear fit (log-log scales)
Plinear = polyfit(times,MSD,1); %linear fit
alpha(i) = Plog(1); % Slope from linear fit (log-log scale)
K(i) = exp(Plog(2))/2; % (e^(intercept))/2 (log-log scale)
D(i) = Plinear(1)/2; %in um^2/s
sigma2(i) = Plinear(2); %intercept in um^2
end

K_alpha_D_sigma2=[K alpha D sigma2];
% dlmwrite('K_alpha_D_sigma2Y.txt',K_alpha_D_sigma2,'\t');

plot(K,alpha,'+')
xlabel('K') % x-axis label
ylabel('alpha') % y-axis label

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots any specific trajectory from an XY file
%
trajectory_file='C:\Users\ryanr\OneDrive\Desktop\Krapf Research\RNA data analysis\Filtered_min64_T001_Tracks.txt';
trajectorynumber=13; % number of the trajectory to plot
%%
trajectories = importdata(trajectory_file, '\t', 1);
trajectories=trajectories.data;
x=trajectories(:,trajectorynumber*2-1);
y=trajectories(:,trajectorynumber*2);
figure(1)
%plot(x/0.193,y/0.193) %%% divided by 0.193 to get original pixels in imageJ %%%
plot(x,y, 'LineWidth',1.5);
xlabel('X (\mum)');
ylabel('Y (\mum)');
%set(gca, 'YDir','reverse') % I'VE ADDED THIS LINE
%hold on
% ------------------------------
% I'VE ADDED PER STEP 3c
% plot(t,Sd)
% xlabel('Time');
% ylabel('Convex Hull');
% ------------------------------
%
%xx=trajectories(:,1:2:length(trajectories)-1);
%yy=trajectories(:,2:2:length(trajectories));
%plot(xx,yy)
daspect([1 1 1])

%figure(2)
%plot(x)
%hold on
%plot(y)


