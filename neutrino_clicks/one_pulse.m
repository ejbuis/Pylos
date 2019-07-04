#!/usr/bin/octave -qf
addpath("/home/ejbuis/research/cosmics/simulations/ACpython/old/Acorne");
args = argv;
filename  = strcat (args{3}, "_", args{1}, "_", args{2}, ".dat")
zpos = str2num(args{1})
rpos = str2num(args{2})

rsc=[.5:9.5 15:10:105]; %radial bin centres (cm)
zsc=10:20:2000;   %Longitudinal Bin Centres (cm)
Eo=1e13; %Primary Energy
Do=[rpos 0 zpos]; % Position of observer 
fs=1e6; %sampling frequency

%t_axis=(-512:511)/fs; %time axis for plot (default 1024 points)
t_axis=(-1024:1023)/fs; % if you redefine it her, also in other scripts...
atten=1; % Learned's attenuation
nmc=1e6; % Number of MC points 
tsmc=ShowerParm(rsc,zsc,Eo,'Sloan');

%as the 10-100cm bins are 10x wider need to scale by a factor of 10 
tsmc=tsmc*diag(kron([1 10],ones(1,10)));

%generate MC points. Note bin EDGES need to be provided 
pointsc=MCGEn(tsmc,[0 zsc+10],[0 rsc+[0.5*ones(1,10) 5*ones(1,10)]],nmc);

%Convert to cartesian
[x,y,z]=pol2cart(rand(nmc,1)*2*pi,pointsc(:,2),pointsc(:,1));

% Convert fom cm to m 
points=[x y z]*1e-2;

p=kernelfr2(points,Do,log10(Eo),atten,10);

% save the data 
%filename = args{2}
tp = [t_axis' p]; % check this weird construction of the transpose columns
save(filename, 'tp');

