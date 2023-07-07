tic
rng(6,'twister');

% loading of displacement
close all

% Loading the dic displacement map
Matdic = importdata('image_dic.csv');
Matdic = Matdic.data;

Xdic = Matdic(:,2);
Ydic = Matdic(:,3);
valuedic = Matdic(:,5);
xyzdic = [Ydic, -Xdic, -valuedic]; 
ptClouddic = pointCloud(xyzdic);
ptCloudOut = pcdenoise(ptClouddic, 'NumNeighbors', 5, 'Threshold', 0.5);
ptCloudOut = pcdenoise(ptCloudOut, 'NumNeighbors', 5, 'Threshold', 0.1);
xyzdic = ptCloudOut.Location;
limXdic = [min(xyzdic(:,1)), max(xyzdic(:,1))];
limYdic = [min(xyzdic(:,2)), max(xyzdic(:,2))];
ptClouddic = pointCloud(xyzdic);

% Showing the dic map
scalev = 1/30; % Z-axis scale
figure
pcshow(ptClouddic, 'MarkerSize', 70)
daspect([max(daspect)*[1 1] scalev]) 
c1 = colorbar;
c1.Label.String = 'uY [mm]';
xlabel('X [mm]');
ylabel('Y [mm]');
title('dic Results')
view(0,90)

% Script that carries out the fitting between the dic results
close all
tic

% Adjusting the initial translation
sX = +0.5447;  
sY = -0.1173;   
sZ = 0.0258;   
xyzdicnew1(:,1) = xyzdic(:,1) + sX;
xyzdicnew1(:,2) = xyzdic(:,2) + sY;
xyzdicnew1(:,3) = xyzdic(:,3) + sZ;

% Adjusting the initial orientation
alphar = -deg2rad(0);
npoints = size(xyzdic,1);
eul = [alphar 9*10^(-4) 0];  
rotm = eul2rotm(eul);
xyzdicnew2 = rotm*xyzdicnew1';
xyzdicnew2 = xyzdicnew2';
xdata = xyzdicnew2(:,1:2);  
ydata = xyzdicnew2(:,3);    
ds = 5; 
xdata = downsample(xdata,ds);
ydata = downsample(ydata,ds);

% Theoretical model
fun = @(par,xdata)matrixcalculator(xdata, par);
nmax = 15; 
           
% Initialization
KIserie = NaN(1,nmax);
KIIserie = NaN(1,nmax);
Tserie = NaN(1,nmax);
sXserie = NaN(1,nmax);
sYserie = NaN(1,nmax);
sZserie = NaN(1,nmax);
alphaserie = NaN(1,nmax);

% STIFs initial values 
KI = 1400;               % [Mpa*mm^0.5]
T =  -22;                % [Mpa]
KII = 0;                 % [Mpa*mm^0.5]

nvet = 2:nmax;

for n = 2:nmax % n is the number of terms
    tic
    L = 2*n + 4;
    x0 = zeros(1,L);  % Initial values initialization
    lb = zeros(1,L);  % Lower bound initialization
    ub = zeros(1,L);  % Upper bound initialization
    % x0, lb, ub are in the form: [a(i)...b(i)...deltaX,deltaY,deltaZ,alpha]; 
    % x0, lb, ub for the terms of type "a" -------------------------------%
    x0(1:2) = [KI/sqrt(2*pi) T/4];      
    lb(1:2) = [0    -300];                 
    ub(1:2) = [3000 +300];               
    for i = 3:(L-4)/2
        x0(i) = 0;        
        lb(i) = -100;             
        ub(i) = +100;       
    end
    % x0, lb, ub for the terms of type "b" -------------------------------%
    Lb = (L-4)/2;
    x0(Lb+1:Lb+2) = [-KII/sqrt(2*pi) 0];
    lb(Lb+1:Lb+2) = [-300 -300];
    ub(Lb+1:Lb+2) = [+300 +300];
    for i = (Lb+3:L-4)
        x0(i) = 0;
        lb(i) = -100;
        ub(i) = +100;
    end
    % x0, lb, ub for the terms of type [deltaX, deltaY, deltaZ, alpha]----%
    delta = 0.1;
    deltaZ = 0.5;
    x0(L-3:L) = [+0;       +0;      +0;       +0   ];  
    lb(L-3:L) = [-delta;   -delta;  -deltaZ;  -0.01];
    ub(L-3:L) = [+delta;   +delta;  +deltaZ;  +0.01];
    
    % Non-linear curve fitting
    [paropt,resnorm] = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
    
    % Saving the values of the parameters on vectors
    KIserie(n) = paropt(1)*sqrt(2*pi);
    KIIserie(n) = -paropt((L-4)/2+1)*sqrt(2*pi);
    sXserie(n) = paropt(L-3);
    sYserie(n) = paropt(L-2);
    sZserie(n) = paropt(L-1);
    alphaserie(n) = paropt(L);
    
    disp(n);
    toc
end

meanKI = (sum(KIserie(8:nmax))/(nmax-7))*10^(-3/2); 


%----- Figures result for different n)----------%


% SIFs
figure 
% KI
plot(nvet,KIserie(2:nmax)*10^(-3/2),'b--o','Linewidth',2);
xlabel('n');
ylabel('SIFs [MPa*m^{0.5}]');
%title('K{I}');
xlim([2 nmax])
ylim([-10 160]) %[-5 100]
xL = xlim;
line(xL, [meanKI meanKI], 'Color','cyan','LineWidth',2);
hold on

% Script that shows the results of the fitting
pardictoTheory = paropt;

% SIFs 
KIfitting = meanKI; 

% Printing the results
fprintf('KI = %4.2f\n',KIfitting);

% Construction of a grid to show the theoretical model. 
limXdic = [min(xyzdicnew1(:,1)), max(xyzdicnew1(:,1))];
limYdic = [min(xyzdicnew1(:,2)), max(xyzdicnew1(:,2))];
xG = limXdic(1):0.05:limXdic(2);
yG = limYdic(1):0.05:limYdic(2);
[XG, YG] = meshgrid(xG,yG);
matxyG = zeros(size(XG,1)*size(XG,2),2);
ind = 0;
for i = 1:size(XG,1)
    for j = 1:size(XG,2)
        ind = ind +1;
        matxyG(ind,:) = [XG(i,j) YG(i,j)];
    end
end
matxyGrot = [matxyG(:,1)*cos(alphar)-matxyG(:,2)*sin(alphar), matxyG(:,1)*sin(alphar)+matxyG(:,2)*cos(alphar)];

% Computing the y-displacement of the points of the grid
uYtheory = matrixcalculator(matxyGrot,pardictoTheory);
xyztheory = [matxyGrot,uYtheory];
ptCloudtheory = pointCloud(xyztheory);

% Computing the difference between dic and the theoretical model
uYtheorydicgrid = matrixcalculator(xyzdicnew2(:,1:2),pardictoTheory);
xyztheorydicgrid = [xyzdicnew2(:,1:2),uYtheorydicgrid];
ptCloudtheorydicgrid = pointCloud(xyztheorydicgrid);
Zdif = (xyzdicnew2(:,3) - uYtheorydicgrid);
cut = 0.010; 
Zdif(Zdif > cut| Zdif < -cut) = NaN; 
xyzdif = [xyzdicnew2(:,1:2), Zdif];
ptClouddiff = pointCloud(xyzdif);



%  -------------RESULTS VISUALIZATION


% Theoretical model
figure
pcshow(ptCloudtheory, 'MarkerSize', 70);
hold on
plot(0,0,'go','LineWidth',3);
daspect([max(daspect)*[1 1] scalev]) 
c1 = colorbar;
c1.Label.String = 'v [mm]';
c1.Label.FontWeight = 'bold';
c1.Limits = [-0.05 0.05];
xlabel('X [mm]');
ylabel('Y [mm]');
title('Theory');
colormap(gca,'jet')
view(0,90)
set(gca,'fontsize',10);
set(gca,'FontWeight','bold');
set(gca,'LineWidth',1.5);
set(gca,'TickLength',[0.018, 0.018]);


% dic results
ptClouddicnew = pointCloud(xyzdicnew2);
figure
pcshow(ptClouddicnew, 'MarkerSize', 70)
hold on
plot(0,0,'go','LineWidth',3);
daspect([max(daspect)*[1 1] scalev]) 
c1 = colorbar;
c1.Label.String = 'v [mm]';
c1.Label.FontWeight = 'bold';
c1.Limits = [-0.05 0.05];
xlabel('X [mm]');
ylabel('Y [mm]');
title('dic')
colormap(gca,'jet')
view(0,90)
set(gca,'fontsize',10);
set(gca,'FontWeight','bold');
set(gca,'LineWidth',1.5);
set(gca,'TickLength',[0.018, 0.018]);


% Difference
figure
pcshow(ptClouddiff, 'MarkerSize', 70)
hold on
daspect([max(daspect)*[1 1] scalev]) 
c1 = colorbar;
c1.Label.String = '\Deltav [mm]';
c1.Label.FontWeight = 'bold';
xlabel('X [mm]');
ylabel('Y [mm]');
title('dic - Theory')
colormap(gca,'jet')
view(0,90)
set(gca,'fontsize',10);
set(gca,'FontWeight','bold');
set(gca,'LineWidth',1.5);
set(gca,'TickLength',[0.018, 0.018]);

toc
