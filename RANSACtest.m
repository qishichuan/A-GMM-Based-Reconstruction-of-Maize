clear
clc

%% ------ load data ------

% M = 10;%the number of views, Nj is the cardinality of view V{j}
% idx = transpose(1:M);
% fname = arrayfun(@(idx) sprintf('C:/Users/hzyll/OneDrive/Desktop/Research/LiuLu/Data/SecondExperiment/%d.0.txt',idx),idx,'uniformoutput',false);
% 
% V = cellfun(@(fname) dlmread(fname,' '),fname,'uniformoutput',false);
% df = 1; % df=1 means no downsampling
% [V,I] = cellfun(@(V) deal(V(1:df:end,1:3)',double(V(1:df:end,4:6))),V,'uniformoutput',false);
% 
% aa = V{1}';
% Ind = aa(:,2)>1.3 & aa(:,2)<2 & aa(:,1)>-0.6 & aa(:,1)<0;
% 
% % bb = pointCloud(aa(Ind,:));%,'Color',iii);
% % pcshow(bb)
% points_ori = aa(Ind,:);
% points = points_ori(:,[2 3]);

%% ------ Generate data ------

% X = mvnrnd([0;0],0.04*eye(2),1000);
% coe = [1;2];
% y = X*coe+normrnd(0,0.5,size(X,1),1);
% outlierIdx = randperm(size(X,1),100);
% y(outlierIdx) = X(outlierIdx,:)*coe+normrnd(2,0.2,numel(outlierIdx),1);
% points = [X,y];
% plot3(points(:,1),points(:,2),points(:,3),'kx')
% hold on
% plot3(points(outlierIdx,1),points(outlierIdx,2),points(outlierIdx,3),'ro')
load single101_200%bb
df = 5;
points_ori = double(ptCloud.Location(1:df:end,:));
points = points_ori;

%% ------ RANSAC ------

sampleSize = 50; % number of points to sample per trial
maxDistance = 0.05; % max allowable distance for inliers

fitLineFcn = @(points) fit3dline(points); % fit function using polyfit
evalLineFcn = ...   % distance evaluation function
  @(model, points) dis2line(model, points);

[modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn, ...
  sampleSize,maxDistance,'MaxNumTrials',1000,'Confidence',50);

color = 1*ones(numel(inlierIdx),3);
color(~inlierIdx,:) = [0 0 1].*color(~inlierIdx,:);
% bb = pointCloud(points(inlierIdx,:),'Color',color);
% cc = pointCloud(points);

% pcshow(bb)
% hold on
figure
pcshow(points,color,'MarkerSize',20)

%%

[coeff,score,latent,tsquared,explained,mu] = pca(points(inlierIdx,:));
pp = points(inlierIdx,:)-mu;
coeff = coeff(:,[2 3 1]);
ppp = pp/coeff;
% ppp = ppp';
scatter3(ppp(:,1),ppp(:,2),ppp(:,3),'rx')

%%

x0 = modelRANSAC.point(1);
y0 = modelRANSAC.point(2);
z0 = modelRANSAC.point(3);
coe = modelRANSAC.coe;
tt = linspace(-1,1,50);
xx = x0+tt*coe(1);
yy = y0+tt*coe(2);
zz = z0+tt*coe(3);
hold on
scatter3(xx,yy,zz,20,'ro')

%% ------ RANSAC ------

% sampleSize = 20; % number of points to sample per trial
% maxDistance = 0.01; % max allowable distance for inliers
% 
% fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
% evalLineFcn = ...   % distance evaluation function
%   @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);
% 
% [modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn, ...
%   sampleSize,maxDistance,'MaxNumTrials',10000,'Confidence',80);
% 
% color = 1*ones(numel(inlierIdx),3);
% color(~inlierIdx,:) = [0 0.2 0.3].*color(~inlierIdx,:);
% % bb = pointCloud(points(inlierIdx,:),'Color',color);
% % cc = pointCloud(points);
% 
% % pcshow(bb)
% % hold on
% pcshow(points_ori,color,'MarkerSize',20)

%%

function mdl = fitmdl(points)

X = points(:,1:2);
y = points(:,3);
mdl.coe = X\y;

end

function y_hat = mdl_eval(mdl,points)

X = points(:,1:2);
y_hat = X*mdl.coe;

end

function mdl = fit3dline(points)

X = points;
N = size(X,1);
% Find line of best fit (in least-squares sense) through X
% -------------------------------------------------------------------------
X_ave=mean(X,1);            % mean; line of best fit will pass through this point  
dX=bsxfun(@minus,X,X_ave);  % residuals
C=(dX'*dX)/(N-1);           % variance-covariance matrix of X
[R,D]=svd(C,0);             % singular value decomposition of C; C=R*D*R'
mdl.coe = R(:,1);
mdl.point = X_ave';

end

function d = dis2line(mdl,points)

X = points-mdl.point';
normX = sqrt(diag(X*X'));
theta = acos((X*mdl.coe)./normX);
d = abs(normX.*sin(theta));

end

