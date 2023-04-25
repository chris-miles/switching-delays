clear all
close all

load('2d_paramsweep.mat');

figure('position',[0,0,500,500]);
pppp=pcolor(X,Y, fftsmix);
pppp.LineStyle='none';

ytickvals = yticks;
yyaxis right
ylim([0.5 5])
yticks(ytickvals)
xlabel('\tau_1')
ylabel('\tau_2','Rotation',0)
colormap(flipud(viridis_white))

colorbar('westoutside')
pbaspect([4,4,1]);

hold on;
load('boundary_sweepvals.mat')


dist_thresh=1;
%scatter(scatter_all(1,:),scatter_all(3,:));
rng(2);

lines=connect_points(scatter_all, dist_thresh);
hold on;

for l = 1:length(lines)
    linel = lines{l};
    pts_on_line = scatter_all(:,linel);
    pts_on_line(:,pts_on_line(3,:)>pts_on_line(1,:))=[];
    plot(pts_on_line(1,:),pts_on_line(3,:),'color','red');
end
% hopf aroun tau1=.861 for these params
tau1hopf = .861;
xline(tau1hopf,'color','red')
yline(tau1hopf,'color','red')

scatter(3,1);

set(gca,'LineWidth',1.5)
set(gcf,'Renderer','painters')