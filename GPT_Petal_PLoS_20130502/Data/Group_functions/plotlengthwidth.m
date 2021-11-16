function plotlengthwidth(varargin)

% Erika 
% plot the leaf canvas length and width over time together with the real
% leaf length and width over time. 
% For this, you need to have the model path in the current working directory. Also,
% currently, the file LeafLengthWidth.xls, which contains the real leaf
% measurements, needs to be in that directory. 

% load in data length-width file
%data = xlsread('\\Nbi-cfs1\coengroup\current lab members\Susana\PetalWidthLength.xlsx');
data = xlsread('D:\ab\Matlab stuff\Growth models\models\PetalWidthLength.xlsx');

if isempty(varargin)
% load in model length-wdith file 
f = uigetfile([pwd, filesep,'Dimensions'], '*.mat');
model = load([pwd,filesep,'Dimensions',filesep,f]);
model = model.data;
else
    model = varargin{1};
end
modeltime = model.time;
modeldata = model.measure;

maxwidth = max(modeldata(:,2:7),[],2); %choose model excel file columns to compare the width


figure
plot(modeltime,1000*modeldata(:,1),'-r','LineWidth',3); hold on; % model length
%plot(modeltime,modeldata(:,2),'-b','LineWidth',3); % model lamina 
plot(modeltime,1000*maxwidth,'-k','LineWidth',3); % model width
plot(data(:,1),data(:,2),'.k','MarkerSize',15); % data width
plot(data(:,1),data(:,3),'.r','MarkerSize',15); % data length
box off;
set(gca,'TickLength',[0 0],'LineWidth',1,'FontSize',14);
xlabel('\fontsize{15}Time (h)');
ylabel('\fontsize{15}Dimensions (\mum)');

legend('Model Length','Model Width','Real Width', 'Real Length');
legend('Location','NorthWest'); legend('boxoff');

% saving the figure
    % set(gcf, 'PaperPositionMode', 'auto');
    % print(gcf,'-r200','-dtiff','Pgraddiff_late3sudden.tif');
