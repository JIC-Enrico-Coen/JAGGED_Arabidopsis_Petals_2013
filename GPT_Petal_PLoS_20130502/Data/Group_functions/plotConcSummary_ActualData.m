function plotConcSummary_ActualData(varargin)

% Erika 30.01.2009
% This tools takes the collected concentration profiles over time from
% leaf_vertex_follow_mgen.m to create a separate graph for each morphogen
% with the Distance from the base on the x-axis and the
% Concentration for morphogens or growth rate in cases of kpar and kper on
% the y-axis. In leaf_vertex_follow_mgen.m an indivisdual starttime does
% not seem settable, so that records are always taken right from the start of the
% simulation.
% By changing the commenting below it is also possible to plot the growth
% rates along the midline over time together with the growth data. Note
% that this part requires quite a bit of user input below.

% leaf_vertex_follow_mgen.m saves a mat-file of the morphogen concentration
% profiles in a subdirectory of the project directory called 'TimeSeries'.
% This function will automatically be called by leaf_vertex_follow_mgen.m,
% but can also open an existing mat-file, when in the project directory.

% loading in the data
if isempty(varargin)
    f = uigetfile([pwd, filesep,'TimeSeries'], '*.mat');
    dat = load([pwd,filesep,'TimeSeries',filesep,f]);
    data = dat.data;
else
    data = varargin{1};
end

%data = data.gradients;


% decide whether to plot the model output individually for each morphogen,
% or plot kml together with the real leaf data (or both).
plotallmorphogens(data) % all morphogens individually along the selected line
%plotmidlinegrowth(data) % plot the midline growth rates
%thresholddistance(data) % Determine distance of a certain growth threshold or morphogen concetration from the petiole-lamina boundary.
end

%% Plot each morphogen on an individual graph with lots of time points.

function plotallmorphogens(data)
% plot each morphogen that was recorded on a seperate figure.

% you may have recorded the morphogen concentrations/ growth rates at many
% time points, but only wish to plot a few. Determine the interval at which
% you want to plot with "plotint".
plotint = 3; %number of group of time points


for num = 1:size(data.CONCACT,2)
    figure
    for i=4:plotint:size(data.CONCACT,1)  % you may wish to adjust the first time point at which things are plotted. For the leaf model, the first few time points are empty.
    %2 first times asre 0, 3 is 12h, 4 is 24h    
        col = rand(1,3); % get a new colour.
        
        % plot the actual canvas distance
        plot(data.POSACT{i,num},data.CONCACT{i,num},'-' ,'Color', col,'LineWidth',2); hold on;
        % plot the relative canvas distance.
        %plot(data.POSACT{i,num}/max(data.POSACT{size(data.CONCACT,1),num}
        %),data.CONCACT{i,num},'-' ,'Color', col,'LineWidth',2); hold on;
        
        title(sprintf('%s', data.name2index{num}),'interpreter','none');
        box off; set(gca,'TickLength',[0 0]);
        xlabel('\fontsize{12}Actual position along canvas length');
        ylabel('\fontsize{12}Growth rate/ Morphogen levels'); % this is likely to differ between plotting growth rate and morphogen levels.
    end
    
end

end

%% Plotting Kml, model and data
function plotmidlinegrowth(data)

%%% this loads in the data and plots it. Note that you need to be in the
%%% directory of Karen's data Excel files (currently this is \\Nbi-fas1a\coengroup\Arabidopsis\Tracking leaf growth\OPT\Tracking_Trichomes_QtVolViewer\CSV_for_Erika)
%%% however, it is likely that the name of the plotting function will change
%%% and also the directory and may need to be updated. Also note that the
%%% subplots etc may change with a new way of plotting the leaf data,
%%% which will need to be adjusted.
hfig = TrichomeMidline_T0T1;

figure(hfig);

%%% requires USER INPUT!!
%%% The data at model can be compared at different time points. data.TIME
%%% contains all the time points at which the model output was measured. From
%%% these the time points to be compared to the model can be selected using
%%% the variable "phase" below. phase contains the indexes of the time points
%%% in data.TIME. It turns out that different model simulation fit better at
%%% slightly different time points.
%%% phase = [3 6 9 12 14 16 19 22 22 25 25 28]; % modelkrns: LATESWITCH model
%%% phase = [3 6 9 12 14 16 18 23 23 26 26 28]; % modelkrns: PGRADBOTH model
%%% phase = [3 6 9 12 14 16 18 23 24 28 28 33]; % modelkrns: LATESWITCHEXPO model (glate 0.009)
%%% phase = [3 6 9 12 14 16 20 23 25 31 30 37]; % modelkrns: LATESWITCHEXPO model (glate 0.01)
phase = [3 6 9 12 15 16 20 23 26 31 31 35]; % modelkrns: PGRADDIFF model

%%% requires USER INPUT!
%%% State the index of Kml in name2index. For instance, in my models kml is
%%% called KML and in the list of name2index: {'ID_PGRAD'  'KAREA'  'KML'}, would be number 3, num = 3;.
num = 3;

%%% requires USER INPUT!
%%% Previously this variable was stored with the mat-file. Since, the real leaf data
%%% is aligned at the petiol-lamina boundary, we want to do the same with
%%% the model data. The parameter "perpos" records how many nodes
%%% along the midline make up the petiole, the position of which can then
%%% be subtracted to get the petiole-lamina boundary. This parameters
%%% changes with the discretisation of the mesh.
petpos = 11;


%%% plotting routine
phase = phase+1; % first entry for data.CONC etc. seems to be empty. Add one to make the same at data.TIME
subplace = [1 4 7 10 2 5 8 11 3 6 9 12]; % order in which the subplots should be populated.
l=0;
%%% for each index in phase plot the corresponding growth profile over canvas
%%% length.
for i=phase
    l = l+1;
    subplot(4,3,subplace(l))
    plot(1000*(data.POSACT{i,num}-data.POSACT{i,num}(petpos)),100*data.CONCACT{i,petpos},'-' ,'Color', [0 0 0],'LineWidth',3); hold on;
    title(sprintf('Model time %d days',ceil(data.TIME(i-1)/24))) % give model time to the nearest day.
end

%%% save the figure if applicable.
%  set(gcf, 'PaperPositionMode', 'auto');
%  print(gcf,'-r300','-dtiff','ModelwithData.tif');

end

%% Function to determine distance of a certain growth threshold or
%% morphogen concetration from the petiole-lamina boundary.

function thresholddistance(data)

fig1 = figure;
fig2 = figure;
fig3 = figure;

% set up thresholds to be measured
pgradthreshlevels = 0.005;
arealgthresholdlevels = 0.005;
kmlthreshlevels = 0.005;

% determine positions of petiole-lamina boundary
% first element in the lamina, along the midline. 
%startFE = data.startFE;
startFE=1;

for n=2:size(data.CONCACT,1)
pgradthresh = max(find(data.CONCACT{n,1}>pgradthreshlevels));
growththresh = max(find(data.CONCACT{n,6}>arealgthresholdlevels));
kmlthresh = max(find(data.CONCACT{n,3}>kmlthreshlevels));
if isempty(pgradthresh)
   pgradthresh = startFE;
end
if isempty(growththresh)
   growththresh = startFE;
end 
if isempty(kmlthresh)
   kmlthresh = startFE;
end 

figure(fig1);
plot(data.TIME(n-1),data.POSACT{n,1}(pgradthresh)-data.POSACT{n,1}(startFE),'+k'); hold on;
plot(data.TIME(n-1),max(data.POSACT{n,1}-data.POSACT{n,1}(startFE)),'+r'); hold on;

figure(fig2);
plot(data.TIME(n-1),data.POSACT{n,2}(growththresh)-data.POSACT{n,2}(startFE),'+k'); hold on;
plot(data.TIME(n-1),max(data.POSACT{n,1}-data.POSACT{n,1}(startFE)),'+r'); hold on;

figure(fig3);
plot(data.TIME(n-1),data.POSACT{n,3}(kmlthresh)-data.POSACT{n,2}(startFE),'+k'); hold on;
plot(data.TIME(n-1),max(data.POSACT{n,1}-data.POSACT{n,1}(startFE)),'+r'); hold on;
end

figure(fig1);
xlabel('Time (hours)'); ylabel(sprintf('PGRAD threshold levels (%1.2f)\n as distance from petiole-lamina boundary (mm)',pgradthreshlevels));
legend('PGRAD threshold position','Canvas lamina length');
%legend('Location','NorthWest'); legend('boxoff');

figure(fig2);
xlabel('Time (hours)'); ylabel(sprintf('Karea threshold (%1.3f /h)\n as distance from petiole-lamina boundary (mm)',arealgthresholdlevels));
legend('Karea threshold position','Canvas lamina length');
%legend('Location','NorthWest'); legend('boxoff');

figure(fig3);
xlabel('Time (hours)'); ylabel(sprintf('Kml threshold (%1.3f /h)\n as distance from petiole-lamina boundary (mm)',kmlthreshlevels));
legend('Location','NorthWest'); 
legend('Kml threshold position','Canvas lamina length');
end