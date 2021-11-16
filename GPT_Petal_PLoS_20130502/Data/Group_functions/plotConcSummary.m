function plotConcSummary(varargin)

% Erika 30.01.2009
% This tools takes the collected concentration profiles over time from
% leaf_vertex_follow_mgen.m to create a separate graph for each morphogen
% with the Distance from the base on the x-axis and the
% Concentration for morphogens or growth rate in cases of kpar and kper on
% the y-axis. In leaf_vertex_follow_mgen.m an indivisdual starttime does
% not seem settable, so that records are always taken right from the start of the
% simulation. 

% leaf_vertex_follow_mgen.m saves a mat file of the morphogen concentration
% profiles in a subdirectory of the project directory called 'TimeSeries'.
% This function will automatically be called by leaf_vertex_follow_mgen.m,
% but can also open an existing .mat file, when in the project directory. 

if isempty(varargin)
f = uigetfile([pwd, filesep,'TimeSeries'], '*.mat');
dat = load([pwd,filesep,'TimeSeries',filesep,f]);
data = dat.data;
else 
  data = varargin{1};  
end

data = data.gradients;
 for num = 1:size(data.CONCACT,2)
      figure
      for i=5:size(data.CONCACT,1)  
          col = rand(1,3);
          plot(data.POSACT{i,num},data.CONCACT{i,num},'o-' ,'Color', col,'LineWidth',2); hold on; % 3
          title(sprintf('%s', data.name2index{num}),'interpreter','none');
          xlabel('Position'); ylabel('Concentration');
      end
  end



