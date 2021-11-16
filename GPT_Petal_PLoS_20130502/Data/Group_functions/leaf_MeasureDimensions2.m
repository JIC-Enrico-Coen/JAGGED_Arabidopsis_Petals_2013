function m = leaf_MeasureDimensions2(m, varargin)
%m = leaf_MeasureDimensions2(m, varargin)
% Find the dimensions of a FEM structure by reading in the positions of the
% element nodes and calculating the distances between them. An excel output
% is generated in the folder with the measurements.
%
% Erika Kuchen February 2011

% The following code has to be inserted in the initialisation section with
% the otions to measure the length of a line, thickness or sector size
% changes over time. 

%% Basic Usage

% m = leaf_MeasureDimensions2(m, 'HEADINGS',{'Line1','Line2'},...
%     'LINES',{vertexset,vertexset},'LENWID',[2 2],...
%     'ENDTIME',100, 'RECORDINTERVAL',1:10:100);

% HEADINGS should be a description of the line along which you are trying
% to measure for each line, e.g. canvas length or canvas width. This will
% appear as a heading in the excel output. 

% LINES this should be a set of vertex indexes that lie on the line to be
% measured.

% LENWID asks if the line could be described as lying parallel to the x-axis or
% y-axis. Enter 1 for the x-axis and 2 for y-axis. 

% Note that all of these entries have to have the same length, i.e. the
% same number of entries. 

% RECORDINTERVAL the times at which the size should be recorded

% ENTIME is the time at which the measurements should be saved as a mat and
% excel file. 

%% with Sectors 

% m.userdata.celledges = 5;
% m = leaf_MeasureDimensions_better(m, 'HEADINGS',{'Line1','Line2'},'LINES',...
%     {m.userdata.longitudinal,m.userdata.longitudinal},'LENWID',[2 2],...
%     'CELLS',1,'ENDTIME',76, 'RECORDINTERVAL', 1:100);

% CELLS either takes the argument 1 or 0.
% If 0 sector sizes will not be calculated. 
% If 1 sector size will be calculated and a positional output image will be
% plotted. In this case you also need to supply the number of edges you use to make a sector by
% defining m.userdata.celledges = ?;


% OR without any line measurements
% m.userdata.celledges = 5;
% m = leaf_MeasureDimensions_better(m, 'CELLS',1,'ENDTIME',76, 'RECORDINTERVAL', 1:100);

%% With Measuring Thickness

% m = leaf_MeasureDimensions_better(m, 'HEADINGS',{'Line1','Thick1','Line2'},'LINES',...
%     {m.userdata.longitudinal,6,m.userdata.longitudinal},'LENWID',[2 0 2],...
%     'ENDTIME',76, 'RECORDINTERVAL', 1:100);

% Entries for thickness work the same as normal line measurements. You will
% also need to define a name under HEADINGS. LINES takes only one vertex
% entry at the point where thickness should be measured. LENWID simply
% takes 0 at the point where thickness should be measured. 



%% Start of the Code


for i=1:2:length(varargin)
    name=upper(varargin{i});
    arg=varargin{i+1};
    switch name
        case 'HEADINGS'
            titles = arg;
        case 'LINES'
            measures = arg;
        case 'LENWID'
            direction = arg;
        case 'CELLS'
            cells = arg;
        case 'ENDTIME'
            endtime = arg;
        case 'RECORDINTERVAL'
            interv = arg;
    end
end

    
    if ~exist('cells','var')
        cells =0;
    end

% record at set intervals
%Stime = m.globalProps.timestep;
%mod(rounddec(m.globalDynamicProps.currenttime, 3),(Stime*interv))==0

if sum(ismember(interv,rounddec(m.globalDynamicProps.currenttime, 3)))
    
   % ~isempty(intersect(rounddec(m.globalDynamicProps.currenttime, 3),interv))
   
    % get the right time of recording
    %m.userdata.clock = m.userdata.clock+1;
    %pos = m.userdata.clock;
    pos = find(interv == rounddec(m.globalDynamicProps.currenttime, 3));
    
    % Record the current time
    m.userdata.dims.time(pos) = m.globalDynamicProps.currenttime;
    
    if exist('titles','var')
    m.userdata.dims.titles = titles;
    
    % measurements ofl lines 
        mindex = find(direction>0);

    for j = mindex
        
        d = m.nodes(measures{j},:);
        %d = sortrows(d1,direction(j));
        % calculate the distance
        for i=1:size(d,1)-1
            x =abs(d(i,1)-d(i+1,1));
            y =abs(d(i,2)-d(i+1,2));
            z=abs(d(i,3)-d(i+1,3));
            len(i)=(x^2+y^2+z^2)^(1/2);
        end
        m.userdata.dims.measure(pos,j) = sum(len);
        clear len d1 d x y z;
    end
    
    % thickness
        for j = find(direction==0)
            % %Thickness
            t = m.prismnodes(measures{j}-1:measures{j},:); % need to save it seperately
            for i=1:2:size(t,1)
                x =abs(t(i,1)-t(i+1,1));
                y =abs(t(i,2)-t(i+1,2));
                z =abs(t(i,3)-t(i+1,3));
                thick(i)=(x^2+y^2+z^2)^(1/2);
            end
            m.userdata.dims.measure(pos,j) = sum(thick);
            clear thick t x y z i;
        end

    
    end
    
    % sectors 
    if cells == 1 && isfield(m.secondlayer,'cell3dcoords')
        
            % sector area sizes at each step
            j = 0;
            for i=1:m.userdata.celledges:size(m.secondlayer.cell3dcoords,1)
                j = j+1;
                
                coords = m.secondlayer.cell3dcoords(i:(i+m.userdata.celledges-1),:);
                % get the coordinates of each sector and calculate the area
                % does not work if the coordinates are not connected,i.e. 1 is also n.
                m.userdata.dims.sectorsize(pos,j) = abs(sum((coords(1,1).*coords(1,2)+coords(1:end,1).*coords([2:end,1],2)) - (coords([2:end,1],1).*coords(1:end,2)+coords(1,1).*coords(1,2)))/2);
            end
    end
end

% save to excel
if rounddec(m.globalDynamicProps.currenttime,3) == rounddec(endtime,3);
    data = m.userdata.dims;
    path = fileparts(which(m.globalProps.modelname));
    %path = pwd;
    if ~exist([path, filesep,'Dimensions'],'dir')
        mkdir(path, 'Dimensions');
    end
    d = date;
    en = num2str(endtime);
    ran = randi(10000, 1);
    % make to morphogens recorded to string to use as file
    % name
    save([path,filesep,'Dimensions',filesep,d,'_end',en,'_',num2str(ran),...
        '.mat'], 'data');
    
    
    sheet{1, 1} = 'Time in Units';
    
    if exist('titles','var')
        for i=1:size(titles,2)
            sheet{1,i+1} = titles{i};
        end
    else
        i = 0;
    end
    
    if cells == 1
        for j=1:size(m.userdata.dims.sectorsize,2)
            sheet{1,i+j+2} = sprintf('Sector %d', j);
        end
    end
    
    for j=1:size(data.time,2)
        sheet{j+1,1} = data.time(j);
        if exist('titles','var')
            for i=1:size(titles,2)
                sheet{j+1, i+1} = data.measure(j,i);
            end
        else
            i = 0;
        end
        
        if cells == 1
            for k=1:size(data.sectorsize,2)
                sheet{j+1,i+2+k} = data.sectorsize(j,k);
            end
        end
    end
    
    
    xlswrite([path,filesep,'Dimensions',filesep,d,'_end',en,'_',num2str(ran),...
        '.xls'],sheet);
    
    if cells == 1
    % plot graph with the sector identities
    leaf_calcsectorarea(m)
    end
end