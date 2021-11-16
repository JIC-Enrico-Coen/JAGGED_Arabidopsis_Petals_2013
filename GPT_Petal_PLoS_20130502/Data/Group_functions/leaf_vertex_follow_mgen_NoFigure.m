function  m=leaf_vertex_follow_mgen_NoFigure(m,varargin) %realtime,RegionLabels,Morphogens,start_figno)
%function  m=leaf_vertex_set_monitor(m,realtime,RegionLabels,Morphogens,start_figno)
%monitor morphogen levels at a set of vertices.
%Several morphogens can be monitored simultaneously.
%The results can be recorded and an output in opne figure given over
%time. The pattern is saved automatically in the model directory.

%m, mesh
%realtime, time starting at 0
%RegionLabels, vertices to be monitored as designated by cell array of strings, i.e. region labels
%Morphogens, cell array of strings, i.e. uppercase morphogen names to
%   be monitored
%start_figno, first figure number ... one figure per morphogen
%
%e.g.
%     leaf_vertex_set_monitor(m,'RealTime',realtime,'ZeroTime',zerotime,...
%         'REGIONLABELS',{'MEDIAL','MEDIAL'},...
%         'MORPHOGENS',{'KPAR','APAR'},...
%         'FigNum',4,...
%         'ENDTIME', 200,'RECORDINTERVAL', 10););
% the  'REGIONLABELS' can be 'SEAM' in which case m.seam is used
% ENDTIME specifies when the last reading should be taken
% RECORDINTERVAL specifies the interval between the readings, this is
% measured in realtime not in number of steps.
%   Topics: Simulation.

if isempty(m), return; end
if length(varargin)<3
    error('leaf_vertex_set_monitor: insufficient arguments');
end
zerotime=0;
profilesflag=true;
allpointsflag=false;
for i=1:2:length(varargin)
    name=upper(varargin{i});
    arg=varargin{i+1};
    switch name
        case 'REALTIME'
            realtime=arg;
        case 'REGIONLABELS'
            if ~iscell(arg)
                error([name,' should be a cell array']);
            end
            RegionLabels=upper(arg);
        case 'MORPHOGENS'
            if ~iscell(arg)
                error([name,' should be a cell array']);
            end
            Morphogens=upper(arg);
        case 'PROFILES'
            profilesflag=arg;
        case 'ALLPOINTS'
            allpointsflag=arg;
        case 'ENDTIME'
            endtime = arg;
        case 'RECORDINTERVAL'
            interv = arg;
    end
end
N=length(RegionLabels);
if N~=length(Morphogens)
    error('monitor called with unbalanced arguments');
end
if realtime<=zerotime*1.00001
    firsttimeflag=true;
else
    firsttimeflag=false;
end

for i=1:N
    % first identify the vertices and their positions along a line
    regionlabel=upper(RegionLabels{i});
    if strcmp(regionlabel,'SEAM')
        if firsttimeflag || allpointsflag
            vertex_set=unique(m.edgeends(m.seams==1,:));
        end
        monitor1_p=zeros(size(m.morphogens,1));
        monitor1_p(vertex_set)=1;
    else
        [monitor1_i,monitor1_p,monitor1_a,monitor1_l] = getMgenLevels( m, regionlabel);
        if firsttimeflag || allpointsflag
            vertex_set=find(monitor1_l>=0.99);
            
        end
    end
    numvert=length(vertex_set);
    if ~allpointsflag && ~firsttimeflag
        list_order=1:numvert;
        error_counter(i)=0;
    else
        if exist('list')
            clear list
        end
        for k=1:numvert
            list{k}=[];
            [ii,jj]=find(m.edgeends==vertex_set(k));
            % ii indexes all edges that include this vertex
            % how many of the other vertices in vertex_sets
            % are in these edges
            for kk=1:length(ii)
                edge=m.edgeends(ii(kk),:); % contains current vertex
                % and neighbour, is this one in the list of other vertices
                nn=intersect(edge,vertex_set);
                list{k}=unique([list{k},nn]);
            end
            startend(k)=length(list{k});
        end
        % pick one end point
        if isempty(startend)
            fprintf(1,'leaf_vertex_set_monitor  %s',Morphogens)
            error('the profile marker either does not have two ends');
        end
        endpoints=find(startend==2);
        if isempty(endpoints)
            fprintf(1,'leaf_vertex_set_monitor ');
            for i=1:length(Morphogens)
                fprintf(1,' %s ',Morphogens{i});
            end
            fprintf(1,'startend = %d ',startend)
            error('cannot find start of morphogen profile')
        end
        endpoint=endpoints(1);
        
        % edge described in list{endpoint} must also exist in one other
        % list
        % for k=1:length(list),disp(list{k}),end
        list_order=endpoint;
        ko=1;
        ed=list{list_order(ko)}; % first current edge, look for next one
        error_counter(i)=0;
        while (length(list_order)<numvert) & (error_counter<length(m.edgeends))
            for k=1:numvert % look through list
                if ~any(k==list_order) % for any edge not already used
                    eds=list{k}; % this might contain a matching edge
                    if length(intersect(ed,eds))==2
                        ko=ko+1; % yes it does
                        list_order(ko)=k; % so add it to the list_order
                        % next edge is this
                        new_node=setdiff(eds,ed);
                        % with one of the others
                        possible1=[new_node,ed(1)];
                        possible2=[new_node,ed(2)];
                        % so look through those remaining in list
                        for kk=1:numvert
                            if ~any(kk==list_order)
                                if length(intersect(list{kk},possible1))==2
                                    ed=possible1;
                                    break % got it
                                end
                                if length(intersect(list{kk},possible2))==2
                                    ed=possible2;
                                    break % got it
                                end
                            end
                        end
                    end
                end
            end
            error_counter(i)=error_counter(i)+1;
        end
    end
    if error_counter(i)<length(m.edgeends) && numvert==length(list_order)
        
        index=list_order(1);
        
        if ~strcmp(regionlabel,'SEAM')
            monitor1_p(vertex_set(index))=2; % mark the origin of the line
            m.morphogens(:,monitor1_i) = monitor1_p;
        end
        % now work along the line estimating distance from origin of line
        origin=m.nodes(vertex_set(list_order(1)),:);
        y=zeros(numvert,1);
        total_so_far=0;
        for k=1:numvert
            point=m.nodes(vertex_set(list_order(k)),:);
            total_so_far=total_so_far+sqrt(sum((point-origin).^2));
            y(k)=total_so_far;
            origin=point;
        end
        Y{i}=y;
        
        % next set up the time axis
        x=repmat(realtime,1,numvert);
        X{i}=x;
        % finally identify the morphogen levels to be plotted on the ordinate
        morphogen=upper(Morphogens{i});
        [morph_i,morph_p,morph_a,morph_l] = getMgenLevels( m, morphogen );
        z=morph_l(vertex_set(list_order)); %data(i).plotz(1:data(i).index,:);
        z(isnan(z))=0;
        Z{i}=z;
        
    end
end
if all(error_counter<length(m.edgeends))&& numvert==length(list_order)
    
    for i=1:length(Z)
        y=Y{i};
        z=Z{i};
        
        Stime = m.globalProps.timestep;
        if mod(rounddec(m.globalDynamicProps.currenttime, 3),(Stime*interv))==0
            % get the right time of recording
            %pos = round(1/Stime/interv*m.globalDynamicProps.currenttime);
            if i == 1
                m.userdata.clock = m.userdata.clock+1;
                pos = m.userdata.clock;
                %m.userdata.gradients.petiole(pos) =  max(m.nodes(m.userdata.petlength,2))-min(m.nodes(m.userdata.petlength,2));
                % time of recording
                m.userdata.gradients.TIME(pos) = rounddec(m.globalDynamicProps.currenttime,2);
            end
            m.userdata.gradients.name2index{i} = Morphogens{i};
            %  morphogen concentration
            m.userdata.gradients.CONC{pos+1,i} = z(:)/max([z(:);1]);
            m.userdata.gradients.CONCACT{pos+1,i} = z(:);
            % relative morphogen position at measured concentration
            m.userdata.gradients.POSREL{pos+1,i} = y(:)/max(y(:));
            % actual morphogen postion at measured concentration
            m.userdata.gradients.POSACT{pos+1,i} = y(:);
        end
    end
    % at the endtime save the recorded measurements in a file3
    if rounddec(m.globalDynamicProps.currenttime,3) == rounddec(endtime,3);
        data = m.userdata.gradients;
        path = fileparts(which(m.globalProps.modelname));
        %path = pwd;
        
         if ~exist([path, filesep,'TimeSeries'],'dir')
             mkdir(path, 'TimeSeries');
         end
        d = date;
        en = num2str(endtime); inte = num2str(interv);
        % make to morphogens recorded to string to use as file
        % name
        namestr = '';
        for k=1:N
            namestr = sprintf('%s_%s',namestr,m.userdata.gradients.name2index{k});
        end
        save([path,filesep,'TimeSeries',filesep,d,'_end',en,'_int',inte,...
            '_',namestr,'.mat'], 'data');
        
        plotConcSummary_ActualData(data);
    end
    
end
