function m = leaf_MeasureDimensions(m, varargin)

% Find the dimensions of a FEM structure by reading in the positions of the
% element nodes and calculating the distances between them. An excel output
% is generated in the folder with the measurements. 
%
% Erika Kuchen January 2010
% 
% The following code has to be inserted during the initialisation phase 
% vertex indecies = have to be the index of two or more verticies that descibe the area/ line to be measured 

%%
% if Steps(m) == 0
%    m.userdata.clock = 0;
%     m.secondlayer.cell3dcoords = [];  % is the cells are not created
%     m.userdata.celledges = 20; % the number of verteces that each sector
%     consists of. 
%     right form the start. 
%    m.userdata.longitudinal = vertex indecies;
%    m.userdata.lateral = vertex indecies; % width
%    m.userdata.petlength = vertex indecies;
%    m.userdata.petiole = vertex indecies;
%    m.userdata.tipdist = vertex indecies;
    % for thickness measurements two vertecies are required 
%     m.userdata.thickness.midmidvein_ind = 2 vertex indecies;
%     m.userdata.thickness.midside_ind =  2 vertex indecies;
%     m.userdata.thickness.basemidvein_ind =  2 vertex indecies;
%     m.userdata.thickness.tip_ind =  2 vertex indecies;
% end
% m = leaf_MeasureDimensions(m, 'ENDTIME', 150,'RECORDINTERVAL', 5);

% recordinterval is the interval at which measurements should be taken,
% this is based on steps rather than realtime, e.g. hours. 
% endtime is the time when the output should be computed
%%
% %2D
% % top rim length before growth in thickness;
% measure = m.nodes(bordernodes,2) > 0.02;
% d = m.nodes(bordernodes(measure),:);
% 
% for i=1:size(d,1)-1
%    x =abs(d(i,1)-d(i+1,1));
%    y =abs(d(i,2)-d(i+1,2));
%    z(i) =(x^2+y^2)^(1/2);
% end
% 
% line = sum(z);
% line2 = sum(z)+0.0078*2;
% 
% % 3D
% %Logitudinal measurment to see if the leaf has the right length
% % Need to do this along a seam, this has been set up unprogramatically. 
% 
%%
% read in the specified endtime and interval time 
for i=1:2:length(varargin)
        name=upper(varargin{i});
        arg=varargin{i+1};
        switch name
            case 'ENDTIME'
                endtime = arg;
            case 'RECORDINTERVAL'
                interv = arg;
        end
end
% record at set intervals - the time step size
Stime = m.globalProps.timestep;

if mod(rounddec(m.globalDynamicProps.currenttime, 3),(Stime*interv))==0                
    
    % get the right time of recording - to place it in the right matrix
    % column 
     m.userdata.clock = m.userdata.clock+1;
     pos = m.userdata.clock;  


% record the time at each interval
m.userdata.dims.time{pos} = m.globalDynamicProps.currenttime;   

% % Get the overall concentration of kaper:
m.userdata.dims.kper{pos} = sum(m.morphogens(:,2))/length(m.morphogens(:,2)); % average kper
m.userdata.dims.kpar{pos} = sum(m.morphogens(:,1))/length(m.morphogens(:,1));
m.userdata.dims.maxkper{pos} = max(m.morphogens(:,2)); % max kper
m.userdata.dims.maxkpar{pos} = max(m.morphogens(:,1));

% % Get overall kmin and kmax
m.userdata.dims.kmin{pos} = sum(m.plotdefaults.tensordata(:,2))/length(m.plotdefaults.tensordata(:,2));
m.userdata.dims.kmax{pos} = sum(m.plotdefaults.tensordata(:,1))/length(m.plotdefaults.tensordata(:,1));
m.userdata.dims.maxkmin{pos} = max(m.plotdefaults.tensordata(:,2));
m.userdata.dims.maxkmax{pos} = max(m.plotdefaults.tensordata(:,1));

% leaf length
d1 = m.nodes(m.userdata.length,:);
d = sortrows(d1,2);
% calculate the distance
for i=1:size(d,1)-1
   x =abs(d(i,1)-d(i+1,1));
   y =abs(d(i,2)-d(i+1,2));
   z=abs(d(i,3)-d(i+1,3));
   len(i)=(x^2+y^2+z^2)^(1/2);   
end
m.userdata.dims.wholelength{pos} = sum(len);
clear len d1 d x y z;

% leaf width
d1 = m.nodes(m.userdata.width,:);
d= sortrows(d1,1);
% calculate the distance
for i=1:size(d,1)-1
   x =abs(d(i,1)-d(i+1,1));
   y =abs(d(i,2)-d(i+1,2));
   z=abs(d(i,3)-d(i+1,3));
   len(i)=(x^2+y^2+z^2)^(1/2);
end
m.userdata.dims.width{pos} = sum(len);
clear len d1 d x y z;

% petiole length
%d1 = m.nodes(m.userdata.petlength,:);
%d = sortrows(d1,2);
% calculate the distance
%for i=1:size(d,1)-1
 %  x =abs(d(i,1)-d(i+1,1));
  % y =abs(d(i,2)-d(i+1,2));
  % z=abs(d(i,3)-d(i+1,3));
  % len(i)=(x^2+y^2+z^2)^(1/2);
%end
%m.userdata.dims.petlength{pos} = sum(len);
%clear len d1 d x y z;

% petiole width
d1 = m.nodes(m.userdata.petwidth,:);
d = sortrows(d1,1);
% calculate the distance
for i=1:size(d,1)-1
   x =abs(d(i,1)-d(i+1,1));
   y =abs(d(i,2)-d(i+1,2));
   z=abs(d(i,3)-d(i+1,3));
   len(i)=(x^2+y^2+z^2)^(1/2);
end
m.userdata.dims.petwid{pos} = sum(len);
clear len d1 d x y z;

% distance from widest point to the tip
%d1 = m.nodes(m.userdata.tipdist,:);
%d = sortrows(d1,2);
% calculate the distance
%for i=1:size(d,1)-1
 %  x =abs(d(i,1)-d(i+1,1));
  % y =abs(d(i,2)-d(i+1,2));
  % z=abs(d(i,3)-d(i+1,3));
   %len(i)=(x^2+y^2+z^2)^(1/2);
%end
%m.userdata.dims.tipdist{pos} = sum(len);
%clear len d1 d x y z;

%Thickness
%t(1:2,:) = m.prismnodes(m.userdata.thickness.midmidvein_ind-1:m.userdata.thickness.midmidvein_ind,:); 
%t(3:4,:)= m.prismnodes(m.userdata.thickness.midside_ind-1:m.userdata.thickness.midside_ind,:); 
%t(5:6,:) = m.prismnodes(m.userdata.thickness.basemidvein_ind-1:m.userdata.thickness.basemidvein_ind,:); 
%t(7:8,:) = m.prismnodes(m.userdata.thickness.tip_ind-1:m.userdata.thickness.tip_ind,:); 
%for i=1:2:size(t,1)
 %  x =abs(t(i,1)-t(i+1,1));
  % y =abs(t(i,2)-t(i+1,2));
  % z =abs(t(i,3)-t(i+1,3));
  % thick(i)=(x^2+y^2+z^2)^(1/2);
%end
%m.userdata.dims.mmv{pos} = thick(1); % middle midvein
%m.userdata.dims.ms{pos} = thick(3); % middle side
%m.userdata.dims.bmv{pos} = thick(5); % base midvein
%m.userdata.dims.tip{pos} = thick(7); % tip
%clear thick t x y z i;

% sector area sizes at each step
j = 0;
for i=1:m.userdata.celledges:size(m.secondlayer.cell3dcoords)
    j = j+1;
    coords = m.secondlayer.cell3dcoords(i:(i+m.userdata.celledges-1),:);
    % get the coordinates of each sector and calculate the area
    % does not work if the coordinates are not connected,i.e. 1 is also n. 
    m.userdata.sectorsize{j,pos} = abs(sum((coords(1,1).*coords(1,2)+coords(1:end,1).*coords([2:end,1],2)) - (coords([2:end,1],1).*coords(1:end,2)+coords(1,1).*coords(1,2)))/2);
end 
     
end

%%
% save to excel and as an mat file
if rounddec(m.globalDynamicProps.currenttime,3) == rounddec(endtime,3); 
   data = m.userdata.dims;               
   path = fileparts(which(m.globalProps.modelname));
   if ~exist([path, filesep,'Dimensions'],'dir')
       mkdir(path, 'Dimensions');
   end
   d = date;
   en = num2str(endtime); inte = num2str(interv);
   h = clock;
   % make to morphogens recorded to string to use as file
   % name
   save([path,filesep,'Dimensions',filesep,d,'_end',en,'_int',inte,
       'time',num2str(h(4)),num2str(h(5)),'.mat'], 'data');
      
s = 1;
t = 1;
u = 2;
sheet{1, t} = 'Location';
sheet{1, t+1} = 'Time in Units';

sheet{s+1, t} = 'Average Kpar';
sheet{s+2, t} = 'Average Kper';
sheet{s+3, t} = 'Average Kmax';
sheet{s+4, t} = 'Average Kmin';
sheet{s+5, t} = 'Whole Length';
sheet{s+6, t} = 'Whole Width';
%sheet{s+7, t} = 'Petiole Length';
sheet{s+8, t} = 'Petiole Width';
%sheet{s+9, t} = 'Tip to Width Point';
%sheet{s+10, t} = 'Thickness Middle Side';
%sheet{s+11, t} = 'Thickness Middle Midvein';
%sheet{s+12, t} = 'Thickness Base Midvein';
%sheet{s+13, t} = 'Thickness Tip';
sheet{s+14, t} = 'Max Kpar';
sheet{s+15, t} = 'Max Kper';
sheet{s+16, t} = 'Max Kmax';
sheet{s+17, t} = 'Max Kmin';

% put all the variable number of sectors in
for j=1:size(m.userdata.sectorsize,1)    
    sheet{s+17+j,t} = sprintf('Sector %d', j);
end

for j=1:size(data.time,2)
sheet{s, u+j} = data.time{j};     
sheet{s+1, u+j} = data.kpar{j}; 
sheet{s+2, u+j} = data.kper{j}; 
sheet{s+3, u+j} = data.kmax{j}; 
sheet{s+4, u+j} = data.kmin{j}; 
sheet{s+5, u+j} = data.wholelength{j}; 
sheet{s+6, u+j} = data.width{j}; 
sheet{s+7, u+j} = data.petlength{j}; 
sheet{s+8, u+j} = data.petwid{j}; 
sheet{s+9, u+j} = data.tipdist{j}; 
sheet{s+10, u+j} = data.ms{j}; 
sheet{s+11, u+j} = data.mmv{j}; 
sheet{s+12, u+j} = data.bmv{j}; 
sheet{s+13, u+j} = data.tip{j}; 
sheet{s+14, u+j} = data.maxkpar{j}; 
sheet{s+15, u+j} = data.maxkper{j}; 
sheet{s+16, u+j} = data.maxkmax{j}; 
sheet{s+17, u+j} = data.maxkmin{j};
for i=1:size(m.userdata.sectorsize,1)    
    sheet{s+17+i,u+j} = m.userdata.sectorsize{i,j};
end
end

h = clock;
xlswrite([path,filesep,'Dimensions',filesep,d,'_end',en,'_int',inte,'_time',num2str(h(4)), num2str(h(5)),...
   '.xls'],sheet);

% plot graph with the sector identities
leaf_calcsectorarea(m)
end