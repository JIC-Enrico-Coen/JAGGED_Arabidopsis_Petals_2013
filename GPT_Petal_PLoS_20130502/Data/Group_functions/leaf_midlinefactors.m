function [Gsmid, Gsnorm] =  leaf_midlinefactors(m, acamount,T,pos)

%pos = m.userdata.longitudinal;
%pos = find(m); % need to change this!!!!!!!!!!!!!!!!!!!!!!!!!

% get the elements that line the seamed midline
if ~isempty(pos)
    
    for i=pos'
        [row,col] = find(m.tricellvxs == i);
        FEpos{find(pos == i)} = row;
    end
    
    % calculate the unit vector used for the rotation: Pick 2 different elements
    % and see if line still in action.
    n1 = m.nodes(pos(4),:);
    n2 = m.nodes(pos(7),:);
    unitmidline = ((n1-n2)/norm(n1-n2));
    
    % for each element
    for i=1:length(pos)
        
        smidline = zeros(1,size(FEpos{i},1));
        smidnorm = smidline;
        
        for j=1:size(FEpos{i},1)
            % calculate the growth rate in the direction of the midline:
            smax(j) = acamount(FEpos{i}(j),1);
            smin(j) = acamount(FEpos{i}(j),2);
            % find the orthogonal plane to calculate the rates normal to the midline
            tripoints = m.tricellvxs(FEpos{i}(j),:);
            p1 = m.nodes(tripoints(1),:);
            p2 = m.nodes(tripoints(2),:);
            p3 = m.nodes(tripoints(3),:);
            b = cross((p2 - p1),(p3-p1));
            b = b/norm(b);
            
            smidline(j) = linear_growth(unitmidline, T{FEpos{i}(j)});
            smidnorm(j) = linear_growth(cross(b,unitmidline),T{FEpos{i}(j)});
            
        end
        
        smid(i) = sum(smidline(j))/length(smidline(j));
        snorm(i) = sum(smidnorm(j))/length(smidnorm(j));
        %smax(i) = sum(smax(j))/length(smax(j));
        %smax(i) = sum(smin(j))/length(smin(j));
        
        % for weighted average:
        %m.cellareas = area of cells
        %m.celldata
    end
    
    Gsmid = zeros(size(m.nodes,1),1);
    Gsnorm = Gsmid;
    Gsmid(pos) = smid;
    Gsnorm(pos) = snorm;
    
else
    Gsmid = 0; Gsnorm =0;
end


end

function [k] = linear_growth(u,Tt)

k = dot(u*Tt',u);
% k = dot(k,u);
end

