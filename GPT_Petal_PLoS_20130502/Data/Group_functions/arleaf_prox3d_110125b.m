function m = arleaf_prox3d_110125b( m )
%m = arleaf_prox3d_110125b( m )
%   Morphogen interaction function.
%   Written at 2011-02-08 15:03:30.
%   GFtbox revision 3381, 2011-02-02 11:45:03.142931.
%   Model last saved to SVN as revision 3172, 2010-10-12 08:38:13.898165.

% The user may edit any part of this function between delimiters
% of the form "USER CODE..." and "END OF USER CODE...".  The
% delimiters themselves must not be moved, edited, deleted, or added.

    if isempty(m), return; end

    fprintf( 1, '%s found in %s\n', mfilename(), which(mfilename()) );

    if exist( 'local_setproperties' )
        m = local_setproperties( m );
    end

    realtime = m.globalDynamicProps.currenttime;

%%% USER CODE: INITIALISATION

if Steps(m) == 0
    % To make measurements of the growing leaf (uses leaf_MeasureDimensions)
    % leaf length
    a1 = find(m.nodes(:,1)>(-0.0008));
    a2 = find(m.nodes(:,1)<0.0007);
    m.userdata.longitudinal = intersect(a1, a2);
    
    % petiole length
    m.userdata.pet_ind=find(m.nodes(:,2)<(min(m.nodes(:,2))+0.012));  
    m.userdata.petlength =intersect(m.userdata.longitudinal,m.userdata.pet_ind);
    
    % width 
    b1= find(m.nodes(:,2) < (-0.01)); 
    b2 = find(m.nodes(:,2)>(-0.013)); 
    m.userdata.lateral = intersect(b1,b2);
    
    % make a straight cut about half way down the leaf (phase C)
    tip1 = find(m.nodes(:,2)>(-0.015));
    tip2 = find(m.nodes(:,2)<-0.012);
    m.userdata.tipoff = intersect(tip1, tip2);
    m.userdata.noncutpart = intersect(m.userdata.longitudinal, tip2);
    
    % seam to measure gradients along the midline
    m = leaf_addseam(m, 'nodes', m.userdata.longitudinal);
    m.userdata.clock = 0;
    
    % to set up the initial conditions
    m = leaf_setproperty( m, 'mingradient', 0, 'timestep', 0.1);
    m = leaf_plotoptions( m, 'hiresdpi', 300, 'highgradcolor',[0,0,1],'lowgradcolor',[0,0,0] );
        
    % Set up names for variant models.
    m.userdata.ranges.modelname.range = { 'LEAF', 'NARROW','MIDLINE','STRAIGHTCUT','MIDPOL'};
    m.userdata.ranges.modelname.index = 5;
    m.userdata.startingshape = {'Flat','Dshape','Primordium'};
    m.userdata.start.index = 1;
end

modelname = m.userdata.ranges.modelname.range{m.userdata.ranges.modelname.index};  % CLUSTER
startname = m.userdata.startingshape{m.userdata.start.index}; 

switch modelname
    case {'LEAF','MIDPOL'}
        if Steps(m) == 0
            m = leaf_setproperty(m, 'userpolarisation',false);
        elseif Steps(m) == 3
            m = leaf_setproperty( m, 'mingradient', 2);
        end
        
    case {'NARROW','MIDLINE'}
        if Steps(m) == 0
            m = leaf_setproperty(m, 'userpolarisation',true);
            for i=1:size(m.gradpolgrowth,1)
                m.gradpolgrowth(i,:) = [0 1 0];
            end
        end
    case 'STRAIGHTCUT'
        if Steps(m) == 0
            m = leaf_setproperty(m, 'userpolarisation',false);
        elseif Steps(m) == 3
            m = leaf_setproperty( m, 'mingradient', 2);
        elseif Steps(m) == 75
            m = leaf_addseam(m, 'nodes', m.userdata.noncutpart);
        end
        
end

if Steps(m) == 0
    switch startname
        case 'Flat'
            m = leaf_setproperty(m,'thicknessMode', 'physical');
            [m,ok] = leaf_setthickness( m, 0.06 );
        case {'Dshape','Primordium'}
            m = leaf_setproperty(m,'thicknessMode', 'direct');
            [m,ok] = leaf_setthickness( m, 0.002 );
    end
end

%m = leaf_MeasureDimensions_short(m, 'ENDTIME',432, 'RECORDINTERVAL', 12);

     m=leaf_vertex_follow_mgen_NoFigure(m,'RealTime',realtime,'ZeroTime',70,...
    'REGIONLABELS',{'SEAM', 'SEAM', 'SEAM', 'SEAM','SEAM', 'SEAM','SEAM'},...
    'MORPHOGENS',{'KAPAR', 'KAPER', 'S_PROX','KMAX','KMIN','KML','KPML'},...
    'PROFILES',false,...
    'ALLPOINTS',true, 'ENDTIME',217, ...
    'RECORDINTERVAL', 5);


% Set priorities for simultaneous plotting of multiple morphogens, if desired.
%m = leaf_plotoptions( m, 'morphogen',{'ID_RIM','ID_MIDVEIN','ID_LAM'});
%m = leaf_mgen_plotpriority( m,
%{'ID_RIM','ID_MIDVEIN','ID_LAM'},[2,1,1],[0.7,0,0]);%,[0.8,1,0.6,1]);
%%% END OF USER CODE: INITIALISATION

%%% SECTION 1: ACCESSING MORPHOGENS AND TIME.
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.

    if isempty(m), return; end

    setGlobals();
    global gNEW_KA_PAR gNEW_KA_PER gNEW_KB_PAR gNEW_KB_PER
    global gNEW_K_NOR gNEW_POLARISER gNEW_STRAINRET gNEW_ARREST
    dt = m.globalProps.timestep;
    polariser_i = gNEW_POLARISER;
    P = m.morphogens(:,polariser_i);
    [kapar_i,kapar_p,kapar_a,kapar_l] = getMgenLevels( m, 'KAPAR' );
    [kaper_i,kaper_p,kaper_a,kaper_l] = getMgenLevels( m, 'KAPER' );
    [kbpar_i,kbpar_p,kbpar_a,kbpar_l] = getMgenLevels( m, 'KBPAR' );
    [kbper_i,kbper_p,kbper_a,kbper_l] = getMgenLevels( m, 'KBPER' );
    [knor_i,knor_p,knor_a,knor_l] = getMgenLevels( m, 'KNOR' );
    [strainret_i,strainret_p,strainret_a,strainret_l] = getMgenLevels( m, 'STRAINRET' );
    [arrest_i,arrest_p,arrest_a,arrest_l] = getMgenLevels( m, 'ARREST' );
    [v_leaf_i,v_leaf_p,v_leaf_a,v_leaf_l] = getMgenLevels( m, 'V_LEAF' );
    [v_petiolelen_i,v_petiolelen_p,v_petiolelen_a,v_petiolelen_l] = getMgenLevels( m, 'V_PETIOLELEN' );
    [v_length_i,v_length_p,v_length_a,v_length_l] = getMgenLevels( m, 'V_LENGTH' );
    [v_petiolewid_i,v_petiolewid_p,v_petiolewid_a,v_petiolewid_l] = getMgenLevels( m, 'V_PETIOLEWID' );
    [v_width_i,v_width_p,v_width_a,v_width_l] = getMgenLevels( m, 'V_WIDTH' );
    [v_thickness_i,v_thickness_p,v_thickness_a,v_thickness_l] = getMgenLevels( m, 'V_THICKNESS' );
    [id_base_i,id_base_p,id_base_a,id_base_l] = getMgenLevels( m, 'ID_BASE' );
    [id_distal_i,id_distal_p,id_distal_a,id_distal_l] = getMgenLevels( m, 'ID_DISTAL' );
    [id_rim_i,id_rim_p,id_rim_a,id_rim_l] = getMgenLevels( m, 'ID_RIM' );
    [id_leaf_i,id_leaf_p,id_leaf_a,id_leaf_l] = getMgenLevels( m, 'ID_LEAF' );
    [s_prox_i,s_prox_p,s_prox_a,s_prox_l] = getMgenLevels( m, 'S_PROX' );
    [id_midvein_i,id_midvein_p,id_midvein_a,id_midvein_l] = getMgenLevels( m, 'ID_MIDVEIN' );
    [v_tip_i,v_tip_p,v_tip_a,v_tip_l] = getMgenLevels( m, 'V_TIP' );
    [kmax_i,kmax_p,kmax_a,kmax_l] = getMgenLevels( m, 'KMAX' );
    [kmin_i,kmin_p,kmin_a,kmin_l] = getMgenLevels( m, 'KMIN' );
    [id_lam_i,id_lam_p,id_lam_a,id_lam_l] = getMgenLevels( m, 'ID_LAM' );
    [id_youth_i,id_youth_p,id_youth_a,id_youth_l] = getMgenLevels( m, 'ID_YOUTH' );
    [kml_i,kml_p,kml_a,kml_l] = getMgenLevels( m, 'KML' );
    [kpml_i,kpml_p,kpml_a,kpml_l] = getMgenLevels( m, 'KPML' );
    [aniso_i,aniso_p,aniso_a,aniso_l] = getMgenLevels( m, 'ANISO' );
    [karea_i,karea_p,karea_a,karea_l] = getMgenLevels( m, 'KAREA' );
    [rot_i,rot_p,rot_a,rot_l] = getMgenLevels( m, 'ROT' );
    [id_cutorg_i,id_cutorg_p,id_cutorg_a,id_cutorg_l] = getMgenLevels( m, 'ID_CUTORG' );
    [id_bend_i,id_bend_p,id_bend_a,id_bend_l] = getMgenLevels( m, 'ID_BEND' );
    [f_seam_i,f_seam_p,f_seam_a,f_seam_l] = getMgenLevels( m, 'F_SEAM' );
    [in_tipinh_i,in_tipinh_p,in_tipinh_a,in_tipinh_l] = getMgenLevels( m, 'IN_TIPINH' );
    [id_borders_i,id_borders_p,id_borders_a,id_borders_l] = getMgenLevels( m, 'ID_BORDERS' );

% Mesh type: lobes
%            base: 0
%        cylinder: 0
%          height: 0.053
%           lobes: 1
%          radius: 0.042
%      randomness: 0
%           rings: 10
%          strips: 18
%         version: 1

%            Morphogen   Diffusion   Decay   Dilution   Mutant
%            -------------------------------------------------
%                KAPAR        ----    ----       ----     ----
%                KAPER        ----    ----       ----     ----
%                KBPAR        ----    ----       ----     ----
%                KBPER        ----    ----       ----     ----
%                 KNOR        ----    ----       ----     ----
%            POLARISER       0.001    0.01       ----     ----
%            STRAINRET        ----    ----       ----     ----
%               ARREST        ----    ----       ----     ----
%               V_LEAF        ----    ----       ----     ----
%         V_PETIOLELEN        ----    ----       ----     ----
%             V_LENGTH        ----    ----       ----     ----
%         V_PETIOLEWID        ----    ----       ----     ----
%              V_WIDTH        ----    ----       ----     ----
%          V_THICKNESS        ----    ----       ----     ----
%              ID_BASE        ----    ----       ----     ----
%            ID_DISTAL        ----    ----       ----     ----
%               ID_RIM        ----    ----       ----     ----
%              ID_LEAF        ----    ----       ----     ----
%               S_PROX        ----    ----       ----     ----
%           ID_MIDVEIN        ----    ----       ----     ----
%                V_TIP        ----    ----       ----     ----
%                 KMAX        ----    ----       ----     ----
%                 KMIN        ----    ----       ----     ----
%               ID_LAM        ----    ----       ----     ----
%             ID_YOUTH        ----    ----       ----     ----
%                  KML        ----    ----       ----     ----
%                 KPML        ----    ----       ----     ----
%                ANISO        ----    ----       ----     ----
%                KAREA        ----    ----       ----     ----
%                  ROT        ----    ----       ----     ----
%            ID_CUTORG        ----    ----       ----     ----
%              ID_BEND        ----    ----       ----     ----
%               F_SEAM        ----    ----       ----     ----
%            IN_TIPINH        ----    ----       ----     ----
%           ID_BORDERS        ----    ----       ----     ----


%%% USER CODE: MORPHOGEN INTERACTIONS

% set up the initial morphogens, diffusioin values etc.
if Steps(m) == 0
    % Set up visualisation morphogens
    v_leaf_p(:) =1;
    v_width_p(m.userdata.lateral) = 1;
    
    %% general identity factors
    
    miny = min(m.nodes(:,2));
    
    % RIM
    borderedges = find( m.edgecells(:,2)==0);
    bordernodes = unique(m.edgeends( borderedges,:));
    margin_nodes = setdiff(bordernodes,find(m.nodes(:,2)<=miny+0.01));
    id_rim_p(margin_nodes) = 1;
    pos = margin_nodes;
    margin = [];
    for i=pos'
        [row,col] = find(m.tricellvxs == i);
        margin =  [margin ;unique(m.tricellvxs(row,:))];
    end
    margin_nds2 = setdiff(margin, margin_nodes);
    id_rim_p(margin_nds2) = 1;
    
    
    % LAM
    m.userdata.lam = find(m.nodes(:,2)>(-0.038)& m.nodes(:,2)<0.013);%0.0085 0.01
    id_lam_p(m.userdata.lam) = 2.05; % 0.036 2.1
    m = leaf_mgen_conductivity(m, id_lam_i, 0.02); %0.011
        
    % MIDVEIN
    id_midvein_p(:) = 0.3;
    m.morphogens(:,id_midvein_i) = id_midvein_p;
    midmax = 1;
    m = leaf_mgen_linear(m,id_midvein_i,midmax,'direction',-90);
    [id_midvein_i,id_midvein_p,id_midvein_a,id_midvein_l] = getMgenLevels( m, 'ID_MIDVEIN' );
    midvein_ind = find(abs(m.nodes(:,1))>0.0059);
    id_midvein_p(m.nodes(:,2)<(-0.022)) = midmax+0.2;
    id_midvein_p([midvein_ind;find(m.nodes(:,2)>0.035)]) = 0;
    m = leaf_mgen_conductivity(m, id_midvein_i, 0.0005);
    
    %% polariser factors
    
    % BASE (+ Organiser)
    base_ind=find(m.nodes(:,2)<=miny+0.001);
    id_base_p(base_ind)=1;
    
    % DISTAL (- ORGANISER)
    id_distal_p(m.nodes(:,2)>0.04)=1.1;
    m = leaf_mgen_conductivity(m, id_distal_i, 0.0001);
    
    %POLARISER
    switch modelname
        case { 'LEAF', 'NARROW','MIDLINE','STRAIGHTCUT'};
            P(base_ind) = 0.1;
            m = leaf_fix_mgen( m, polariser_i,'vertex',base_ind,'fix', 1);
        case 'MIDPOL'
            midpol = find(abs(m.nodes(:,1))<0.0025); %0.0059
            lowerpol = find(m.nodes(:,2)<0.03);
            midpol_ind = intersect(midpol, lowerpol);
            id_leaf_p(midpol_ind) = 1;
            P(midpol_ind) = 0.1;
            m = leaf_fix_mgen( m, polariser_i,'vertex',midpol_ind,'fix', 1);
    end
    m = leaf_mgen_conductivity(m, polariser_i, 0.001);
    m = leaf_mgen_absorption(m, polariser_i, 0.01);
    
    %% Growth Factors
    
    % Growth Promotion proximal-distally
    s_prox_p(:) = 1;
    m.morphogens(:,s_prox_i) = s_prox_p;
    m = leaf_mgen_linear(m,s_prox_i,0.7,'direction',-90);
    [s_prox_i,s_prox_p,s_prox_a,s_prox_l] = getMgenLevels( m, 'S_PROX' );

        %% Thickness Factors for initilisation
    
    switch startname
        case {'Dshape','Primordium'}
            % knor inhibitor from tip
            in_tipinh_p(m.nodes(:,2)>0.04) = 1;
            m = leaf_mgen_conductivity(m, in_tipinh_i, 0.0005);
            
            % know inhibitor around margins
            borders = setdiff(bordernodes,find(m.nodes(:,2)<=miny)); % leave out to make curvature at bottom
            id_borders_p(borders) = 3;
            m = leaf_mgen_conductivity(m, id_borders_i, 0.00002);
    end
    %% Cutting
    
    switch modelname
        case 'STRAIGHTCUT'
            f_seam_p(m.userdata.tipoff) = 1;
            id_leaf_p(:) = 1;
            id_cutorg_p(:) = 1;
            id_cutorg_p(m.nodes(:,2)>(-0.017)) = 0;
    end
elseif Steps(m) == 1   %% Continue modifiying the identity factors and sprox
    

    m = leaf_mgen_conductivity(m, id_distal_i, 0);
    m = leaf_mgen_conductivity(m, id_midvein_i, 0);
    id_midvein_p(id_rim_l>0) = 0;
    
    m = leaf_mgen_conductivity(m, id_lam_i, 0);
    id_lam_p(m.nodes(:,2)<(min(m.nodes(:,2))+0.012)) = 0.6;
    
    m = leaf_setproperty(m, 'timestep', 1);
    m.userdata.time = 72;
    
    % Induce sectors
%     m = leaf_makesecondlayer( m, ...  % This function adds biological cells.
%         'mode', 'each', ...  % Make biological cells randomly scattered over the leaf.
%         'relarea', 1/80,...%(10*8)/(84*95), ...  % Each cell has area was 1/2000 of the initial area of the flower.
%         'probpervx', 'V_LEAF', ...
%         'numcells',50,...%number of clones
%         'sides', 20, ...  % Each cell is approximated as a 6-sided regular polygon.
%         'colorvariation', 0.5, ...  % Each cell is has a different colour [0.8 0.8 0.8];
%         'add', true,...
%         'allowoverlap', false);  % These cells are added to any cells existing already.
%     m = leaf_setproperty( m,'colorvariation', 0, 'bioAsplitcells', false, 'allowSplitBio',false);%'bioAsplitcells', false);
%     
   
    
    %% Biological Growth
else    
      %@@GRN
      % logistic function of kpar
      age = 0.018/(1+exp(0.018*(m.userdata.time-240))); % 319 240 for phase B and C, 255 good all phases 250
      id_youth_p(:) = age;
      id_youth_p(id_lam_l==0.6) = age*1.1; % 1.15
      
      % logistic function of kper
      agelam = 0.018/(1+exp(0.018*(m.userdata.time-280)));% 304  290   
      m.userdata.time =  m.userdata.time+1;
      
      % increase bending of lamina
      id_bend_p(id_lam_l>0.6) = 0.001;
     
     if realtime >= 144  
      id_leaf_p = id_cutorg_l;
     end
    
    switch modelname
        case {'LEAF','STRAIGHTCUT','MIDLINE'}
            % @@ KRN
            kapar_p(:) = s_prox_l.*id_youth_l;%.*id_leaf_l;
            kbpar_p(:) = kapar_p;%+id_bend_l;%0.001; %id_bend_l;
            kaper_p(:) = id_lam_l.*agelam....*id_leaf_l
                .*inh(0.7, id_midvein_l)...
                .*inh(0.7, id_rim_l); 
            kbper_p(:) = kaper_p(:);
        case 'MIDPOL'
            kapar_p(:) = id_lam_l.*inh(0.7, id_midvein_l)...
                .*inh(0.7, id_rim_l).*agelam;
            kbpar_p(:) = kapar_p;
            kaper_p(:) = s_prox_l.*id_youth_l;
            kbper_p(:) = kaper_p(:);
            
        case 'NARROW'
            
            % @@ KRN
            kapar_p(:) = s_prox_l.*id_youth_l;
            kbpar_p(:) = kapar_p;
            kaper_p(:) = 0;
            kbper_p(:) = kaper_p(:);
            
    end
    
end

%%     Other Stuff
  
 switch modelname
        case {'LEAF','STRAIGHTCUT'}
        %@@PRN
        m.mgen_production(:,polariser_i) = 0.001+id_base_l-P.*id_distal_l*0.1;
     case 'MIDPOL'
         m.mgen_production(:,polariser_i) = 0.001+id_leaf_l-P.*id_rim_l*0.1;  %leaf = midvein
 end
 
 switch startname
    case 'Dshape'
        if Steps(m) < 6
            knor_p(:) = 1*inh(2,in_tipinh_l).*inh(2,id_borders_l);
        else
            knor_p(:) = 0;
            m = leaf_setproperty( m,'thicknessMode', 'physical');
        end
    case 'Primordium'
        if Steps(m) < 6
            knor_p(:) = 1*inh(2,in_tipinh_l).*inh(2,id_borders_l);
            if Steps(m) == 5 
           knor_p(:) = 0;
           underside = m.prismnodes( 1:2:end, 3 ); 
           m.prismnodes(:,3) = m.prismnodes(:,3) - reshape( [(-0.3*underside');underside'], [], 1 );
           m.nodes(:,3) = m.nodes(:,3) - underside;
           m = recalc3d(m);
           m = leaf_setproperty( m,'thicknessMode', 'physical');
            end            
        end        
end
 
         m = calculateOutputs(m);
        [acA,frames,T1] = tensorsToComponentsTmatrix( m.outputs.actualstrain.A, m.cellFrames );
        kmax_p(:) = perFEtoperVertex( m, max(acA(:,[1 2]),[],2));
        kmin_p(:) = perFEtoperVertex( m, min(acA(:,[1 2]),[],2));
        
        karea_p(:) = kmax_p+kmin_p;
        aniso_p(:) = abs(kmax_p-kmin_p)./max(kmax_p,kmin_p);
        
        [smid, smidnorm] = leaf_midlinegrowth(m,acA,T1); 
        kml_p(:) = smid;
        kpml_p(:) = smidnorm;
        
        [inplane,outofplane] = splitVector( m.outputs.rotations,m.unitcellnormals );
        rot_p(:) = perFEtoperVertex( m, inplane);
 
        
%path = fileparts(which(m.globalProps.modelname));
%[m,ok] = leaf_snapshot( m,[path,filesep,'snapshots',filesep,'wtPhaseF_idbend0p001upto144h.png'], 'resolution',[]);
  

%%
%%% END OF USER CODE: MORPHOGEN INTERACTIONS

%%% SECTION 3: INSTALLING MODIFIED VALUES BACK INTO MESH STRUCTURE
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.
    m.morphogens(:,polariser_i) = P;
    m.morphogens(:,kapar_i) = kapar_p;
    m.morphogens(:,kaper_i) = kaper_p;
    m.morphogens(:,kbpar_i) = kbpar_p;
    m.morphogens(:,kbper_i) = kbper_p;
    m.morphogens(:,knor_i) = knor_p;
    m.morphogens(:,strainret_i) = strainret_p;
    m.morphogens(:,arrest_i) = arrest_p;
    m.morphogens(:,v_leaf_i) = v_leaf_p;
    m.morphogens(:,v_petiolelen_i) = v_petiolelen_p;
    m.morphogens(:,v_length_i) = v_length_p;
    m.morphogens(:,v_petiolewid_i) = v_petiolewid_p;
    m.morphogens(:,v_width_i) = v_width_p;
    m.morphogens(:,v_thickness_i) = v_thickness_p;
    m.morphogens(:,id_base_i) = id_base_p;
    m.morphogens(:,id_distal_i) = id_distal_p;
    m.morphogens(:,id_rim_i) = id_rim_p;
    m.morphogens(:,id_leaf_i) = id_leaf_p;
    m.morphogens(:,s_prox_i) = s_prox_p;
    m.morphogens(:,id_midvein_i) = id_midvein_p;
    m.morphogens(:,v_tip_i) = v_tip_p;
    m.morphogens(:,kmax_i) = kmax_p;
    m.morphogens(:,kmin_i) = kmin_p;
    m.morphogens(:,id_lam_i) = id_lam_p;
    m.morphogens(:,id_youth_i) = id_youth_p;
    m.morphogens(:,kml_i) = kml_p;
    m.morphogens(:,kpml_i) = kpml_p;
    m.morphogens(:,aniso_i) = aniso_p;
    m.morphogens(:,karea_i) = karea_p;
    m.morphogens(:,rot_i) = rot_p;
    m.morphogens(:,id_cutorg_i) = id_cutorg_p;
    m.morphogens(:,id_bend_i) = id_bend_p;
    m.morphogens(:,f_seam_i) = f_seam_p;
    m.morphogens(:,in_tipinh_i) = in_tipinh_p;
    m.morphogens(:,id_borders_i) = id_borders_p;

%%% USER CODE: FINALISATION
 switch modelname
        case 'STRAIGHTCUT'
if Steps(m) == 74
    m=leaf_set_seams(m,f_seam_p);
    m = leaf_dissect(m);
    m = leaf_explode(m,3);  
end
 end
% In this section you may modify the mesh in any way whatsoever.
%%% END OF USER CODE: FINALISATION

end


%%% USER CODE: SUBFUNCTIONS


function m = initproperties( m )
% This function is called at time zero in the INITIALISATION section of the
% interaction function.  It provides commands to set each of the properties
% that are contained in m.globalProps.  Uncomment whichever ones you would
% like to set yourself, and put in whatever value you want.
%
% Some of these properties are for internal use only and should never be
% set by the user.  At some point these will be moved into a different
% component of m, but for the present, just don't change anything unless
% you know what it is you're changing.

%    m = leaf_setproperty( m, 'trinodesvalid', true );
%    m = leaf_setproperty( m, 'prismnodesvalid', true );
%    m = leaf_setproperty( m, 'thicknessRelative', 0.010700 );
%    m = leaf_setproperty( m, 'thicknessArea', 1.000000 );
%    m = leaf_setproperty( m, 'activeGrowth', 1.000000 );
%    m = leaf_setproperty( m, 'displayedGrowth', 17.000000 );
%    m = leaf_setproperty( m, 'displayedMulti', [] );
%    m = leaf_setproperty( m, 'allowNegativeGrowth', true );
%    m = leaf_setproperty( m, 'usePrevDispAsEstimate', true );
%    m = leaf_setproperty( m, 'mingradient', 1.000000 );
%    m = leaf_setproperty( m, 'relativepolgrad', false );
%    m = leaf_setproperty( m, 'userpolarisation', false );
%    m = leaf_setproperty( m, 'thresholdsq', 0.000061 );
%    m = leaf_setproperty( m, 'splitmargin', 1.400000 );
%    m = leaf_setproperty( m, 'splitmorphogen', '' );
%    m = leaf_setproperty( m, 'thresholdmgen', 0.500000 );
%    m = leaf_setproperty( m, 'bulkmodulus', 1.000000 );
%    m = leaf_setproperty( m, 'unitbulkmodulus', true );
%    m = leaf_setproperty( m, 'poissonsRatio', 0.300000 );
%    m = leaf_setproperty( m, 'starttime', 68.000000 );
%    m = leaf_setproperty( m, 'timestep', 1.000000 );
%    m = leaf_setproperty( m, 'timeunitname', 'hours' );
%    m = leaf_setproperty( m, 'distunitname', 'mm' );
%    m = leaf_setproperty( m, 'scalebarvalue', 0.000000 );
%    m = leaf_setproperty( m, 'validateMesh', true );
%    m = leaf_setproperty( m, 'rectifyverticals', false );
%    m = leaf_setproperty( m, 'allowSplitLongFEM', false );
%    m = leaf_setproperty( m, 'longSplitThresholdPower', 0.000000 );
%    m = leaf_setproperty( m, 'allowSplitBentFEM', false );
%    m = leaf_setproperty( m, 'allowSplitBio', false );
%    m = leaf_setproperty( m, 'allowFlipEdges', false );
%    m = leaf_setproperty( m, 'allowElideEdges', false );
%    m = leaf_setproperty( m, 'mincellangle', 0.200000 );
%    m = leaf_setproperty( m, 'alwaysFlat', 1.000000 );
%    m = leaf_setproperty( m, 'flattenforceconvex', true );
%    m = leaf_setproperty( m, 'flatten', false );
%    m = leaf_setproperty( m, 'flattenratio', 1.000000 );
%    m = leaf_setproperty( m, 'useGrowthTensors', false );
%    m = leaf_setproperty( m, 'plasticGrowth', false );
%    m = leaf_setproperty( m, 'totalinternalrotation', 0.000000 );
%    m = leaf_setproperty( m, 'stepinternalrotation', 2.000000 );
%    m = leaf_setproperty( m, 'showinternalrotation', false );
%    m = leaf_setproperty( m, 'performinternalrotation', false );
%    m = leaf_setproperty( m, 'internallyrotated', false );
%    m = leaf_setproperty( m, 'maxFEcells', 0.000000 );
%    m = leaf_setproperty( m, 'inittotalcells', 0.000000 );
%    m = leaf_setproperty( m, 'maxBioAcells', 0.000000 );
%    m = leaf_setproperty( m, 'maxBioBcells', 0.000000 );
%    m = leaf_setproperty( m, 'colors', (3 values) );
%    m = leaf_setproperty( m, 'colorvariation', 0.000000 );
%    m = leaf_setproperty( m, 'colorparams', (6 values) );
%    m = leaf_setproperty( m, 'freezing', 0.000000 );
%    m = leaf_setproperty( m, 'canceldrift', 0.000000 );
%    m = leaf_setproperty( m, 'mgen_interaction', (unknown type ''function_handle'') );
%    m = leaf_setproperty( m, 'mgen_interactionName', 'arleaf_ver3187_ad0p4goodunfixrim_101021' );
%    m = leaf_setproperty( m, 'allowInteraction', true );
%    m = leaf_setproperty( m, 'interactionValid', true );
%    m = leaf_setproperty( m, 'gaussInfo', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'stitchDFs', [] );
%    m = leaf_setproperty( m, 'D', (36 values) );
%    m = leaf_setproperty( m, 'C', (36 values) );
%    m = leaf_setproperty( m, 'G', (6 values) );
%    m = leaf_setproperty( m, 'solver', 'cgs' );
%    m = leaf_setproperty( m, 'solvertolerance', 0.000100 );
%    m = leaf_setproperty( m, 'diffusiontolerance', 0.000001 );
%    m = leaf_setproperty( m, 'allowsparse', true );
%    m = leaf_setproperty( m, 'maxIters', 166.000000 );
%    m = leaf_setproperty( m, 'maxsolvetime', 1000.000000 );
%    m = leaf_setproperty( m, 'cgiters', 232.000000 );
%    m = leaf_setproperty( m, 'simsteps', 0.000000 );
%    m = leaf_setproperty( m, 'stepsperrender', 0.000000 );
%    m = leaf_setproperty( m, 'growthEnabled', true );
%    m = leaf_setproperty( m, 'diffusionEnabled', true );
%    m = leaf_setproperty( m, 'makemovie', 0.000000 );
%    m = leaf_setproperty( m, 'moviefile', '' );
%    m = leaf_setproperty( m, 'codec', 'None' );
%    m = leaf_setproperty( m, 'autonamemovie', true );
%    m = leaf_setproperty( m, 'overwritemovie', false );
%    m = leaf_setproperty( m, 'framesize', [] );
%    m = leaf_setproperty( m, 'mov', [] );
%    m = leaf_setproperty( m, 'jiggleProportion', 1.000000 );
%    m = leaf_setproperty( m, 'cvtperiter', 0.200000 );
%    m = leaf_setproperty( m, 'boingNeeded', false );
%    m = leaf_setproperty( m, 'initialArea', 0.008228 );
%    m = leaf_setproperty( m, 'bendunitlength', 0.090709 );
%    m = leaf_setproperty( m, 'targetRelArea', 1.000000 );
%    m = leaf_setproperty( m, 'defaultinterp', 'min' );
%    m = leaf_setproperty( m, 'readonly', false );
%    m = leaf_setproperty( m, 'projectdir', 'D:\DArT_Models\trunk\Growth Arrest\PaperWorthy\UnfixedBase' );
%    m = leaf_setproperty( m, 'modelname', 'ArLeaf_ver3187_Ad0p4GoodUnfixRim_101021' );
%    m = leaf_setproperty( m, 'allowsave', 0.000000 );
%    m = leaf_setproperty( m, 'addedToPath', true );
%    m = leaf_setproperty( m, 'bendsplit', 0.300000 );
%    m = leaf_setproperty( m, 'usepolfreezebc', false );
%    m = leaf_setproperty( m, 'dorsaltop', true );
%    m = leaf_setproperty( m, 'defaultazimuth', -45.000000 );
%    m = leaf_setproperty( m, 'defaultelevation', 33.750000 );
%    m = leaf_setproperty( m, 'defaultroll', 0.000000 );
%    m = leaf_setproperty( m, 'defaultViewParams', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'comment', '' );
%    m = leaf_setproperty( m, 'legendTemplate', '%T: %q\n%m' );
%    m = leaf_setproperty( m, 'bioAsplitcells', false );
%    m = leaf_setproperty( m, 'bioApullin', 0.142857 );
%    m = leaf_setproperty( m, 'bioAfakepull', 0.202073 );
%    m = leaf_setproperty( m, 'interactive', false );
%    m = leaf_setproperty( m, 'coderevision', 3187 );
%    m = leaf_setproperty( m, 'coderevisiondate', '2010-10-15 08:03:44.157047' );
%    m = leaf_setproperty( m, 'modelrevision', 3172.000000 );
%    m = leaf_setproperty( m, 'modelrevisiondate', '2010-10-12 08:38:13.898165' );
%    m = leaf_setproperty( m, 'vxgrad', (108 values) );
%    m = leaf_setproperty( m, 'lengthscale', 0.107000 );
%    m = leaf_setproperty( m, 'targetAbsArea', 0.007388 );
%    m = leaf_setproperty( m, 'RecordMeshes', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'usefrozengradient', true );
%    m = leaf_setproperty( m, 'solverprecision', 'double' );
%    m = leaf_setproperty( m, 'thicknessMode', 'physical' );
%    m = leaf_setproperty( m, 'savedrunname', '' );
%    m = leaf_setproperty( m, 'savedrundesc', '' );
%    m = leaf_setproperty( m, 'solvertolerancemethod', 'norm' );
%    m = leaf_setproperty( m, 'physicalThickness', 0.000000 );
end





% Here you may add any functions of your own, that you want to call from
% the interaction function, but never need to call from outside it.
% Whichever section they are called from, they must respect the same
% restrictions on what modifications they are allowed to make to the mesh.
% This comment can be deleted.of sink