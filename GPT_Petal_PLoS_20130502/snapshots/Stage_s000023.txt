function m = gpt_petal_plos_20130502( m )
%m = gpt_petal_plos_20130502( m )
%   Morphogen interaction function.
%   Written at 2013-05-02 13:51:36.
%   GFtbox revision 4575, .

% The user may edit any part of this function between delimiters
% of the form "USER CODE..." and "END OF USER CODE...".  The
% delimiters themselves must not be moved, edited, deleted, or added.

    if isempty(m), return; end

    fprintf( 1, '%s found in %s\n', mfilename(), which(mfilename()) );

    try
        m = local_setproperties( m );
    catch
    end

    realtime = m.globalDynamicProps.currenttime;

%%% USER CODE: INITIALISATION
% Model combinations explored in the paper: All the models are Organiser-based Dynamic Models. 
% They differ in the (-)organiser presence or positon and in the distribution and strength of the DGRAD factor promoting Kper.
% Convergent model:                                                                modelname='PROXORG',  (1)               submodel='G_DGRAD_topdist',(2)   
% Divergent model:                                                                 modelname='PROXDISTbroad',(3)           submodel='G_DGRAD_lobe',(3)     
% Divergent model with reduced or null DGRAD activity:                             modelname='PROXDISTbroad',(3)           submodel='G_DGRAD_reduced',(4) 
% Divergent model with reduced or null DGRAD activity and discontinuous DISTORG:   modelname='PROXDISTbroad_discont',(4)   submodel='G_DGRAD_reduced',(4)   
% Divergent model with ectopic DGRAD activity and extended DISTORG:                modelname='PROXDISTbroad_extended',(5)  submodel='G_DGRAD_extended',(5)   

if Steps(m) == 0
    
    % Model Visualisation and Setup Parameters
    m = leaf_setproperty( m, 'mingradient', 0, 'timestep', 1, 'userpolarisation',false, 'bioAsplitcells', false, 'allowSplitBio',false);
    m = leaf_plotoptions( m, 'hiresdpi', 900,'FEthinlinesize', 3, 'layeroffset', 0.01,'arrowthickness', 2.8,'highgradcolor',[0,0,0],'lowgradcolor',[0,0,1] );
       
    % If wanting to fix which morphogens to plot and in which layers
    %  m = leaf_plotoptions( m,'morphogen',{'ID_PROXORG','ID_DISTORG','s_gdist', 's_gpet'}); %to define which morphogens to plot
    %  m = leaf_mgen_plotpriority( m, {'ID_PROXORG','ID_DISTORG','s_gdist', 's_gpet'},[2,2,1,1]);% to specify the layers in which they appear
    %  m = leaf_mgen_plotpriority( m, {'s_gdist', 's_gtopmar'},[2,1]); %same

    % to specify certain vertexs to measure width at different positions along the proximodistal axis of the petal// the maximum lenght of the petal // and the lenght of petal claw
    m.userdata.length1 = [17 1190];
    m.userdata.width1 = [1201 1179];
    m.userdata.width2 = [1206 1174];
    m.userdata.width3 = [1211 1169];
    m.userdata.width4 = [1216 1164];
    m.userdata.width5 = [1221 1159];
    m.userdata.width6 = [307 645];
    m.userdata.petwidth = [1 373];
    m.userdata.lengthpet = [17 357];
    
    % to mark the position of the region where DGRAG activity is 1 for Figure preparation porpouses.
    m.userdata.dgrad1 = [1200 1180]; % convergent model
    m.userdata.dgrad2 = [1213 1167]; % divergent model
    m.userdata.dgrad3 = [239 597]; % extended model. 
    % clock time step 1h
    m.userdata.clock = 1;
    
    %stablish array for the sectors and number of edges in the sectors
    m.secondlayer.cell3dcoords = [];
    m.userdata.celledges = 20;
    %% Models
    % Organiser-based models with a (+) organiser
    m.userdata.ranges.modelname.range{1}='PROXORG'; % no distal organiser
    m.userdata.ranges.modelname.range{2}='Y_axis'; %POL based on XY
    m.userdata.ranges.modelname.index = 1;
    
    % Organiser-based model adding a (-) organiser
    m.userdata.ranges.submodel.range{1}='NONE'; %only proxorg
    m.userdata.ranges.submodel.range{2}='PROXDISTORG'; % narrow distal organiser
    m.userdata.ranges.submodel.range{3}='PROXDISTbroad'; % broad distal organiser  THIS ONE (Elasticity is 0.001)
    m.userdata.ranges.submodel.range{4}='PROXDISTbroad_extended'; %extended broad distal organiser
    m.userdata.ranges.submodel.range{5}='PROXDISTbroad_discont'; % broad distal organiser discontinuos// need to copy this from previous model
    m.userdata.ranges.submodel.index =3;
    
    % Organiser-based models with different organiser positions
    m.userdata.ranges.modelversion.range{1}='NO_G';
    m.userdata.ranges.modelversion.range{2}='G_DGRAD_topdist'; %convergent model
    m.userdata.ranges.modelversion.range{3}='G_DGRAD_lobe'; % divergent model
    m.userdata.ranges.modelversion.range{4}='G_DGRAD_reduced'; % jag mutant 30% AND THIS ONE (Elasticity is 0.001)
    m.userdata.ranges.modelversion.range{5}='G_DGRAD_extended'; %JAG ectopic expression, broader domain of JAG
     m.userdata.ranges.modelversion.range{6}='G_DGRAD_nul'; %jag mutant null AND THIS ONE (Elasticity is 0.001)
    m.userdata.ranges.modelversion.range{7}='G_DGRAD_elas'; %to try with reduce elasticity (instead of 0.001 used in models so far)
    m.userdata.ranges.modelversion.index =4;
    
    % choosing the time to induce the virtual clones   //I NEED TO LEVE ONLY
    % THE ONES I USE
    m.userdata.ranges.modelcells.range{1} = 'D20';
    m.userdata.ranges.modelcells.range{2} = 'E170';
    m.userdata.ranges.modelcells.range{3} = 'nocells';
    m.userdata.ranges.modelcells.index =3;
    %% Parameters
    % fixed levels of POL at bottom
    m.userdata.ranges.bpol.range = [0.08 0.8 0.12]; % levels of POL at the PROXORG region are fixed to this value
    m.userdata.ranges.bpol.index = 2;
    %   Distorg_id region,  +/- 20% of original canvas lenght (0.03mm), ie +/-0.006
    m.userdata.ranges.distorgpoint.range = [0.0207 0.0147 0.0087];
    m.userdata.ranges.distorgpoint.index = 2;
    m.userdata.ranges.distorgbroad.range = [0.006 0 -0.006];
    m.userdata.ranges.distorgbroad.index = 2;
    m.userdata.ranges.distorgextended.range = [-0.001 -0.007 -0.013];
    m.userdata.ranges.distorgextended.index = 2;
    %   Parameters for plateau position for DGRAD (region with levels=1),  +/- 20% of original canvas lenght (0.03mm), ie +/-0.006
    m.userdata.ranges.pldgradconv.range = [0.019 0.013 0.007];
    m.userdata.ranges.pldgradconv.index = 2;
    m.userdata.ranges.pldgraddiv.range = [0.012 0.006 0];
    m.userdata.ranges.pldgraddiv.index = 2;
    m.userdata.ranges.pldgradext.range = [0.001 -0.005 -0.011];
    m.userdata.ranges.pldgradext.index = 2;
    %   Parameters for basic growth rates Kpar and Kper, +/- 20%
    m.userdata.ranges.bkparconv.range = [0.0144 0.018 0.0216];
    m.userdata.ranges.bkparconv.index = 2;
    m.userdata.ranges.bkpardiv.range = [0.0152 0.019 0.0228];
    m.userdata.ranges.bkpardiv.index = 2;
    m.userdata.ranges.bkper.range = [0.0078 0.0065 0.0052];
    m.userdata.ranges.bkper.index = 2;
    %  Parameters for promotion of Kper by DGRAD   +/- 20%
    m.userdata.ranges.pdgradconv.range = [2.16 1.8 1.44];
    m.userdata.ranges.pdgradconv.index = 2;
    % parameters for the DRAG model and the DGRAG_reduced model to 0% or 30% activity
    m.userdata.ranges.pdgraddiv.range = [0.33 1.1 0];
    m.userdata.ranges.pdgraddiv.index = 2;
end

modelname = m.userdata.ranges.modelname.range{m.userdata.ranges.modelname.index};
submodel = m.userdata.ranges.submodel.range{m.userdata.ranges.submodel.index};
modelversion = m.userdata.ranges.modelversion.range{m.userdata.ranges.modelversion.index};
modelcells = m.userdata.ranges.modelcells.range{m.userdata.ranges.modelcells.index};

disp(sprintf('\n------> modelname %s submodel %s modelversion %s modelcells %s\n',modelname,submodel,modelversion,modelcells));

bpol = m.userdata.ranges.bpol.range(m.userdata.ranges.bpol.index);
distorgpoint = m.userdata.ranges.distorgpoint.range(m.userdata.ranges.distorgpoint.index);
distorgbroad = m.userdata.ranges.distorgbroad.range(m.userdata.ranges.distorgbroad.index);
distorgextended = m.userdata.ranges.distorgextended.range(m.userdata.ranges.distorgextended.index);
pldgradconv = m.userdata.ranges.pldgradconv.range(m.userdata.ranges.pldgradconv.index);  
pldgraddiv = m.userdata.ranges.pldgraddiv.range(m.userdata.ranges.pldgraddiv.index); 
pldgradext = m.userdata.ranges.pldgradext.range(m.userdata.ranges.pldgradext.index); 
bkper = m.userdata.ranges.bkper.range(m.userdata.ranges.bkper.index); 
bkparconv = m.userdata.ranges.bkparconv.range(m.userdata.ranges.bkparconv.index); 
bkpardiv = m.userdata.ranges.bkpardiv.range(m.userdata.ranges.bkpardiv.index); 
pdgradconv = m.userdata.ranges.pdgradconv.range(m.userdata.ranges.pdgradconv.index); 
pdgraddiv = m.userdata.ranges.pdgraddiv.range(m.userdata.ranges.pdgraddiv.index); 

%function m=leaf_vertex_set_monitor(m,realtime,RegionLabels,Morphogens,start_figno) monitor morphogen levels at a set of vertices.
%Several morphogens can be monitored simultaneously.
%The results can be recorded and an output in opne figure given over time. The pattern is saved automatically in the model directory. 
% made by Erika

% if Steps(m)>1 ;% PUT it back!!!!
%     m=leaf_vertex_follow_mgen_NoFigure(m,'RealTime',realtime,'ZeroTime',0,...
%     'REGIONLABELS',{'V_SEAM', 'V_SEAM', 'V_SEAM', 'V_SEAM','V_SEAM','V_SEAM'},...
%     'MORPHOGENS',{'KAPAR', 'KAPER', 'Kmax', 'Kmin', 'KAREA','KML'},... %, 'KAREA','KML'
%     'PROFILES',false,...
%     'ALLPOINTS',true, 'ENDTIME',12, ...
%     'RECORDINTERVAL', 12);
% 
% end   

% FOR THIS EXPERIMENTAL STUFF PLEASE CONTACT THE AUTHORS
% m = leaf_MeasureDimensions2(m, 'HEADINGS',{'L','W1','W2','W3','W4','W5','W6','petW','Lpet'},...
%     'LINES',{m.userdata.length1,m.userdata.width1,m.userdata.width2,m.userdata.width3,m.userdata.width4,m.userdata.width5,m.userdata.width6,...
%     m.userdata.petwidth, m.userdata.lengthpet},...
%     'LENWID',[2 1 1 1 1 1 1 1 2],'ENDTIME',312, 'RECORDINTERVAL',0:12:312);
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
    [id_rim_i,id_rim_p,id_rim_a,id_rim_l] = getMgenLevels( m, 'ID_RIM' );
    [v_length_i,v_length_p,v_length_a,v_length_l] = getMgenLevels( m, 'V_LENGTH' );
    [v_seam_i,v_seam_p,v_seam_a,v_seam_l] = getMgenLevels( m, 'V_SEAM' );
    [kmax_i,kmax_p,kmax_a,kmax_l] = getMgenLevels( m, 'KMAX' );
    [kmin_i,kmin_p,kmin_a,kmin_l] = getMgenLevels( m, 'KMIN' );
    [karea_i,karea_p,karea_a,karea_l] = getMgenLevels( m, 'KAREA' );
    [kaniso_i,kaniso_p,kaniso_a,kaniso_l] = getMgenLevels( m, 'KANISO' );
    [v_subdivision_i,v_subdivision_p,v_subdivision_a,v_subdivision_l] = getMgenLevels( m, 'V_SUBDIVISION' );
    [id_proxorg_i,id_proxorg_p,id_proxorg_a,id_proxorg_l] = getMgenLevels( m, 'ID_PROXORG' );
    [id_distorg_i,id_distorg_p,id_distorg_a,id_distorg_l] = getMgenLevels( m, 'ID_DISTORG' );
    [id_gpet_i,id_gpet_p,id_gpet_a,id_gpet_l] = getMgenLevels( m, 'ID_GPET' );
    [s_gpet_i,s_gpet_p,s_gpet_a,s_gpet_l] = getMgenLevels( m, 'S_GPET' );
    [kml_i,kml_p,kml_a,kml_l] = getMgenLevels( m, 'KML' );
    [id_gdistgrad_i,id_gdistgrad_p,id_gdistgrad_a,id_gdistgrad_l] = getMgenLevels( m, 'ID_GDISTGRAD' );
    [s_gdistgrad_i,s_gdistgrad_p,s_gdistgrad_a,s_gdistgrad_l] = getMgenLevels( m, 'S_GDISTGRAD' );
    [v_dgrad_i,v_dgrad_p,v_dgrad_a,v_dgrad_l] = getMgenLevels( m, 'V_DGRAD' );
    [v_sect_i,v_sect_p,v_sect_a,v_sect_l] = getMgenLevels( m, 'V_SECT' );

% Mesh type: lobes
%            base: 0
%        cylinder: 0
%          height: 0.015
%           lobes: 1
%          radius: 0.015
%      randomness: 0
%           rings: 24
%          strips: 30
%         version: 1

%            Morphogen    Diffusion   Decay   Dilution   Mutant
%            --------------------------------------------------
%                KAPAR         ----    ----       ----     ----
%                KAPER         ----    ----       ----     ----
%                KBPAR         ----    ----       ----     ----
%                KBPER         ----    ----       ----     ----
%                 KNOR         ----    ----       ----     ----
%            POLARISER       0.0005   0.003       ----     ----
%            STRAINRET         ----    ----       ----     ----
%               ARREST         ----    ----       ----     ----
%               ID_RIM         ----    ----       ----     ----
%             V_LENGTH         ----    ----       ----     ----
%               V_SEAM         ----    ----       ----     ----
%                 KMAX         ----    ----       ----     ----
%                 KMIN         ----    ----       ----     ----
%                KAREA         ----    ----       ----     ----
%               KANISO         ----    ----       ----     ----
%        V_SUBDIVISION         ----    ----       ----     ----
%           ID_PROXORG         ----    ----       ----     ----
%           ID_DISTORG         ----    ----       ----     ----
%              ID_GPET         ----    ----       ----     ----
%               S_GPET         ----    ----       ----     ----
%                  KML         ----    ----       ----     ----
%         ID_GDISTGRAD         ----    ----       ----     ----
%          S_GDISTGRAD         ----    ----       ----     ----
%              V_DGRAD         ----    ----       ----     ----
%               V_SECT         ----    ----       ----     ----


%%% USER CODE: MORPHOGEN INTERACTIONS


if Steps(m) == 1 % a bit later to give the cluster the chance to catch up.
    
    % id_v-length will include all the points for length and width set up at the begining
     v_length_p([m.userdata.length1;m.userdata.width1;m.userdata.width2;m.userdata.width3;m.userdata.width4;m.userdata.width5;m.userdata.width6;...
         m.userdata.petwidth;m.userdata.lengthpet]) =1;  
        
     % creates the path for the morphogen gradient measuring function
     v_seam_p(abs(m.nodes(:,1))<0.0005)=1;
     
     borderedges = find( m.edgecells(:,2)==0 );
     bordernodes = unique( m.edgeends( borderedges, : ) );
     id_rim_p( bordernodes ) = 1;
     
     %am I using it?
    v_subdivision_p(m.nodes(:,2)> -0.002)=1;  % 65
end

if Steps(m)==2
    %induce the virtual clones.
    switch modelcells
        case 'D20'
            m = leaf_makesecondlayer( m, ...  % This function adds biological cells.
                'mode', 'each', ...  % Make biological cells randomly scattered over the leaf.
                'relarea', 1/20,...% 1/30. try with bigger cells.not sure about the total area. beforer 1/35
                'probpervx', 'V_SECT', ...
                'numcells',70,...%number of clones
                'sides', 20, ...  % Each cell is approximated as a 6-sided regular polygon.
                'colors', [0,0.6,0] , ...  % [1,0.3,0.2]   Each cell is has a different colour;
                'colorvariation', 0.25,...
                'add', true,...
                'allowoverlap', false);
            % 'allowoveredge', false);  % use to prevent cells to go over the edge of the canvas. These cells are added to any cells existing already.
            m = leaf_setproperty( m,'colorvariation', 0, 'bioAsplitcells', false, 'allowSplitBio',false);
        case 'nocells'
    end
end

if Steps(m)>2
    switch modelcells
        case {'D20'}
            %morphogen that will represent the sectors on
            v_sect_p(:) =1;
    end
end

if Steps(m) == 4
   %% set up organisers and Polariser
   switch modelname
       
%        % in this way you fix the polairty to the initial Y position of each
%        % node but then it does not change with the organ deformation
        case 'Y_axis'
%            % setup polarity parallel to the external y-axis throughout the
%            % simulation
%            P = max(m.nodes(:,2)) - m.nodes(:,2);
%            P = P/max(P)*bpol;
%            
%            m.mgen_production(:,polariser_i) = 0;
%            m = leaf_mgen_conductivity( m, 'polariser', 0 );
%            m = leaf_mgen_absorption( m, 'polariser', 0 );
%            
           
           %set up the (+) organiser region
       case 'PROXORG'
           proxorg_ind= m.nodes(:,2)<=min(m.nodes(:,2))+0.0006; % this value could be a variable
           id_proxorg_p(proxorg_ind) = 1;
           % id_proxorg_l = id_proxorg_p * id_proxorg_a;
           
           %fix polariser levels at the PROXORg region
           P(proxorg_ind) = bpol; %0.8
           m = leaf_fix_mgen( m, polariser_i,'vertex',find(proxorg_ind),'fix', 1); %this is another way of fixing, instead of doing the clamp thing
           
           
           % add the distal organiser region
           switch submodel
               case 'PROXDISTORG' % narrow distal
                   distorg_ind= id_rim_p.*m.nodes(:,2)>distorgpoint; % distorgpoint 0.0147
                   id_distorg_p(distorg_ind) = 1;
                   
               case 'PROXDISTbroad' % broad distal. WT petal
                   distorg_ind= m.nodes(:,2)> distorgbroad; % distorgbroad 0
                   id_distorg_p(distorg_ind) = 1;
                   id_distorg_p = id_distorg_p.*id_rim_p;
                   id_distorg_l = id_distorg_p * id_distorg_a;
                   
               case 'PROXDISTbroad_extended' %extended broad distal organiser
                   distorg_ind= m.nodes(:,2)> distorgextended; % distorgextended -0.007
                   id_distorg_p(distorg_ind) = 1;
                   id_distorg_p = id_distorg_p.*id_rim_p;
                   id_distorg_l = id_distorg_p * id_distorg_a;
                   
               case 'PROXDISTbroad_discont' %broad distal organiser discontinuos
                   id_distorg_p(id_rim_p.*m.nodes(:,2)>0 & m.nodes(:,2)<0.0085)=1 ; % 
                   id_distorg_p(id_rim_p.*m.nodes(:,2)>0.011 & m.nodes(:,2)<0.0149 ) =1;% 
                   %id_distorg_p(distorg_ind) = 1;
                   id_distorg_p = id_distorg_p.*id_rim_p;
                   id_distorg_l = id_distorg_p * id_distorg_a;
                   
               case 'NONE'
               otherwise
                   error('You have to specify an existing submodel');
           end    
                   % Polariser parameters for all models
                   m = leaf_mgen_conductivity(m, polariser_i, 0.0005); % difussion constant
                   m = leaf_mgen_absorption(m, polariser_i, 0.003); % decay constant
                   m = leaf_mgen_dilution(m, polariser_i, 0);
           
           
               otherwise
                   error('You have to specify an existing modelname');
   end  
    
        %% FACTORS
    % Factors influencing growth
    switch modelversion
        case 'G_DGRAD_topdist' %DGRAD domain used for convergent model
            %domain of DGRAD
            id_gdistgrad_p(m.nodes(:,2)> pldgradconv)=1; %pldgradconv 0.013
            id_gdistgrad_l=id_gdistgrad_p*id_gdistgrad_a;
            %visualise DGRAG plateau region
            v_dgrad_p([m.userdata.dgrad1]) =1;
            
        case {'G_DGRAD_lobe', 'G_DGRAD_reduced','G_DGRAD_elas'} %DGRAD domain used for divergent model
            %  domain of DGRAD
            id_gdistgrad_p(m.nodes(:,2)> pldgraddiv)=1; %pldgraddiv 0.006
            id_gdistgrad_l=id_gdistgrad_p*id_gdistgrad_a;
            %visualise DGRAG plateau region
            v_dgrad_p([m.userdata.dgrad2]) =1;
            
        case  'G_DGRAD_extended'
            %  domain of DGRAD %DGRAD domain used for ectopic JAG model
            id_gdistgrad_p(m.nodes(:,2)> pldgradext)=1; %pldgradext -0.005
            id_gdistgrad_l=id_gdistgrad_p*id_gdistgrad_a;
            %visualise DGRAG plateau region
            v_dgrad_p([m.userdata.dgrad3]) =1;
        case {'NO_G','G_DGRAD_nul'} %look how to put best te Dgrad reduced
            
        otherwise
    end
            % DRAG diffusion parameters
            m = leaf_mgen_conductivity(m, s_gdistgrad_i, 0.01);
            m = leaf_mgen_absorption(m, s_gdistgrad_i, 0.01);
            m = leaf_mgen_dilution(m, s_gdistgrad_i, 0);
            m.mgen_production(:, s_gdistgrad_i) = id_gdistgrad_p;
                        
            % fix DRAG levels to 1 in DRAG domain
            m.morphogenclamp( id_gdistgrad_p==1, s_gdistgrad_i ) = 1;
            s_gdistgrad_p = 1 * id_gdistgrad_p;
            
            %define petiol region
            gpet_ind= m.nodes(:,2)<(-0.014);
            id_gpet_p(gpet_ind) = 1;
            %define levels of DRAG in petiole region to be 0
            s_gdistgrad_p(gpet_ind) = 0;
            m = leaf_fix_mgen( m, s_gdistgrad_i,'vertex',find(gpet_ind),'fix', 1);
 
elseif Steps(m)==5
 %% fix DRAG
     switch modelversion
       case {'G_DGRAD_topdist', 'G_DGRAD_lobe','G_DGRAD_reduced','G_DGRAD_extended','G_DGRAD_elas' }
       m = leaf_mgen_conductivity(m, s_gdistgrad_i, 0);
       m = leaf_mgen_absorption(m, s_gdistgrad_i, 0);
       m = leaf_mgen_dilution(m, s_gdistgrad_i, 0);
       m.mgen_production(:, s_gdistgrad_i) = 0;
       otherwise
     end    
elseif Steps(m)==23  
    switch modelname
        case 'Y_axis_fixed'
             m = leaf_setproperty( m, 'mingradient', 200); %2
             
             m = leaf_mgen_conductivity(m, polariser_i, 0);
             m = leaf_mgen_absorption(m, polariser_i, 0);
        otherwise
    end
 elseif Steps(m)>23
 %%  KRN
 switch modelversion
     case 'NO_G'% gives a bit short petal but wider as well, so i can not increase L wo increasing further W
         if Steps(m)>=24 && Steps(m)<312
             kapar_p= bkparconv;
             kbpar_p=kapar_p;
             
             kaper_p= bkper;
             kbper_p=kaper_p;
         end
         %
     case 'G_DGRAD_topdist'
         if Steps(m)>=24 && Steps(m)<312
             kapar_p= bkparconv; % bkparconv 0.018
             kbpar_p=kapar_p;
             
             kaper_p= bkper.*pro(pdgradconv,s_gdistgrad_l); % bkper 0.0065// pdgradconv 1.8 // change to 0.5 just to try issue of heigh
             kbper_p=kaper_p;
         end
         
     case {'G_DGRAD_lobe','G_DGRAD_extended' }% 
         if Steps(m)>=24 && Steps(m)<312
             kapar_p= bkpardiv; %bkpardiv 0.019
             kbpar_p=kapar_p;
             
             kaper_p= bkper.*pro(pdgraddiv,s_gdistgrad_l); %bkper// pdgraddiv 1.1
             kbper_p=kaper_p;
         end
     case 'G_DGRAD_reduced'
         if Steps(m)>=24 && Steps(m)<312
             kapar_p= bkpardiv; %bkpardiv 0.019
             kbpar_p=kapar_p;
             
             kaper_p= bkper.*pro(0.33,s_gdistgrad_l); %bkper// pdgraddiv 1.1// 30% is 0.33, changed to 75% just to chekc the length and now 50%
             kbper_p=kaper_p;
         end
     case 'G_DGRAD_nul'
         if Steps(m)>=24 && Steps(m)<312
             kapar_p= 0.018; %bkpardiv 0.019 just to try now with 0.0001, changed to 0.018
             kbpar_p=kapar_p;
             
             kaper_p= bkper.*pro(0,s_gdistgrad_l); %bkper// pdgraddiv 1.1
             kbper_p=kaper_p;
         end
         
     case'G_DGRAD_elas';
         if Steps(m)>=24 && Steps(m)<312
             kapar_p= 0.018; % instead of bkpardiv 0.019
             kbpar_p=kapar_p;
             
             kaper_p= bkper.*pro(1.1,s_gdistgrad_l); %bkper// pdgraddiv 1.1 
             kbper_p=kaper_p;
         end
         
     otherwise
         error('You have to specify an existing modelversion');
         
 end
end
%%induce virtual clones later in petal development
if Steps(m)==71
    switch modelcells
        case 'E170'
            m = leaf_makesecondlayer( m, ...  % This function adds biological cells.
                'mode', 'each', ...  % Make biological cells randomly scattered over the leaf.
                'relarea', 1/170,...% 1/30. try with bigger cells.not sure about the total area. beforer 1/35
                'probpervx', 'V_SECT', ...
                'numcells',120,...%number of clones
                'sides', 20, ...  % Each cell is approximated as a 6-sided regular polygon.
                'colors', [0,0.7,0] , ...  % [1,0.3,0.2]   Each cell is has a different colour;
                'colorvariation', 0.3,...
                'add', true,...
                'allowoverlap', false);
            % 'allowoveredge', false);  % These cells are added to any cells existing already.
            m = leaf_setproperty( m,'colorvariation', 0, 'bioAsplitcells', false, 'allowSplitBio',false);
        case 'nocells'
    end
elseif Steps(m)>71
    switch modelcells
        case 'E170'
            %morphogen that will represent the sectors on
            v_sect_p(:) =1;
    end
end
%% Polariser import at (-) organiser for all models
switch modelname
    case 'Y_axis'% now it is a Y_axis dynamic
        P = max(m.nodes(:,2)) - m.nodes(:,2);
        P = P/max(P)*bpol;
        
           m.mgen_production(:,polariser_i) = 0;
           m = leaf_mgen_conductivity( m, 'polariser', 0 );
           m = leaf_mgen_absorption( m, 'polariser', 0 );
    otherwise
        if Steps(m)>=2  && Steps(m)<=120 %
            m.mgen_production(:, polariser_i) =  -id_distorg_l*5.*P; %
        elseif Steps(m)>120   && Steps(m)<=312 %need to go down as levels go down
            m.mgen_production(:, polariser_i) =  -id_distorg_l*1.*P;
        end
end

%%calucalte model outputs
% m = calculateOutputs(m);
% [acA,frames,T1] = tensorsToComponentsTmatrix( m.outputs.actualstrain.A, m.cellFrames );
% kmax_p(:) = perFEtoperVertex( m, max(acA(:,[1 2]),[],2));
% kmin_p(:) = perFEtoperVertex( m, min(acA(:,[1 2]),[],2));
% 
% karea_p(:) = kmax_p+kmin_p;
% kaniso_p(:) = abs(kmax_p-kmin_p)./max(kmax_p,kmin_p);
% 
% 
% % determine the K along the midline.
% pos = find(v_seam_p>0.5);
% [smid, smidnorm] = leaf_midlinefactors(m,acA,T1,pos);
% kml_p(:) = smid; %parallel to the seam line // to y axis
%klat_p(:) = smidnorm; %perpendicular to v-seam line




%    %at the moment is doing it every time     
%         num = size(m.secondlayer.cellcolor,1);
%         col = zeros(num,3);
%         col(:,2) = 0.6; %which of the 3 columns in RBG you replace
% % %        %col(:,3)=1;% or whatever
%         m.secondlayer.cellcolor = col;%replace color with new array of colors.
% %         
%         [inplane,outofplane] = splitVector( m.outputs.rotations,m.unitcellnormals );
%         rot_p(:) = perFEtoperVertex( m, inplane);
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
    m.morphogens(:,id_rim_i) = id_rim_p;
    m.morphogens(:,v_length_i) = v_length_p;
    m.morphogens(:,v_seam_i) = v_seam_p;
    m.morphogens(:,kmax_i) = kmax_p;
    m.morphogens(:,kmin_i) = kmin_p;
    m.morphogens(:,karea_i) = karea_p;
    m.morphogens(:,kaniso_i) = kaniso_p;
    m.morphogens(:,v_subdivision_i) = v_subdivision_p;
    m.morphogens(:,id_proxorg_i) = id_proxorg_p;
    m.morphogens(:,id_distorg_i) = id_distorg_p;
    m.morphogens(:,id_gpet_i) = id_gpet_p;
    m.morphogens(:,s_gpet_i) = s_gpet_p;
    m.morphogens(:,kml_i) = kml_p;
    m.morphogens(:,id_gdistgrad_i) = id_gdistgrad_p;
    m.morphogens(:,s_gdistgrad_i) = s_gdistgrad_p;
    m.morphogens(:,v_dgrad_i) = v_dgrad_p;
    m.morphogens(:,v_sect_i) = v_sect_p;

%%% USER CODE: FINALISATION
% m = leaf_mgen_plotpriority( m, {'id_proxorg','ID_DISTORG','s_gdistgrad'},[2,2,1]);
% 
%   % change offset of cells layers
 %  m = leaf_plotoptions( m, 'hiresdpi', 900,'FEthinlinesize',3, 'layeroffset', 0.04,'arrowthickness', 2.8,'highgradcolor',[0,0,0],'lowgradcolor',[0,0,1] );
% 
%   %this one to chnage cell colours   
%         num = size(m.secondlayer.cellcolor,1);
%         col = zeros(num,3);
%         col(:,2) = 0.7; %which of the 3 columns in RBG you replace
% % %        %col(:,3)=1;% or whatever
%         m.secondlayer.cellcolor = col;%replace color with new array of colors.
  
%  if Steps(m) == 1
%          m = leaf_subdivide(m,'morphogen','V_SUBDIVISION','min',0.1,'max',2,'mode','mid');  
%   end
%  
% In this section you may modify the mesh in any way whatsoever.
%%% END OF USER CODE: FINALISATION

end


%%% USER CODE: SUBFUNCTIONS

function m = local_setproperties( m )
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
%    m = leaf_setproperty( m, 'thicknessRelative', 0.003000 );
%    m = leaf_setproperty( m, 'thicknessArea', 0.000000 );
%    m = leaf_setproperty( m, 'thicknessMode', 'physical' );
%    m = leaf_setproperty( m, 'activeGrowth', 1.000000 );
%    m = leaf_setproperty( m, 'displayedGrowth', 1.000000 );
%    m = leaf_setproperty( m, 'displayedMulti', [] );
%    m = leaf_setproperty( m, 'allowNegativeGrowth', true );
%    m = leaf_setproperty( m, 'usePrevDispAsEstimate', true );
%    m = leaf_setproperty( m, 'mingradient', 0.000000 );
%    m = leaf_setproperty( m, 'relativepolgrad', false );
%    m = leaf_setproperty( m, 'usefrozengradient', true );
%    m = leaf_setproperty( m, 'userpolarisation', false );
%    m = leaf_setproperty( m, 'thresholdsq', 0.000009 );
%    m = leaf_setproperty( m, 'splitmargin', 1.400000 );
%    m = leaf_setproperty( m, 'splitmorphogen', '' );
%    m = leaf_setproperty( m, 'thresholdmgen', 0.500000 );
%    m = leaf_setproperty( m, 'bulkmodulus', 1.000000 );
%    m = leaf_setproperty( m, 'unitbulkmodulus', true );
%    m = leaf_setproperty( m, 'poissonsRatio', 0.300000 );
%    m = leaf_setproperty( m, 'starttime', 0.000000 );
%    m = leaf_setproperty( m, 'timestep', 0.010000 );
%    m = leaf_setproperty( m, 'timeunitname', '' );
%    m = leaf_setproperty( m, 'distunitname', 'mm' );
%    m = leaf_setproperty( m, 'scalebarvalue', 0.000000 );
%    m = leaf_setproperty( m, 'validateMesh', true );
%    m = leaf_setproperty( m, 'rectifyverticals', false );
%    m = leaf_setproperty( m, 'allowSplitLongFEM', true );
%    m = leaf_setproperty( m, 'longSplitThresholdPower', 0.000000 );
%    m = leaf_setproperty( m, 'allowSplitBentFEM', false );
%    m = leaf_setproperty( m, 'allowSplitBio', true );
%    m = leaf_setproperty( m, 'allowFlipEdges', false );
%    m = leaf_setproperty( m, 'allowElideEdges', true );
%    m = leaf_setproperty( m, 'mincellangle', 0.200000 );
%    m = leaf_setproperty( m, 'alwaysFlat', true );
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
%    m = leaf_setproperty( m, 'maxFEcells', 0 );
%    m = leaf_setproperty( m, 'inittotalcells', 0 );
%    m = leaf_setproperty( m, 'maxBioAcells', 0 );
%    m = leaf_setproperty( m, 'maxBioBcells', 0 );
%    m = leaf_setproperty( m, 'colors', (6 values) );
%    m = leaf_setproperty( m, 'colorvariation', 0.050000 );
%    m = leaf_setproperty( m, 'colorparams', (12 values) );
%    m = leaf_setproperty( m, 'freezing', 0.000000 );
%    m = leaf_setproperty( m, 'canceldrift', false );
%    m = leaf_setproperty( m, 'mgen_interaction', '' );
%    m = leaf_setproperty( m, 'mgen_interactionName', 'gpt_110202_kw1' );
%    m = leaf_setproperty( m, 'allowInteraction', true );
%    m = leaf_setproperty( m, 'interactionValid', true );
%    m = leaf_setproperty( m, 'gaussInfo', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'stitchDFs', [] );
%    m = leaf_setproperty( m, 'D', (36 values) );
%    m = leaf_setproperty( m, 'C', (36 values) );
%    m = leaf_setproperty( m, 'G', (6 values) );
%    m = leaf_setproperty( m, 'solver', 'cgs' );
%    m = leaf_setproperty( m, 'solverprecision', 'double' );
%    m = leaf_setproperty( m, 'solvertolerance', 0.001000 );
%    m = leaf_setproperty( m, 'solvertolerancemethod', 'norm' );
%    m = leaf_setproperty( m, 'diffusiontolerance', 0.000010 );
%    m = leaf_setproperty( m, 'allowsparse', true );
%    m = leaf_setproperty( m, 'maxIters', 0 );
%    m = leaf_setproperty( m, 'maxsolvetime', 1000.000000 );
%    m = leaf_setproperty( m, 'cgiters', 0 );
%    m = leaf_setproperty( m, 'simsteps', 0 );
%    m = leaf_setproperty( m, 'stepsperrender', 0 );
%    m = leaf_setproperty( m, 'growthEnabled', true );
%    m = leaf_setproperty( m, 'diffusionEnabled', true );
%    m = leaf_setproperty( m, 'makemovie', false );
%    m = leaf_setproperty( m, 'moviefile', '' );
%    m = leaf_setproperty( m, 'codec', 'None' );
%    m = leaf_setproperty( m, 'autonamemovie', true );
%    m = leaf_setproperty( m, 'overwritemovie', false );
%    m = leaf_setproperty( m, 'framesize', [] );
%    m = leaf_setproperty( m, 'mov', [] );
%    m = leaf_setproperty( m, 'jiggleProportion', 1.000000 );
%    m = leaf_setproperty( m, 'cvtperiter', 0.200000 );
%    m = leaf_setproperty( m, 'boingNeeded', false );
%    m = leaf_setproperty( m, 'initialArea', 0.000653 );
%    m = leaf_setproperty( m, 'bendunitlength', 0.025553 );
%    m = leaf_setproperty( m, 'targetRelArea', 1.000000 );
%    m = leaf_setproperty( m, 'defaultinterp', 'min' );
%    m = leaf_setproperty( m, 'readonly', false );
%    m = leaf_setproperty( m, 'projectdir', 'D:\MATLAB\Petal-models' );
%    m = leaf_setproperty( m, 'modelname', 'GPT_110202_Kw1' );
%    m = leaf_setproperty( m, 'allowsave', true );
%    m = leaf_setproperty( m, 'addedToPath', false );
%    m = leaf_setproperty( m, 'bendsplit', 0.300000 );
%    m = leaf_setproperty( m, 'usepolfreezebc', false );
%    m = leaf_setproperty( m, 'dorsaltop', true );
%    m = leaf_setproperty( m, 'defaultazimuth', -45.000000 );
%    m = leaf_setproperty( m, 'defaultelevation', 33.750000 );
%    m = leaf_setproperty( m, 'defaultroll', 0.000000 );
%    m = leaf_setproperty( m, 'defaultViewParams', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'comment', '' );
%    m = leaf_setproperty( m, 'legendTemplate', '%T: %q\n%m' );
%    m = leaf_setproperty( m, 'bioAsplitcells', true );
%    m = leaf_setproperty( m, 'bioApullin', 0.142857 );
%    m = leaf_setproperty( m, 'bioAfakepull', 0.202073 );
%    m = leaf_setproperty( m, 'interactive', false );
%    m = leaf_setproperty( m, 'coderevision', 0 );
%    m = leaf_setproperty( m, 'coderevisiondate', '' );
%    m = leaf_setproperty( m, 'modelrevision', 0 );
%    m = leaf_setproperty( m, 'modelrevisiondate', '' );
%    m = leaf_setproperty( m, 'savedrunname', '' );
%    m = leaf_setproperty( m, 'savedrundesc', '' );
%    m = leaf_setproperty( m, 'vxgrad', (108 values) );
%    m = leaf_setproperty( m, 'lengthscale', 0.030000 );
end

% Here you may write any functions of your own, that you want to call from
% the interaction function, but never need to call from outside it.
% Remember that they do not have access to any variables except those
% that you pass as parameters, and cannot change anything except by
% returning new values as results.
% Whichever section they are called from, they must respect the same
% restrictions on what modifications they are allowed to make to the mesh.

% For example:

% function m = do_something( m )
%   % Change m in some way.
% end

% Call it from the main body of the interaction function like this:
%       m = do_something( m );
