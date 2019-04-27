%% Plot Brain

% Generating individual hemispheres
[cortex.vert_lh,cortex.tri_lh]= read_surf('./surf/lh.pial'); % Reading left side pial surface
[cortex.vert_rh,cortex.tri_rh]= read_surf('./surf/rh.pial'); % Reading right side pial surface
% Generating entire cortex
cortex.vert = [cortex.vert_lh; cortex.vert_rh]; % Combining both hemispheres
cortex.tri = [cortex.tri_lh; (cortex.tri_rh + length(cortex.vert_lh))]; % Combining faces (Have to add to number of faces)
cortex.tri=cortex.tri+1; % freesurfer starts at 0 for indexing


%Reading Annotations
[v,l,ct] = read_annotation(fullfile(cd,'label','lh.aparc.a2009s.annot'));
[v1,l1,~] = read_annotation(fullfile(cd,'label','rh.aparc.a2009s.annot'));
v = [v;v1+max(v)]; l = [l;l1];
Filter = {'G_and_S_cingul-Ant','G_and_S_cingul-Mid-Ant','G_and_S_cingul-Mid-Post'};
%corresponds to 7,8,9 in the table
Filter1 = ct.table(7:9,5);
[v_filt] = find(ismember(l,Filter1));


% Reading in MRI parameters
f = MRIread(fullfile(cd,'mri','orig.mgz'));

% Translating into the appropriate space
for k=1:size(cortex.vert,1)
    a=f.vox2ras/f.tkrvox2ras*[cortex.vert(k,:) 1]';
    cortex.vert(k,:)=a(1:3)';
end

%Apply filter to cortex
cortex1 = cortex;
adata = ones(size(cortex1.vert,1),1)*1;
adata(v_filt) = 1;
cdata = repmat([.65,.65,.65],size(cortex1.vert,1),1);
cdata(v_filt,:) = repmat([1,.25,.25],length(v_filt),1);


% Displaying cortex
clf
Hp = patch('vertices',cortex1.vert,'faces',cortex1.tri(:,[1 3 2]),...
    'facecolor','flat', 'facevertexcdata',cdata,...
    'facealpha','flat','facevertexalphadata',adata,...
    'edgecolor','none','facelighting', 'gouraud', 'specularstrength', .25);
camlight('headlight','infinite');
fh(1)=gcf;
vertnormals = get(Hp,'vertexnormals');
axis off; axis equal

%% Plot electrodes
sbj = 'Sbj2';
load(fullfile(cd,'Electrode_Locations', sprintf('electrodes_projected_%s',sbj)))

for i_plot = 1:numel(eleccell)
    electrodeH{i_plot} = el_add(eleccell{i_plot},'b',15); % Plots the electrodes
    
    % Sets the handles for each of the electrodes and enables them to be
    % visible or not
    hcb{i_plot} = ['a=get(electrodeH(' num2str(i_plot) '),''visible''); if strcmp(a,''on'')==1, set(electrodeH(' num2str(i_plot) '),''visible'',''off''), else set(electrodeH(' num2str(i_plot) '),''visible'',''on''), end'];
    
    %cm_label{i_plot} = uimenu(hcmenu,'Label',electrode_names{i_plot},'Callback',hcb{i_plot});
end

cameratoolbar
%% Adjust Parameters

cortex_color = [.65,.65,.65];
target_color = [1,.25,.25];
target_alpha = 1.5;

%Apply filter to cortex
cortex1 = cortex;
adata = ones(size(cortex1.vert,1),1)*1;
adata(v_filt) = target_alpha;
cdata = repmat(cortex_color,size(cortex1.vert,1),1);
cdata(v_filt,:) = repmat(target_color,length(v_filt),1);

set(Hp,'facevertexalphadata',adata,'facevertexcdata',cdata);