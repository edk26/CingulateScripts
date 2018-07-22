%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electrode Localization Instructions                                    %
% This program takes the output of Freesurfer and post-op CT scan to     %
% perform electrode localization                                         %
%                                                                        %
% Michael Randazzo 8/6/14                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set path to MRI and CT image (If any variables are not set, GUI
% will allow for input)
Subject_ID = 'AA0716';
cortical_disp = 'freesurfer';
Subject_Path = '/Users/Witek/Documents/Data';
Grid_Names = {'LG','LPST'};
Grid_Lower = [1,81];
Grid_Upper = [64,96];
Strip_Names = {'LAPH','LPIH', 'LFP', 'LTP'};
Strip_Lower = [65,69,73,97];
Strip_Upper = [68,72,80,100];
Depth_Names = {'LHIPP','LPH'};
Depth_Lower = [101,109];
Depth_Upper = [108,116];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

% If variables are not set above then GUI will be used

% Input Subject ID and Number of grids and strips
if isempty(Subject_ID) && isempty(Grid_Names) && isempty(Strip_Names)
    prompt_subject = {'Enter the Subject ID','Enter the Number of Grids','Enter the Number of Strips'};
    dlg_title_subject = 'Subject Information';
    num_lines = 1;
    defaults_subject = {Subject_ID,num2str(length(Grid_Names)),num2str(length(Strip_Names))};
    Subject_Info = inputdlg(prompt_subject,dlg_title_subject,num_lines,defaults_subject);
    Subject_ID = Subject_Info{1};
    
    % Define strips and grids
    prompt_electrode = {'Enter Electrode Name','Enter Minimum Electrode Number','Enter Maximum Electrode Number'};
    dlg_title_grid = 'Grid Information';
    defaults_electrode = {'LP','0','0'};
    % Grids
    for i_grids = 1:str2double(Subject_Info{2})
        Grid_Info = inputdlg(prompt_electrode,dlg_title_grid,num_lines,defaults_electrode);
        Grid_Names{i_grids} = Grid_Info{1};
        Grid_Lower(i_grids) = str2double(Grid_Info{2});
        Grid_Upper(i_grids) = str2double(Grid_Info{3});
    end
    % Strips
    dlg_title_strips = 'Strip Informarion';
    for i_strips = 1:str2double(Subject_Info{3})
        Strip_Info = inputdlg(prompt_electrode,dlg_title_strips,num_lines,defaults_electrode);
        Strip_Names{i_strips} = Strip_Info{1};
        Strip_Lower(i_strips) = str2double(Strip_Info{2});
        Strip_Upper(i_strips) = str2double(Strip_Info{3});
    end
end

% Get Subject Directory
if isempty(Subject_Path)
    Subject_Path = uigetdir('/Users/richardsonlab/Subjects/');
    Subject_Path = strcat(Subject_Path,'/');
end

%%

% Step 1 
% Coregister and reslice the CT image (source) to the MRI image (reference)
cd(fullfile(Subject_Path,[Subject_ID,'_FS']));

% SPM Commands
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[Subject_Path,'/',Subject_ID,'_FS/mri/T1.nii,1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[Subject_Path,'/',Subject_ID,'_FS/ct/ct_',Subject_ID,'.nii,1']};
%matlabbatch{1}.spm.spatial.coreg.estwrite.source = fullfile(Subject_Path,[Subject_ID,'_FS'],'ct',['ct_',Subject_ID,'.nii,1']);
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% Executes the SPM commands
spm('defaults', 'FMRI');
spm_jobman('serial',matlabbatch);

%%

% Step 2
% Localize the electrodes
mkdir Electrode_Locations;
ctmr % Runs the electrode localization program
% Need to load co-registered and realigned CT
% Identfiy the electrodes

%%

% Renaming the electrode output
movefile(fullfile(Subject_Path,[Subject_ID,'_FS'],'Electrode_Locations','electrodes1.hdr'),  fullfile(Subject_Path,[Subject_ID,'_FS'],'Electrode_Locations',['electrodes_',Subject_ID,'.hdr']) );
movefile(fullfile(Subject_Path,[Subject_ID,'_FS'],'Electrode_Locations','electrodes1.img'),  fullfile(Subject_Path,[Subject_ID,'_FS'],'Electrode_Locations',['electrodes_',Subject_ID,'.img']) );

% Step 3
% Sorting the electrodes
sortElectrodes(Subject_ID);

% Renaming the electrode locations .mat file
% Renaming the electrode output
movefile([Subject_Path,Subject_ID,'_FS/Electrode_Locations/electrodes_loc1.mat'],[Subject_Path,Subject_ID,'_FS/Electrode_Locations/electrodes_loc_',Subject_ID,'.mat']);

%%

% Step 4
% Project the electrodes onto the surface

% Creating a smooth brain surface to project to
get_mask_V3(Subject_ID,[Subject_Path,'/',Subject_ID,'_FS/mri/t1_class.nii'],'./Electrode_Locations/','l',13,0.3);
get_mask_V3(Subject_ID,[Subject_Path,'/',Subject_ID,'_FS/mri/t1_class.nii'],'./Electrode_Locations/','r',13,0.3);

pre_projected_electrode_loc = load(['./Electrode_Locations/electrodes_loc_',Subject_ID,'.mat']);


% Single array of electrodes (projected to the closest point of the
% cortical surface)
for i_strips = 1:length(Strip_Names)
    if pre_projected_electrode_loc.elecmatrix(Strip_Lower(i_strips),1)>0
        [out_els] = electrodes2surf(Strip_Names{i_strips},0,2,[Strip_Lower(i_strips):Strip_Upper(i_strips)],['./Electrode_Locations/electrodes_loc_',Subject_ID,'.mat'],['./Electrode_Locations/',Subject_ID,'_surface1_13_03_r.img'],[Subject_Path,'/',Subject_ID,'_FS/mri/T1.nii']);
    else
        [out_els] = electrodes2surf(Strip_Names{i_strips},0,2,[Strip_Lower(i_strips):Strip_Upper(i_strips)],['./Electrode_Locations/electrodes_loc_',Subject_ID,'.mat'],['./Electrode_Locations/',Subject_ID,'_surface1_13_03_l.img'],[Subject_Path,'/',Subject_ID,'_FS/mri/T1.nii']);
    end
    elecmatrix(Strip_Lower(i_strips):Strip_Upper(i_strips),:) = out_els;
    eleccell{i_strips} = out_els;
end

% Grids of electrodes (projected to the vector perpendicular to defined by
% principal component analysis
for i_grids = 1:length(Grid_Names)
    if pre_projected_electrode_loc.elecmatrix(Grid_Lower(i_grids),1)>0
        [out_els] = electrodes2surf(Grid_Names{i_grids},5,1,[Grid_Lower(i_grids):Grid_Upper(i_grids)],['./Electrode_Locations/electrodes_loc_',Subject_ID,'.mat'],['./Electrode_Locations/',Subject_ID,'_surface1_13_03_r.img'],[Subject_Path,'/',Subject_ID,'_FS/mri/T1.nii']);
    else
        [out_els] = electrodes2surf(Grid_Names{i_grids},5,1,[Grid_Lower(i_grids):Grid_Upper(i_grids)],['./Electrode_Locations/electrodes_loc_',Subject_ID,'.mat'],['./Electrode_Locations/',Subject_ID,'_surface1_13_03_l.img'],[Subject_Path,'/',Subject_ID,'_FS/mri/T1.nii']);
    end
    eleccell{i_grids+i_strips} = out_els;
    elecmatrix(Grid_Lower(i_grids):Grid_Upper(i_grids),:) = out_els;
end

% Saving the electrode positions
save(['./Electrode_Locations/electrodes_projected_',Subject_ID,'.mat'],'elecmatrix','eleccell');

%%
figure;
switch cortical_disp
    
% Step 5a
% Displays the cortical surface using 1 of 2 methods
    case 'orig'
        gen_cortex_click_V3([Subject_ID,'_R'],0.4,[15 3],'r'); % Generates right cortex
        gen_cortex_click_V3([Subject_ID,'_L'],0.4,[15 3],'l'); % Generates left cortex
        
        % Need to load the cortex
        figure;
        load(['./Electrode_Locations/',Subject_ID,'_R_cortex.mat']);
        
        % Outputs the surface rendering
        ctmr_gauss_plot(cortex,[0 0 0],0); % right side
        
        load(['./Electrode_Locations/',Subject_ID,'_L_cortex.mat']);
        
        % Outputs the surface rendering
        ctmr_gauss_plot(cortex,[0 0 0],0); % left side

    case 'freesurfer'
        % Generating individual hemispheres
        [cortex.vert_lh,cortex.tri_lh]= read_surf('./surf/lh.pial'); % Reading left side pial surface
        [cortex.vert_rh,cortex.tri_rh]= read_surf('./surf/rh.pial'); % Reading right side pial surface
        
        % Generating entire cortex
        cortex.vert = [cortex.vert_lh; cortex.vert_rh]; % Combining both hemispheres
        cortex.tri = [cortex.tri_lh; (cortex.tri_rh + length(cortex.vert_lh))]; % Combining faces (Have to add to number of faces)
        
        cortex.tri=cortex.tri+1; % freesurfer starts at 0 for indexing
        
        % Reading in MRI parameters
        f=MRIread('./mri/orig.mgz');
        
        % Translating into the appropriate space
        for k=1:size(cortex.vert,1)
            a=f.vox2ras/f.tkrvox2ras*[cortex.vert(k,:) 1]';
            cortex.vert(k,:)=a(1:3)';
        end
        
        % Displaying cortex
        Hp = patch('vertices',cortex.vert,'faces',cortex.tri(:,[1 3 2]),...
            'facecolor',[.65 .65 .65],'edgecolor','none',...
            'facelighting', 'gouraud', 'specularstrength', .25);
        camlight('headlight','infinite');
        fh(1)=gcf;
        vertnormals = get(Hp,'vertexnormals');
        axis off; axis equal
end

%%

% Step 6
% Add the electrodes onto the surface
electrode_names = [Strip_Names,Grid_Names];

% Allows the user to right click and display only selected electrodes
hcmenu = uicontextmenu;

% Loops through all of the electrodes and plots them
for i_plot = 1:numel(eleccell)
    electrodeH{i_plot} = el_add(eleccell{i_plot},'b',80); % Plots the electrodes
    
    % Sets the handles for each of the electrodes and enables them to be
    % visible or not
    hcb{i_plot} = ['a=get(electrodeH(' num2str(i_plot) '),''visible''); if strcmp(a,''on'')==1, set(electrodeH(' num2str(i_plot) '),''visible'',''off''), else set(electrodeH(' num2str(i_plot) '),''visible'',''on''), end'];
    
    cm_label{i_plot} = uimenu(hcmenu,'Label',electrode_names{i_plot},'Callback',hcb{i_plot});
end


% set([fh(1); allchild(fh(1)); Hp; electrodeH'],'uicontextmenu',hcmenu);


% Adding labels to electrodes
labels_choice = questdlg('Would you like to add labels?','Yes','No');
switch labels_choice
    case 'Yes'
        % Plotting text labels
        for k_labels=1:length(elecmatrix(:,1))
            th(k_labels)=text(elecmatrix(k_labels,1)*1.01,elecmatrix(k_labels,2)*1.01,elecmatrix(k_labels,3)*1.01,num2str(k_labels),'FontSize',13,'HorizontalAlignment','center','VerticalAlignment','middle');
        end
end

 