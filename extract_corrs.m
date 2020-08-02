% EXTRACT_CORRS Convert deformed source model to correspondence and view 
% correspondences of shape pairs.
% Filename format of deformed source mesh: sourcename_targetname.obj
%
% Author Roberto M. Dyke

addpath('Tools;Tools/projection');

models_dir = 'original_models'; % original dataset models
deformed_models_dir = 'mesh_files'; % deformed source mesh input
submission_dir = 'mat_files'; % correspondences

draw_figures = false; % view debug tools

matches = dir(fullfile(deformed_models_dir,'*.obj'));

assert(numel(matches) > 0, ['No deformed models found in directory: ', deformed_models_dir]);

errors = cell(1,numel(matches));
for fi=1:numel(matches)
    % get names of file pairs
    expression = '([a-zA-Z]*(_a|_b)(_|\.)|[a-zA-Z]*)';
    chars = regexp(matches(fi).name,expression,'match');
    fn1 = strip(chars{1},'right','_');
    fn2 = strip(chars{2},'right','.');
    fprintf('%s,%s\n',fn1,fn2);
    
    % load mesh data
    M = load_obj(fullfile(models_dir,fn1));
    N = load_obj(fullfile(models_dir,fn2));
    O = load_obj(fullfile(deformed_models_dir,[fn1,'_',fn2]));
    
    % compute correspondence
    corr = [(1:O.n)',knnsearch(N.VERT,O.VERT)]; % vertex to vertex
    baryc_corr = points2barycentric(N.VERT,N.TRIV,O.VERT); % vertex to point on face
    
    save(fullfile(submission_dir,[fn1,'_',fn2]),'M','N','O','corr','baryc_corr');
    
    if draw_figures
        figure;
        subplot(1,3,1);
        trisurf(O.TRIV,O.VERT(:,1),O.VERT(:,2),O.VERT(:,3),'EdgeAlpha',0,'FaceColor','g');
        hold on;
        trisurf(N.TRIV,N.VERT(:,1),N.VERT(:,2),N.VERT(:,3),'EdgeAlpha',.5,'FaceAlpha',0.5,'FaceColor','r');
        camlight; lighting gouraud; axis equal; axis off; 
        hold off;
        title('Deformed source and target mesh');

        subplot(1,3,2);
        trisurf(M.TRIV,M.VERT(:,1),M.VERT(:,2),M.VERT(:,3),'EdgeAlpha',0,'FaceColor','g');
        hold on;
        trisurf(N.TRIV,N.VERT(:,1),N.VERT(:,2),N.VERT(:,3),'EdgeAlpha',0,'FaceColor','r');
        camlight; lighting gouraud; axis equal; axis off;
        for it=1:25:size(corr,1)
            pts = [M.VERT(corr(it,1),:);N.VERT(corr(it,2),:)];
            plot3(pts(:,1),pts(:,2),pts(:,3),'b-');
        end
        hold off;
        title('Vertex-to-vertex correspondences');

        M_point = M.VERT(baryc_corr(:,1) > 0,:);
        N_point = N.VERT(N.TRIV(baryc_corr(baryc_corr(:,1)>0,1),1),:).*baryc_corr(baryc_corr(:,1)>0,2)+N.VERT(N.TRIV(baryc_corr(baryc_corr(:,1)>0,1),2),:).*baryc_corr(baryc_corr(:,1)>0,3)+N.VERT(N.TRIV(baryc_corr(baryc_corr(:,1)>0,1),3),:).*baryc_corr(baryc_corr(:,1)>0,4);

        assert(size(M_point,1) == size(N_point,1),'Unexpected error.');

        subplot(1,3,3);
        trisurf(M.TRIV,M.VERT(:,1),M.VERT(:,2),M.VERT(:,3),'EdgeAlpha',0,'FaceColor','g');
        hold on;
        trisurf(N.TRIV,N.VERT(:,1),N.VERT(:,2),N.VERT(:,3),'EdgeAlpha',0,'FaceColor','r');
        camlight; lighting gouraud; axis equal; axis off;
        for it=1:25:size(M_point,1)
            pts = [M_point(it,:);N_point(it,:)];
            plot3(pts(:,1),pts(:,2),pts(:,3),'b-');
        end
        hold off;
        title('Vertex-to-point correspondences');
        drawnow;
    end
end