% COMPUTE_COVERAGE Compute coverage measure for a set of correspondences.
%
% Notes:
% Results should be stored in a distinc folder with each file named
% sourcename_targetname.mat with the following structure(s):
%
%   data.corr - vertex-to-vertex correspondences
%   or
%   data.baryc_corr - vertex-to-barycentric correspondences
%
% For registration methods use extract_corrs.m to extract
% correspondences from deformed .obj files.
% N.B. please run toolbox_fast_marching/compile_mex.m first.
%
% Author Roberto M. Dyke

addpath('toolbox_fast_marching;toolbox_fast_marching/toolbox;'); % Gabriel Peyre's fast matching implementation
addpath('FMM'); % Ron Kimmel's fast marching implementation (much faster) available here: https://github.com/orlitany/3D_shapes_tools/
addpath('Tools;Tools/projection');

models_dir = 'original_models'; % original dataset models
submission_dir = 'user_results/mat_files/'; % a correspondence/registration method's predicted correspondences
cache_dir = 'cache'; % cache directory for model quantities
cache_geodesics = true; % caching geodeics is very costly

resolution = 1000; % The number of uniformly spaced steps taken to compute Voronoi cells from 1 to the no. of vertices on a mesh

models = dir(fullfile(models_dir,'*.obj'));

assert(numel(models) > 0, ['No models found in directory: ',models_dir])

if ~exist(cache_dir, 'dir')
   mkdir(cache_dir)
end

%% Precompute quantities (vertex weights, geodesics)
for fi=1:numel(models)
    fname = models(fi).name;
    [~,name,~] = fileparts(models(fi).name);
    
    fprintf('%i/%i %s\n',fi,numel(models),fname);
    
    N = load_obj(fullfile(models_dir,fname));
    
    % compute vertex weights
    N.tri_areas = face_areas(N.VERT,N.TRIV);
    
    assert(~any((1:N.n)'-unique(N.TRIV)),'unreferenced vertices exist in mesh');
    N.vert_areas = accumarray(N.TRIV(:),repmat(N.tri_areas/3,[3,1])); % assumes all vertices are references in face list
    
    % check if distances have already been computed
    dist_cached = false;
    if isfile(fullfile(cache_dir,[name,'.mat'])) && cache_geodesics
        data = load(fullfile(cache_dir,[name,'.mat']));
        if all([data.N.VERT;N.VERT] == [N.VERT;data.N.VERT],'all') && ...
                all([data.N.TRIV;N.TRIV]==[N.TRIV;data.N.TRIV],'all') && ...
                isfield(data,'D')
            dist_cached = true;
            D = data.D;
        end
    end
    % compute geodesics
    if ~dist_cached
        D = compute_distances(N);
        %D = compute_distancesFMM(N);
        if cache_geodesics
            save(fullfile(cache_dir,[name,'.mat']),'N','D','-v7.3');
        end
    end
    
    % uniformly sample the no. of Voronoi cells between 1 and the no. of vertices
    ns = int32(linspace(1,N.n,resolution));
    
    % check if Voronoi cells have already been computed
    voronoi_cached = false;
    if isfile(fullfile(cache_dir,[name,'_voronoi.mat']))
        data = load(fullfile(cache_dir,[name,'_voronoi.mat']));
        if all([data.N.VERT;N.VERT] == [N.VERT;data.N.VERT],'all') && ...
                all([data.N.TRIV;N.TRIV]==[N.TRIV;data.N.TRIV],'all') && ...
                all([data.ns,ns] == [ns,data.ns]) && isfield(data,'C')
            voronoi_cached = true;
            C = data.C;
        end
    end
    % compute Voronoi cells
    if ~voronoi_cached
        upd = textprogressbar(numel(ns),...
            'startmsg', 'Computing voronoi cells... ',...
            'endmsg', ' done.');
        C = zeros(numel(ns),N.n); % cells
        for jt=1:numel(ns)
            upd(jt);
            n=ns(jt);
            C(jt,:) = compute_voronoi(N,D,n);
        end
        save(fullfile(cache_dir,[name,'_voronoi']),'C','N','ns','-v7.3');
    end
end


matches = dir(fullfile(submission_dir,'*.mat'));

assert(numel(matches) > 0, ['No matches found in directory: ', submission_dir])

errors = cell(1,numel(matches));
for fi=1:numel(matches)
    % get names of file pairs
    expression = '([a-zA-Z]*(_a|_b)(_|\.)|[a-zA-Z]*)';
    chars = regexp(matches(fi).name,expression,'match');
    fn1 = strip(chars{1},'right','_');
    fn2 = strip(chars{2},'right','.');
    
    % load match data
    match_data = load(fullfile(matches(fi).folder,matches(fi).name));
    corr = match_data.corr;
    
    % when using barycentric coodinates, uncomment this code
    %baryc_corr = match_data.baryc_corr;
    % barycentric correspondence to vertex correspondence
    %[~,idx] = max(baryc_corr(:,2:4),[],2);
    %corr = [(1:size(baryc_corr,1))',N.TRIV(sub2ind(size(N.TRIV),baryc_corr(:,1),idx))];
    
    % load Voronoi data
    voronoi_data = load(fullfile(cache_dir,[fn2,'_voronoi']));
    C = voronoi_data.C;
    ns = voronoi_data.ns;
    assert(match_data.N.n == voronoi_data.N.n,'Meshes differ in size');
    N = voronoi_data.N;
    
    % compute error for each Voronoi segmentation
    R = zeros(numel(ns),1);
    for it=1:numel(ns)
        R(it) = compute_coverage_measure(N,corr(:,2),C(it,:)');
    end
    errors{fi} = R;
end

error_curve = compute_error_curve(errors,resolution);

%% Visualise results
figure, plot(linspace(1,100,resolution),error_curve);
set(gca,'XGrid','on','YGrid','on');
set(gca,'XLim',[0 100],'YLim',[0 100]);
xlabel('% No. of Voronoi cells');
ylabel('% Coverage');
title('Total corrrespondence coverage');

function D = compute_distances(N)
    % COMPUTE_DISTANCES computes geodesic distances using fast matching
    % method implemented by Gabriel Peyre
    %
    % Inputs:
    %   N.VERT, N.TRIV - target surface N
    %
    % Outputs:
    %   D - dense distance matrix
    
    upd = textprogressbar(N.n,...
        'startmsg', 'Computing geodesics... ',...
        'endmsg', ' done.');
    
    D = inf(N.n,N.n);
    for it=1:N.n
        [d,~,~] = perform_fast_marching_mesh(N.VERT, N.TRIV, it);
        D(:,it) = d;
        upd(it);
    end
end

function D = compute_distancesFMM(N)
    % COMPUTE_DISTANCESFMM computes geodesic distances using the fast
    % marching method implemented by Ronny Kimmel
    %
    % Inputs:
    %   N.VERT, N.TRIV - target surface N
    %
    % Outputs:
    %   D - dense distance matrix
    
    upd = textprogressbar(N.n,...
        'startmsg', 'Computing geodesics... ',...
        'endmsg', ' done.');
    
    geo = fastmarchmex('init', int32(N.TRIV-1), double(N.VERT(:,1)), double(N.VERT(:,2)), double(N.VERT(:,3)));

    D = inf(N.n,N.n);
    for it=1:N.n
        upd(it);
        d = fastmarchmex('march', geo, it);
        D(it,:) = d;
    end
    fastmarchmex('deinit', geo);
end

function C = compute_voronoi(N,D,n)
    % COMPUTE_VORONOI segment vertices into voronoi cells
    %
    % Inputs:
    %   N.VERT, N.TRIV - target surface N
    %   D - distance matrix
    %   n - no. of voronoi cells
    %
    % Outputs:
    %   C - cell classification of each vertex
    
    rng(1);
    seed = randi([1,N.n]); % select a (pseudo-)random start point
    S = fps_geodesic(D, n, seed);
    assert(numel(unique(S))==numel(S),'same point selected twice by fps_geodesic');
    
    [~,C] = min(D(S,:));
end


function ratio = compute_coverage_measure(N,idx,Q)
    % COMPUTE_COVERAGE_MEASURE counts the number of
    %
    % Inputs:
    %   N.VERT, N.TRIV - target surface N
    %   idx - verex indices of points on N
    %   Q - per vertex voronoi cell classification
    %   n - no. of voronoi cells
    %
    % Outputs:
    %   ratio - the area covered by correspondences normalised by the total
    %   shape area
    
    matches = false(N.n,1); % indices where there is a correspondence
    matches(idx) = true;
    
    ratio = sum(ismember(Q,unique(Q(matches))) .* N.vert_areas) / sum(N.vert_areas);
end

function curve = compute_error_curve(errors,resolution)
    % COMPUTE_ERROR_CURVE combines the error curves of each shape into one
    %
    % Inputs:
    %   errors - cell array of errors
    

    curve = zeros(resolution,1);
    for it=1:numel(errors)
        curve = curve + errors{it};
    end

    curve = curve / numel(errors)*100;

end