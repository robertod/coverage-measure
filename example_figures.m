% EXAMPLE_FIGURES code used to reproduce Fig. 7 & 8 in the associated
% report [1]. Explanation of these figures can be found in Section 5.
%
% Dyke RM, et al. Shape Correspondence with Non-Isometric Deformations. C&G 2020.
%
% Author Roberto M. Dyke

addpath('Tools');

name = 'leopard'; % examplar shape name
output_dir = 'Output';
if ~isfolder(output_dir)
   mkdir(output_dir)
end

N = load_obj(fullfile('Data',[name,'.obj']));
N.tri_areas = face_areas(N.VERT,N.TRIV);
N.vert_areas = accumarray(N.TRIV(:),repmat(N.tri_areas/3,[3,1])); % assumes all vertices are references in face list

%% Generage correspondences with different characteristics
% 1. bijective
ids = 1:N.n;

fid = fopen(fullfile(output_dir,'bijective.txt'),'w');
fprintf(fid,'%f\n',ids);
fclose(fid);

% 2. part1
part = ones(N.n,1);
part(knnsearch(N.VERT,N.VERT(1,:),'k',floor(N.n*0.1878))) = 0;
ids = find(part);

fid = fopen(fullfile(output_dir,'part1.txt'),'w');
fprintf(fid,'%f\n',ids);
fclose(fid);

% 3. part2
part = ones(N.n,1);
part(knnsearch(N.VERT,N.VERT(1,:),'k',floor(N.n*0.4062))) = 0;
ids = find(part);

fid = fopen(fullfile(output_dir,'part2.txt'),'w');
fprintf(fid,'%f\n',ids);
fclose(fid);

% 4. sparse1
ids = fps_euclidean(N.VERT, floor(N.n*0.37), 1);

fid = fopen(fullfile(output_dir,'seed1.txt'),'w');
fprintf(fid,'%f\n',ids);
fclose(fid);

% 5. sparse2
ids = fps_euclidean(N.VERT, floor(N.n*0.1777), 1);

fid = fopen(fullfile(output_dir,'seed2.txt'),'w');
fprintf(fid,'%f\n',ids);
fclose(fid);

% 6. sparse3
ids = fps_euclidean(N.VERT, floor(N.n*0.0692), 1);

fid = fopen(fullfile(output_dir,'seed3.txt'),'w');
fprintf(fid,'%f\n',ids);
fclose(fid);

clear fid ids N part;

%% Compute coverage measure
T = ["bijective";"part1";"part2";"seed1";"seed2";"seed3"];

%addpath('toolbox_fast_marching');
%D = compute_distances(N); % geodesics
%ns = int32(linspace(1,N.n,1000));
%C = zeros(numel(ns),N.n); % Voronoi cells
%for jt=1:numel(ns)
%    C(jt,:) = compute_voronoi(N,D,ns(jt));
%end

% load Voronoi data
voronoi_data = load([name,'_voronoi']);
N = voronoi_data.N;
C = voronoi_data.C;
ns = voronoi_data.ns;
resolution = numel(ns);

% compute coverage
curves = cell(numel(T),1);
fig_corrs = figure;
fig_results = figure;
for it=1:numel(T)
    fn1 = char(T(it));
    ids = load(fullfile(output_dir,[fn1,'.txt']));
    
    R = zeros(numel(ns),1);
    for jt=1:numel(ns)
        R(jt) = compute_coverage_measure(N,ids,C(jt,:)');
    end
    
    curves{it} = R*100;
    error_curve = curves{it};
    
    % plot error curve
    figure(fig_results);
    hold on;
    plot(linspace(1,100,resolution),error_curve);
    hold off;
    
    % draw correspondences
    figure(fig_corrs);
    subplot(2,ceil(numel(T)/2),it);
    trisurf(N.TRIV,N.VERT(:,1),N.VERT(:,2),N.VERT(:,3),'FaceColor',[0.8,0.8,0.8],'EdgeAlpha',0.0);
    axis equal; axis off; view(90,0); camroll(90);
    hold on;
    scatter3(N.VERT(ids,1),N.VERT(ids,2),N.VERT(ids,3),10,'filled');
    hold off;
    title(T(it))
end

% format figure
figure(fig_results);
set(gca,'XGrid','on','YGrid','on');
set(gca,'XLim',[0 100],'YLim',[0 100]);
xlabel('% No. of Voronoi cells');
ylabel('% Coverage');
title('Total corrrespondence coverage');
legend(T);

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