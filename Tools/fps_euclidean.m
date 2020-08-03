function S = fps_euclidean(V, n, seed)
    % FPS_EUCLIDEAN Sample n points using a farthest point sampling strategy w.r.t.
    % the distance between vertices by using farthest point sampling.
    %
    % Inputs:
    %   V - #V by 3 list of vertices
    %   n - no. of points sampled
    %   seed - the index of the initial starting point
    % Outputs:
    %   S - #n by 1 indices of sampled points
    %

    if n > size(V, 1)
        n = size(V, 1);
        warning('The no. of points to sample is greater than the no. of points on the mesh.');
    end

    S = zeros(n,1);
    S(1) = seed;
    d = pdist2(V,V(seed,:));

    for i=2:n
        [~,m] = max(d);
        S(i) = m(1);
        d = min(pdist2(V,V(S(i),:)),d);
    end
end