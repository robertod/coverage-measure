function S = fps_geodesic(D, n, seed)
    % FPS_GEODESIC Sample n points using a farthest point sampling strategy w.r.t.
    % distance matrix D by using farthest point sampling.
    %
    % Inputs:
    %   D - #v by #v dense distance matrix of geodesics (or any unit distance)
    %   n - no. of points sampled
    %   seed - the index of the initial starting point
    % Outputs:
    %   S - #n by 1 indices of sampled points
    %

    if n > size(D, 1)
        n = size(D, 1);
        warning('The no. of points to sample is greater than the no. of points on the mesh.');
    end

    S = zeros(n,1);
    S(1) = seed;
    d = D(seed,:);

    for it=2:n
        [~,m] = max(d);
        S(it) = m(1);
        d = min(D(S(it),:),d);
    end
end