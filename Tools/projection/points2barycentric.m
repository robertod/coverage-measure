function b = points2barycentric(v,f,p,alpha)
    % POINTS2BARYCENTRIC project a set of 3D points onto a triangular
    % mesh and compute barycentric coordinates
    %
    % b = points_to_barycentric(v,f,p)
    % b = points_to_barycentric(v,f,p,alpha)
    %
    % Inputs:
    %   v - #v by 3 vertices
    %   f - #f by 3 faces
    %   p - #p by 3 points
    %   alpha - max projection error distance threshold
    % Output:
    %   b - #p by 4 barycentric coordinates
    %
    % Author Roberto M. Dyke
    
    if nargin < 3 || nargin > 4
        error('Invalid input');
    end
    if size(v,2) ~= 3
        v=v';
    end
    if size(f,2) ~= 3
        f=f';
    end
    if size(p,2) ~= 3
        p=p';
    end
    if nargin <= 3
        alpha = 1.1; % max ratio difference in area between point and face
    end
    
    % estimate face normals
    nf = estimate_face_normals(v,f);
    
    % compute face centroids
    cf = face_centroids(v,f);
    
    % compute face edge vectors
    e1 = v(f(:,2),:) - v(f(:,1),:);
    e2 = v(f(:,3),:) - v(f(:,2),:);
    e3 = v(f(:,1),:) - v(f(:,3),:);
    
    P = sqrt(sum(e1.*e1,2)); % compute edge length
    Q = sqrt(sum(e2.*e2,2));
    R = sqrt(sum(e3.*e3,2));
    S = (P+Q+R)./2;
    area = sqrt(S.*(S-P).*(S-Q).*(S-R)); % Heron's method for area
    
    % project points to surface
    [~,p] = point2trimesh(struct('vertices',v,'faces',f), 'QueryPoints', p);
    
    distances = zeros(size(p,1),1);
    b = -1*ones(size(p,1),4);
    for it=1:size(p,1)
        %fprintf('it: %i\tp: %f, %f, %f\n',it,p(it,:));
        p_i = p(it,:);
        
        %e1 = v(f(:,2),:) - v(f(:,1),:); % precomuted
        %e2 = v(f(:,3),:) - v(f(:,2),:);
        %e3 = v(f(:,1),:) - v(f(:,3),:);
        e4 = p_i - v(f(:,1),:);
        e5 = p_i - v(f(:,2),:);
        e6 = p_i - v(f(:,3),:);
        
        [area1,area2,area3] = point_in_face_area(e1,e2,e3,e4,e5,e6);
        
        [value,idx] = max(min((area1+area2+area3)./area, area./(area1+area2+area3))-1);
        
        
        if value > alpha
            warning('point too distant, skipping point %i\n',it);
            continue
        end
        
        d = abs(dot(p_i-cf(idx,:),nf(idx,:),2)); % distance to centre
        p1 = p_i - d.*nf(idx,:); % project point
        
        % check point is projected onto face
        e1_ = e1(idx,:);
        e2_ = e2(idx,:);
        e3_ = e3(idx,:);
        e4_ = p1 - v(f(idx,1),:);
        e5_ = p1 - v(f(idx,2),:);
        e6_ = p1 - v(f(idx,3),:);
        
        [area1,area2,area3] = point_in_face_area(e1_,e2_,e3_,e4_,e5_,e6_);
        value = min((area1+area2+area3)./area(idx), area(idx)./(area1+area2+area3));
        if value > alpha
            warning('point too distant, skipping point %i\n',it);
            continue
        end
        
        % compute barycentric coordinate
        b(it,:) = [idx, area2/area(idx), 1-(area1+area2)/area(idx), area1/area(idx)];
        
        % recover orginal point
        p2 = v(f(b(it,1),1),:)*b(it,2) + v(f(b(it,1),2),:)*b(it,3) + v(f(b(it,1),3),:)*b(it,4);
        
        distances(it) = vecnorm(p1-p2);
    end
end

function [area1,area2,area3] = point_in_face_area(e1,e2,e3,e4,e5,e6)
    % sub-triangle 1 area
    P1 = sqrt(sum(e1.*e1,2));
    Q1 = sqrt(sum(e5.*e5,2));
    R1 = sqrt(sum(e4.*e4,2));
    S1 = (P1+Q1+R1)./2;
    area1 = sqrt(S1.*(S1-P1).*(S1-Q1).*(S1-R1)); % Heron's method for area

    % sub-triangle 2 area
    P2 = sqrt(sum(e2.*e2,2));
    Q2 = sqrt(sum(e6.*e6,2));
    R2 = sqrt(sum(e5.*e5,2));
    S2 = (P2+Q2+R2)./2;
    area2 = sqrt(S2.*(S2-P2).*(S2-Q2).*(S2-R2)); % Heron's method for area

    % sub-triangle 1 area
    P3 = sqrt(sum(e3.*e3,2));
    Q3 = sqrt(sum(e4.*e4,2));
    R3 = sqrt(sum(e6.*e6,2));
    S3 = (P3+Q3+R3)./2;
    area3 = sqrt(S3.*(S3-P3).*(S3-Q3).*(S3-R3)); % Heron's method for area
    
    % fixes rounding error where S3-P/Q/R3 is slightly less than 0
    area1(imag(area1) ~= 0) = 0;
    area2(imag(area2) ~= 0) = 0;
    area3(imag(area3) ~= 0) = 0;
end

function c = face_centroids(v,f)
    % compute face centroids
    if size(v,2) ~= 3
        v=v';
    end
    if size(f,2) ~= 3
        f=f';
    end
    
    c = (v(f(:,1),:) + v(f(:,2),:) + v(f(:,3),:)) ./ 3;
end

function n = estimate_face_normals(v,f)
    % estimate face normals
    if nargin ~= 2
        error('Invalid input');
    end    
    if size(v,2) ~= 3
        v=v';
    end
    if size(f,2) ~= 3
        f=f';
    end
    
    v0 = v(f(:,1),:) - v(f(:,3),:); % u
    v1 = v(f(:,2),:) - v(f(:,3),:); % v
    
    n = cross(v0,v1);
    
    n = n ./ sqrt(sum(n.*n,2));
end