function A = face_areas(v,f)
    % FACE_AREAS calculate the area of each face
    %
    % Inputs:
    %   v - #v by 3 mesh vertices
    %   f - #f by 3 mesh faces
    %
    % Outputs:
    %   A - #f by 1 face areas
    %
    % Author Roberto M. Dyke
    
    if nargin ~= 2
        error('Invalid input');
    end
    if size(v,2) ~= 3
        v=v';
    end
    if size(f,2) ~= 3
        f=f';
    end
    
    e1 = v(f(:,3),:) - v(f(:,1),:);
    e2 = v(f(:,2),:) - v(f(:,1),:);
    A = 0.5*sqrt(sum(cross(e1,e2).^2,2));
end