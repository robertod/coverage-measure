function shape = load_obj(filename)
% LOAD_OBJ load an .obj file.
%
%   shape = load_obj(filename);
%
% Inputs:
%   filename - location of file
%
% Outputs:
%   shape.VERT - node vertexinatates
%   shape.TRIV - list of facesangle elements
%   shape.n - no. of vertices
%   shape.m - no. of faces
%
%   Copyright (c) 2008 Gabriel Peyre

fid = fopen(filename);
if fid<0
    fid = fopen([char(filename),'.obj']);
    if fid<0
        error(['Cannot open ' filename '.']);
    end
end

shape = struct('VERT',[],'TRIV',[],'n',0,'m',0);

frewind(fid);

shape.VERT = [];
shape.TRIV = [];
while 1
    s = fgetl(fid);
    if ~ischar(s)
        break;
    end
    if ~isempty(s) && strcmp(s(1), 'f')
        % face
        shape.TRIV(end+1,:) = sscanf(s(3:end), '%d %d %d');
    end
    if ~isempty(s) && strcmp(s(1), 'v')
        % vertex
        shape.VERT(end+1,:) = sscanf(s(3:end), '%f %f %f');
    end
end
fclose(fid);

shape.n = size(shape.VERT,1);
shape.m = size(shape.TRIV,1);