function [success] = ckmkdir(dir_path, call_mkdir)
%CKMKDIR creates a new directory or simply returns true if it already exists.
% returns true if dir exists or if it was created sucessfully.
if(nargin < 2)
    call_mkdir = true;
end
success = exist(dir_path, 'dir');
if(~success && call_mkdir)
    success = mkdir(dir_path);
end
end

