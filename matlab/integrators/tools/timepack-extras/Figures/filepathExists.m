function flag = filepathExists(full_file_name, call_mkdir) % checks if directory of filename exists. If call_mkdir = true, then path is created if it does not exist.
[file_path] = fileparts(full_file_name);
if(isempty(file_path) || ckmkdir(file_path, call_mkdir))
    flag = true;
else
    flag = false;
end
end

