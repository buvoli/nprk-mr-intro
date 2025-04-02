function exportFigure(figure, save_options)
%SAVEPLOT Prints file as specified by save_options
% Useful Documentation
% http://www.mathworks.com/help/matlab/ref/figure_props.html
% http://www.mathworks.com/help/matlab/ref/print.html

% -- Set Options -------------------------------------------------------------------------------------------------------
field_value_pairs = { ...
    {'Renderer',        'zbuffer'} ...
    {'PaperUnits',      'centimeters'} ...
    {'PaperPosition',   [0 0 15 15]} ...
    {'dpi',             500} ...
    {'SavePath',        'figure_output'} ...
    {'Format',          'png'} ...
    {'CallMakeDir',     true} % will create save directory structure if it does not exist
    };
if(nargin == 1)
    save_options = struct();
end
save_options = setDefaultOptions(save_options, field_value_pairs);
% ----------------------------------------------------------------------------------------------------------------------

if(filepathExists(save_options.SavePath, save_options.CallMakeDir))    
    if(strcmp(save_options.Format, 'fig'))
        saveas(figure, strcat(save_options.SavePath, '.fig'));
    else
        set(figure, 'renderer', save_options.Renderer);
        set(figure, 'PaperUnits', save_options.PaperUnits);
        set(figure, 'PaperPosition', save_options.PaperPosition);
        set(figure, 'PaperPositionMode', 'manual'); % force size to PaperPosition specification
        if(strcmp(save_options.Format, 'eps'))
            print(figure, '-depsc2', sprintf('-r%d', save_options.dpi), strcat(save_options.SavePath, '.eps'));
        elseif(strcmp(save_options.Format, 'jpg'))
            print(figure, '-djpeg', sprintf('-r%d', save_options.dpi), strcat(save_options.SavePath, '.jpg'));
        elseif(strcmp(save_options.Format, 'png'))
            print(figure, '-dpng', sprintf('-r%d', save_options.dpi), strcat(save_options.SavePath, '.png'));
        elseif(strcmp(save_options.Format, 'pdf')) % Matlab PDF save works poorly. If possible use epstopdf ------------
            appendPathForEPS();
            [status_eps, ~] = system('epstopdf --version');
            [status_gs, ~] = system('gs --version');
            if(status_eps == 0 && status_gs == 0)
            	temp_file = strcat(save_options.SavePath, '.tmp.eps');
                if(isfile(temp_file))
                    error('export figure cannot create temp file since it already exists');
                else
                    print(figure, '-depsc2', sprintf('-r%d', save_options.dpi), temp_file);
                    str_command = ['epstopdf "', temp_file, '" --outfile "', strcat(save_options.SavePath, '.pdf'), '"'];
                    [status, ~] = system(str_command); % call epstopdf                    
                    delete(temp_file); % remove temp file
                    if(status ~= 0)
                        error('error creating pdf using epstopdf');
                    end                    
                end                
            else
                warning('eps2pdf or gs command not found, Using default MATLAB pdf renderer. To resolve, update local appendPathForEPS() function with correct path.')
                print(figure, '-dpdf', sprintf('-r%d', save_options.dpi), strcat(save_options.SavePath, '.pdf'));
            end
        end
        
    end
    
else
    warning('Could not save figure! Could not create directory:%s', fileparts(save_options.SavePath));
end

end

function appendPathForEPS() % modifies bash path variable so that it system() can call epstopdf and gs commands

 epstopdf_path = '/Library/TeX/texbin';  % path to command epstopdf
 gs_path       = '/usr/local/bin';       % path to ghostscript required by epstopdf
 cur_path_str  = getenv('PATH');        % current path 
 cur_paths     = split(cur_path_str, ':');
 
 if(~ismember(epstopdf_path,cur_paths))
    setenv('PATH', [getenv('PATH'), ':', epstopdf_path]);
 end
 
 if(~ismember(gs_path,cur_paths))
    setenv('PATH', [getenv('PATH'), ':', gs_path]);
 end

end