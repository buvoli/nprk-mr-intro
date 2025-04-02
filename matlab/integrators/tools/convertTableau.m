function [outputArg1, outputArg2, name] = convertTableau(coeffHandle, format)
%CONVERTTABLEAU returns NP-RK coefficients in desired format
%   [Al, Ad] = convertTableau(coeffHandle, 'sparse')
%   [A, b] = convertTableau(coeffHandle, 'full')

    [outputArg1, outputArg2, name] = coeffHandle();
        
    switch(format)
        case 'sparse'
           [outputArg1, outputArg2] = toSparse(outputArg1, outputArg2);
        case 'full'
           [outputArg1, outputArg2] = toFull(outputArg1, outputArg2);
    end

end

function [arg1, arg2] = toSparse(arg1, arg2)
            
    if(numel(size(arg1)) == 3 && numel(size(arg2)) == 2) % arguments describe full format
        [arg1, arg2] = diAtoSparse(arg1, arg2);
    elseif(iscell(arg1) && iscell(arg2)) % arguments describe sparse format
        return;
    else
        error('invalid coefficient format');
    end

end

function [arg1, arg2] = toFull(arg1, arg2)
            
    if(numel(size(arg1)) == 3 && numel(size(arg2)) == 2) % arguments describe full format
        return
    elseif(iscell(arg1) && iscell(arg2)) 
        [arg1, arg2] = sparsetodiA(arg1, arg2); % arguments describe sparse format
    else
        error('invalid coefficient format');
    end

end
