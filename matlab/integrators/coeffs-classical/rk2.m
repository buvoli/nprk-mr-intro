function [A, b, c] = rk2(type)
%RK2 2 stage optimal ssp2 method

if(nargin < 1)
    type = 'ralston';
end

switch(type)

    case 'ralston' % optimal error

        A = [
            0   0;
            2/3 0;
        ];

        b = [1/4 3/4];
    
    case 'huen' % ssp method
        
        A = [
                0 0;
                1 0;
            ];

        b = [1 1] / 2;

    case 'midpoint'

        A = [
                0   0;
                1/2 0;
            ];

        b = [0 1];

end

c = sum(A, 2);

end
