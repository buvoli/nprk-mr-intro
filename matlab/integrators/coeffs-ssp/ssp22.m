function [A, b, c] = ssp22()
%SSP22 Tableau for second-order, 2-stage SSP2 method

    A = [
        0 0;
        1 0;
    ];
    
    b = [1 1] / 2;
    c = sum(A, 2);

end

