function [A, b, c] = SDIRK3_Sa()
%SDIRK3_Sa Stiffly accurate, L-stable three stage SDIRK method

lm = 0.43586652150845899941601945119356;

A = [
    lm                                  0                               0;
    (1 - lm) / 2                        lm                              0;
    (-6 * lm^2 + 16 * lm - 1) / 4       (6 * lm^2 - 20 * lm + 5) / 4    lm;
    ];

b = A(3,:);

c = sum(A, 2);

end

