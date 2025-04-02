function [A, b, c] = rk4()
%RK4 Clasical fourth-order RK

  A = [
      0     0       0   0;
      1 / 2 0       0   0;
      0     1 / 2   0   0;
      0     0       1   0
  ];

  b = [1 2 2 1] / 6;

  c = sum(A, 2);

end

