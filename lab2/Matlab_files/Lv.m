function pixels = Lv(inpic, shape)

    if (nargin < 2)
        shape = 'same';
    end

    % Difference operators
    dxmask = [1 0 -1; 2 0 -2; 1 0 -1];
    dymask = dxmask';

    % Discrete derivation approximations
    Lx = conv2(inpic, dxmask, shape);
    Ly = conv2(inpic, dymask, shape);

    % Approximation of the gradient magnitude
    pixels = sqrt(Lx.^2 + Ly.^2);