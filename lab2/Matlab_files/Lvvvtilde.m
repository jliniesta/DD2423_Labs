function Lvvv_tilde = Lvvvtilde(inpic, shape)

    if (nargin < 2)
        shape = 'same';
    end

    % First order derivatives
    delta_x = [0,   0, 0,    0, 0; 
               0,   0, 0,    0, 0; 
               0, 1/2, 0, -1/2, 0;
               0,   0, 0,    0, 0;
               0,   0, 0,    0, 0];
    
    delta_y = delta_x';

    % Second order derivatives
    delta_xx = [0, 0,  0, 0, 0; 
                0, 0,  0, 0, 0; 
                0, 1, -2, 1, 0;
                0, 0,  0, 0, 0;
                0, 0,  0, 0, 0];
    
    delta_yy = delta_xx';
    
    delta_xxx = conv2(delta_x, delta_xx, 'same');
    delta_yyy = conv2(delta_y, delta_yy, 'same');
    delta_xxy = conv2(delta_xx, delta_y, 'same');
    delta_xyy = conv2(delta_x, delta_yy, 'same');

    Lx   = conv2(inpic, delta_x,   shape);
    Lxxx = conv2(inpic, delta_xxx, shape);
    Ly   = conv2(inpic, delta_y,   shape);
    Lyyy = conv2(inpic, delta_yyy, shape);
    Lxxy = conv2(inpic, delta_xxy, shape);
    Lxyy = conv2(inpic, delta_xyy, shape);

    Lvvv_tilde = Lx.^3.*Lxxx + 3.*Lx.^2.*Ly.*Lxxy + 3.*Lx.*Ly.^2.*Lxyy + Ly.^3 .* Lyyy;
    
end