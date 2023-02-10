function Lvv_tilde = Lvvtilde(inpic, shape)

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
    
    delta_xy = conv2(delta_x, delta_y, 'same');
    
    Lx = conv2(inpic, delta_x, 'same');
    Lxx = conv2(inpic, delta_xx, 'same');
    Ly = conv2(inpic, delta_y, 'same');
    Lyy = conv2(inpic, delta_yy, 'same');
    Lxy = conv2(inpic, delta_xy, 'same');
    
    Lvv_tilde = (Lx.^2 .* Lxx)+(2.*Lx.*Ly.*Lxy) + (Ly.^2 .* Lyy);

end