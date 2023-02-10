function [linepar, acc] = houghline(curves, magnitude, nrho, ntheta, threshold, nlines, verbose)

% Check if input appear to be valid

% Allocate accumulator space
acc = zeros(nrho,ntheta);

% Define a coordinate system in the accumulator space
x_max = size(magnitude,1);
y_max = size(magnitude,2);

R = sqrt(x_max.^2 + y_max.^2);
r_space = linspace(-R,R,nrho);
r_step = r_space(2)-r_space(1);
theta_space = linspace(-pi/2,pi/2, ntheta);

% Loop over all the input curves (cf. pixelplotcurves)
insize = size(curves,2);
trypointer=1;
numcurves=0;

% For each point on each curve
while trypointer <= insize
    polylength = curves(2, trypointer);
    numcurves = numcurves + 1;
    trypointer = trypointer + 1;
    
    for polyidx = 1:polylength
        x = curves(2, trypointer);
        y = curves(1, trypointer);
        
        trypointer = trypointer + 1;
        % Check if valid point with respect to threshold
        % Optionally, keep value from magnitude image
        if abs(magnitude(round(x),round(y))) > threshold
            % Loop over a set of theta values
            for theta_idx = 1:ntheta
                % Compute rho for each theta value
                r = x*cos(theta_space(theta_idx))+ y*sin(theta_space(theta_idx));
                % Compute index values in the accumulator space
                r_idx = find(r_space<r+r_step/2,1,'last');

                % Update the accumulator
                acc(r_idx,theta_idx) = acc(r_idx,theta_idx) + 1;

                % Question 10: Choice of accumulator incrementation function
%                 acc(r_idx,theta_idx) = acc(r_idx,theta_idx) + log(abs(magnitude(round(x),round(y))));

            end
        end
    end
end

% Extract local maxima from the accumulator
acctmp = acc;
[pos, value] = locmax8(acctmp);
[dummy, indexvector] = sort(value);
nmaxima = size(value, 1);

% Delimit the number of responses if necessary
if nmaxima < nlines
    nlines = nmaxima;
end
fprintf("The number of lines will be: %1.0f \n", nlines);

% Compute a line for each one of the strongest responses in the accumulator
outcurves = zeros(2, 4*nlines);

% Overlay these curves on the gradient magnitude image
for idx = 1:nlines

    % Extract the index values in the accumulator that correspond to the nlines strongest responses
    rhoidxacc = pos(indexvector(nmaxima - idx + 1), 1);
    thetaidxacc = pos(indexvector(nmaxima - idx + 1), 2);
    r = r_space(rhoidxacc);
    theta = theta_space(thetaidxacc);
    
    x0 = 0;
    y0 = r/sin(theta);
    dx = R^2;
    dy = (-cos(theta)/sin(theta)) * dx + r/sin(theta);  
    
    outcurves(1, 4*(idx-1) + 1) = 0; % level, not significant
    outcurves(2, 4*(idx-1) + 1) = 3; % number of points in the curve
    outcurves(2, 4*(idx-1) + 2) = x0-dx;
    outcurves(1, 4*(idx-1) + 2) = y0-dy;
    outcurves(2, 4*(idx-1) + 3) = x0;
    outcurves(1, 4*(idx-1) + 3) = y0;
    outcurves(2, 4*(idx-1) + 4) = x0+dx;
    outcurves(1, 4*(idx-1) + 4) = y0+dy;
end

% Return the output data
linepar = outcurves;

end