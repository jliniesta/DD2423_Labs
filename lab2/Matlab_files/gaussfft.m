function f_out = gaussfft(f_in, t)
% INPUTS:
% f_in: image
% t: variance
% 
% OUTPUT: 
% f_out: filtered image with a gaussian filter

[x_dim, y_dim] = size(f_in);
[x, y] = meshgrid(-((x_dim/2)-1):(x_dim/2), -((y_dim/2)-1):(y_dim/2));

% 1. Generate a filter based on a sampled version of the Gaussian function.
g=(1/(2*pi*t))*exp(-(x.^2+y.^2)/(2*t));

% 2. Fourier transform the original image and the Gaussian filter
Ghat = fft2(fftshift(g));
F_in_hat = fft2(f_in);

% 3. Multiply the Fourier transforms.
F_out_hat = Ghat.*F_in_hat;

% 4. Invert the resulting Fourier transform.
f_out = real(ifft2(F_out_hat));

end