%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          LAB 1 - SCRIPT 
%
% Javier Lopez Iniesta Diaz del Campo
%         Mathias Näreaho
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close('all')

%% Basis functions

sz = 128;

q1 = [[5,9];[9,5];[17,9];[17,121];[5,1];[62,1]];

for i = 1:6
    pq = q1(i,:);
    figure('Name','Basic functions');
    sgtitle("Plot " + i + ": (p,q) = (" + pq(1) + ", " + pq(2) + ")");
    fftwave(pq(1),pq(2),sz)
end

%% Linearity

close('all')

% Rectangular test images
F = [zeros(56, 128); ones(16, 128); zeros(56, 128)];
G = F';
H = F + 2*G;

% Spectra
Fhat = fft2(F);
Ghat = fft2(G);
Hhat = fft2(H);

figure('Name','Linearity','NumberTitle','off');
subplot(3, 3, 1);
    showgrey(F);
    title('F')
subplot(3, 3, 2);
    showgrey(log(1 + abs(Fhat)));
    title('Fhat with log')
subplot(3, 3, 3);
    showgrey(log(1 + abs(fftshift(Fhat))));
    title('Fhat with log and ffftshift')
subplot(3, 3, 4);
    showgrey(G);
    title('G')
subplot(3, 3, 5);
    showgrey(log(1 + abs(Ghat)));
    title('Ghat with log')
subplot(3, 3, 6);
    showgrey(log(1 + abs(fftshift(Ghat))));
    title('Ghat with log and ffftshift')
subplot(3, 3, 7);
    showgrey(H);
    title('H')
subplot(3, 3, 8);
    showgrey(log(1 + abs(Hhat)));
    title('Hhat with log')
subplot(3, 3, 9);
    showgrey(log(1 + abs(fftshift(Hhat))));
    title('Hhat with log and ffftshift')

figure('Name','Linearity - logarithm','NumberTitle','off');
subplot(1, 2, 1);
    showgrey(abs(Fhat(1:30,1:20)));
    title('Fhat')
subplot(1, 2, 2);
    showgrey(log(1 + abs(Fhat(1:30,1:20))));
    title('Fhat with log')

%% Multiplication

F = [zeros(56, 128); ones(16, 128); zeros(56, 128)];
G = F';

% Zhat with multiplication in the spatial domain
Zhat_1 = fft2(F.* G);

% Zhat with convolution in the Fourier domain
Fhat = fft2(F);
Ghat = fft2(G);
Zhat_2 = conv2(Fhat,Ghat);

M = length(Fhat); 
N = length(Ghat); 
Zhat_2_norm = Zhat_2(1:M, 1:N)/(M*N);

figure('Name','Multiplication','NumberTitle','off');
subplot(1, 2, 1);
    showfs(Zhat_1);
    title('Zhat with multiplication')
subplot(1, 2, 2);
    showfs(Zhat_2_norm);
    title('Zhat with convolution')
 
%% Scaling

% Test image
F_scale = [zeros(60, 128); ones(8, 128); zeros(60, 128)] .* [zeros(128, 48) ones(128, 32) zeros(128, 48)];
F_scale_hat = fft2(F_scale);

figure('Name','Scaling','NumberTitle','off');
subplot(2, 2, 1);
    showgrey(F.* G);
    title('Original test image')
subplot(2, 2, 2);
    showfs(Zhat_1);
    title('Fourier transform of the original test image')
subplot(2, 2, 3);
    showgrey(F_scale);
    title('Test image scaled')
subplot(2, 2, 4);
    showfs(F_scale_hat);
    title('Fourier transform of the test image scaled')

%% Rotation

alpha = [30, 45, 60, 90];

nplots = length(alpha) + 1;
% nplots = 2;

figure('Name','Rotation','NumberTitle','off');
subplot(nplots, 2, 1);
    showfs(F_scale_hat);
    title('Spectrum of the original image')
subplot(nplots, 2, 2);
    H_scale_hat = rot(fftshift(F_scale_hat),-0);
    showgrey(log(1 + abs(H_scale_hat)));
    title('Original spectrum rotated back 0º')

for i = 1:4

    G = rot(F_scale, alpha(i));
    G_hat = fft2(G);
    H_hat = rot(fftshift(G_hat), -alpha(i));

    subplot(nplots, 2, 1+2*i);
    showfs(G_hat);
    title(sprintf('Spectrum rotated with alpha = %d', alpha(i)))
    subplot(nplots, 2, 2+2*i);
    showgrey(log(1 + abs(H_hat)));
    title(sprintf('Spectrum back with alpha = %d', alpha(i)))  
end

%% Information in Fourier phase and magnitude

image1 = phonecalc128;
% image2 = few128;
% image3 = nallo128;
a = exp(-10);

image1_function1 = pow2image(image1, a);
image1_function2 = randphaseimage(image1);

figure('Name','Information in Fourier phase and magnitude','NumberTitle','off');
subplot(1, 3, 1);
    showgrey(image1);
    title('Original image')
subplot(1, 3, 2);
    showgrey(image1_function1);
    title('Image with pow2image function')
subplot(1, 3, 3);
    showgrey(image1_function2);
    title('Image with randphaseimage function')

%% Gaussian 

clear
clc
close('all')

t = [0.1, 0.3, 1.0, 10.0, 100.0];
t_img = [1.0, 4.0, 16.0, 64.0, 256.0];

img = phonecalc128;

figure('Name','Impulse response','NumberTitle','off');
for i=1:length(t)
    psf = gaussfft(deltafcn(128, 128), t(i));
    subplot(1, 5, (i-1)+1);
    showgrey(psf);
    cov = variance(psf);
    var = cov(1,1);
    title(sprintf('t = %1.1f \n variance = %1.3f', t(i), var));
    fprintf("t = %1.1f ; variance = %1.3f\n", t(i), var);
    fprintf("Covariance: \n");
    disp(cov)
end 

figure('Name','Image filtered with a Gaussian filter','NumberTitle','off');
% sgtitle("phonecalc128 image with different values of t");
for i=1:length(t_img)
    img_filtered = gaussfft(img, t_img(i));
    subplot(1, 5, (i-1)+1);
    showgrey(img_filtered);
    title(sprintf('t = %1.1f', t_img(i)));
end 

%% Smoothing
close('all')
office = office256;
office_hat = real(fft2(office));

% Filter parameters
variance = 1;
ws = 3; % Window size
fc = 0.25; % Cut-off frequency

% Add noise
add = gaussnoise(office, 16);

figure('Name','Gaussian noise','NumberTitle','off');
subplot(1, 5, 1);
    showgrey(office);
    title("Original image");
subplot(1, 5, 2);
    showgrey(add);
    title("Image with Gaussian noise");
subplot(1, 5, 3);
    % Gaussian smoothing
    add_filtered_gaussfft = gaussfft(add, variance);
    showgrey(add_filtered_gaussfft);
    title(sprintf("Gaussian filter with t=%1.1f", variance));
subplot(1, 5, 4);
    % Median filtering
    add_filtered_median = medfilt(add, ws);
    showgrey(add_filtered_median);
    title(sprintf("Median filter with ws=%1.0f", ws));
subplot(1, 5, 5);
    % Ideal low-pass filtering
    add_filtered_ideal = ideal(add, fc);
    showgrey(add_filtered_ideal);
    title(sprintf("Ideal low pass filter with fc=%1.1f", fc));

variance = 6;
ws = 4; % Window size
fc = 0.1; % Cut-off frequency

sap = sapnoise(office, 0.1, 255);

figure('Name','Sap noise','NumberTitle','off');
subplot(1, 5, 1);
    showgrey(office);
    title("Original image");
subplot(1, 5, 2);
    showgrey(sap);
    title("Sap noise");
subplot(1, 5, 3);
    % Gaussian smoothing
    sap_filtered_gaussfft = gaussfft(sap, variance);
    showgrey(sap_filtered_gaussfft);
    title(sprintf("Gaussian filter with t=%1.1f", variance));
subplot(1, 5, 4);
    % Median filtering
    sap_filtered_median = medfilt(sap, ws);
    showgrey(sap_filtered_median);
    title(sprintf("Median filter with ws=%1.0f", ws));
subplot(1, 5, 5);
    % Ideal low-pass filtering
    sap_filtered_ideal = ideal(sap, fc);
    showgrey(sap_filtered_ideal);
    title(sprintf("Ideal low pass filter with fc=%1.1f", fc));


%% Smoothing and subsampling

close('all')

img = phonecalc256;

smoothimg_gauss = img;
smoothimg_low = img;

% Parameters
N=5;
variance = 0.6;
fc = 0.3;

figure('Name','Smoothing and subsampling','NumberTitle','off');
for i=1:N
    if i>1 % generate subsampled versions
        img = rawsubsample(img);
        smoothimg_gauss = gaussfft(smoothimg_gauss, variance);
        smoothimg_gauss = rawsubsample(smoothimg_gauss);
        smoothimg_low = ideal(smoothimg_low, fc);
        smoothimg_low = rawsubsample(smoothimg_low);
    end

    subplot(3, N, i)
        showgrey(img)
        title(sprintf("Iteration %1.0f without filter", i));
    subplot(3, N, i+N)
        showgrey(smoothimg_gauss)
        title(sprintf("Iteration %1.0f with Gaussian filter", i));
    subplot(3, N, i+2*N)
        showgrey(smoothimg_low)
        title(sprintf("Iteration %1.0f without low-pass filter", i));

%     if i==4
%         figure('Name','Smoothing and subsampling: Iteration 4','NumberTitle','off');
%         subplot(1, 3, 1)
%             showgrey(img)
%             title("Without filter");
%         subplot(1, 3, 2)
%             showgrey(smoothimg_gauss)
%             title("Gaussian filter");
%         subplot(1, 3, 3)
%             showgrey(smoothimg_low)
%             title("Low-pass filter");
%     end

end
