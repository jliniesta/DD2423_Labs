function [segmentation, centers] = kmeans_segm(image, K, L, seed)

% INPUTS:
% image
% K: number of cluster centers
% L: number of iterations
% seed: for initializing randomization

% OUTPUT:
% segmentation: with a colour index per pixel
% centers: centers of all clusters in 3D colour space.

rng(seed);

width = size(image,1);
height = size(image,2);

image_vec = double(reshape(image, width*height, 3));
n_pixels = size(image_vec,1);
n_colors = size(image_vec,2);

% Randomly initialize the K cluster centers
centers = zeros(K,3);
dist = zeros(K,n_pixels);
pixel_k = zeros(n_pixels,1);
for k = 1:K
    % Form 1: Assign a random color
%     centers(k,:) = randi(255,[1,n_colors]);
    
    % Form 2: Pick the colors of random pixels of the image
    centers_coord = randi(n_pixels,1,1);
    centers(k,:) = image_vec(centers_coord, :);
    
end

% Compute all distances between pixels and cluster centers
for i = 1:n_pixels
    for k = 1:K
        dist(k,i) = pdist2(image_vec(i,:),centers(k,:));
    end
end

old_pixel = zeros(n_pixels,1);

% Iterate L times
for l = 1:L

    disp("Iteration: " + l)

    % Assign each pixel to the cluster center for which the distance is minimum
    for i = 1:n_pixels
        [~,argmin] = min(dist(:,i));
        pixel_k(i) = argmin;
    end

    if pixel_k == old_pixel
        disp("Covergence at L = " + l)
        break
    else
        old_pixel = pixel_k;
    end

    % Recompute each cluster center by taking the mean of all pixels assigned to it
    for k = 1:K
        rgb = zeros(1,n_colors);
        t = 0;
        for i = 1:n_pixels
            if pixel_k(i) == k
                t = t + 1; 
                rgb = rgb + image_vec(i,:);
            end
        end
        
        if t > 0 % if t is still 0 then no pixels picked this cluster
            centers(k,:) = rgb/t;
        end
    end

    % Recompute all distances between pixels and cluster centers
    for i = 1:n_pixels
        for k = 1:K
            dist(k,i) = pdist2(image_vec(i,:),centers(k,:));
        end
    end
end

segmentation = reshape(pixel_k,[width,height]);

end