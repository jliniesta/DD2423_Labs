function [ segmentation, centers ] = kmeans_segm(image, K, L, seed);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
rng(seed);

n_pixels = size(image, 1);



image_vec = image;


% Randomly initialize the K cluster centers

centers = zeros(K,3);
dist = zeros(K,n_pixels);
pixel_k = zeros(n_pixels,1);
for k = 1:K
    %centers(k,:) = randi(255,[1,3]);
    centers_coord = randi(n_pixels,1,1);
    centers(k,:) = image_vec(centers_coord, :);
end
% Compute all distances between pixels and cluster centers
for i = 1:n_pixels
    
    for k = 1:K
        dist(k,i) = pdist2(image_vec(i,:),centers(k,:));

    end
end
% Iterate L times
old_pixel = zeros(n_pixels,1);
for l = 1:L
    % Assign each pixel to the cluster center for which the distance is minimum
    for i = 1:n_pixels
        [~,argmin] = min(dist(:,i));
        pixel_k(i) = argmin;
        
    end
    if pixel_k == old_pixel
        disp("Covergence at L = "+l)
        break
    else
        old_pixel = pixel_k;
    end

    % Recompute each cluster center by taking the mean of all pixels assigned to it
    for k = 1:K
        rgb = zeros(1,3);
        t = 0;
        for i = 1:n_pixels
            if pixel_k(i) == k
                t = t+ 1; 
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

segmentation = pixel_k;