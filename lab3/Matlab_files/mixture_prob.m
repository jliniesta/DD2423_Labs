function prob = mixture_prob(image, K, L, mask)
% Let I be a set of pixels and V be a set of K Gaussian components in 3D (R,G,B).
width = size(image,1);
height = size(image,2);

n_pixels = width * height;
I_vec = im2double(reshape(image,n_pixels, 3));
M_vec = reshape(mask,n_pixels,1);

% Store all pixels for which mask=1 in a Nx3 matrix
I_vec_masked = I_vec(M_vec==1,:);

seed = 1;

% Randomly initialize the K components using masked pixels

[segmentation, centers] = kmeans_segm_noreshape(I_vec_masked, K, L, seed);
g = zeros(size(I_vec_masked,1),K);
w = zeros(1,K);

for k=1:K
    w(k) = sum(segmentation==k)/size(segmentation,1);
end

cov = cell(K,1);
cov(:) = {rand*eye(3)};



% Iterate L times
for i=1:L
    for k=1:K
    % Expectation: Compute probabilities P_ik using masked pixels
        mean = centers(k,:);
        cov_k = cov{k};
        diff = bsxfun(@minus, I_vec_masked, mean);
        g(:,k) = 1 / sqrt((2*pi)^3*det(cov{k})) * exp(-0.5*sum((diff/cov_k.*diff),2));
       
    end

    p = bsxfun(@times, g, w); % p is NxK*1xK = NxK
    p_sum = sum(p,2);
    p = bsxfun(@rdivide, p, p_sum);

    % Maximization: Update weights, means and covariances using masked pixels
    w = sum(p,1) / size(p,1);

    for k=1:K
        %w(k) = sum(p(:,k))/n_pixels_m; % p(:,k) probabilities of all pixels belonging to cluster k
        centers(k,:) = p(:,k)'*I_vec_masked / sum(p(:,k),1); % p(:,k) is N_pixels x 1, I_vec_masked is N_pixels x 3
        diff = bsxfun(@minus, I_vec_masked, centers(k,:));
        cov{k} = (diff'*bsxfun(@times, p(:,k), diff)) / sum(p(:,k),1);
    end
end
    
    
% Compute probabilities p(c_i) in Eq.(3) for all pixels I.
g = zeros(n_pixels,K);
for k = 1:K
    mean = centers(k,:);
    cov_k = cov{k};
    diff = bsxfun(@minus, I_vec, mean);
    g(:,k) = 1 / sqrt((2*pi)^3*det(cov_k)) * exp(-0.5*sum((diff/cov_k.*diff),2));
end
prob = sum(bsxfun(@times, g, w),2);
prob = reshape(prob,width,height,1);

end