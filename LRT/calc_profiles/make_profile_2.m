load('profile_128_64_5T2prep_50undersampled_60TSE.mat');


profile_reshaped = reshape(profile_order, 2, size(profile_order, 2) * size(profile_order, 3))';

nsa = zeros(size(profile_reshaped, 1), 1);
for k = 1 : size(profile_order, 3)
        profile_this_dynamic = profile_order(:,:,k)';
    for p = 1 : size(profile_order, 2)
        
        ky_match = (profile_this_dynamic(1:p-1, 1) == profile_this_dynamic(p, 1));
        kz_match = (profile_this_dynamic(1:p-1, 2) == profile_this_dynamic(p, 2));
        
        if(~isempty(ky_match)&&~isempty(kz_match))
            match_profiles_idx = find( ky_match .* kz_match );
            if(~isempty(match_profiles_idx))
                nsa(p + (k-1) * size(profile_order, 2)) = nsa(match_profiles_idx(end) + (k-1) * size(profile_order, 2)) + 1;
            end
        end
        
    end
end

kspa = zeros(64, 128);
for m = 1:9840
    x_idx = profile_this_dynamic(m, 1) + 32;
    y_idx = profile_this_dynamic(m, 2) + 64;
    kspa(x_idx, y_idx) = kspa(x_idx, y_idx) + 1;
end
figure; imshow(kspa,[]); colormap jet



dlmwrite('profile_128_64_5T2prep_50undersampled_60TSE.dat',profile_final,'delimiter','\t');


