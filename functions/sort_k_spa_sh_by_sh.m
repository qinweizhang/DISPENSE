%Sort K space data by kx, ky, kz, ch, shot

function kspa_sortted = sort_k_spa_sh_by_sh(ima_k_spa_data, shot_range, TSE, pars)

kx_dim = TSE.kx_dim;
ky_dim = TSE.ky_dim;  %max_ky * 2 + 1;
kz_dim = TSE.kz_dim;
ch_dim = length(pars.enabled_ch); %TSE.ch_dim;

%% Sortting; should be quick
%step1: sort all dims except sh_dim
kspa_temp = zeros(size(ima_k_spa_data,1), ky_dim, kz_dim, ch_dim);
shot_all_temp = zeros(1, ky_dim, kz_dim);

used_profile_nr = 0;
total_profile = size(ima_k_spa_data, 2);

h=waitbar(0, '[1/2] Sortting all dims except shot dim + correct checkerboard pattern');

for prof_idx = 1:total_profile
  
    if(mod(prof_idx,1000)==0)
        waitbar(prof_idx./total_profile);
    end
    
    ch_nr = mod(prof_idx-1, TSE.ch_dim)+1;
    ch_idx = find(pars.enabled_ch==ch_nr);
    sh_nr = TSE.shot_matched(prof_idx);
    sh_idx = find(shot_range==sh_nr);
    if(~isempty(ch_idx)&&~isempty(sh_idx))
        ky_idx = TSE.ky_matched(prof_idx) + floor(ky_dim/2)+1;
        kz_idx = TSE.kz_matched(prof_idx) + floor(kz_dim/2)+1;
        used_profile_nr = used_profile_nr+1;
        checker_board_factor = double((-1).^(ky_idx+kz_idx+1));  %-1 top left corner
        kspa_temp(:,ky_idx,kz_idx,ch_idx) =  ima_k_spa_data(:,prof_idx) .* checker_board_factor;
        shot_all_temp(1,ky_idx,kz_idx,ch_idx) = sh_idx;
    end
end
close(h)

temp = kspa_temp(round(kx_dim/2),:,:,:);
assert(used_profile_nr == sum(abs(temp(:))>0)); clear temp
size(kspa_temp)

%% step2: srot shot dim
sh_dim = length(shot_range);
kspa_sortted = zeros(kx_dim, ky_dim, kz_dim, ch_dim, sh_dim);
padding_size = kx_dim-size(kspa_temp,1);
padding_left = ceil(padding_size/2);

kx_sort_range = padding_left+1:padding_left+size(kspa_temp,1);

h=waitbar(0, '[2/2] Sortting shot dim  + kx size matching');
for sh=1:sh_dim
    waitbar(sh/sh_dim);
    kspa_sortted(kx_sort_range,:,:,:,sh) = bsxfun(@times, kspa_temp, shot_all_temp==sh);
end
close(h)

temp = kspa_sortted(round(kx_dim/2),:,:,:,:);
assert(used_profile_nr == sum(abs(temp(:))>0));
assert(kx_dim==size(kspa_sortted,1));
clear temp  kspa_b0_temp shot_b0_all_temp


end