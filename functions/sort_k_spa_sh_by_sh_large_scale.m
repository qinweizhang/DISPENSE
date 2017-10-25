%Optimized sort_k_spa_sh_by_sh for large scale recon called by DPsti_TSE_phase_error_cor_large_scale
%x is preselected. ima_k_spa_data in size of [1, profiles] instead of [kx, profiles]
%
%Sort K space data by kx, ky, kz, ch, shot

function kspa_sortted = sort_k_spa_sh_by_sh_large_scale(ima_k_spa_data, shot_range, TSE, pars)
ky_dim = TSE.ky_dim;  %max_ky * 2 + 1;
kz_dim = TSE.kz_dim;
ch_dim = length(pars.enabled_ch); %TSE.ch_dim;

%% Sortting; should be quick
%step1: sort all dims except sh_dim
kspa_temp = zeros(1, ky_dim, kz_dim, ch_dim);
shot_all_temp = zeros(1, ky_dim, kz_dim);

used_profile_nr = 0;
total_profile = size(ima_k_spa_data, 2);

for prof_idx = 1:total_profile
      
    ch_nr = mod(prof_idx-1, TSE.ch_dim)+1;
    ch_idx = find(pars.enabled_ch==ch_nr);
    sh_nr = TSE.shot_matched(prof_idx);
    sh_idx = find(shot_range==sh_nr);
    if(~isempty(ch_idx)&&~isempty(sh_idx))
        
        ky_idx = TSE.ky_matched(prof_idx) + floor(ky_dim/2)+1;  
        kz_idx = TSE.kz_matched(prof_idx) + floor(kz_dim/2)+1;  
        
        ky_idx_sense = TSE.ky_matched(prof_idx).*TSE.SENSE_ky + floor(ky_dim/2)+1;  %temporary solution only for SENSE = int
        kz_idx_sense = TSE.kz_matched(prof_idx).*TSE.SENSE_kz + floor(kz_dim/2)+1;  %temporary solution only for SENSE = int
        
        used_profile_nr = used_profile_nr+1;
        checker_board_factor = double((-1).^(ky_idx+kz_idx+1));  %-1 top left corner
        
        kspa_temp(:,ky_idx_sense,kz_idx_sense,ch_idx) =  ima_k_spa_data(:,prof_idx) .* checker_board_factor;
        shot_all_temp(1,ky_idx_sense,kz_idx_sense,ch_idx) = sh_idx;
    end
end

temp = kspa_temp(round(size(kspa_temp,1)/2),:,:,:);
assert(used_profile_nr == sum(abs(temp(:))>0)); clear temp
size(kspa_temp)

%% step2: srot shot dim
sh_dim = length(shot_range);
kspa_sortted = zeros(1, ky_dim, kz_dim, ch_dim, sh_dim);


for sh=1:sh_dim
    kspa_sortted(:,:,:,:,sh) = bsxfun(@times, kspa_temp, shot_all_temp==sh);
end

temp = kspa_sortted(1,:,:,:,:);
assert(used_profile_nr == sum(abs(temp(:))>0));
clear temp  kspa_b0_temp shot_b0_all_temp


end