%Sort K space data by kx, ky, kz, ch, shot

function kspa_sortted = sort_k_spa_sh_by_sh_3(ima_k_spa_data, shot_range, TSE, pars)

ky_dim = TSE.ky_dim;  %max_ky * 2 + 1;
kz_dim = TSE.kz_dim;
ch_dim = length(pars.enabled_ch); %TSE.ch_dim;
sh_dim = length(shot_range);

%% Sortting; all dim including "shot" in one go

kspa_temp = zeros(size(ima_k_spa_data,1), ky_dim, kz_dim, ch_dim, sh_dim);
shot_all_temp = zeros(1, ky_dim, kz_dim);

used_profile_nr = 0;
total_profile = size(ima_k_spa_data, 2);

h=waitbar(0, 'Sortting all dims (INCLUDING shot dim) + correct checkerboard pattern');

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
        
        ky_idx_sense = TSE.ky_matched(prof_idx).*TSE.SENSE_ky + floor(ky_dim/2)+1;  %temporary solution only for SENSE = int
        kz_idx_sense = TSE.kz_matched(prof_idx).*TSE.SENSE_kz + floor(kz_dim/2)+1;  %temporary solution only for SENSE = int
        
        used_profile_nr = used_profile_nr+1;
        checker_board_factor = double((-1).^(ky_idx+kz_idx+1));  %-1 top left corner
        
        kspa_temp(:,ky_idx_sense,kz_idx_sense,ch_idx, sh_idx) =  ima_k_spa_data(:,prof_idx) .* checker_board_factor;
    end
end
close(h)

temp = kspa_temp(round(size(kspa_temp,1)/2),:,:,:);
assert(used_profile_nr == sum(abs(temp(:))>0)); clear temp
size(kspa_temp)

kspa_sortted = kspa_temp;

clear temp  kspa_b0_temp shot_b0_all_temp


end