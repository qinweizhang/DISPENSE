% Function of calculating temporal interpolator for time segmented NUFFT
% 
% INPUT
% 
% tSeg:                 Time-segment structure contains nr of segments and delta_t
% samples_idx:          Sample index. sampling time = sample_idx * delta_t
% b0_maps_per_ms:       Offresonance map with correct size in mHz
% 
% OUT
% 
% tSeg_nuuft_interpolator:    temporal interpolator in size of [samplled_points time_segments]. Intuition: The contribution of segments l to every sampled point S(t).
%                             S(l,t) = interpolator * NUFFT * omega(l) * x      and    S(t) = sum(S(l,t)) over l
% 
% (C) q.zhang @AMC 2017 ref:  "Fast, Iterative Imaging Reconstruction for MRI in the Presence of Field Inhomogeneities" Bradley Sutton et al. 

function  [tSeg_nuuft_interpolator, tSeg_updated] = tSeg_nuuft_interpolator_init(tSeg, samples_idx,  b0_maps_per_ms)

%% preparations

samples_t = samples_idx * tSeg.aq_interval;  %in ms
tau = range(samples_t) / tSeg.nr_segments;   %intervals for every segment; in ms

tSeg_updated = tSeg;
tSeg_updated.tau = tau;
%% calculations

h=waitbar(0, 'Calculating temporal interpolator for time segmented NUFFT');

clear GG Gb_t
for l = 0: tSeg.nr_segments
    phi_t = bsxfun(@times, b0_maps_per_ms(:),  (samples_t - tau * l));
    Gb_t(l+1,:) =  mean(exp(-1i * phi_t),1);

    for l_p = 0: tSeg.nr_segments     
        temp =  exp(-1i * b0_maps_per_ms * tau * (l_p - l));
        GG(l+1, l_p+1)= mean(temp(:));
    end
    waitbar(l./tSeg.nr_segments);
end
close(h)

tSeg_nuuft_interpolator_raw = GG \ Gb_t; % same as: inv(GG) * Gb_t;

% normalize interpolator
tSeg_nuuft_interpolator_mag = colnorm2(tSeg_nuuft_interpolator_raw); % 2nd norm along first dim
tSeg_nuuft_interpolator = bsxfun(@rdivide, tSeg_nuuft_interpolator_raw,  tSeg_nuuft_interpolator_mag);
tSeg_nuuft_interpolator(find(isnan(tSeg_nuuft_interpolator)+isinf(tSeg_nuuft_interpolator))) = 0;


figure(31); plot((abs(tSeg_nuuft_interpolator))'); title('abs(tSeg nuuft interpolator)'); 


end





