function sens = gen_sense_map(raw_data_fn, OutputSizeSensitivity)
% GENERATE SENSE MAP
% 
% Input:    raw_data_fn: the targeted image filename
%           OutputSizeSensitivity: recon image size [x y z]
% 
% 
% Output:   sens: sense map in the size of [x y z channel]

    sense_ref = 'dp_17052017_1403271_1000_11_wipsenserefscanexperiment1clearV4.raw';
    MR2 = raw_data_fn;
    coil_survey = 'dp_17052017_1401228_1000_8_wipcoilsurveyscanV4.raw';
    MR_sense = MRsense(sense_ref, MR2, coil_survey);

    MR_sense.Smooth = 1;
    MR_sense.Mask = 1;
    MR_sense.Rotate = 0;
    MR_sense.Extrapolate = 1;
    MR_sense.MatchTargetSize = 1;
    % MR_sense.RemoveMOversampling = 1;
    MR_sense.OutputSizeReformated = OutputSizeSensitivity;
    MR_sense.OutputSizeSensitivity = OutputSizeSensitivity;
    disp('MRSense Perform...');
    MR_sense.Perform;
%     sens = conj(flipdim(MR_sense.Sensitivity, 1));
%     sens = conj((MR_sense.Sensitivity));
    sens = ((MR_sense.Sensitivity));
    sens = double(sens); %./max(sens(:))
end