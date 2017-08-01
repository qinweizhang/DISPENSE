 function sens_map = get_sense_map_ecalib(fn, recon_dim )
       
    MR_TSEDPima_data = MRecon(fn);

    MR_TSEDPima_data_recon1 = MR_TSEDPima_data.Copy;

    MR_TSEDPima_data_recon1.Parameter.Parameter2Read.typ = 1;
    MR_TSEDPima_data_recon1.Parameter.Parameter2Read.mix = 0;  %for DPnav Spirals
    MR_TSEDPima_data_recon1.ReadData;
    MR_TSEDPima_data_recon1.RandomPhaseCorrection;
    MR_TSEDPima_data_recon1.RemoveOversampling;
    MR_TSEDPima_data_recon1.PDACorrection;
    MR_TSEDPima_data_recon1.DcOffsetCorrection;
    MR_TSEDPima_data_recon1.MeasPhaseCorrection;
    MR_TSEDPima_data_recon1.SortData;
    kspa_sorted = double(squeeze(MR_TSEDPima_data_recon1.Data));

    bart_command_1 = sprintf('resize -c 0 %d 1 %d 2 %d',recon_dim(1),recon_dim(2),recon_dim(3) )
    kspa_sorted = kspa_sorted(:,:,:,:,1);
    kspa_TSE_resize = bart(bart_command_1, kspa_sorted);
    che=create_checkerboard([1,size(kspa_TSE_resize,2),size(kspa_TSE_resize,3)]);
    kspa_TSE_resize=bsxfun(@times,kspa_TSE_resize,che);
    sens=bart('ecalib -m1 -r 20 -c0.1',kspa_TSE_resize);
    sens_map =  sens;
    figure(701); montage(abs(sens_map(:,:,11,:)),'displayrange',[])
        
 end