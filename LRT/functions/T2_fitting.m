% INPUT
%
% T2PREP_data: 4D image data in size of [x y z T2PREP]
% T2PREP:      T2PREP valus that matches T2PREP_data
%
% OUTPUT
%
%(c) q.zhang . 2017 . AMC

function [T2_mono_all_allslice, T2_mono_all_rgb_allslice, rsquare_mono_all_allslice] = T2_fitting(T2PREP_data, T2PREP,varargin)

%% options
drawROI=0; %1 means manually draw ROI, otherwize ROI determined by cut_off
auto_loop_all_slice = 1; % 1 means loop over all slices. 0 means manually choose slice
slice_id = 0;
if(nargin==3)
    intensity_cut_off = varargin{1};
else
    intensity_cut_off = 100;
end
rsquare_cutoff = 0.5;
T2_lower_bound=0;
T2_upper_bound=800;
%%
datasource = double(abs(T2PREP_data));

[xdim, ydim, slices,b_num] = size(datasource);
button = 1;


xls_line_position=1;
VOI_pixel = 0;
VOI_selected_pixel = 0;
while (button ~= 27) %27 corresponds to key 'Esc'
    %%%%%%%%%%%%%%%%%%%%%%%%%%% ROI Bypass by Abs %%%%%%%%%%%%%%%%%%%
    Bypass = 27;
    while (Bypass == 27)
        if (~auto_loop_all_slice)
            slice_id = input(['Please input the slice index. 1~',num2str(slices),'\n']);
        else
            slice_id = slice_id + 1;
        end
        
        disp(['Slice:', num2str(slice_id),'    /',num2str(slices)]);
        
        while((slice_id > slices)||(slice_id < 0))
            fprintf('Error! The maximum slice id is %d, input again.\n', slices);
            slice_id = input('Please input the slice index.\n');
        end
        
        log_sldata(:,:,slice_id,:) = log(datasource(:,:,slice_id, 1:b_num));
        
        
        img_b0 = datasource(:, :, slice_id, 1);   %b=400 for better lesion detection
        layer = datasource(:, :, slice_id, 1)+1;
        
        
        figure(41); h_im = imagesc(img_b0); axis off; axis square; colormap gray
        title(['Image of slice ' num2str(slice_id) ', b=' num2str(T2PREP(1)) ]);
        
        if (~auto_loop_all_slice)
            [x2, y2, Bypass] = ginput(1);
        else
            Bypass = 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    button = 1;
    
    while ((button ~= 27)&&(button ~= 3))   % 3 means right mouse button click
        
        if drawROI==1
            flag = 27;
            while (flag == 27)

                fprintf('Please select  ROI:\n');
                button = 1;
                figure(42); h_im = imagesc(img_b0); axis off; axis square; colormap gray
                title(['Image of slice ' num2str(slice_id) ', b=' num2str(T2PREP(1))]);
                e = imfreehand(gca);
                [x4, y4, flag] = ginput(1);
            end
            
            position = wait(e);
            ROI = createMask(e,h_im);
        else
            ROI = img_b0>intensity_cut_off;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        %% T2 mapping (T2 all T2 high withb0 T2 high nob0 T2 few)
        
        % %%%%%%%%%%%%%%%% T2 mapping ALL%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cfun1 = fittype('poly1');
        rsquare_mono_all = zeros(xdim, ydim);
        T2_mono_all = zeros(xdim, ydim);
        sse_mono_all = zeros(xdim, ydim);
        fprintf('running mono-exponential T2_all mapping\n');
        valid_pixel_for_average=zeros(xdim, ydim);
        for(i1=1:xdim)
            if(mod(i1, 20)==0)
                fprintf('running line [%d]\n', i1);
            end
            for (j1=1:ydim)
                if(ROI(i1,j1)==1)
                    if(all(isfinite(log_sldata(i1, j1, slice_id, :))))
                        if (mean(datasource(i1, j1, slice_id, 1))>intensity_cut_off)
                            valid_pixel_for_average(i1,j1)=1;
                            [c_100,gof_100] = fit(col(T2PREP), col(squeeze(log_sldata(i1, j1, slice_id, :))), cfun1);
                            rsquare_mono_all(i1,j1) = gof_100.rsquare;
                            sse_mono_all(i1,j1) = gof_100.sse;
                            
                            if(gof_100.rsquare>=rsquare_cutoff)
                                values_pixel = coeffvalues(c_100);
                                T2_mono_all(i1,j1)= -1./values_pixel(1);
                                %                 m_rmse(i1, j1) = 1-(-values_pixel(1)/gof_100.rmse);
                                if ((T2_mono_all(i1,j1)<0)||(T2_mono_all(i1,j1)>1000))  %1000 cratiria is temperary
                                    T2_mono_all(i1,j1)=0;
                                    rsquare_mono_all(i1,j1)=0;
                                    ROI(i1,j1)=0;
                                    %                     m_rmse(i1,j1)=0;
                                end
                            else
                                T2_mono_all(i1,j1)= 0;
                                ROI(i1,j1)=0;
                            end
                        end
                    else
                        T2_mono_all(i1,j1)=0;
                        ROI(i1,j1)=0;
                    end
                end
            end
        end
        
        ROI(find(T2_mono_all(:)==0))=0;
        
        %-----------------------------------------------------------
        T2_mono_all_rgb = ind2rgb(round(layer), gray(max(layer(:))));
        %     rescale = 256/max(T2_mono_all(:));
        rescale = 256/T2_upper_bound;
        impose = ind2rgb(round(T2_mono_all*rescale), hot(256));
        for(i1=1:xdim)
            for (i2=1:ydim)
                if(ROI(i1, i2)==1)
                    T2_mono_all_rgb(i1, i2, 1) = impose(i1, i2, 1);
                    T2_mono_all_rgb(i1, i2, 2) = impose(i1, i2, 2);
                    T2_mono_all_rgb(i1, i2, 3) = impose(i1, i2, 3);
                end
            end
        end
        figure(43);
        imagesc(T2_mono_all); colormap hot;  axis equal;   axis off;  title 'T2 all'; colorbar;
        drawnow;
        figure(44);
        subplot(1,3,1); imagesc(T2_mono_all_rgb); axis equal; axis off; title 'T2 all';
        subplot(1,3,2); imagesc(T2_mono_all, [T2_lower_bound T2_upper_bound]); colormap hot;  axis equal;   axis off;  title 'T2 all'; colorbar;
        subplot(1,3,3); imagesc(rsquare_mono_all);  colormap hot; axis equal; axis off; title 'T2 all R^2 map'; colorbar;
        drawnow;
        
        figure(45);
        imagesc(T2_mono_all_rgb); axis equal; axis off; title(['T2 slice:',num2str(slice_id)])
        T2_mono_all_rgb_allslice(:,:,:,slice_id)=T2_mono_all_rgb;
        drawnow;
        
        figure(46);
        imagesc(T2_mono_all,[T2_lower_bound T2_upper_bound]); axis equal; axis off; title(['T2 slice:',num2str(slice_id)]); colorbar
        T2_mono_all_allslice(:,:,slice_id)=T2_mono_all;
        drawnow;
   
        figure(47);
        imagesc(rsquare_mono_all, [rsquare_cutoff 1]); axis equal; axis off; title(['R^2 slice: ',num2str(slice_id)]); colorbar
        rsquare_mono_all_allslice(:,:,slice_id)=rsquare_mono_all;
        drawnow;

        %%%%%%%%%%%%%%% Mono-exponential T2 all mapping end %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        figure(48);
        pd=histfit_qw(T2_mono_all(T2_mono_all>0),50); title(['fitted mean: ',  num2str(pd.mu), ';  fitted std: ', num2str(pd.sigma)]);
        display(['Raw data: Mean: ',num2str(mean(T2_mono_all(T2_mono_all>0))),'  Median: ',num2str(median(T2_mono_all(T2_mono_all>0))),'  Std: ',num2str(std(T2_mono_all(T2_mono_all>0)))]);
        display(['Raw data: Mean-Std=', num2str(mean(T2_mono_all(T2_mono_all>0))-std(T2_mono_all(T2_mono_all>0))), ';  Mean+Std=', num2str(mean(T2_mono_all(T2_mono_all>0))+std(T2_mono_all(T2_mono_all>0))) ]);
        display(['Fitted data: Mean: ', num2str(pd.mu),'  Std: ', num2str(pd.sigma)]);
        display(['Fitted data: Mean-Std=', num2str(pd.mu-pd.sigma), ';  Mean+Std=', num2str(pd.mu+pd.sigma) ]);
        display(['Fitted data: Mean-3*Std=', num2str(pd.mu-3*pd.sigma), ';  Mean+3*Std=', num2str(pd.mu+3*pd.sigma) ]);
        %KERRY: histfit_qw is a revision of histfit. Add the function of export
        %fit resutls (min and sigma of normal distribution)
        
        
        
        %     T2_lower_bound=input('Input the lower bound of tumor core: ');
        %     T2_upper_bound=input('Input the upper bound of tumor core: ');
        
        T2_mono_selected=T2_mono_all;
        T2_mono_selected((T2_mono_selected<T2_lower_bound))=0;
        T2_mono_selected((T2_mono_selected>T2_upper_bound))=0;
        rsquare_mono_selected=rsquare_mono_all;
        rsquare_mono_selected((T2_mono_selected<T2_lower_bound))=0;
        rsquare_mono_selected((T2_mono_selected>T2_upper_bound))=0;
        
        %-----------------------------------------------------------
        T2_mono_selected_rgb = ind2rgb(round(layer), gray(max(layer(:))));
        rescale = 256/max(T2_mono_selected(:));
        impose = ind2rgb(round(T2_mono_selected*rescale), hot(256));
        for(i1=1:xdim)
            for (i2=1:ydim)
                if(ROI(i1, i2)==1)
                    T2_mono_selected_rgb(i1, i2, 1) = impose(i1, i2, 1);
                    T2_mono_selected_rgb(i1, i2, 2) = impose(i1, i2, 2);
                    T2_mono_selected_rgb(i1, i2, 3) = impose(i1, i2, 3);
                end
            end
        end
        figure(49);
        imagesc(T2_mono_selected, [T2_lower_bound, T2_upper_bound]); colormap hot;  axis equal;   axis off;  title 'T2 selected'; colorbar;
        
        figure(50);
        subplot(1,3,1); imagesc(T2_mono_selected_rgb); axis equal; axis off; title 'T2 selected';
        subplot(1,3,2); imagesc(T2_mono_selected, [T2_lower_bound, T2_upper_bound]); colormap hot;  axis equal;   axis off;  title 'T2 selected'; colorbar;
        subplot(1,3,3); imagesc(rsquare_mono_selected);  colormap hot; axis equal; axis off; title 'T2 selected R^2 map'; colorbar;
        
        
        if (~auto_loop_all_slice)
            [x1, y1, button] = ginput(1);
        else
            if(slice_id == slices)
                button = 27;
            else
                button = 3;
            end                    
        end
        drawnow;
        
        
    end
    
end     % for while loop of button detection

end