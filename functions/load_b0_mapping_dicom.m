function data = load_b0_mapping_dicom(path)

current_path = pwd;
cd(path)

filelist=dir;
% fid=fopen('dicominfo.txt','w+');
%fprintf(fid,'%s\t\t\t%s\t\t\t%s\t\t\t\n', 'Protocol Name', 'Rescale Slop', 'Rescale Slop Intercept');
for i=3:length(filelist)
    filename=filelist(i).name;
    x=dicominfo(filename);
    %     tmp=x.ProtocolName;
    tmp=x.SeriesDescription;
    instanceNo=x.InstanceNumber;             % Slice No.
    %     RS=1;                    % Rescales Slop
    %     RSI=0;                %Rescale intercept
    %     SS=1;
    RS=x.RescaleSlope;                    % Rescales Slop
    RSI=x.RescaleIntercept;                %Rescale intercept
    SS=x.ScaleSlope;                      %philips private parameter tag: (2005,100e)
    %       scannumber=x.AcquisitionNumber;
    %     fprintf(fid,'%s\t\t\t%2.3f\t\t\t%2.2f\t\t\t\n', tmp, RS, RSI);
    matname= regexprep(tmp,'[\s'']','_'); % Eliminate white space - not tested
    matname= regexprep(matname,'-','NEGATIVE'); % Eliminate white space - not tested
    matname= regexprep(matname,'\.','dot'); % Eliminate white space - not tested
    %matname= strcat(matname,num2str(sliceNo));
    matname= strcat(matname,num2str(i-2));
    %         dicomdata=(double(dicomread(x)).*RS+RSI)
    %         dicomdata=(double(dicomread(x)).*RS+RSI)/SS;
    dicomdata=(double(dicomread(x)).*RS+RSI)/(RS*SS);
        %       dicomdata=double(dicomread(x));
    while exist(matname)
        matname=strcat(matname,'_D');
    end
    eval([matname,'=dicomdata;']);
    
    %     dyn_no = mod(instanceNo-1,2) + 1;
    %     slice_no = ceil(instanceNo/2);
    %     data(:,:,slice_no,dyn_no) = dicomdata;
    data(:,:,instanceNo) = dicomdata;

        
    if mod(i,10)==0;
        disp([num2str(i), ' of ',num2str(length(filelist)),'...'] );
    end
    
end
dim = size(data);
data = reshape(data, dim(1), dim(2), dim(3)/2, 2 );
data = double(data);
cd(current_path);


%%