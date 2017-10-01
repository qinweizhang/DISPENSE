 % ===== Set the parameters ======
 function Mxy_at_echo_time = CUBE_Mxy_calculation_lite(T1, T2, TEes, ETL)
    

%==================== Get the flipangles for CUBE

    angles=CUBEangles(1000,50,TEes,0.4,ETL); %fix FA for all
    angles=angles*180/pi;


%====================OVERRULL the FA here. 
    
%     angles = 180 * ones(1,ETL);


%     angles = [134.5000    79.0000   58.1300   51.2500   51.2500   55.2900   62.1400   72.1400   86.7000  110.2800  120.0000  122.1100  124.2100  126.3200  128.4200  130.5300  132.6300  134.7400  136.8400  138.9500  141.0500  143.1600  145.2600  147.3700  149.4700  151.5800  153.6800  155.7900  157.8900  160.0000];
    % angg=[160 140 145 150 155 160:-4:20 20:4:160]
    % angg=[160 150 130 110 130 145 150 155 160:-4:30 28 27 25 24 22 21 20 21 22 22 23 24 25 26 28 30:4:120 126:6:160]
    % angg = [160 147 156 164 168 173 177 180*ones(1,121)];
    % plot(1:length(angg),angg)

    
    
    assert(length(angles)==ETL); %want same number of bvals as T1vals (= number of phantoms)

    % ===== Simulate the Decay ======
    
    [Mxy_at_echo_time,m_long,FZall] = epg_forward(angles,ETL,TEes*ones(1,ETL),T1,T2,1); 
    
    figure(400)
    hold on;
    plot(TEes.*[1:ETL], Mxy_at_echo_time);
    xlabel('Echo (ms)'); ylabel('Mxy');
    
    figure(401);
    hold on;
    plot(angles);
    xlabel('Echo (#)'); ylabel('degree');

    

 end

    
    