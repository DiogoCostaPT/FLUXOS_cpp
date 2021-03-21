
% Posprocessing the results (calculating the connected area to the outlet
% over time)

clear all
close all

PROCESS_data_3event_databased_flag = 0;
PLOT_data_3event_databased_flag = 0;

PROCESS_data_SCENARIOS_SENS = 0;
PLOT_data_SCENARIOS_SENS = 1;

% PROCESSING THE DATA
if PROCESS_data_3event_databased_flag
    
    delete('batch_1_select_paper.mat');
    
    bathdir = '/media/dcosta/DATADRIVE1/fluxos_tests/SIMULATIONS_sync/batch_1_select_paper/';
    outlet_xy = [658,892];

    sims_folders_info = dir(bathdir);
    sims_folders_nam = {sims_folders_info.name};
    sims_folders_nam = sims_folders_nam(3:end);

    simtsteps = numel(sims_folders_nam);

    Qmelt_allsim = {};
    connect_area_outlet_timser_allsim = {}; 
    totalwet_area_timser_allsim = {}; 
    activecellhours_timser_allim = {};
    timeOutletConnect_timeser_allsim = {};
    connectbio_timser_allsim = {};
    disp('Processing the simulations...')
    for s = 1:simtsteps

        % Snowmelt load
        try
            modset = importdata([bathdir,sims_folders_nam{s},'/modset.fluxos']);
            Qmeltfile = modset.colheaders;
            Qmelt = importdata([bathdir,sims_folders_nam{s},'/',Qmeltfile{1}]);
        catch
            continue;
        end
        
        % Model Results
        sim_dir = [bathdir,sims_folders_nam{s},'/Results/'];
        steps_i = dir(sim_dir);
        steps_i_nam = {steps_i.name};
        disp(['  >',sims_folders_nam{s}]);
        steps_i_nam_num = sort(gettimesteps(steps_i_nam));
        wetarea_timser = steps_i_nam_num * NaN;
        connectedarea_timser = steps_i_nam_num * NaN;
        activecellhours_timser = steps_i_nam_num * NaN;
        timeOutletConnect_timeser = steps_i_nam_num * NaN;
        connectbio_timser = steps_i_nam_num * NaN;

        hbar = parfor_progressbar(numel(steps_i_nam_num), 'Calculating...');
        %h = waitbar(0,'Calculating...');
        
        for t = 1:numel(steps_i_nam_num)    
        %%for t = 1:numel(steps_i_nam_num)
        
            % Snowmelt forcing
            
            
            % Results data
            dataraw = importdata([sim_dir,num2str(steps_i_nam_num(t)),'.txt']);
                        
            % total wet area
            h_res_2D = get2Dmatrixresults(dataraw,steps_i_nam_num(t),sims_folders_nam{s},4);
            wetcells = binnary4desiredcells(h_res_2D);              
            wetarea_timser(t) = nansum(nansum(wetcells)) * 9 / 1000000; % km2
            
            % Total connected area to the outlet
            connect_cells_outlet = processdata(h_res_2D,outlet_xy);
            connectedarea_timser(t) = nansum(nansum(connect_cells_outlet)) * 9 / 1000000; % km2
            
            % Total active cell hours
            timeactive_res_2D = get2Dmatrixresults(dataraw,steps_i_nam_num(t),sims_folders_nam{s},15);           
            activecellhours_timser(t) = nansum(nansum(timeactive_res_2D)); % km2
            
            % Total connection hours to the outlet
            loc = connect_cells_outlet;
            loc(isnan(loc))=0;
            if sum(sum(loc))~=0
                timeOutletConnect_2D = timeactive_res_2D(find(loc==1));           
                timeOutletConnect_timeser(t) = nansum(nansum(timeOutletConnect_2D)); % km2
            else
                timeOutletConnect_timeser(t) = 0;
            end
            
            % biogeochemical connectivity
            concsoil_res_2D = get2Dmatrixresults(dataraw,steps_i_nam_num(t),sims_folders_nam{s},12);
            concsoil_res_2D(concsoil_res_2D==0) = NaN;
            connectbio_timser(t) = 1-nanmean(nanmean(concsoil_res_2D));
            
            %figure
            %surf(wetcells)
            %view(90,90)
            %axis tight

            %waitbar(t / numel(steps_i_nam_num))
            hbar.iterate(1)
        end
        close(hbar)
        
        % Saving total values in array
        Qmelt_allsim.(genvarname(sims_folders_nam{s})) = Qmelt;
        totalwet_area_timser_allsim.(genvarname(sims_folders_nam{s})) = wetarea_timser';
        connect_area_outlet_timser_allsim.(genvarname(sims_folders_nam{s})) = connectedarea_timser';
        activecellhours_timser_allim.(genvarname(sims_folders_nam{s})) = activecellhours_timser';
        timeOutletConnect_timeser_allsim.(genvarname(sims_folders_nam{s})) = timeOutletConnect_timeser';
        connectbio_timser_allsim.(genvarname(sims_folders_nam{s})) = connectbio_timser';
        
    end

    save('batch_1_select_paper','Qmelt_allsim')
    save('batch_1_select_paper','totalwet_area_timser_allsim','-append')
    save('batch_1_select_paper','connect_area_outlet_timser_allsim','-append')
    save('batch_1_select_paper','activecellhours_timser_allim','-append')
    save('batch_1_select_paper','timeOutletConnect_timeser_allsim','-append')
    save('batch_1_select_paper','connectbio_timser_allsim','-append')

end

if PROCESS_data_SCENARIOS_SENS
    
    delete('batch_CRHM_scenarios.mat');
    
    bathdir = '/media/dcosta/DATADRIVE1/fluxos_tests/SIMULATIONS_sync/batch_CRHM_scenarios/';
    outlet_xy = [658,892];

    sims_folders_info = dir(bathdir);
    sims_folders_nam = {sims_folders_info.name};
    sims_folders_nam = sims_folders_nam(3:end);

    simtsteps = numel(sims_folders_nam);

    Qmelt_allsim = {};
    connect_area_outlet_timser_allsim = {}; 
    totalwet_area_timser_allsim = {}; 
    activecellhours_timser_allim = {};
    timeOutletConnect_timeser_allsim = {};
    connectbio_timser_allsim = {};
    ic_allsim = {};
    disp('Processing the simulations...')
    for s = 1:simtsteps

        % Snowmelt load
        try
            modset = importdata([bathdir,sims_folders_nam{s},'/modset.fluxos']);
            IC_val = modset.data(2);
            Qmeltfile = modset.colheaders;
            Qmelt = importdata([bathdir,sims_folders_nam{s},'/',Qmeltfile{1}]);
        catch
            continue;
        end
        
        % Model Results
        sim_dir = [bathdir,sims_folders_nam{s},'/Results/'];
        steps_i = dir(sim_dir);
        steps_i_nam = {steps_i.name};
        disp(['  >',sims_folders_nam{s}]);
        steps_i_nam_num = sort(gettimesteps(steps_i_nam));
        wetarea_timser = steps_i_nam_num * NaN;
        connectedarea_timser = steps_i_nam_num * NaN;
        activecellhours_timser = steps_i_nam_num * NaN;
        timeOutletConnect_timeser = steps_i_nam_num * NaN;
        connectbio_timser = steps_i_nam_num * NaN;

        hbar = parfor_progressbar(numel(steps_i_nam_num), 'Calculating...');
        %h = waitbar(0,'Calculating...');
        
        for t = 1:numel(steps_i_nam_num)    
        %%for t = 1:numel(steps_i_nam_num)
        
            % Snowmelt forcing
            
            
            % Results data
            dataraw = importdata([sim_dir,num2str(steps_i_nam_num(t)),'.txt']);
                        
            % total wet area
            h_res_2D = get2Dmatrixresults(dataraw,steps_i_nam_num(t),sims_folders_nam{s},4);
            wetcells = binnary4desiredcells(h_res_2D);              
            wetarea_timser(t) = nansum(nansum(wetcells)) * 9 / 1000000; % km2
            
            % Total connected area to the outlet
            connect_cells_outlet = processdata(h_res_2D,outlet_xy);
            connectedarea_timser(t) = nansum(nansum(connect_cells_outlet)) * 9 / 1000000; % km2
            
            % Total active cell hours
            timeactive_res_2D = get2Dmatrixresults(dataraw,steps_i_nam_num(t),sims_folders_nam{s},15);           
            activecellhours_timser(t) = nansum(nansum(timeactive_res_2D)); % km2
            
            % Total connection hours to the outlet
            loc = connect_cells_outlet;
            loc(isnan(loc))=0;
            if sum(sum(loc))~=0
                timeOutletConnect_2D = timeactive_res_2D(find(loc==1));           
                timeOutletConnect_timeser(t) = nansum(nansum(timeOutletConnect_2D)); % km2
            else
                timeOutletConnect_timeser(t) = 0;
            end
            
            % biogeochemical connectivity
            concsoil_res_2D = get2Dmatrixresults(dataraw,steps_i_nam_num(t),sims_folders_nam{s},12);
            concsoil_res_2D(concsoil_res_2D==0) = NaN;
            connectbio_timser(t) = 1-nanmean(nanmean(concsoil_res_2D));
            
            %figure
            %surf(wetcells)
            %view(90,90)
            %axis tight

            %waitbar(t / numel(steps_i_nam_num))
            hbar.iterate(1)
        end
        close(hbar)
        
        % Saving total values in array
        Qmelt_allsim.(genvarname(sims_folders_nam{s})) = Qmelt;
        totalwet_area_timser_allsim.(genvarname(sims_folders_nam{s})) = wetarea_timser';
        connect_area_outlet_timser_allsim.(genvarname(sims_folders_nam{s})) = connectedarea_timser';
        activecellhours_timser_allim.(genvarname(sims_folders_nam{s})) = activecellhours_timser';
        timeOutletConnect_timeser_allsim.(genvarname(sims_folders_nam{s})) = timeOutletConnect_timeser';
        connectbio_timser_allsim.(genvarname(sims_folders_nam{s})) = connectbio_timser';
        ic_allsim.(genvarname(sims_folders_nam{s})) = IC_val;
        
    end

    save('batch_CRHM_scenarios','Qmelt_allsim')
    save('batch_CRHM_scenarios','totalwet_area_timser_allsim','-append')
    save('batch_CRHM_scenarios','connect_area_outlet_timser_allsim','-append')
    save('batch_CRHM_scenarios','activecellhours_timser_allim','-append')
    save('batch_CRHM_scenarios','timeOutletConnect_timeser_allsim','-append')
    save('batch_CRHM_scenarios','connectbio_timser_allsim','-append')
    save('batch_CRHM_scenarios','ic_allsim','-append')

end


% PLOTTING THE DATA
if PLOT_data_3event_databased_flag
   load('batch_1_select_paper.mat')
   
   structnam = fieldnames(totalwet_area_timser_allsim);
   structnam = {structnam{2},structnam{3},structnam{1}};
   legendstr = [];
   
   mean_area_connected = zeros(3,1);
   snowmelt_duraction = zeros(3,1);
   peak_snowmelt_rate = zeros(3,1);
   total_snowmelt = zeros(3,1);
   total_active_connect_outlet_m2hours = zeros(3,1);   
   
   
   stylelist = {'-','--'};
   colorlist = {[0 0 0],[0.7 0.7 0.7],[0.4 0.4 0.4]};
    % Hydrology
   fig1 = figure('position',[100 100 800 550])
   set(fig1,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
   for i=1:numel(structnam)
       subplot(3,2,1)
       yyaxis left
       Qmelt = Qmelt_allsim.(genvarname(structnam{i}))(:,2)/24; % mm/day to mm/h
       Qmelt_time = Qmelt_allsim.(genvarname(structnam{i}))(:,1)/3600; % sec to h
       h1(i) = plot(Qmelt_time,Qmelt,'color',colorlist{i},'linestyle',stylelist{1},'linewidth',1.5)
       hold on
       yyaxis right
       h11(i) = plot(Qmelt_time,cumsum(Qmelt),'color',colorlist{i},'linestyle',stylelist{2},'linewidth',1.5)
       hold on
       subplot(3,2,3)
       h2(i) = plot(totalwet_area_timser_allsim.(genvarname(structnam{i})),'color',colorlist{i},'linestyle',stylelist{1},'linewidth',1.5)
       hold on
       subplot(3,2,3)
       area_con_outlet = connect_area_outlet_timser_allsim.(genvarname(structnam{i}));
       h3(i) = plot(area_con_outlet,'color',colorlist{i},'linestyle',stylelist{2},'linewidth',1.5)
       hold on
       subplot(3,2,5)
       activem2hours_all = activecellhours_timser_allim.(genvarname(structnam{i}))*3*3; % cell -> m2 (area of cells = 3 x 3 m2)
       h4(i) = plot(activem2hours_all,'color',colorlist{i},'linestyle',stylelist{1},'linewidth',1.5)
       hold on
       subplot(3,2,5)
       activem2hours_onlycon2outlet = timeOutletConnect_timeser_allsim.(genvarname(structnam{i}))*3*3; % cell -> m2 (area of cells = 3 x 3 m2)
       h5(i) = plot(activem2hours_onlycon2outlet,'color',colorlist{i},'linestyle',stylelist{2},'linewidth',1.5)
       hold on
       
       mean_area_connected(i) = mean(area_con_outlet);
       snowmelt_duraction(i) = Qmelt_time(end);
       peak_snowmelt_rate(i) = max(Qmelt);
       total_snowmelt(i) = sum(Qmelt); % mm
       total_active_connect_outlet_m2hours(i) = sum(activem2hours_onlycon2outlet);
   
   end 
   
   subplot(3,2,1)
   grid on
   yyaxis left
   ylabel({'Snowmelt Input','[mm/h]'})
   yyaxis right
   ylabel({'Cumulative Snowmelt',' Input [mm]'})
   legend([h1(1),h1(3),h1(2),h11(1),h11(3),h11(2)],'2009','2011','2010','2009 (cumulative)','2011 (cumulative)','2010 (cumulative)')
   xlabel('Time [h]')
   axis tight
   subplot(3,2,3)
   grid on
   ylabel({'Active (flowing)','area [km^2]'})
   xlabel('Time [h]')
    axis tight
   subplot(3,2,5)
   grid on
    axis tight
   ylabel({'Cumulative active','m^2-hours [h]'})
   legend([h2(1),h3(1),h2(3),h3(3),h2(2),h3(2)],'2009 (total)','2009 (contributing area to outlet)',...
         '2011 (total)','2011 (contributing area to outlet)',...
            '2010 (total)','2010 (contributing area to outlet)')
   xlabel('Time [h]')

   subplot(3,2,4)
   yyaxis left
   plot(mean_area_connected,peak_snowmelt_rate,'ko')
   ylabel({'Peak snowmelt','rate [mm/h]'})
   grid on
   xlabel('Average contributing area to outlet [km^2]')
   yyaxis right
   plot(mean_area_connected,total_snowmelt,'k*')
   ylabel({'Total snowmelt','[mm]'})
   grid on
   legend('left axis','right axis')
   
   
   subplot(3,2,6)
   yyaxis left
   plot(total_active_connect_outlet_m2hours,snowmelt_duraction,'ko')
   ylabel({'Total snowmelt','duraction [h]'})
   grid on
   xlabel('Total active m^2-hours [h]')
   yyaxis right
   plot(total_active_connect_outlet_m2hours,total_snowmelt,'k*')
   ylabel({'Total snowmelt,','[mm]'})
   grid on
   legend('left axis','right axis')
   
   
     
   % Biogeochemistry
   final_connectbio_timser_allsim = zeros(3,1);
   fig2 = figure('position',[100 100 800 550])
   set(fig2,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    for i=1:numel(structnam)
       subplot(2,2,1)
       Qmelt = Qmelt_allsim.(genvarname(structnam{i}))(:,2)/24; % mm/day to mm/h
       Qmelt_time = Qmelt_allsim.(genvarname(structnam{i}))(:,1)/3600; % sec to h
       yyaxis left
       h1(i) = plot(Qmelt_time,Qmelt,'color',colorlist{i},'linestyle',stylelist{1},'linewidth',1.5)
       hold on
       yyaxis right
       h11(i) = plot(Qmelt_time,cumsum(Qmelt),'color',colorlist{i},'linestyle',stylelist{2},'linewidth',1.5)
       hold on
       grid on
       subplot(2,2,3)
       val = connectbio_timser_allsim.(genvarname(structnam{i}));
        h2 (i) = plot(val,'color',colorlist{i},'linestyle',stylelist{1},'linewidth',1.5)
        hold on
        grid on
        final_connectbio_timser_allsim(i) = val(end);
    end
    subplot(2,2,1)
    xlabel('Time [h]')
   yyaxis left
   ylabel({'Snowmelt Input','[mm/h]'})
   yyaxis right
   ylabel({'Cumulative Snowmelt',' Input [mm]'})
   legend([h1(1),h1(3),h1(2),h11(1),h11(3),h11(2)],'2009','2011','2010','2009 (cumulative)','2011 (cumulative)','2010 (cumulative)')
    subplot(2,2,3)
   legend([h2(1),h2(3),h2(2)],'2009','2011','2010')
   xlabel('Time [h]')
   ylabel({'Biogeochemical','connectivity [-]'})
       
      
       subplot(2,2,4)
       yyaxis left
       plot(final_connectbio_timser_allsim,snowmelt_duraction,'ko')
       ylabel({'Total snowmelt','duration [h]'})
       grid on
       xlabel('Total active m^2-hours [h]')
       yyaxis right
       plot(final_connectbio_timser_allsim,total_snowmelt,'k*')
       ylabel({'Total snowmelt,','[mm]'})
       grid on
       legend('left axis','right axis')
           
end


%%
if PLOT_data_SCENARIOS_SENS

    load('batch_CRHM_scenarios.mat')
   

   structnam = fieldnames(totalwet_area_timser_allsim);
   %structnam = {structnam{2},structnam{3},structnam{1}};
   legendstr = [];
   
   mean_area_connected = zeros(3,1);
   snowmelt_duraction = zeros(3,1);
   peak_snowmelt_rate = zeros(3,1);
   total_snowmelt = zeros(3,1);
   total_active_connect_outlet_m2hours = zeros(3,1);   
   
   
   stylelist = {'-','--'};
   colorlist = {[0 0 0],[0.7 0.7 0.7],[0.4 0.4 0.4]};
    % Hydrology
   fig1 = figure('position',[100 100 800 550]);
   set(fig1,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    for i=1:numel(structnam)
       subplot(3,2,1)
       %yyaxis left
       Qmelt = Qmelt_allsim.(genvarname(structnam{i}))(:,2)/24; % mm/day to mm/h
       Qmelt_time = Qmelt_allsim.(genvarname(structnam{i}))(:,1)/3600; % sec to h
       %h1(i) = plot(Qmelt_time,Qmelt,'color',colorlist{i},'linestyle',stylelist{1},'linewidth',1.5)
       %hold on
       %yyaxis right
       h11(i) = plot(Qmelt_time,cumsum(Qmelt),'color',colorlist{1},'linestyle',stylelist{1},'linewidth',1.5);
       hold on
       subplot(3,2,3)
       h2(i) = plot(totalwet_area_timser_allsim.(genvarname(structnam{i})),'color',colorlist{1},'linestyle',stylelist{1},'linewidth',1.5);
       hold on
       subplot(3,2,3)
       area_con_outlet = connect_area_outlet_timser_allsim.(genvarname(structnam{i}));
       h3(i) = plot(area_con_outlet,'color',colorlist{1},'linestyle',stylelist{1},'linewidth',1.5);
       hold on
       %subplot(3,2,5)
       %activem2hours_all = activecellhours_timser_allim.(genvarname(structnam{i}))*3*3; % cell -> m2 (area of cells = 3 x 3 m2)
       %h4(i) = plot(activem2hours_all,'color',colorlist{3},'linestyle',stylelist{1},'linewidth',1.5)
       %hold on
       subplot(3,2,5)
       activem2hours_onlycon2outlet = timeOutletConnect_timeser_allsim.(genvarname(structnam{i}))*3*3; % cell -> m2 (area of cells = 3 x 3 m2)
       h5(i) = plot(activem2hours_onlycon2outlet,'color',colorlist{1},'linestyle',stylelist{1},'linewidth',1.5);
       hold on
       
       
       snowmelt_duraction(i) = Qmelt_time(end);
       mean_area_connected(i) = median(area_con_outlet);%/(Qmelt_time(end));
       peak_snowmelt_rate(i) = max(Qmelt);
       total_snowmelt(i) = sum(Qmelt); % mm
       total_active_connect_outlet_m2hours(i) = sum(activem2hours_onlycon2outlet);
   
   end 
   
   subplot(3,2,1)
   grid on
   %yyaxis left
   %ylabel({'Snowmelt Input','[mm/h]'})
   yyaxis right
   ylabel({'Cumulative Snowmelt',' Input [mm]'})
   %legend([h1(3),h1(2),h11(1),h11(3),h11(2)],'2009','2011','2010','2009 (cumulative)','2011 (cumulative)','2010 (cumulative)')
   xlabel('Time [h]')
   axis tight
   subplot(3,2,3)
   grid on
   ylabel({'Active (flowing)','area [km^2]'})
   xlabel('Time [h]')
    axis tight
   subplot(3,2,5)
   grid on
    axis tight
   ylabel({'Cumulative active','m^2-hours [h]'})
   %legend([h2(1),h3(1),h2(3),h3(3),h2(2),h3(2)],'2009 (total)','2009 (contributing area to outlet)',...
   %      '2011 (total)','2011 (contributing area to outlet)',...
   %         '2010 (total)','2010 (contributing area to outlet)')
   xlabel('Time [h]')

   subplot(3,2,4)
   yyaxis left
   plot(mean_area_connected,peak_snowmelt_rate,'ko')
   ylabel({'Peak snowmelt','rate [mm/h]'})
   grid on
   xlabel('Average contributing area to outlet [km^2]')
   yyaxis right
   plot(mean_area_connected,total_snowmelt,'k*')
   ylabel({'Total snowmelt','[mm]'})
   grid on
   legend('left axis','right axis')
   
   
   subplot(3,2,6)
   yyaxis left
   plot(total_active_connect_outlet_m2hours,snowmelt_duraction,'ko')
   ylabel({'Total snowmelt','duraction [h]'})
   grid on
   xlabel('Total active m^2-hours [h]')
   yyaxis right
   plot(total_active_connect_outlet_m2hours,total_snowmelt,'k*')
   ylabel({'Total snowmelt,','[mm]'})
   grid on
   legend('left axis','right axis')
   
    % Biogeochemistry
   final_connectbio_timser_allsim = zeros(3,1);
   fig2 = figure('position',[100 100 800 550])
   set(fig2,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
   ic_allsim_vals = [];
    for i=1:numel(structnam)
       subplot(2,2,1)
       ic_allsim_vals = [ic_allsim_vals,ic_allsim.(genvarname(structnam{i}))];
       
       Qmelt = Qmelt_allsim.(genvarname(structnam{i}))(:,2)/24; % mm/day to mm/h
       Qmelt_time = Qmelt_allsim.(genvarname(structnam{i}))(:,1)/3600; % sec to h
       %yyaxis left
       %h1(i) = plot(Qmelt_time,Qmelt,'color',colorlist{1},'linestyle',stylelist{1},'linewidth',1.5)
       %hold on
       %yyaxis right
       h11(i) = plot(Qmelt_time,cumsum(Qmelt),'color',colorlist{1},'linestyle',stylelist{1},'linewidth',1.5);
       hold on
       grid on
       subplot(2,2,3)
       val = connectbio_timser_allsim.(genvarname(structnam{i}));
        h2 (i) = plot(val,'color',colorlist{1},'linestyle',stylelist{1},'linewidth',1.5);
        hold on
        grid on
        final_connectbio_timser_allsim(i) = val(end);
    end
    subplot(2,2,1)
    xlabel('Time [h]')
   %yyaxis left
   %ylabel({'Snowmelt Input','[mm/h]'})
   yyaxis right
   ylabel({'Cumulative Snowmelt',' Input [mm]'})
   %legend([h1(1),h1(3),h1(2),h11(1),h11(3),h11(2)],'2009','2011','2010','2009 (cumulative)','2011 (cumulative)','2010 (cumulative)')
    subplot(2,2,3)
   %legend([h2(1),h2(3),h2(2)],'2009','2011','2010')
   xlabel('Time [h]')
   ylabel({'Biogeochemical','connectivity [-]'})

    subplot(2,2,4)
       yyaxis left
       plot(final_connectbio_timser_allsim,snowmelt_duraction,'ko')
       ylabel({'Total snowmelt','duration [h]'})
       grid on
       xlabel('Biogeochemical connectivity [-]')
       yyaxis right
       plot(final_connectbio_timser_allsim,total_snowmelt,'k*')
       ylabel({'Total snowmelt,','[mm]'})
       grid on
       legend('left axis','right axis')
         
       

   figure('name','Contributing area','position',[100 100 800 600])
   for i = 1:4
       subplot(2,2,i)
       rangeval = 1+(i-1)*23:23+(i-1)*23;
       peak_snowmelt_rate_i = peak_snowmelt_rate(rangeval);
       mean_area_connected_i = mean_area_connected(rangeval);
       total_snowmelt_i = total_snowmelt(rangeval);
       snowmelt_duraction_i = snowmelt_duraction(rangeval);
        F = scatteredInterpolant(peak_snowmelt_rate_i,total_snowmelt_i,mean_area_connected_i,'linear','none');
        Xvec = linspace(min(peak_snowmelt_rate_i),max(peak_snowmelt_rate_i),500);
        Yvec = linspace(min(total_snowmelt_i),max(total_snowmelt_i),500);
        [X,Y] = meshgrid(Xvec,Yvec);
        Zc = F(X,Y);
        surf(X,Y,Zc,Zc)
        intes = 1-((snowmelt_duraction_i-min(snowmelt_duraction_i))/(max(snowmelt_duraction_i)-min(snowmelt_duraction_i)));
        hold on
        max_loc_all = find(snowmelt_duraction_i == max(snowmelt_duraction_i));
        max_loc = max_loc_all(1);
        min_loc_all = find(snowmelt_duraction_i == min(snowmelt_duraction_i));
        min_loc = max_loc_all(1);
        s1 = scatter3(peak_snowmelt_rate_i(max_loc),...
            total_snowmelt_i(max_loc),...
            mean_area_connected_i(min_loc)-10000,...
            40*ones(numel(mean_area_connected_i(max_loc)),1),...
            [0,0,0],'filled','o','MarkerEdgeColor','k');
        hold on
        s2 = scatter3(peak_snowmelt_rate_i(min_loc),...
            total_snowmelt_i(min_loc),...
            mean_area_connected_i(min_loc)-10000,...
            40*ones(numel(mean_area_connected_i(min_loc)),1),...
            [1,1,1],'filled','o','MarkerEdgeColor','k');
         scatter3(peak_snowmelt_rate_i,...
            total_snowmelt_i,...
            mean_area_connected_i-10000,...
            40*ones(numel(mean_area_connected_i),1),...
            [intes,intes,intes],'filled','o','MarkerEdgeColor','k');
        hold on
        xlabel('Peak snowmelt rate [mm/hour]')
        ylabel('Total snowmelt [mm]')
        zlabel('Contributing area [km^2]')
        view(90,-90)
        hc = colorbar;
        title(hc,'[km^2]')
        caxis([0.8 1.6])
        shading interp
        alpha 0.9
   end
    legend([s1,s2],'Max snowmelt duration simulated: 12 days','Min snowmelt duration simulated: 3 days','Location','best')  
  
   
  figure('name','Conection time','position',[100 100 800 600]) 
   for i = 1:4
       subplot(2,2,i)
       rangeval = 1+(i-1)*23:23+(i-1)*23;
       peak_snowmelt_rate_i = peak_snowmelt_rate(rangeval);
       total_active_connect_outlet_m2hours_i = total_active_connect_outlet_m2hours(rangeval)/(24*365);
       total_snowmelt_i = total_snowmelt(rangeval)
        F = scatteredInterpolant(peak_snowmelt_rate_i,total_snowmelt_i,total_active_connect_outlet_m2hours_i,'linear','none');
        Xvec = linspace(min(peak_snowmelt_rate_i),max(peak_snowmelt_rate_i),500);
        Yvec = linspace(min(total_snowmelt_i),max(total_snowmelt_i),500);
        [X,Y] = meshgrid(Xvec,Yvec);
        Zc = F(X,Y);
        surf(X,Y,Zc,Zc)
        intes = 1-((snowmelt_duraction_i-min(snowmelt_duraction_i))/(max(snowmelt_duraction_i)-min(snowmelt_duraction_i)));
        hold on
        max_loc_all = find(snowmelt_duraction_i == max(snowmelt_duraction_i));
        max_loc = max_loc_all(1);
        min_loc_all = find(snowmelt_duraction_i == min(snowmelt_duraction_i));
        min_loc = max_loc_all(1);
        s1 = scatter3(peak_snowmelt_rate_i(max_loc),...
            total_snowmelt_i(max_loc),...
            mean_area_connected_i(min_loc)-10000,...
            40*ones(numel(mean_area_connected_i(max_loc)),1),...
            [0,0,0],'filled','o','MarkerEdgeColor','k');
        hold on
        s2 = scatter3(peak_snowmelt_rate_i(min_loc),...
            total_snowmelt_i(min_loc),...
            mean_area_connected_i(min_loc)-10000,...
            40*ones(numel(mean_area_connected_i(min_loc)),1),...
            [1,1,1],'filled','o','MarkerEdgeColor','k');
         scatter3(peak_snowmelt_rate_i,...
            total_snowmelt_i,...
            mean_area_connected_i-10000,...
            40*ones(numel(mean_area_connected_i),1),...
            [intes,intes,intes],'filled','o','MarkerEdgeColor','k');
        hold on
        xlabel('Peak snowmelt rate [mm/hour]')
        ylabel('Total snowmelt [mm]')
        zlabel('Connectivity years [m^2 \cdot years]')
        view(90,-90)
        caxis([0 (9.5*10^9)/24/365])
        hc = colorbar;
        ylabel(hc,'[m^2 \cdot years]')
        shading interp
        alpha 0.9
 end
 legend([s1,s2],'Max snowmelt duration simulated: 12 days','Min snowmelt duration simulated: 3 days','Location','best')  
  
   

figure
stylelist = {'k*','ro','g^','bs'};
hp = [];
x1_all = [];
y1_all = [];
x2_all = [];
y2_all = [];
   for i = 1:4
        rangeval = 1+(i-1)*23:23+(i-1)*23;
        
       snowmelt_duraction_i = snowmelt_duraction(rangeval);%/max(snowmelt_duraction(rangeval));
       total_snowmelt_i = total_snowmelt(rangeval);%/max(peak_snowmelt_rate(rangeval));
       total_active_connect_outlet_m2hours_i = total_active_connect_outlet_m2hours(rangeval);      
       ic_allsim_vals_i = ic_allsim_vals(rangeval);%/max(ic_allsim_vals(rangeval));
       peak_snowmelt_rate_i = peak_snowmelt_rate(rangeval);%/max(peak_snowmelt_rate(rangeval));
       mean_area_connected_i = mean_area_connected(rangeval);
       
       subplot(1,2,1)
       x1 = (peak_snowmelt_rate_i.*ic_allsim_vals_i')./snowmelt_duraction_i;
       y1 = mean_area_connected_i;
       x1_all = [x1_all;x1];
       y1_all = [y1_all;y1];
       hp_i = plot(x1,y1,stylelist{i});
       hold on
              
       subplot(1,2,2)
       x2 = snowmelt_duraction_i.*total_snowmelt_i;
       y2 = total_active_connect_outlet_m2hours_i;
       x2_all = [x2_all;x2];
       y2_all = [y2_all;y2];
       hp_i = plot(x2,y2,stylelist{i});
       hp = [hp,hp_i(1)];
       hold on
    
   end
subplot(1,2,1)
fit1 = fit(x1_all,y1_all,'power1');
plot(fit1,'k');
xlabel('(Q_p \times S_{init}) / D_m')
ylabel('Average contributing area [km^2]')
grid on
legend('IC = 0.001','IC = 0.005','IC = 0.010','IC = 0.015','power law fit')
dim = [.2 .6 .3 .3];
eq1 = '$\big( \frac{3.65 \cdot Q_p \cdot S_{init}}{D_m} \big)^{1/10}$';
annotation('textbox',dim,'String',eq1,'FitBoxToText','on','Interpreter','latex','Fontsize',14);

subplot(1,2,2)
fit2 = fit(x2_all,y2_all,'power1');
plot(fit2,'k');
hold on
plot(x2,y_pd2);
xlabel('D_m \times V_m')
ylabel('Connection hours [m^2 \cdot years]')
grid on
legend('IC = 0.001','IC = 0.005','IC = 0.010','IC = 0.015','power law fit')
dim = [.2 .6 .3 .3];
eq2 = '$\big(2.7 \cdot D_m \cdot V_{m} \big)^{0.57}$';
annotation('textbox',dim,'String',eq2,'FitBoxToText','on','Interpreter','latex','Fontsize',13);

   
   
end


% Binary 2D array to identify the desired cells
function outputvar = binnary4desiredcells(inputvar)

    outputvar = inputvar * NaN;
    outputvar(inputvar>0) = 1;
            
end


function connect_cells_outlet = processdata(h_res_2D,outlet_xy)

    eval_cells = h_res_2D * 0;
    eval_cells(h_res_2D==0) = NaN;
    eval_cells(outlet_xy(1),outlet_xy(2)) = 1;
    connect_cells_outlet = h_res_2D * NaN;
    intpoints_loc = outlet_xy;
    contsearch = 1;
   
    while contsearch % iteracte until covering the entire domain
        
        contsearch = 0;
        num_aciveIntpoints = numel(intpoints_loc(:,1));
        
        locpairs_total = [];
        
        for p = 1:num_aciveIntpoints
            ie = intpoints_loc(p,1) + 1;
            iw = intpoints_loc(p,1) - 1;
            in = intpoints_loc(p,2) + 1;
            is = intpoints_loc(p,2) - 1;
            
            locpairs = [intpoints_loc(p,1),in;
                        intpoints_loc(p,1),is;
                        ie,intpoints_loc(p,2);
                        iw,intpoints_loc(p,2)];   
                    
            locpairs_rev = [];
            
            for i = 1:numel(locpairs(:,1))
                try
                    if eval_cells(locpairs(i,1),locpairs(i,2))==0
                        hcell = h_res_2D(locpairs(i,1),locpairs(i,2));
                        if hcell > 0
                            eval_cells(locpairs(i,1),locpairs(i,2)) = 1;
                            connect_cells_outlet(locpairs(i,1),locpairs(i,2)) = 1;
                            contsearch = 1;
                        else
                            eval_cells(locpairs(i,1),locpairs(i,2)) = NaN;
                        end
                        locpairs_rev = [locpairs_rev;locpairs(i,:)];
                    end
                catch
                    
                end
            end     
            
            locpairs_total = [locpairs_total;locpairs_rev];
        
        end
        
        intpoints_loc = locpairs_total;
        
    end   
end



function h_res_2D = get2Dmatrixresults(dataraw_rel,timestep,simname,coloc)

h_res_2D = zeros(max(dataraw_rel(:,1)),max(dataraw_rel(:,2)));

for i = 1:numel(dataraw_rel(:,1))
    ix = dataraw_rel(i,1);
    iy = dataraw_rel(i,2);
    try
        h_res_2D(ix,iy) = dataraw_rel(i,coloc);
    catch
        disp(['Invalid entry (=',num2str(dataraw_rel(i,coloc)) ,') at sim = ', simname, ', timestep = ', num2str(timestep),', line ',num2str(i),', ix = ',num2str(ix),' and iy = ',num2str(iy)])
    end
end

end


function steps_i_nam_num = gettimesteps(steps_i_nam)

    steps_i_nam_num = [];

    for i = 1:numel(steps_i_nam)

        steps_i_nam_str = steps_i_nam{i};

        if numel(steps_i_nam_str) > 3
            steps_i_nam_str = steps_i_nam_str(1:end-4);

            try
                steps_i_nam_num_i = str2num(steps_i_nam_str);
                steps_i_nam_num = [steps_i_nam_num,steps_i_nam_num_i];
            catch
            end

    end

    end
    
end

