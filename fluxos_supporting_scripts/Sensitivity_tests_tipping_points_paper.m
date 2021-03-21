
% PREPARE THE SENSITIVITY ANALYSIS FOR FLUXOS (looking for tipping points)

clc
close all
clear all

format long g

Qmelt_CRHM_analysis_STC_flag = 1 % 2) STC Diogo
Qmelt_CRHM_analysis_SD_flag = 0; % 1) VB Kevin (reanalysis)

generate_FLUXOS_files_flag = 0;

gen_Qmelt_files_sens_gammadistrib_flag = 1; % obsolete

%% CRHM-based
if Qmelt_CRHM_analysis_STC_flag
    
    CRHM_data_hourly = importdata('/media/dcosta/data/megasync/ec_main/models/crhm/support/PROJECTS/CRHM_for_FLUXOS_STC_Diogo//CRHM_output_1.txt'); % my STC basin (Diogo)
    snomelt_col_start = 1+42+42+1;
    snomelt_D_col_start = 1+42+42+42+1;
    snowmelt_runoff_col_start = 1+42+42+42+42+1;
    snowmewlt = CRHM_data_hourly.data(:,snomelt_col_start:snomelt_col_start+41);
    snowmewlt_D = CRHM_data_hourly.data(:,snomelt_D_col_start:snomelt_D_col_start+41);
    snowmelt_runoff = CRHM_data_hourly.data(:,snowmelt_runoff_col_start:snowmelt_runoff_col_start+41);
    infilt_col_start = 2;
    infiltration = CRHM_data_hourly.data(:,infilt_col_start:infilt_col_start+41);
    snowmelt_time = datenum(CRHM_data_hourly.data(:,1) + 693960);
    snowmelt_time_vec = datevec(snowmelt_time);

    snowmelt_time_vec_yeardyn = snowmelt_time_vec;
    snowmelt_time_vec_yeardyn(:,1) = 0;

    % separate the events
    Snowmelt_events = {};
    Snowmelt_events_runoff = {};
    event_i = 1;
    snowmewlt_iprev_D = 0;
    Snowmelt_events_i = [];
    Snowmelt_events_runoff_i = [];
    event_duration_all = [];
    snowmelt_peak_all = [];
    snowmelt_vol_all = [];

    Snowmelt_events_infilt = {};
    Snowmelt_events_infilt_i = [];
    infilt_rate_i = [];
   
    snowmewlt_avrgHRU = mean(snowmewlt');
    snowmewlt_avrgHRU_D = mean(snowmewlt_D');
    snowmelt_runoff_avrgHRU = mean(snowmelt_runoff');
    
    figure
    for i=1:numel(snowmelt_time)
        
        
        snowmewlt_i = snowmewlt_avrgHRU(i);
        snowmewlt_i_D = snowmewlt_avrgHRU_D(i);
        meltrunoff_i = snowmelt_runoff_avrgHRU(i);
        infilt_rate_i = infiltration(i);
        
        if snowmewlt_i_D ~= 0
            Snowmelt_events_i = [Snowmelt_events_i,snowmewlt_i];  
            Snowmelt_events_runoff_i = [Snowmelt_events_runoff_i,meltrunoff_i];
             Snowmelt_events_infilt_i = [Snowmelt_events_infilt_i,infilt_rate_i];
        end
        
        if snowmewlt_i_D == 0 && snowmewlt_iprev_D~=0
            plot(Snowmelt_events_i)
            hold on
            grid on
            Snowmelt_events = [Snowmelt_events;{Snowmelt_events_i}];
            Snowmelt_events_runoff = [Snowmelt_events_runoff;{Snowmelt_events_runoff_i}];
            event_duration_all = [event_duration_all,numel(Snowmelt_events_i)];
            snowmelt_peak_all = [snowmelt_peak_all,max(Snowmelt_events_i)];
            snowmelt_vol_all = [snowmelt_vol_all,sum(Snowmelt_events_i)];
             Snowmelt_events_infilt = [Snowmelt_events_infilt;{Snowmelt_events_infilt_i}];
             Snowmelt_events_infilt_i = [];
             Snowmelt_events_runoff_i = [];
            Snowmelt_events_i = [];
            
        end
        
        snowmewlt_iprev_D = snowmewlt_i_D;
    end
    xlabel('Time [days]')
    ylabel('Snowmelt + rainfall [mm/hour]')
    
   % Statistical analysis
   figure
   subplot(121)
   histogram(event_duration_all)
   xlabel('Snowmelt duration [days]')
   ylabel('# of events')
   grid on
   hold on
   subplot(122)
   pd = fitdist(event_duration_all','Lognormal');
   x = [1:.1:max(event_duration_all)];
   y = lognpdf(x,pd.mu,pd.sigma);
   plot(x,y)
   grid on
   xlabel('Snowmelt duration [days]')
   ylabel('lognormal pdf []')
   
   figure
   subplot(121)
   histogram(snowmelt_vol_all)
   ylabel('Total Snowmelt Volume  [mm]')
   ylabel('# of events')
   grid on
   hold on
   subplot(122)
   %pd = fitdist(snowmelt_vol_all','Lognormal');
   x = [1:.1:max(snowmelt_vol_all)];
   y = lognpdf(x,pd.mu,pd.sigma);
   plot(x,y)
   grid on
   xlabel('Total Snowmelt Volume  [mm]')
   ylabel('lognormal pdf []')
   
   figure
   subplot(121)
   histogram(snowmelt_peak_all)
   xlabel('Peak snowmelt rate [mm/hour]')
   ylabel('# of events')
   grid on
   hold on
   subplot(122)
%   pd = fitdist(snowmelt_peak_all','Lognormal');
   x = [1:.1:max(snowmelt_peak_all)];
   y = lognpdf(x,pd.mu,pd.sigma);
   plot(x,y)
   grid on
   xlabel('Peak snowmelt rate [mm/hour]')
   ylabel('lognormal pdf []')
   
   % linear comparison and looking for the events
   figure
   subplot(1,3,1)
   plot(event_duration_all,snowmelt_vol_all,'ok');
   grid on
   xlabel('Snowmelt duration [days]')
   ylabel('Total Snowmelt Volume  [mm]')
   %sim_scenar = [2,1.54; 3,4.88;9,3.10;8,7.99;13,14.23;13,2.93;19,15.86;21,7.33;10,19.42;6,0.66];
   sim_scenar = [[48,7.17840294664286];[240,39.3202350735714];[240,8.07470955280952];[72,34.2433452952381];[48,1.18305000000000];[264,73.2577423294048];[144,27.2008422538095];[72,11.4236627904762]];
   sim_scenar = round(sim_scenar,2);
   col_find_all_1 = [];
   for i = 1:numel(sim_scenar(:,1))
    col_find_i = find(event_duration_all==sim_scenar(i,1) & round(snowmelt_vol_all,2) == sim_scenar(i,2));
    col_find_all_1 = [col_find_all_1,col_find_i];
   end
   hold on
   scatter(event_duration_all(col_find_all_1),snowmelt_vol_all(col_find_all_1),'or','filled');

        subplot(1,3,2)
   plot(event_duration_all,snowmelt_peak_all,'ok');
   grid on
   xlabel('Snowmelt duration [days]')
   ylabel('Peak snowmelt rate [mm/hour]')
   %sim_scenar = [3,0.49;10,0.48;1,1.32;13,1.67;4,2.80;17,2.40;5,4.41;8,1.24;21,1.44;8,2.11];
   sim_scenar = [[168,1.95368761904762];[264,3.71773142857143];[72,3.59332833333334];[240,2.55173333333333];[48,2.18039000000000];[48,1.18305000000000];[48,0.332653023809524];[240,0.606458642857143]];
   sim_scenar = round(sim_scenar,2); 
   col_find_all_2 = [];
   for i = 1:numel(sim_scenar(:,1))
    col_find_i = find(event_duration_all==sim_scenar(i,1) & round(snowmelt_peak_all,2) == sim_scenar(i,2));
    col_find_all_2 = [col_find_all_2,col_find_i];
   end
    hold on
   scatter(event_duration_all(col_find_all_2),snowmelt_peak_all(col_find_all_2),'or','filled');
    
   subplot(1,3,3)
   plot(snowmelt_vol_all,snowmelt_peak_all,'ok');
   grid on
   xlabel('Snowmelt volume / SWE [mm]')
   ylabel('Peak snowmelt rate [mm/hour]')
   % sim_scenar = [0.54,0.29;2.69,1.05; 6.11,0.95;2.19,2.19;10.01,2.00;10.25,4.41;19.28,2.40;5.51,1.65;15.49,1.44;4.88,2.97;7.33,0.91];
   sim_scenar = [[9.04691000000000,2.18039000000000];[14.1793183002381,1.48434619047619];[73.2577423294048,3.71773142857143];[34.2433452952381,3.59332833333334];[4.13012348000000,2.80013000000000];[33.3545104907143,2.39682500000000];[8.07470955280952,0.606458642857143]];
    sim_scenar = round(sim_scenar,2)
   col_find_all_3 = [];
   for i = 1:numel(sim_scenar(:,1))
    col_find_i = find(round(snowmelt_vol_all,2)==sim_scenar(i,1) & round(snowmelt_peak_all,2) == sim_scenar(i,2));
    col_find_all_3 = [col_find_all_3,col_find_i];
   end
    hold on
   scatter(snowmelt_vol_all(col_find_all_3),snowmelt_peak_all(col_find_all_3),'or','filled');
   
   col_find_all = [];
   col_find_all = [col_find_all_1,col_find_all_2,col_find_all_3];
   
   % Plotting the snowmelt events to simulate
   figure
   for i = 1:numel(col_find_all)
    subplot(2,1,1)
    plot(Snowmelt_events{col_find_all(i)})
    hold on
    subplot(2,1,2)
    plot(Snowmelt_events_infilt{col_find_all(i)}/24)
    hold on
   end
   subplot(2,1,1)
   grid on
   xlabel('Time [days]')
   ylabel('Snowmelt rate [mm/hour]')
   subplot(2,1,2)
   grid on
   xlabel('Time [days]')
   ylabel('Infiltration rate [mm/hour]')
   
   % Generate the FLUXOS files
   if generate_FLUXOS_files_flag
        generate_FLUXOS_files(Snowmelt_events(col_find_all)* 24); %mm/h -> mm/day;
   end
   
end



%% Gamma distribution-based
if gen_Qmelt_files_sens_gammadistrib_flag
%% Parameter combinations

    SnowmeltVol_SWE_1 = [50,150,300];
    PeakSnowmelt_rate_2 = [0.1,1,3,5];
    Snowmelt_duration_3 = [2,15,30] * 24; % days -> hour
    Initial_conditions_4 = [0.005,0.01,0.015];


    combinations_all = [];

    for i1 = 1:numel(SnowmeltVol_SWE_1)  
        for i2 = 1:numel(PeakSnowmelt_rate_2)   
            for i3 = 1:numel(Snowmelt_duration_3)
                for i4 = 1:numel(Initial_conditions_4)
                    combinations_i = [SnowmeltVol_SWE_1(i1),PeakSnowmelt_rate_2(i2),Snowmelt_duration_3(i3),Initial_conditions_4(i4)];
                    combinations_all = [combinations_all;combinations_i];
                end 
            end  
        end 
    end


    num_comb = numel(combinations_all(:,1));
    storage_needed_TB = num_comb * 20 / 1000;

    %% Determining the gamma curve parameter A and B
    parA_range = 1:0.5:100;
    parB_range = 1:2:500;
    par_numel = numel(parA_range) * numel(parB_range);
    par_A_B_all = [];
    parA_selected = zeros(num_comb,1) * NaN;
    parB_selected = zeros(num_comb,1) * NaN;

    % get all parA and parB combinations
    for i = 1:numel(parA_range)
        for j = 1:numel(parB_range)
            par_A_B_all = [par_A_B_all;[parA_range(i),parB_range(j)]];
        end
    end

    hbar = parfor_progressbar(num_comb, 'Calculating...');
    fig1 = figure;
    fig2 = figure;
    for sim_i = 1:num_comb

        SnowmeltVol_SWE_1_i = combinations_all(sim_i,1);
        PeakSnowmelt_rate_2_i = combinations_all(sim_i,2);
        Snowmelt_duration_3_i = combinations_all(sim_i,3); % days -> hour
        Initial_conditions_4_i = combinations_all(sim_i,4);

        Qmelt_time = 0:1:Snowmelt_duration_3_i; % time (this is also fixed)
        Qmelt = 0;
        % look for best combination
        i = 1;
        figure(fig1)
        plot(Qmelt_time([1 end]),[SnowmeltVol_SWE_1_i  SnowmeltVol_SWE_1_i],'r-','linewidth',3)
        ylim([0 SnowmeltVol_SWE_1_i*2])
        for i = 1:par_numel  % match the peak snowmelt rate
            parA_i = par_A_B_all(i,1);
            parB_i = par_A_B_all(i,2);
            Y = gampdf(Qmelt_time,parA_i,parB_i);
            Qmelt = Y * SnowmeltVol_SWE_1_i; % max volume is guaranteed because the cdf is = 1 (so multiplying by SWE will give that SWE   
            figure(fig1)
            plot(Qmelt_time,Qmelt)
            ylim([0 SnowmeltVol_SWE_1_i*2])
            hold on

            Ycum = gamcdf(Qmelt_time,parA_i,parB_i);
            if round(max(Qmelt),2) == PeakSnowmelt_rate_2_i && Ycum(end) * SnowmeltVol_SWE_1_i > 0.95 * SnowmeltVol_SWE_1_i
                break
            end

        end
        if round(max(Qmelt),2) == SnowmeltVol_SWE_1_i
           parA_selected(sim_i) = parA_i;
           parB_selected(sim_i) = parB_i;
           figure(fig2)
           plot(Qmelt_time,Qmelt)
           ylim([0 combinations_all(sim_i,2)*2])
           hold on
        end


        hbar.iterate(1)

    end
    close(h)

end


%% St Denis
if Qmelt_CRHM_analysis_SD_flag
    
    %CRHM_output_all = importdata('Scenarios_FLUXOS/CRHM_output_1901-2000_final.obs');
    runoff_CRHM_daily = importdata('/media/dcosta/DATADRIVE1/MegaSync/2_CRHM_WQ/1_Model/PROJECTS/CRHM_VB_Kevin/SD_runoff_1901-2000.dat'); % VB project (Kevin)
    runoff = runoff_CRHM_daily.data;
    runoff_time = datenum(runoff_CRHM_daily.textdata(2:end,1));

    figure;
    plot(runoff_time,runoff)
    datetick('x','yyyy')
    
end



%% Generate FLUXOS files
function generate_FLUXOS_files(meltrunoff)
    %% Get snowmelt rate to print in "Qmelt.fluxos
   
    IC = [0.001,0.005,0.01,0.015];
    numIC = numel(IC);
    
    numSnowmeltScen = numel(meltrunoff);
    
    dirfold = '/media/dcosta/DATADRIVE1/fluxos_tests/SIMULATIONS_sync/batch_CRHM_scenarios/';
    zerofile = importdata([dirfold,'0.txt']);
    
    for i = 1:numIC
        h = waitbar(0,'Calculating...')
        for s = 1:numSnowmeltScen
            
            IC_i = IC(i);
            S_i = meltrunoff{s};
            S_i = [S_i,zeros(24,1)']; % extend the simulation by 1 day so that the runoff stabilizes
            
            simname = ['IC',num2str(i),'_S',num2str(s)]; % folder name
            simfolder = [dirfold,simname];
            
            % Make directory for simulation
            mkdir(simfolder);
            mkdir([dirfold,simname,'/Results']);
            
            % Prepare job.sh file
            fid = fopen([simfolder,'/job.sh'],'W');
            jobshfile= {'#!/bin/bash',...
                '#SBATCH --account=rpp-hwheater',...
                '#SBATCH --time=0-240:0',...
                '#SBATCH --cpus-per-task=8',...
                'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK',...
                './fluxos_cpp'};
            for l=1:numel(jobshfile)
                jobshfile_i = jobshfile{l};
                if l<numel(jobshfile)
                   styleformat= ['%',num2str(numel(jobshfile_i)),'s\n'];
                else
                    styleformat= ['%',num2str(numel(jobshfile_i)),'s'];
                end
                 fprintf(fid,styleformat,jobshfile_i);
            end
            fclose(fid)
            
            % Copy ersi_grid_dem_alldomain,ersi_grid_dem_basin and fluxos_cpp 
            copyfile([dirfold,'ersi_grid_dem_alldomain'],simfolder);
            copyfile([dirfold,'ersi_grid_dem_basin'],simfolder);
            copyfile([dirfold,'fluxos_cpp'],simfolder);
            
            % Create modset.fluxos file
            fid = fopen([simfolder,'/modset.fluxos'],'W');
            modsetfile= {'### Batch_CRHM_scenarios ###';...
                        'ersi_grid_dem_alldomain';...
                        'ersi_grid_dem_basin';...
                        'Qmelt_info.fluxos';...
                        '3600';...
                        num2str(IC_i);...
                        '3';...
                        '0.08';...
                        '1';...
                        '9.5';...
                        '9.5'};
            for l=1:numel(modsetfile)
                modsetfile_i = modsetfile{l};
                if l<numel(modsetfile)
                   styleformat= ['%',num2str(numel(modsetfile_i)),'s\n'];
                else
                    styleformat= ['%',num2str(numel(modsetfile_i)),'s'];
                end
                 fprintf(fid,styleformat,modsetfile_i);
            end
            fclose(fid)
            
            % Print Qmelt_info
            Qmelt_i = S_i;
            
            time_print = [0:3600:numel(Qmelt_i)*3600]';
            Qmelt_i_print = [0,Qmelt_i]';           
            array_2print = [time_print,Qmelt_i_print];
            
            fid = fopen([simfolder,'/Qmelt_info.fluxos'],'W');
            formatstyle = '%1.0f,%f\n';
            for l = 1:numel(array_2print(:,1))
                fprintf(fid,formatstyle,array_2print(l,1),array_2print(l,2));
            end
            fclose(fid)
            
            % Prepare and print the 0.txt file
            zerofile(:,4) = IC_i;
            csvwrite([simfolder,'/Results/0.txt'],zerofile);
            waitbar(s/numSnowmeltScen)
        end
        close(h)
    end
   
end
















