
CrossSections_outFiles_flag = 1;
MedianMax_velocity_flag = 0;


FLUXOS_res_dir = '/media/dcosta/data/megasync/my_server/fluxos/';
batch_dir = 'batch_1_selected_paper_additional_previousFLUXOSversion';
%batch_dir = 'batch_1_selected_paper_additional_previousFLUXOSversion_graham';
%batch_dir = 'batch_1_select_paper';
%batch_dir = 'batch_1_selected_paper_additional';
%batch_dir = 'batch_1_selected_paper_additional_graham';
%batch_dir = 'batch_1_selected_paper_additional_noIC';
%batch_dir = 'batch_1_selected_paper_additional_noIC_graham';

if CrossSections_outFiles_flag
     %%%%%%%%%% 
    yearselect = 2012;
    ResType = 1; %1-flow, 2-WQ, 3-SQ
    Obs_col = 2;      % FLOW: Obs_col = 2 
                      % NH4: Obs: Obs_col = 2 
                      % DOC: Obs_col = 3
                      % SUSPC: Obs_col = 4 	
                      % TOC: Obs_col = 5
                      % NO3NO2: Obs_col = 6 
                      % SUSPN: Obs_col = 7 
                      % TDN: Obs_col = 8
                      % TN: Obs_col = 9
                      % SUSPP: Obs_col = 10 	
                      % TP: Obs_col = 11
                      % TDP: Obs_col = 12 
                      % SRP: Obs_col = 13
                      % TSS: Obs_col = 14
    lag = 8; % in hours

    if ResType == 1
        outfilenam = 'f.out';
    elseif ResType == 2
        outfilenam = 'wq.out';
    elseif ResType == 3
        outfilenam = 'sq.out';
    end

    if (yearselect==2005)
        fluxos_timestart = 38440.8437499995 + 695422 - lag/24; 
    elseif (yearselect==2006)
        fluxos_timestart = 38809.1874999995 + 695422 - lag/24; 
    elseif (yearselect==2009)
       fluxos_timestart = 39913.01042 + 695422 - lag/24; 
    elseif (yearselect==2010)
       fluxos_timestart = 40252.03125 + 695422 - lag/24;       
    elseif (yearselect==2011)
        fluxos_timestart = 40633 + 695422 - lag/24;
    elseif (yearselect==2012)
        fluxos_timestart = 40979.3333333333 + 695422 - lag/24;  
    elseif (yearselect==2013)
        fluxos_timestart = 41389.90625 + 695422 - lag/24;  
    elseif (yearselect==2014)
        fluxos_timestart = 41737.1354166667 + 695422 - lag/24;  
    elseif (yearselect==2015)
        fluxos_timestart = 42156.0729166667 + 695422 - lag/24;  
    end

    [resultdir_list, obsPath] = get_resultdir_list(FLUXOS_res_dir,batch_dir,yearselect,ResType);

    % Load Obs
    if ~isempty(obsPath)
        obsdata = importdata(obsPath);
        %time_obs = cumsum([0; diff(obsdata.data(:,1))])*24; % days -> hour
        time_obs = obsdata.data(:,1) +  695422;
        data_obs = obsdata.data(:,Obs_col);
    end

    %%%%%%%%% MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure
    resultdir_legend = {};
    if ~isempty(obsPath)
        p(1) = plot(time_obs,data_obs,'or');
    end
    hold on
    timmod_min = [];
    timmod_max = [];
    hw = waitbar(0,'Calculating...');
    numloop = numel(resultdir_list);
    for i = 1:numloop

    % load FLUXOS cs results
    try
        res = importdata([resultdir_list{i},'cs/',outfilenam]);
        resultdir_legend = [resultdir_legend,resultdir_list{i}];
    catch
        disp(['Result (f.out) not found for: "',resultdir_list{i},'" (SKIPPED)'])
        continue
    end

    % Results
    time_mod = fluxos_timestart + res(:,1)/(3600*24); % sec -> hour
    data_mod = res(:,2:end);
    timmod_min = min([timmod_min; time_mod]);
    timmod_max = max([timmod_max; time_mod]);

    %% Model results to plot
    fluxos_time  = time_mod;
    fluxos_results = (sum(data_mod'))'; % summing up flow on all cells of the cross-section
    fluxos_results(fluxos_results<0.2) = max( fluxos_results(fluxos_results<0.2) - 0.1,0);
    
    %% Calculate NSE, RMSE and Bias
     if ~isempty(obsPath)
        Model_data_interc = interp1(fluxos_time,fluxos_results,time_obs);

        %[r2 rmse] = rsquare(Z_obs_mesh,Model_data_interc);

         % remove NaNs
         Model_data_interc_use = Model_data_interc;
         WQ_use = data_obs;


        nanloc = find(isnan(Model_data_interc_use));
        Model_data_interc_use(nanloc) = [];
        WQ_use(nanloc) = [];

        nanloc = find(isnan(WQ_use));
        Model_data_interc_use(nanloc) = [];
        WQ_use(nanloc) = [];

        % Nash
        numerator=(WQ_use-Model_data_interc_use).^2;
        denominator=(WQ_use-mean(WQ_use)).^2;
        Nash =1-(sum(numerator)/sum(denominator));

        % RMSE
        Sumcal = (Model_data_interc_use-WQ_use).^2;
        numerator = sum(Sumcal);
        n=numel(WQ_use);
        RMSE=(numerator/n)^(1/2);

        % BIAS
        numerator = sum(WQ_use);
        denominator = sum(Model_data_interc_use);
        Bias = numerator/denominator-1;
        
    
     end
    %% Plot
    p(i+1) = plot(fluxos_time,fluxos_results,'linewidth',1.5,'Color',[0.65 0.65 0.65])%[0 0 0]+1/(numel(resultdir_list)+5) * i)
    hold on
    waitbar(i / numloop)
    end
    %xlabel('Time [hour]')
    axis tight

    if ResType == 1
        ylabel('Flow [m3/s]')
        ylim([0 1.8])
    elseif ResType == 2
        ylabel('Conc [mg/l]')
        ylim([0 1.8])
    end


    alpha 0.7
    if ~isempty(obsPath)
        %xlim([timmod_min timmod_max])
    end
    grid on
    legend(['obs',resultdir_legend],'interpreter','none')
    title(yearselect)
    datetick('x','dd-mmm','keeplimits','keepticks')
end


if MedianMax_velocity_flag
    
    newdir = [FLUXOS_res_dir,batch_dir];
    
     resultdir_list = {
                    [newdir,'/t_49_paper/Results/'], % 2009
                    [newdir,'/t_65_paper/Results/'], % 2010
                    [newdir,'/t_36_paper/Results/'], % 2011
                    };               
                
     fn_1_col = 13;
     fe_1_col = 14;
     h_col = 4;
     
    qmag_median_allsim = {};
    qmag_median_nozero_allsim = {};
    qmag_std_allsim = {};
    qmag_std_nozero_allsim = {};
    qmag_mean_allsim = {};
    qmag_mean_nozero_allsim = {};
    qmag_max_allsim = {};

    vmag_median_allsim = {};
    vmag_median_nozero_allsim = {};
    vmag_std_allsim = {};
    vmag_std_nozero_allsim = {};
    vmag_mean_allsim = {};
    vmag_mean_nozero_allsim = {};
    vmag_max_allsim = {};
    
     for i = 1:numel(resultdir_list)
        
        resultdir = resultdir_list{i};
        resfiles_raw = dir(resultdir);
        resfiles = {resfiles_raw.name};
        resfiles_loc = contains(resfiles,'.txt');
        resfiles = resfiles(resfiles_loc);
        resfiles_loc2 = ~contains(resfiles,'lock');
        resfiles = resfiles(resfiles_loc2);
        simname = ['Sim ',num2str(i),' (',resultdir,')'];
        
        qmag_median = [];
        qmag_median_nozero = [];
        qmag_std = [];
        qmag_std_nozero = [];
        qmag_mean = [];
        qmag_mean_nozero = [];
        qmag_max = [];
        
        vmag_median = [];
        vmag_median_nozero = [];
        vmag_std = [];
        vmag_std_nozero = [];
        vmag_mean = [];
        vmag_mean_nozero = [];
        vmag_max = [];
        
        hw = waitbar(0,simname)
        for s = 1:numel(resfiles)
            resfiles_i = resfiles{s};
            dataraw = importdata([resultdir,resfiles_i]);
            
            qx = dataraw(:,fn_1_col);
            qy = dataraw(:,fe_1_col);
            h = dataraw(:,h_col);
            qmag = (qx.^2+qy.^2).^0.5;
            %iloc = qmag>0;
            %qmag = qmag(iloc);
            %h = h(iloc);
            
            vmag = qmag./h;
            
            qmag_median = [qmag_median,median(qmag)];
            qmag_median_nozero = [qmag_median_nozero,median(qmag(qmag>0))];
            qmag_std = [qmag_std,std(qmag)];
            qmag_std_nozero = [qmag_std_nozero,std(qmag(qmag>0))];
            qmag_mean = [qmag_mean,mean(qmag)];
            qmag_mean_nozero = [qmag_mean_nozero,mean(qmag(qmag>0))]; 
            qmag_max = [qmag_max,max(qmag)];
            
            vmag_median = [vmag_median,median(vmag)];
            vmag_median_nozero = [vmag_median_nozero,median(vmag(vmag>0))];
            vmag_std = [vmag_std,std(vmag)];
            vmag_std_nozero = [vmag_std_nozero,std(vmag(vmag>0))];
            vmag_mean = [vmag_mean,mean(vmag)];
            vmag_mean_nozero = [vmag_mean_nozero,mean(vmag(vmag>0))]; 
            vmag_max = [vmag_max,mean(vmag(qmag == max(qmag)))];
            
            waitbar(s / numel(resfiles))
        end
        close(hw)
        
        qmag_median_allsim = [qmag_median_allsim;qmag_median];
        qmag_median_nozero_allsim = [qmag_median_nozero_allsim;qmag_median_nozero];
        qmag_std_allsim = [qmag_std_allsim;qmag_std];
        qmag_std_nozero_allsim = [qmag_std_nozero_allsim;qmag_std_nozero];
        qmag_mean_allsim = [qmag_mean_allsim;qmag_mean];
        qmag_mean_nozero_allsim = [qmag_mean_nozero_allsim;qmag_mean_nozero];
        qmag_max_allsim = [qmag_max_allsim;qmag_max];

        vmag_median_allsim = [vmag_median_allsim;vmag_median];
        vmag_median_nozero_allsim = [vmag_median_nozero_allsim;vmag_median_nozero];
        vmag_std_allsim = [vmag_std_allsim;vmag_std];
        vmag_std_nozero_allsim = [vmag_std_nozero_allsim;vmag_std_nozero];
        vmag_mean_allsim = [vmag_mean_allsim;vmag_mean];
        vmag_mean_nozero_allsim = [vmag_mean_nozero_allsim;vmag_mean_nozero];
        vmag_max_allsim = [vmag_max_allsim;vmag_max];
     end
   
   results_all = {qmag_median_allsim,...
                qmag_median_nozero_allsim,...
                qmag_std_allsim,...
                qmag_std_nozero_allsim,...
                qmag_mean_allsim,...
                qmag_mean_nozero_allsim,...
                qmag_max_allsim,...
                vmag_median_allsim,...
                vmag_median_nozero_allsim,...
                vmag_std_allsim,...
                vmag_std_nozero_allsim,...
                vmag_mean_allsim,...
                vmag_mean_nozero_allsim,...
                vmag_max_allsim};
   
    title_all = {'qmag_median_allsim',... %1
                'qmag_median_nozero_allsim',... %2
                'qmag_std_allsim',... %3
                'qmag_std_nozero_allsim',... %4
                'qmag_mean_allsim',... %5
                'qmag_mean_nozero_allsim',... %6
                'qmag_max_allsim',... %7
                'vmag_median_allsim',... %8
                'vmag_median_nozero_allsim',... %9
                'vmag_std_allsim',... %10
                'vmag_std_nozero_allsim',... %11
                'vmag_mean_allsim',... %12
                'vmag_mean_nozero_allsim',... %13
                'vmag_max_allsim'}; %14        
            
   years = [2009,2010,2011];
            
   lag = 8;
   fluxos_timestart(1) = 39913.01042; %2009
   fluxos_timestart(2) = 40252.03125; %2010
   fluxos_timestart(3) = 40633;       %2011
    
    % All Results
    for v = 1:14
    
        results = results_all{v};    

        figure('name',title_all{v})
        for i = 1:3
            subplot(3,1,i)
            t_init = fluxos_timestart(i) + + 695422 - lag/24;
            time = [0:3600:numel(results{i})*3600-3600]/(3600*24) + t_init;
            plot(time,results{i},'k-');
            title(num2str(years(i)))
            datetick('x','dd-mmm','keeplimits','keepticks')
            grid on
        end   
    end
    
    % Costumized plots      
    v_selected_all = {[9,12,2,5],...
                   [14,7]};
    
    color_list = [0 0 0;...
                 0 0 0;...
                 1 0 0;...
                 1 0 0];
    linewidth_list = [2,...
                      1,...
                      2,...
                      1];
    style_list = {'-',...
                  '--',...
                  '-',...
                  '--'};
    
     ymax_all = [0.1,0.006;...
                50,1];         
              
    for f = 1:2          
   
        v_selected = v_selected_all{f};
        
        fig(f) = figure;
        set(fig(f),'defaultAxesColorOrder',[0 0 0; 0 0 0]);

        for vi = 1:numel(v_selected)
            results = results_all{v_selected(vi)};  

            for i = 1:3
                subplot(3,1,i)
                t_init = fluxos_timestart(i) + + 695422 - lag/24;
                time = [0:3600:numel(results{i})*3600-3600]/(3600*24) + t_init;

                if v_selected(vi) >= 8
                    yyaxis left
                    ymax = ymax_all(f,1);
                else
                    yyaxis right
                    ymax = ymax_all(f,2);
                end
                ressplie = csaps(time,results{i},time);
                plot(time,ressplie,'color',color_list(vi,:),'linewidth',linewidth_list(vi),'linestyle',style_list{vi});
                hold on
                title(num2str(years(i)))
                datetick('x','dd-mmm','keeplimits','keepticks')
                ylim([0 ymax]);
                grid on
                if v_selected(vi) <= 8
                    ylabel('Runoff / Flow [m^2/s]')
                else
                   ylabel('Velocity [m/s]')
                end

            end

        end
        if f == 1
            legend('Median velocity [m/s]','Mean velocity [m/s]','Median runoff/flow [m^2/s]','Mean runoff/flow [m^2/s]','orientation','vertical')
        elseif f ==2
            legend('Max velocity [m/s]','Max runoff/flow [m^2/s]','orientation','vertical')
        end
    end
end

