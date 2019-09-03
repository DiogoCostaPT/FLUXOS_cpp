
%%%%%%%%%%% OBSERVATIONS %%%%%%%%%%%%%%%%%% 
yearselect = 2011;

% Load Obs
if (yearselect==2009)
	obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/1_2009_Compiled/Streamflow_MS9C_2009_trimmed_for_simulation.csv';
elseif (yearselect==2010)
    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/2_2010_Compiled/Streamflow_MS9C_2010_trimmed_to_simulation.csv';
elseif (yearselect==2011)
    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/0_Obs_used_in_first_2011_tests/Snowmelt_Runoff_MS9_2011_justflow.csv';
end
obsdata = importdata(obsPath);

% Obs
time_obs = cumsum([0; diff(obsdata.data(:,1))])*24; % days -> hour
data_obs = obsdata.data(:,2:end);

%%%%%%%%% MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (yearselect==2009)
	resultdir_list = {...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_55/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_52/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_40/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_45/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_57/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_46/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_59/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_49/Results/',...
                };
elseif (yearselect==2010)
    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/2_2010_Compiled/Streamflow_MS9C_2010_trimmed_to_simulation.csv';
    resultdir_list = {...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_56/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_53/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_41/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_44/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_58/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_47/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_60/Results/',...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_50/Results/',...
                };
elseif (yearselect==2011)
    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/0_Obs_used_in_first_2011_tests/Snowmelt_Runoff_MS9_2011_justflow.csv';
    resultdir_list = {...
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_39/Results/',...
               % '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_54/Results/',...
                %'/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_36/Results/',...
               % '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_43/Results/',...
                %'/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_35/Results/',...
               % '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_48/Results/',...
               % '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_37/Results/',...
               % '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_51/Results/',...
                };
end


figure
resultdir_legend = {};
plot(time_obs,data_obs,'ok')
hold on
for i = 1:numel(resultdir_list)

% load FLUXOS cs results
try
    res = importdata([resultdir_list{i},'/cs/f.out']);
    resultdir_legend = [resultdir_legend,resultdir_list{i}];
catch
    continue
end

% Results
time_mod = res(:,1)/3600; % sec -> hour
data_mod = res(:,2:end);


%subplot(211)
%time_surf = repmat(time_mod',numel(data_mod(1,:)),1)';
%crosec_surf = repmat((1:1:numel(data_mod(1,:)))',1,numel(time_mod))';
%surf(time_surf,crosec_surf,data_mod)
%axis tight
%shading interp
%alpha 0.8
%view(0,90)
%subplot(212)
plot(time_mod,(sum(data_mod'))','linewidth',2)
hold on

end
xlabel('Time [hour]')
ylabel('Flow [m3/s]')
axis tight
grid on
legend(['obs',resultdir_list],'interpreter','none')
title(yearselect)

