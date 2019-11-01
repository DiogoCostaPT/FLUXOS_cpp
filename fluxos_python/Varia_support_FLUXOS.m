

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete Results folder inside simulation folders and add 0.txt if desired
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simfolder_dir_path = '/media/dcosta/DATADRIVE1/fluxos_tests/SIMULATIONS_sync/batch_2/';

% Tool 1: remove the Results files 
simfolder_names_cell = dir(simfolder_dir_path);
sim_folder_names = {simfolder_names_cell.name};
sim_folder_names = sim_folder_names(3:end)';

for i=1:numel(sim_folder_names)
    sim_i = sim_folder_names{i}
    path = [simfolder_dir_path,sim_i,'/'];
    path_results = [path,'Results'];
    try
        rmdir(path_results,'s')
        continue
    catch
        disp(['problem at: ',path_results])
    end
end

% Tool 2: create Results folder
simfolder_names_cell = dir(simfolder_dir_path);
sim_folder_names = {simfolder_names_cell.name};
sim_folder_names = sim_folder_names(3:end)';

for i=1:numel(sim_folder_names)
    sim_i = sim_folder_names{i}
    path_results = [simfolder_dir_path,sim_i,'/Results'];
    try
        mkdir(path_results)
        continue
    catch
        disp(['problem at: ',path_results])
    end
end

% Tool 3: copy 0.txt file from Seagate hardisk (backup) - looks for same simulation
% names and copies the respective 0.txt files

source_folder = '/media/dcosta/Seagate Backup Plus Drive/FLUXOS/SIMULATIONS_sync/batch_2_plato_incomplete/';

simfolder_names_cell = dir(source_folder);
filenam_src = {simfolder_names_cell.name};
filenam_src = sim_folder_names(3:end)';
%files = dir(dest_folder);
%filenam_dst = {files.name};
%filenam_dst = filenam(3:end)';

for i=1:numel(filenam_src)
    sim_i = filenam_src{i}
    path_0txt_file_src = [source_folder,sim_i,'/Results/0.txt'];
    path_dst_Resfolder = [simfolder_dir_path,sim_i,'/Results'];
    try
        copyfile(path_0txt_file_src,path_dst_Resfolder,'f');
        continue
    catch
        disp(['problem at: ',path_results])
    end
end


% Tool 5: copy fluxo_cpp executable (or any other file) from one location to all simulation folders
fluxos_cpp_folder_path = '/media/dcosta/DATADRIVE1/fluxos_tests/SIMULATIONS_sync/0_fluxos_graham/build/fluxos_cpp';
%fluxos_cpp_folder_path = '/media/dcosta/DATADRIVE1/fluxos_tests/SIMULATIONS_sync/job.sh';

simfolder_names_cell = dir(simfolder_dir_path);
sim_folder_names = {simfolder_names_cell.name};
sim_folder_names = sim_folder_names(3:end)';

for i=1:numel(sim_folder_names)
    sim_i = sim_folder_names{i}
    path_results = [simfolder_dir_path,sim_i,'/'];
    try
        copyfile(fluxos_cpp_folder_path,path_results,'f');
        continue
    catch
        disp(['problem at: ',path_results])
    end
end

% Tool 6: delete certain files in all simulation folders if needed

%files_2_delete = {'*.out','fluxos_run.log'};
files_2_delete = {'fluxos_cpp'};

simfolder_names_cell = dir(simfolder_dir_path);
sim_folder_names = {simfolder_names_cell.name};
sim_folder_names = sim_folder_names(3:end)';

for f=1:numel(files_2_delete)
    for i=1:numel(sim_folder_names)
        sim_i = sim_folder_names{i}
        path_file2del = [simfolder_dir_path,sim_i,'/',files_2_delete{f}];
        try
            delete(path_file2del)
            continue
        catch
            disp(['problem at: ',path_results, '(file or folder dont exist'])
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEM manipulations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirfolder = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/FLUXOS_cpp/'

dem_denoise = importdata([dirfolder, 'ersi_grid_dem_alldomain_denoised']);

dem_basin_walls = importdata([dirfolder, 'ersi_grid_dem_alldomain_basinwalls']);


dem_denoise_walls = dem_denoise;
dem_denoise_walls(dem_basin_walls==99999)=99999;

figure
surf(dem_denoise)
hold on
surf(dem_denoise_walls)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OTHERS - 2 delete %%%%%%%%%%%

folder_batch2 = '/media/dcosta/DATADRIVE1/fluxos_tests/SIMULATIONS_sync/batch_2/';
year_i = 2011

if year_i==2009
    year_file = 67:88;
elseif  year_i==2010
    year_file = 93:114;
elseif  year_i==2011
    year_file = 119:140;
end

for i = 1:numel(year_file)
    modset_file = [folder_batch2,'t_',mat2str(year_file(i)),'/modset.fluxos']
    mod_file = importdata(modset_file);
    
    fid = fopen(modset_file, 'w');
    fprintf(fid,'%s\n', mod_file.textdata{1});
    fprintf(fid,'%s\n', mod_file.textdata{2});
    fprintf(fid,'%s\n', mod_file.textdata{3});
    %fprintf(fid,'%s\n', ['Qmelt_info_',mat2str(year_i),'_T-index.fluxos']);
    fprintf(fid,'%s\n', ['Qmelt_info_',mat2str(year_i),'.fluxos']);
    fprintf(fid,'%s\n', mat2str(mod_file.data(1)));
    fprintf(fid,'%s\n', mat2str(mod_file.data(2)));
    fprintf(fid,'%s\n', mat2str(mod_file.data(3)));
    fprintf(fid,'%s\n', mat2str(mod_file.data(4)));
    fprintf(fid,'%s\n', mat2str(mod_file.data(5)));
    fprintf(fid,'%s\n', mat2str(mod_file.data(6)));
    fprintf(fid,'%s\n', mat2str(mod_file.data(7)));
end
fclose(fid);

timeinit = 40252.03125 + 695422 - 8/24;
timepeak = 7.356767812500000e+05;

t_diff_secs = etime(datevec(timepeak),datevec(timeinit));


