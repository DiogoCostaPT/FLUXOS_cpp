
dirfolder = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/FLUXOS_cpp/'

dem_denoise = importdata([dirfolder, 'ersi_grid_dem_alldomain_denoised']);

dem_basin_walls = importdata([dirfolder, 'ersi_grid_dem_alldomain_basinwalls']);



dem_denoise_walls = dem_denoise;
dem_denoise_walls(dem_basin_walls==99999)=99999;

figure
surf(dem_denoise)
hold on
surf(dem_denoise_walls)

%%%%%%%%%%%%%%

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

%%%%%%%%% delete Results folder
files = dir(folder_batch2);
filenam = {files.name};
filenam = filenam(3:end)';

for i=1:numel(filenam)
    sim_i = filenam{i}
    path = [folder_batch2,sim_i,'/'];
    path_results = [path,'Results/'];
    try
        rmdir(path_results)
        continue
    catch
        disp(['problem at: ',path_results])
    end
end

% create results folder and new 0.txt file
file_0_path = '/media/dcosta/DATADRIVE1/fluxos_tests/SIMULATIONS_sync/0.txt';
file_0_data = importdata(file_0_path);

files = dir(folder_batch2);
filenam = {files.name};
filenam = filenam(3:end)';

for i=1:numel(filenam)
    sim_i = filenam{i}
    path = [folder_batch2,sim_i,'/'];
    path_results = [path,'Results/'];
    save('MyMatrix.txt', 'A', '-double', ',')
end