
clear all

sites = {'LON','ESS','RIS2','FW4','CMCDC2'};

site_i = 4;

plot_options = {'dem_elevation','water_depth'};
plot_i = 2;

%% Input data

site = sites{site_i};

if (strcmp(site,'LON'))
    % LON
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_2/LON/ersi_grid_dem_LON_basin';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_2/LON/Results/151200.txt';
    delimiterIn = '\t';
    utm_zone = '17 T';
    ix = 1; % for simulations
    iy = 2; % for simulations
    x_inverse = false;
    row_i_north_utm = 3; % for DEM
    row_i_east_utm = 4; % for DEM
elseif (strcmp(site,'ESS'))
    % Essex
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_2/ESS/Essex_DEM_ascii_original';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_2/ESS/Results/151200.txt';
    delimiterIn = ' ';
    utm_zone = '17 T';
    ix = 2; % for simulations
    iy = 1; % for simulations
    x_inverse = false;
    row_i_north_utm = 3; % for DEM
    row_i_east_utm = 4; % for DEM
elseif (strcmp(site,'RIS2'))
    % Essex
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/RIS2/RIS2.asc';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/RIS2/Results/208800.txt';
    delimiterIn = ' ';
    utm_zone = '14 T';
    ix = 2; % for simulations
    iy = 1; % for simulations
    x_inverse = true;
    row_i_north_utm = 4; % for DEM
    row_i_east_utm = 3; % for DEM
elseif (strcmp(site,'FW4'))
    % Essex
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/FW4/FW4.asc';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/FW4/Results/259200.txt';
    delimiterIn = ' ';
    utm_zone = '14 U';
    ix = 2; % for simulations
    iy = 1; % for simulations
    x_inverse = true;
    row_i_north_utm = 4; % for DEM
    row_i_east_utm = 3; % for DEM
elseif (strcmp(site,'CMCDC2'))
    % Essex
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/CMCDC2/CMCDC2.asc';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/CMCDC2/Results/90000.txt';
    delimiterIn = ' ';
    utm_zone = '14 U';
    ix = 2; % for simulations
    iy = 1; % for simulations
    x_inverse = true;
    row_i_north_utm = 4; % for DEM
    row_i_east_utm = 3; % for DEM
end

cmax_waterdepth = 10; % in cm

%% Manipulations
headerlinesIn = 6; % in dem

% DEM
[dem_raw,delimiterOut,headerlinesOut] = importdata(demfile,delimiterIn,headerlinesIn);

dem = dem_raw.data;
[nx,ny] = size(dem);
dem_info = dem_raw.textdata;
b=regexp(dem_info,'\d+(\.)?(\d+)?','match');
out=str2double([b{:}]);
nrows = out(1);
ncols = out(2);
utm_north = out(row_i_north_utm);
utm_east = out(row_i_east_utm);
gridsize = out(5);
NODATA_value = out(6);
dem(dem == NODATA_value) = NaN;
dem(dem == -NODATA_value) = NaN;

% for DEM
if plot_i == 1
    dem_resh = reshape(dem,numel(dem),1);
    dem_array = zeros(numel(dem_resh),3)*NaN;
    dem_array(:,3) = dem_resh;
    for rowi = 1:nrows
        for coli = 1:ncols
            linei = (rowi-1)*ncols + coli;
            ncols_r = ncols - coli + 1;
            utm_east_i = utm_east + (rowi - 1)* gridsize;
            utm_north_i = utm_north + (ncols_r - 1)* gridsize;
            [lon, lat] = utm2deg(utm_east_i,utm_north_i,utm_zone);
            
            dem_array(linei,1) = lon;
            dem_array(linei,2) = lat;
        end
    end
% for water depth
elseif plot_i == 2

    results_raw = importdata(results_file);
    results_array = results_raw(:,1:4);
    row_num = numel(results_array(:,1));
    results_geoscatter = results_raw(:,[1:2,4]);
    

    for rowi = 1:row_num
        x = results_raw(rowi,ix);
        y = results_raw(rowi,iy);
        
        if (x_inverse); y = ncols - y; end

        utm_east_i = utm_east + (x - 1)* gridsize;
        utm_north_i = utm_north + (y+1 - 1)* gridsize;

        [lon, lat] = utm2deg(utm_east_i,utm_north_i,utm_zone);
        results_geoscatter(rowi,1) = lon;
        results_geoscatter(rowi,2) = lat;

    end
end

%% Plotting
figure('name',site)
if plot_i == 1
    gs1 = geoscatter(dem_array(:,1),dem_array(:,2),gridsize*1.5,dem_array(:,3));
    hcb = colorbar;
    colormap(othercolor('BrBG10'))
    set(get(hcb,'Title'),'String','Elevation [m]')   
elseif plot_i == 2
    gs1 = geoscatter(results_geoscatter(:,1),results_geoscatter(:,2),gridsize,results_geoscatter(:,3)*100)   
    caxis([0 cmax_waterdepth])
    hcb = colorbar;
    colormap('jet')
    set(get(hcb,'Title'),'String','Maximum water depth [cm]')   
end
geobasemap satellite
alpha(gs1,0.5)
title(['Site: ',site])

