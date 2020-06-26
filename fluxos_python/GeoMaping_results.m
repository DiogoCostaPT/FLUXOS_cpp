
clear all

sites = {'LON','ESS','RIS2','FW4','CMCDC2'};

site_i = 1;

plot_options = {'dem_elevation','water_depth'};
plot_i = 1;

%% Input data

site = sites{site_i};

if (strcmp(site,'LON'))
    % LON
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_2/LON/ersi_grid_dem_LON_basin';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_2/LON/Results/151200.txt';
    delimiterIn = '\t';
    utm_zone = '17 T';
    ix = 1;
    iy = 2;
elseif (strcmp(site,'ESS'))
    % Essex
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_2/ESS/Essex_DEM_ascii_original';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_2/ESS/Results/151200.txt';
    delimiterIn = ' ';
    utm_zone = '17 T';
    ix = 2;
    iy = 1;
elseif (strcmp(site,'RIS2'))
    % Essex
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/RIS2/RIS2.asc';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/RIS2/Results/208800.txt';
    delimiterIn = ' ';
    utm_zone = '14 U';
    ix = 2;
    iy = 1;
elseif (strcmp(site,'FW4'))
    % Essex
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/FW4/FW4.asc';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/FW4/Results/259200.txt';
    delimiterIn = ' ';
    utm_zone = '14 U';
    ix = 2;
    iy = 1;
elseif (strcmp(site,'CMCDC2'))
    % Essex
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/CMCDC2/CMCDC2.asc';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/CMCDC2/Results/90000.txt';
    delimiterIn = ' ';
    utm_zone = '14 U';
    ix = 2;
    iy = 1;
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
utm_north = out(3);
utm_east = out(4);
gridsize = out(5);
NODATA_value = out(6);
dem(dem == NODATA_value) = NaN;
dem(dem == -NODATA_value) = NaN;

% if plotting water depth
if plot_i == 2

    results_raw = importdata(results_file);
    results_array = results_raw(:,1:4);
    row_num = numel(results_array(:,1));
    results_geoscatter = results_raw(:,[1:2,4]);

    for rowi = 1:row_num
        x = results_raw(rowi,ix);
        y = results_raw(rowi,iy);

        utm_east_i = utm_east + (x - 1)* gridsize;
        utm_north_i = utm_north + (y - 1)* gridsize;

        [lon, lat] = utm2deg(utm_east_i,utm_north_i,utm_zone);
        results_geoscatter(rowi,1) = lon;
        results_geoscatter(rowi,2) = lat;

    end
elseif plot_i == 1
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
alpha(gs1,1)
title(['Site: ',site])

