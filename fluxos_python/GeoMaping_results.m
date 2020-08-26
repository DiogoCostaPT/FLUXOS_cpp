
clear all

sites = {'LON','ESS','RIS2','FW4','CMCDC2','STC'};

site_i = 5;

plot_options = {'dem_elevation','water_depth'};
plot_i = 2;

cmax_waterdepth = 305; % in cm
add_sites = 1;
MarkerFaceAlpha = 0.05;

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
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_3/RIS2/Results/259200.txt';
    delimiterIn = ' ';
    utm_zone = '14 T';
    ix = 2; % for simulations
    iy = 1; % for simulations
    x_inverse = true;
    row_i_north_utm = 4; % for DEM
    row_i_east_utm = 3; % for DEM
    sitesshapefile_up = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/Janina_Henry_Merrin/RIS2_up.shp';
     sitesshapefile_down = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/Janina_Henry_Merrin/RIS2_down.shp';
elseif (strcmp(site,'FW4'))
    % Essex
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_4/FW4/FW4_extended_2.asc';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_4/FW4/Results/259200.txt';
    delimiterIn = ' ';
    utm_zone = '14 U';
    ix = 2; % for simulations
    iy = 1; % for simulations
    x_inverse = true;
    row_i_north_utm = 4; % for DEM
    row_i_east_utm = 3; % for DEM
    sitesshapefile_up = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/Janina_Henry_Merrin/FW4_up.shp';
    sitesshapefile_down = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/Janina_Henry_Merrin/FW4_down.shp';
elseif (strcmp(site,'CMCDC2'))
    % Essex
    demfile = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_4/CMCDC2/CMCDC2_extended_2.asc';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/Henry_Janina_4/CMCDC2/Results/259200.txt';
    delimiterIn = ' ';
    utm_zone = '14 U';
    ix = 2; % for simulations
    iy = 1; % for simulations
    x_inverse = true;
    row_i_north_utm = 4; % for DEM
    row_i_east_utm = 3; % for DEM
    sitesshapefile_up = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/Janina_Henry_Merrin/CA_up.shp';
    sitesshapefile_down = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/Janina_Henry_Merrin/CA_down.shp';
elseif (strcmp(site,'STC'))
    demfile = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/STC_data_pre-processing/DEM_ASCII/steppler_dem_ascii_3m.txt';
    results_file = '/media/dcosta/data/megasync/my_server/fluxos/batch_1_select_paper/t_65_paper/Results/723600.txt';
    delimiterIn = ' ';
    utm_zone = '14 U';
    ix = 2; % for simulations
    iy = 1; % for simulations
    x_inverse = false;
    row_i_north_utm = 4; % for DEM
    row_i_east_utm = 3; % for DEM
    %sitesshapefile_up = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/Janina_Henry_Merrin/CA_up.shp';
    %sitesshapefile_down = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/Janina_Henry_Merrin/CA_down.shp';
end


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
    results_array = results_raw(:,1:end);
    row_num = numel(results_array(:,1));
    results_geoscatter = results_raw(:,[1:2,end]);
    small_loc = find(results_geoscatter(:,end)< 1);
    results_geoscatter(small_loc,:) = NaN;
    
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
alpha(gs1,MarkerFaceAlpha)
title(['Site: ',site])

% adding sites
if add_sites
    sitesshapefile_all = {sitesshapefile_up,sitesshapefile_down};
    styleform = {'wo','ws'};
   for f = 1:2
       sitesshapefile_i = sitesshapefile_all{f};
       site_data = shaperead(sitesshapefile_i);
       sites_XY = [[site_data.X]',[site_data.Y]'];
       sites_lonlat = zeros(size(sites_XY));
       for i = 1:numel(sites_lonlat(:,1))
           [lon, lat] = utm2deg(sites_XY(i,1),sites_XY(i,2),utm_zone);
           sites_lonlat(i,1) = lon;
           sites_lonlat(i,2) = lat;
       end
       hold on
       geoscatter(sites_lonlat(:,1),sites_lonlat(:,2),70,styleform{f},'linewidth',1.5)
   end
end

% To delete
KVA = [0579921,5512270];
KVB = [0579891,5513032];

STM = [43+19/60,-81+07/60];
lon = STM(1);
lat = STM(2);

ILD = [43+30/60,-74+23/60];
lon = ILD(1);
lat = ILD(2);

[lon, lat] = utm2deg(KVB(1),KVB(2),'14 U');
hold on
geoscatter(lon,lat,70,styleform{f},'linewidth',1.5)
