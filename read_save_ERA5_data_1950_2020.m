close all
clear all

addpath(genpath('/Users/patrickbrown/Documents/MATLAB/'))

%% readme

% This folder contains two classes of files: 
% 
% 1) files containing variables that cover the land area of the WECC, and 
% 2) files containing variables that cover all longitudes between the latitudes of 15 N and 65 N ? I am happy to expand these if you want to include the tropics.
% 3) files containing variables that cover all longitudes between the latitudes of 20 S and 60 N ? I am happy to expand these
% 
% 1) Variables to define the events of interest over the WECC:
% These variables are at 0.5 by 0.5 degrees and hourly resolutions.
% These are already subsetted to include only land cells. 
% (a) WECC_WP.nc ? estimated wind power (capacity factor derived from 100m wind speed data)
% (b) WECC_ssrd.nc ? surface solar radiation downwards
% (c) WECC_2m_temp.nc ? 2m surface temperature
% (a) WECC_WP_W_per_sq_m.nc ? estimated wind power density -- Watts per square meter derived from 100m wind speed data and a 90 m turbine rotor diameter and assumed spacing of 4 times rotor diameter in one direction and 7 times rotor diameter in second direction. IN this case, it means that each 2 MW turbine occupies 4 * 90 m * 7 * 90 m = 226,800 m^2.
% 
% 
% The file lat_lon_index_key.csv will help you map the ?lat_lon_index? variable in the .nc files ? I had to use this mapping because the grid is not rectangular (remember I am only saving land cells to save space).
% 
% 
% 2) Large scale variables to use for event diagnostics:
% These variables are at 2.5 by 2.5 degrees and daily resolution.
% (a) NH_V_x_HPA.nc ? v wind component at x pressure level
% (b) NH_U_x_HPA.nc ? u wind component at x pressure level
% (c) NH_GEO_x_HPA.nc ? geopotential at x pressure level
% (d) NH_TEMP_850_hPA.nc ? temperature at 850 hPa
% (e) NH_VOR_500_hPa.nc ? vorticity at 500 hPa
% (f) NH_MSL.nc ? mean sea level pressure
% 
% 
% 3) Large scale variables to use for event diagnostics:
% These variables are at 0.25 by 0.25 degrees and monthly resolution. (The spatial resolution is the default for downloading ERA5)
% (a) GLOBAL_SST_monthly ? SSTs
% 
% 
% All data was originally downloaded using an API from:
% https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form and 
% https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=form
% 
% Please reach out to me at dfarnham@carnegiescience.edu if anything is unclear, if you think there are mistakes, if you have other questions or if you would like to me download other data with different spatial/temporal domains and/or resolutions.
%         

%% load files

filepath = strcat('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/David_ERA5_data/');

NH_lat_era5 = ncread(strcat(filepath,'NH_GEO_250_hPa_1950_2020','.nc'),'lat');
NH_lon_era5 = ncread(strcat(filepath,'NH_GEO_250_hPa_1950_2020','.nc'),'lon');
%NH_time_era5 = ncread(strcat(filepath,'NH_GEO_250_hPa','.nc'),'time');

% descriptive time array

t1 = datetime(1950,1,1,12,0,0);
t2 = datetime(2020,12,31,12,0,0);
datetime_desc_time_daily = t1:t2;

% gridded

NH_GEO_250_hPa = ncread(strcat(filepath,'NH_GEO_250_hPa_1950_2020','.nc'),'z_250');
NH_GEO_500_hPa = ncread(strcat(filepath,'NH_GEO_500_hPa_1950_2020','.nc'),'z_500');
NH_GEO_700_hPa = ncread(strcat(filepath,'NH_GEO_700_hPa_1950_2020','.nc'),'z_700');
%NH_GEO_850_hPa = ncread(strcat(filepath,'NH_GEO_850_hPa_1950_2020','.nc'),'z_850');

%NH_TEMP_850_hPa = ncread(strcat(filepath,'NH_TEMP_850_hPa_1950_2020','.nc'),'t_850');
NH_TEMP_700_hPa = ncread(strcat(filepath,'NH_TEMP_700_hPa_1950_2020','.nc'),'t_700');

NH_U_250_hPa = ncread(strcat(filepath,'NH_U_250_hPa_1950_2020','.nc'),'u_250');
NH_U_500_hPa = ncread(strcat(filepath,'NH_U_500_hPa_1950_2020','.nc'),'u_500');
NH_U_700_hPa = ncread(strcat(filepath,'NH_U_700_hPa_1950_2020','.nc'),'u_700');
%NH_U_850_hPa = ncread(strcat(filepath,'NH_U_850_hPa_1950_2020','.nc'),'u_850');

NH_V_250_hPa = ncread(strcat(filepath,'NH_V_250_hPa_1950_2020','.nc'),'v_250');
NH_V_500_hPa = ncread(strcat(filepath,'NH_V_500_hPa_1950_2020','.nc'),'v_500');
NH_V_700_hPa = ncread(strcat(filepath,'NH_V_700_hPa_1950_2020','.nc'),'v_700');
%NH_V_850_hPa = ncread(strcat(filepath,'NH_V_850_hPa_1950_2020','.nc'),'v_850');

NH_VOR_500_hPa = ncread(strcat(filepath,'NH_VOR_500_hPa_1950_2020','.nc'),'vo_500');
NH_SLP = ncread(strcat(filepath,'NH_MSL_1950_2020','.nc'),'msl');
NH_VV_700 = ncread(strcat(filepath,'NH_VV_700_hPa_1950_2020','.nc'),'w_700');

all_gridded_NH_data = cat(4,NH_GEO_250_hPa,...
                            NH_GEO_500_hPa,...
                            NH_GEO_700_hPa,...
                            NH_U_250_hPa,...
                            NH_U_500_hPa,...
                            NH_U_700_hPa,...
                            NH_V_250_hPa,...
                            NH_V_500_hPa,...
                            NH_V_700_hPa,...
                            NH_VOR_500_hPa,...
                            NH_TEMP_700_hPa,...
                            NH_VV_700,...
                            NH_SLP);
                        
all_gridded_NH_data_varnames = {'NH_GEO_250_hPa',...
                                'NH_GEO_500_hPa',...
                                'NH_GEO_700_hPa',...
                                'NH_U_250_hPa',...
                                'NH_U_500_hPa',...
                                'NH_U_700_hPa',...
                                'NH_V_250_hPa',...
                                'NH_V_500_hPa',...
                                'NH_V_700_hPa',...
                                'NH_VOR_500_hPa',...
                                'NH_TEMP_700_hPa',...
                                'NH_VV_700',...
                                'NH_SLP'};                            
clear NH_GEO_250_hPa
clear NH_GEO_500_hPa
clear NH_GEO_700_hPa
clear NH_U_250_hPa
clear NH_U_500_hPa
clear NH_U_700_hPa
clear NH_V_250_hPa
clear NH_V_500_hPa
clear NH_V_700_hPa
clear NH_VOR_500_hPa
clear NH_TEMP_700_hPa
clear NH_VV_700
clear NH_SLP
                            
%% WECC data

WECC_lat_lon_index_key = readmatrix(strcat(filepath,'WECC_lat_lon_index_key','.csv'));

WECC_2m_temp_daily = ncread(strcat(filepath,'WECC_2m_temp_daily_1950_2020','.nc'),'tmp');
WECC_ssrd_daily = ncread(strcat(filepath,'WECC_ssrd_adj_daily_1950_2020','.nc'),'ssrd_adj');
WECC_WP_daily = ncread(strcat(filepath,'WECC_WP_W_per_sq_m_daily_1950_2020','.nc'),'WP');
                    
WECC_daily_var_names = {'WECC_2m_temp_daily',...
                        'WECC_ssrd_daily',...
                        'WECC_WP_daily'};
                    
WECC_2m_all_var_daily = cat(3,WECC_2m_temp_daily,...
                          WECC_ssrd_daily,...
                          WECC_WP_daily);

% WECC_2m_temp_hourly = ncread(strcat(filepath,'WECC_2m_temp','.nc'),'tmp');
% WECC_ssrd_hourly = ncread(strcat(filepath,'WECC_ssrd','.nc'),'ssrd');
% WECC_WP_hourly = ncread(strcat(filepath,'WECC_WP','.nc'),'WP_CF');

WECC_lat_min = min(WECC_lat_lon_index_key(:,2));
WECC_lat_max = max(WECC_lat_lon_index_key(:,2));
WECC_lon_min = min(WECC_lat_lon_index_key(:,1));
WECC_lon_max = max(WECC_lat_lon_index_key(:,1));

WECC_lat = WECC_lat_min:WECC_lat_max;
WECC_lon = WECC_lon_min:WECC_lon_max;

gridded_WECC_all_var_daily = NaN(length(WECC_lat),length(WECC_lon),length(datetime_desc_time_daily),length(WECC_daily_var_names));

for var_i = 1:length(WECC_daily_var_names)
    for day_i = 1:length(datetime_desc_time_daily)
        for lin_loc_i = 1:size(WECC_2m_temp_daily,1)

            lin_lat_now = WECC_lat_lon_index_key(lin_loc_i,2);
            lin_lon_now = WECC_lat_lon_index_key(lin_loc_i,1);

            gridded_lat_ind = find(WECC_lat == lin_lat_now);
            gridded_lon_ind = find(WECC_lon == lin_lon_now);

            gridded_WECC_all_var_daily(gridded_lat_ind,gridded_lon_ind,day_i,var_i) = WECC_2m_all_var_daily(lin_loc_i,day_i,var_i);

        end
    end
end

%make 1d time series

%t_series_WECC_all_var_daily
    
%% save

save('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/David_ERA5_data/read_save_ERA5_data_1950_2020.mat',...
    'all_gridded_NH_data',...
    'gridded_WECC_all_var_daily',...
    'NH_lat_era5',...
    'NH_lon_era5',...
    'all_gridded_NH_data_varnames',...
    'WECC_daily_var_names',...
    'WECC_lat',...
    'WECC_lon',...
    'datetime_desc_time_daily',...
    '-v7.3')





%% load files

%%

% filepath = strcat('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/clim_expl_data/');
% 
% %lat/lon
%     global_lat_era5 = ncread(strcat(filepath,'era5_t2m_daily_2_2','.nc'),'lat');
%     global_lon_era5 = ncread(strcat(filepath,'era5_t2m_daily_2_2','.nc'),'lon');
%     global_lat_erai = ncread(strcat(filepath,'erai_v850_daily_2_2','.nc'),'lat');
%     global_lon_erai = ncread(strcat(filepath,'erai_v850_daily_2_2','.nc'),'lon');
% 
% %big global files
%     %ERA5
%         era5_t2m_daily_2_2 = ncread(strcat(filepath,'era5_t2m_daily_2_2','.nc'),'t2m');
%         era5_prcp_daily_2_2 = ncread(strcat(filepath,'era5_prcp_daily_2_2','.nc'),'tp');
%         era5_z500_daily_2_2 = ncread(strcat(filepath,'era5_z500_daily_2_2','.nc'),'z500');
%         era5_slp_daily_2_2 = ncread(strcat(filepath,'era5_slp_daily_2_2','.nc'),'msl');
%     %ERA-I
%         erai_u850_daily_2_2 = ncread(strcat(filepath,'erai_u850_daily_2_2','.nc'),'u850');
%         erai_v850_daily_2_2 = ncread(strcat(filepath,'erai_v850_daily_2_2','.nc'),'v850');
%         
% %% interpolate ERA-I to ERA5
% 
% [era5_lat_grid,era5_lon_grid] = meshgrid(global_lat_era5,global_lon_era5);
% 
% era5_lat_grid = double(era5_lat_grid);
% era5_lon_grid = double(era5_lon_grid);
% global_lat_era5 = double(global_lat_era5);
% global_lon_era5 = double(global_lon_era5);
% 
% erai_u850_daily_2_2_regrid = NaN(length(global_lon_era5),length(global_lat_era5),size(erai_u850_daily_2_2,3));
% for time_i = 1:size(erai_u850_daily_2_2,3)
%     erai_u850_daily_2_2_regrid(:,:,time_i) = interp2(global_lat_erai,global_lon_erai,erai_u850_daily_2_2(:,:,time_i),era5_lat_grid,era5_lon_grid);
% end
% 
% erai_v850_daily_2_2_regrid = NaN(length(global_lon_era5),length(global_lat_era5),size(erai_v850_daily_2_2,3));
% for time_i = 1:size(erai_u850_daily_2_2,3)
%     erai_v850_daily_2_2_regrid(:,:,time_i) = interp2(global_lat_erai,global_lon_erai,erai_v850_daily_2_2(:,:,time_i),era5_lat_grid,era5_lon_grid);
% end
% 
% %cut off SH
% good_lat_inds = find(global_lat_era5 >= -10);
% 
% era5_t2m_daily_2_2 = era5_t2m_daily_2_2(:,good_lat_inds,:);
% era5_prcp_daily_2_2 = era5_prcp_daily_2_2(:,good_lat_inds,:);
% era5_z500_daily_2_2 = era5_z500_daily_2_2(:,good_lat_inds,:);
% era5_slp_daily_2_2 = era5_slp_daily_2_2(:,good_lat_inds,:);
% erai_u850_daily_2_2_regrid = erai_u850_daily_2_2_regrid(:,good_lat_inds,:);
% erai_v850_daily_2_2_regrid = erai_v850_daily_2_2_regrid(:,good_lat_inds,:);
% 
% NH_lat_era5 = global_lat_era5(good_lat_inds);
% 
% all_gridded_global_data = NaN(length(global_lon_era5),length(NH_lat_era5),size(era5_t2m_daily_2_2,3),6);
% 
% all_gridded_global_data(:,:,1:size(era5_t2m_daily_2_2,3)) = era5_t2m_daily_2_2;
% all_gridded_global_data(:,:,1:size(era5_prcp_daily_2_2,3)) = era5_prcp_daily_2_2;
% all_gridded_global_data(:,:,1:size(era5_z500_daily_2_2,3)) = era5_z500_daily_2_2;
% all_gridded_global_data(:,:,1:size(era5_slp_daily_2_2,3)) = era5_slp_daily_2_2;
% all_gridded_global_data(:,:,1:size(erai_u850_daily_2_2,3)) = erai_u850_daily_2_2_regrid;
% all_gridded_global_data(:,:,1:size(erai_v850_daily_2_2,3)) = erai_v850_daily_2_2_regrid;
% 
% var_names_all_gridded_global_data = {'t2m',...
%                                      'tp',...
%                                      'z500',...
%                                      'msl',...
%                                      'u850',...
%                                      'v850'};
%                                 
% clear era5_t2m_daily_2_2
% clear era5_prcp_daily_2_2
% clear era5_z500_daily_2_2
% clear era5_slp_daily_2_2
% clear erai_u850_daily_2_2
% clear erai_v850_daily_2_2
% 
% %% read the domain-gridded variables
% 
% % era5_t2m_daily_220-255E_30-60N
% % era5_sfcWind_daily_220-255E_30-60N
% % erai_rsds_daily_220-255E_30-60N
% 
% %lat/lon
%     domain_lat_era5 = ncread(strcat(filepath,'era5_t2m_daily_220-255E_30-60N','.nc'),'lat');
%     domain_lon_era5 = ncread(strcat(filepath,'era5_t2m_daily_220-255E_30-60N','.nc'),'lon');
%     domain_lat_erai = ncread(strcat(filepath,'erai_rsds_daily_220-255E_30-60N','.nc'),'lat');
%     domain_lon_erai = ncread(strcat(filepath,'erai_rsds_daily_220-255E_30-60N','.nc'),'lon');
% 
% %big global files
%     %ERA5
%         era5_t2m_daily_220_255E_30_60N = ncread(strcat(filepath,'era5_t2m_daily_220-255E_30-60N','.nc'),'t2m');
%         era5_sfcWind_daily_220_255E_30_60N = ncread(strcat(filepath,'era5_sfcWind_daily_220-255E_30-60N','.nc'),'sfcWind');
%     %ERA-I
%         erai_rsds_daily_220_255E_30_60N = ncread(strcat(filepath,'erai_rsds_daily_220-255E_30-60N','.nc'),'rsds');
%         
% % interpolate ERA-I to ERA5
% 
% [domain_era5_lat_grid,domain_era5_lon_grid] = meshgrid(domain_lat_era5,domain_lon_era5);
% 
% domain_era5_lat_grid = double(domain_era5_lat_grid);
% domain_era5_lon_grid = double(domain_era5_lon_grid);
% domain_lat_era5 = double(domain_lat_era5);
% domain_lon_era5 = double(domain_lon_era5);
% 
% erai_rsds_daily_220_255E_30_60N_regrid = NaN(length(domain_lon_era5),length(domain_lat_era5),size(erai_rsds_daily_220_255E_30_60N,3));
% for time_i = 1:size(erai_rsds_daily_220_255E_30_60N,3)
%     erai_rsds_daily_220_255E_30_60N_regrid(:,:,time_i) = interp2(domain_lat_erai,domain_lon_erai,erai_rsds_daily_220_255E_30_60N(:,:,time_i),domain_era5_lat_grid,domain_era5_lon_grid);
% end
% 
% all_gridded_domain_data = NaN(length(domain_lon_era5),length(domain_lat_era5),size(era5_t2m_daily_220_255E_30_60N,3),6);
% 
% all_gridded_domain_data(:,:,1:size(era5_t2m_daily_220_255E_30_60N,3)) = era5_t2m_daily_220_255E_30_60N;
% all_gridded_domain_data(:,:,1:size(era5_sfcWind_daily_220_255E_30_60N,3)) = era5_sfcWind_daily_220_255E_30_60N;
% all_gridded_domain_data(:,:,1:size(erai_rsds_daily_220_255E_30_60N_regrid,3)) = erai_rsds_daily_220_255E_30_60N_regrid;
% 
% var_names_all_gridded_domain_data = {'t2m',...
%                                      'sfcWind',...
%                                      'rsds'};
%                                 
% clear era5_t2m_daily_220_255E_30_60N
% clear era5_sfcWind_daily_220_255E_30_60N
% clear erai_rsds_daily_220_255E_30_60N
% 
% 
% 
% %% load netcdf files
% 
% % 'iera5_t2m_daily_220-255E_30-60N_n_5lan__yr'
% % 'iera5_t2m_daily_220-255E_30-60N_n_5lan'
% % 'iera5_wspd_daily_220-255E_30-60N_n_5lan__yr'
% % 'iera5_wspd_daily_220-255E_30-60N_n_5lan'
% % 'ierai_rsds_daily_220-255E_30-60N_n_5lan__yr'
% % 'ierai_rsds_daily_220-255E_30-60N_n_5lan'
% 
% %load the data with the seasonal cycle
% domain_mean_iera5_t2m_daily_220_255E_30_60N = ncread(strcat(filepath,'iera5_t2m_daily_220-255E_30-60N_n_5lan','.nc'),'t2m');
% domain_mean_iera5_sfcWind_daily_220_255E_30_60N = ncread(strcat(filepath,'iera5_wspd_daily_220-255E_30-60N_n_5lan','.nc'),'sfcWind');
% domain_mean_ierai_rsds_daily_220_255E_30_60N = ncread(strcat(filepath,'ierai_rsds_daily_220-255E_30-60N_n_5lan','.nc'),'rsds');
% 
% domain_mean_with_seas = cat(2,domain_mean_iera5_t2m_daily_220_255E_30_60N,...
%                               domain_mean_iera5_sfcWind_daily_220_255E_30_60N,...    
%                               domain_mean_ierai_rsds_daily_220_255E_30_60N);
% 
% %load the data without the seasonal cycle
% domain_mean_deseas_iera5_t2m_daily_220_255E_30_60N = ncread(strcat(filepath,'iera5_t2m_daily_220-255E_30-60N_n_5lan_a.txt','.nc'),'t2m');
% domain_mean_deseas_iera5_sfcWind_daily_220_255E_30_60N = ncread(strcat(filepath,'iera5_wspd_daily_220-255E_30-60N_n_5lan_a.txt','.nc'),'sfcWind');
% domain_mean_deseas_ierai_rsds_daily_220_255E_30_60N = ncread(strcat(filepath,'ierai_rsds_daily_220-255E_30-60N_n_5lan_a.txt','.nc'),'rsds');
% 
% domain_mean_without_seas = cat(2,domain_mean_deseas_iera5_t2m_daily_220_255E_30_60N,...
%                                  domain_mean_deseas_iera5_sfcWind_daily_220_255E_30_60N,...    
%                                  domain_mean_deseas_ierai_rsds_daily_220_255E_30_60N);
%                              
% domain_means_w_and_wo_seas = cat(3,domain_mean_with_seas,domain_mean_without_seas);
%                              
% var_names_domain_means = {'t2m',...
%                           'sfcWind',...
%                           'rsds'};
% 
% %load the seasonal ycle
% 
% domain_mean_seascyc_iera5_t2m_daily_220_255E_30_60N = readmatrix(strcat(filepath,'iera5_t2m_daily_220-255E_30-60N_n_5lan__yr_clean','.txt'));
% domain_mean_seascyc_iera5_sfcWind_daily_220_255E_30_60N = readmatrix(strcat(filepath,'iera5_wspd_daily_220-255E_30-60N_n_5lan__yr_clean','.txt'));
% domain_mean_seascyc_ierai_rsds_daily_220_255E_30_60N = readmatrix(strcat(filepath,'ierai_rsds_daily_220-255E_30-60N_n_5lan__yr_clean','.txt'));
% 
% domain_means_annual_seascyc = cat(3,domain_mean_seascyc_iera5_t2m_daily_220_255E_30_60N(:,1:7),...
%                                     domain_mean_seascyc_iera5_sfcWind_daily_220_255E_30_60N(:,1:7),...
%                                     domain_mean_seascyc_ierai_rsds_daily_220_255E_30_60N(:,1:7));
% 
% annual_var_names = {'day of year',...
%                     'mean',...
%                     '2.5%',...
%                     '17%',...
%                     '50%',...
%                     '83%',...
%                     '97.5%'};
%                 
% %% save
% 
% %calculate time
% 
% t1 = datetime(1979,1,1,12,0,0);
% t2 = datetime(2019,7,31,12,0,0);
% datetime_desc_time_daily = t1:t2;
% 
% %it says that 01Jan1979 to 31Jul2019 is 14853 
% % July 31st is the 212th day of the year
% desc_time_daily = linspace(1979.0,2019+212/365,14853);
% 
% % But datetime says that there are 14822 days between 01Jan1979 and
% % 31Jul2019 - must figure out this descrep
% 
% %14730 is the shortest time series of all the data so cut everythign to
% %that discrepancy
% 
% save('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/clim_expl_data/read_save_clim_expl_data.mat',...
%     'all_gridded_global_data',...
%     'var_names_all_gridded_global_data',...
%     'NH_lat_era5',...
%     'global_lon_era5',...
%     'all_gridded_domain_data',...
%     'var_names_all_gridded_domain_data',...
%     'domain_lat_era5',...
%     'domain_lon_era5',...
%     'domain_means_w_and_wo_seas',...
%     'var_names_domain_means',...
%     'domain_means_annual_seascyc',...
%     'annual_var_names',...
%     'desc_time_daily',...
%     'datetime_desc_time_daily',...
%     '-v7.3')