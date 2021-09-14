close all
clear all

addpath(genpath('/Users/patrickbrown/Documents/MATLAB/'))
addpath(genpath('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/'))

%% load david data

load('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/David_ERA5_data/read_save_ERA5_data_1950_2020.mat',...
    'all_gridded_NH_data',...
    'gridded_WECC_all_var_daily',...
    'NH_lat_era5',...
    'NH_lon_era5',...
    'all_gridded_NH_data_varnames',...
    'WECC_daily_var_names',...
    'WECC_lat',...
    'WECC_lon',...
    'datetime_desc_time_daily');

%% make a distance array to get correct units on divergence, vorticty advection, temperature advection

meridional_distances = NaN(length(NH_lat_era5),length(NH_lon_era5));
zonal_distances = NaN(length(NH_lat_era5),length(NH_lon_era5));

for lat_i = 1:length(NH_lat_era5)-1
    for lon_i = 1:length(NH_lon_era5)-1
        
            %   [KM, NMI, MI] = HAVERSINE(LOC1, LOC2) returns the computed distance in
            %   kilometers (KM)
            %
            %   Examples
            %       haversine('53 08 50N, 001 50 58W', '52 12 16N, 000 08 26E') returns
            %           170.2547
            %       haversine([53.1472 -1.8494], '52 12.16N, 000 08.26E') returns
            %           170.2508
            %       haversine([53.1472 -1.8494], [52.2044 0.1406]) returns 170.2563
            
            %   Inputs
            %       LOC must be either a string specifying the location in degrees,
            %       minutes and seconds, or a 2-valued numeric array specifying the
            %       location in decimal degrees.  If providing a string, the latitude
            %       and longitude must be separated by a comma.
            %
            %       The first element indicates the latitude while the second is the
            %       longitude.
            
            LOC1 = [NH_lat_era5(lat_i) NH_lon_era5(lon_i)];
            LOC2_meridional = [NH_lat_era5(lat_i+1) NH_lon_era5(lon_i)];
            LOC2_zonal = [NH_lat_era5(lat_i) NH_lon_era5(lon_i+1)];
            
            distance_meridional = haversine(LOC1,LOC2_meridional);
            distance_zonal = haversine(LOC1,LOC2_zonal);
            
            meridional_distances(lat_i,lon_i) = 1000*distance_meridional;
            zonal_distances(lat_i,lon_i) = 1000*distance_zonal;
        
    end
end

meridional_distances(end,:) = meridional_distances(end-1,:);
meridional_distances(:,end) = meridional_distances(:,end-1);
zonal_distances(end,:) = zonal_distances(end-1,:);
zonal_distances(:,end) = zonal_distances(:,end-1);

%make cartesian coords

zonal_x_grid = NaN(size(zonal_distances));
meridional_y_grid = NaN(size(meridional_distances));

for lat_i = 1:length(NH_lat_era5)
    zonal_x_grid(lat_i,:) = cumsum(squeeze(meridional_distances(lat_i,:)));
end
for lon_i = 1:length(NH_lon_era5)
    meridional_y_grid(:,lon_i) = cumsum(squeeze(zonal_distances(:,lon_i)));
end

%% make some derived variables

%     {'NH_GEO_250_hPa' } 1
%     {'NH_GEO_500_hPa' } 2
%     {'NH_GEO_700_hPa' } 3
%     {'NH_U_250_hPa'   } 4
%     {'NH_U_500_hPa'   } 5
%     {'NH_U_700_hPa'   } 6
%     {'NH_V_250_hPa'   } 7
%     {'NH_V_500_hPa'   } 8
%     {'NH_V_700_hPa'   } 9
%     {'NH_VOR_500_hPa' } 10
%     {'NH_TEMP_700_hPa'} 11
%     {'NH_VV_700'      } 12
%     {'NH_SLP'         } 13

wind_speed_250 = sqrt(all_gridded_NH_data(:,:,:,4).^2 + all_gridded_NH_data(:,:,:,7).^2);
wind_speed_500 = sqrt(all_gridded_NH_data(:,:,:,5).^2 + all_gridded_NH_data(:,:,:,8).^2);
wind_speed_850 = sqrt(all_gridded_NH_data(:,:,:,6).^2 + all_gridded_NH_data(:,:,:,9).^2);
wind_speeds_gridded_NH_data = cat(4,wind_speed_250,wind_speed_500,wind_speed_850);

clear wind_speed_250
clear wind_speed_500
clear wind_speed_850

divergences = NaN(length(NH_lat_era5),length(NH_lon_era5),length(datetime_desc_time_daily),3);

for time_i = 1:length(datetime_desc_time_daily)

    divergences(:,:,time_i,1) = divergence(zonal_x_grid,meridional_y_grid,squeeze(all_gridded_NH_data(:,:,time_i,4)),squeeze(all_gridded_NH_data(:,:,time_i,7)));
    divergences(:,:,time_i,2) = divergence(zonal_x_grid,meridional_y_grid,squeeze(all_gridded_NH_data(:,:,time_i,5)),squeeze(all_gridded_NH_data(:,:,time_i,8)));
    divergences(:,:,time_i,3) = divergence(zonal_x_grid,meridional_y_grid,squeeze(all_gridded_NH_data(:,:,time_i,6)),squeeze(all_gridded_NH_data(:,:,time_i,9)));

end

vort_advection = NaN(length(NH_lat_era5),length(NH_lon_era5),length(datetime_desc_time_daily));

for time_i = 1:length(datetime_desc_time_daily)
    
    vort_advection(1:end-1,1:end-1,time_i) = -squeeze(all_gridded_NH_data(1:end-1,1:end-1,time_i,5)).*(diff(squeeze(all_gridded_NH_data(1:end-1,:,time_i,10)),1,2)./zonal_distances(1:end-1,1:end-1)) - ...
                                              squeeze(all_gridded_NH_data(1:end-1,1:end-1,time_i,8)).*(diff(squeeze(all_gridded_NH_data(:,1:end-1,time_i,10)),1,1)./meridional_distances(1:end-1,1:end-1));

end

temp_advection = NaN(length(NH_lat_era5),length(NH_lon_era5),length(datetime_desc_time_daily));

for time_i = 1:length(datetime_desc_time_daily)
    
     temp_advection(1:end-1,1:end-1,time_i) = -squeeze(all_gridded_NH_data(1:end-1,1:end-1,time_i,6)).*(diff(squeeze(all_gridded_NH_data(1:end-1,:,time_i,11)),1,2)./zonal_distances(1:end-1,1:end-1)) - ...
                                               squeeze(all_gridded_NH_data(1:end-1,1:end-1,time_i,9)).*(diff(squeeze(all_gridded_NH_data(:,1:end-1,time_i,11)),1,1)./meridional_distances(1:end-1,1:end-1));
 
end

%put SLP at end to maintain previous index order
all_gridded_NH_data = cat(4,all_gridded_NH_data,...
                            wind_speeds_gridded_NH_data,...
                            divergences,...
                            vort_advection,...
                            temp_advection);
                        
clear wind_speeds_gridded_NH_data
clear divergences
clear vort_advection
clear temp_advection

all_gridded_NH_data_varnames_desc = {'250 hPa heights',...
                                     '500 hPa heights',...
                                     '700 hPa heights',...
                                     '250 hPa U wind',...
                                     '500 hpa U wind',...
                                     '700 hpa U wind',...
                                     '250 hPa V wind',...
                                     '500 hpa V wind',...
                                     '700 hpa V wind',...
                                     '500 hpa vorticity',...
                                     '700 hpa temperature',...
                                     '700 Omega',...
                                     'sea level pressure',...
                                     '250 hPa wind speed',...
                                     '500 hpa wind speed',...
                                     '700 hpa wind speed',...
                                     '250 hpa divergence',...
                                     '500 hpa divergence',...
                                     '700 hpa divergence',...
                                     '500 hpa vorticity advection',...
                                     '700 hpa temperature advection'};
                                 
%% make HDD and CDD
 
degree_day_thresh = 291.483;
 
temp_degree_departure = squeeze(gridded_WECC_all_var_daily(:,:,:,1)) - degree_day_thresh;

gridded_WECC_CDD_daily = zeros(size(temp_degree_departure));
gridded_WECC_HDD_daily = zeros(size(temp_degree_departure));

gridded_WECC_CDD_daily(temp_degree_departure >= 0) = temp_degree_departure(temp_degree_departure >= 0);
gridded_WECC_HDD_daily(temp_degree_departure < 0) = abs(temp_degree_departure(temp_degree_departure < 0));

%replce zeros over ocean with NaNs

for lat_i = 1:length(WECC_lat)
    for lon_i = 1:length(WECC_lon)
        
        if isnan(gridded_WECC_all_var_daily(lat_i,lon_i,100,3)) == 1
            
            gridded_WECC_CDD_daily(lat_i,lon_i,:) = NaN;
            gridded_WECC_HDD_daily(lat_i,lon_i,:) = NaN;
            
        end
        
    end
end


gridded_WECC_all_var_daily_2 = cat(4,gridded_WECC_all_var_daily(:,:,:,3),...
                                     gridded_WECC_all_var_daily(:,:,:,2),...
                                     gridded_WECC_CDD_daily,...
                                     gridded_WECC_HDD_daily,...
                                     temp_degree_departure);
                                 
WECC_daily_var_names = {'wind power (W/m^2)',...
                        'surface solar power (W/m^2)',...
                        'cooling degree-days (C*days)',...
                        'heating degree-days (C*days)',...
                        'temperature departure (C)'};
                                 
clear gridded_WECC_all_var_daily
clear temp_degree_departure
clear gridded_WECC_HDD_daily
clear gridded_WECC_CDD_daily

%% plot clim

gridded_WECC_all_var_daily_ann_clim = squeeze(mean(gridded_WECC_all_var_daily_2,3));
 
min_color_ranges = [0 100 0 0 0 ];
max_color_ranges = [5 270 10 30 30];
 
row_num = 3;
col_num = 2;
 
coast = load('coastlines');
 
WECC_lat = double(WECC_lat);
WECC_lat = double(WECC_lat);
 
FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',12);
set(0,'defaultAxesFontName', 'helvetica')
 
for var_i = 1:length(WECC_daily_var_names)

    h = subplot(row_num,col_num,var_i);
    hold on

    axesm('robinson',...
    'Frame', 'on',...
    'Grid', 'off',...
    'maplatlim',[min(WECC_lat) max(WECC_lat)],...
    'maplonlim',[min(WECC_lon) max(WECC_lon)])
    tightmap

    %if var_i == 4; colormap(redblue); end
%      min_color_range = 0.0;
%      max_color_range = 0.35;

    num_contours = 40;

    %contourfm(WECC_lat,WECC_lon,squeeze(gridded_WECC_all_var_daily_ann_clim(:,:,var_i)), 'LineStyle','none');
    pcolorm(WECC_lat,WECC_lon,squeeze(gridded_WECC_all_var_daily_ann_clim(:,:,var_i)), 'LineStyle','none');        
    caxis([min_color_ranges(var_i) max_color_ranges(var_i)])
    
    %if var_i == 4; colormap(h,redblue); end
    %if var_i ~= 4; colormap(h,parula); end
    
    geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')

    t=contourcbar;
    %set(get(t,'ylabel'),'String', 'Fraction of days in season');
    
    title(WECC_daily_var_names{var_i})
    
end

%% reorganize matricies to have these dimensions:

% 1) lat
% 2) lon
% 3) week of year
% 5) year
% 6) variable

datetime_desc_time_daily_numeric = linspace(1950.0,2020+364/365,length(datetime_desc_time_daily));

% ----- first put into day-year ---------

    % ---- WECC ------

        ann_desc_time_daily_numeric = 1950:2021; %code below needs the handle of 2021 but data goes through end of 2020

        gridded_WECC_day_year = NaN(length(WECC_lat),length(WECC_lon),366,length(ann_desc_time_daily_numeric),length(WECC_daily_var_names));

        datetime_desc_time_daily_numeric_day_year = datetime;

        for var_i = 1:length(WECC_daily_var_names)
            for year_i = 1:length(ann_desc_time_daily_numeric)-1

                day_segment_for_this_year_inds = find(isbetween(datetime_desc_time_daily,datetime(ann_desc_time_daily_numeric(year_i),01,01),datetime(ann_desc_time_daily_numeric(year_i+1),01,01)) == 1);

                gridded_WECC_day_year(:,:,1:length(day_segment_for_this_year_inds),year_i,var_i) = gridded_WECC_all_var_daily_2(:,:,day_segment_for_this_year_inds,var_i);

                datetime_desc_time_daily_numeric_day_year(1:length(day_segment_for_this_year_inds),year_i) = datetime_desc_time_daily(day_segment_for_this_year_inds);

            end
        end

    % ---- NH ------

        all_gridded_NH_data_day_year = NaN(length(NH_lat_era5),length(NH_lon_era5),366,length(ann_desc_time_daily_numeric),length(all_gridded_NH_data_varnames_desc));

        for var_i = 1:length(all_gridded_NH_data_varnames_desc)
            for year_i = 1:length(ann_desc_time_daily_numeric)-1

                day_segment_for_this_year_inds = find(isbetween(datetime_desc_time_daily,datetime(ann_desc_time_daily_numeric(year_i),01,01),datetime(ann_desc_time_daily_numeric(year_i+1),01,01)) == 1);

                all_gridded_NH_data_day_year(:,:,1:length(day_segment_for_this_year_inds),year_i,var_i) = all_gridded_NH_data(:,:,day_segment_for_this_year_inds,var_i);

            end
        end

% put into week discretization

disc_timescale = 7;

    % ---- WECC ------

        gridded_WECC_week_year = NaN(length(WECC_lat),length(WECC_lon),52,length(ann_desc_time_daily_numeric),length(WECC_daily_var_names));

        for var_i = 1:length(WECC_daily_var_names)
            for year_i = 1:length(ann_desc_time_daily_numeric)-1

                week_i = 1;

                for end_day_i = disc_timescale:disc_timescale:366

                    gridded_WECC_week_year(:,:,week_i,year_i,var_i) = mean(gridded_WECC_day_year(:,:,end_day_i-disc_timescale+1:end_day_i,year_i,var_i),3);

                    week_i = week_i + 1;

                end
            end
        end

    % ---- NH ------

        gridded_NH_week_year = NaN(length(NH_lat_era5),length(NH_lon_era5),52,length(ann_desc_time_daily_numeric),length(all_gridded_NH_data_varnames_desc));

        for var_i = 1:length(all_gridded_NH_data_varnames_desc)
            for year_i = 1:length(ann_desc_time_daily_numeric)-1

                week_i = 1;

                for end_day_i = disc_timescale:disc_timescale:366

                    gridded_NH_week_year(:,:,week_i,year_i,var_i) = mean(all_gridded_NH_data_day_year(:,:,end_day_i-disc_timescale+1:end_day_i,year_i,var_i),3);

                    week_i = week_i + 1;

                end
            end
        end
        
%% make anomaly arrays

years = 1950:2020;
weeks = 1:52; %with exactly 52 weeks, at least 1 day per year (2 on leap years) are thrown out and its the last days of each year)

% size(gridded_WECC_week_year)
% ans =
%     31    36    52    40     4

gridded_WECC_week_year_anoms = NaN(size(gridded_WECC_week_year));

for var_i = 1:length(WECC_daily_var_names)
    var_i/length(WECC_daily_var_names)
    for year_i = 1:length(years)
        for week_i = 1:length(weeks)
            for lat_i = 1:length(WECC_lat)
                for lon_i = 1:length(WECC_lon)
                    
                    gridded_WECC_week_year_anoms(lat_i,lon_i,week_i,year_i,var_i) = gridded_WECC_week_year(lat_i,lon_i,week_i,year_i,var_i) - nanmean(squeeze(gridded_WECC_week_year(lat_i,lon_i,week_i,:,var_i)));
                    
                end
            end
        end
    end
end

% size(gridded_NH_week_year)
% ans =
%     21   144    52    40    19
    
gridded_NH_week_year_anoms = NaN(size(gridded_NH_week_year));

for var_i = 1:length(all_gridded_NH_data_varnames_desc)
    var_i/length(all_gridded_NH_data_varnames_desc)
    for year_i = 1:length(years)
        for week_i = 1:length(weeks)
            for lat_i = 1:length(NH_lat_era5)
                for lon_i = 1:length(NH_lon_era5)
                    
                    gridded_NH_week_year_anoms(lat_i,lon_i,week_i,year_i,var_i) = gridded_NH_week_year(lat_i,lon_i,week_i,year_i,var_i) - nanmean(squeeze(gridded_NH_week_year(lat_i,lon_i,week_i,:,var_i)));
                    
                end
            end
        end
    end
end

%% make anom deviations

gridded_WECC_week_year_frac_devs_wrt_ann_mean = NaN(size(gridded_WECC_week_year_anoms));
gridded_WECC_week_year_frac_devs_wrt_woy_mean = NaN(size(gridded_WECC_week_year_anoms));


for var_i = 1:length(WECC_daily_var_names)
    for lat_i =  1:length(WECC_lat)
        for lon_i = 1:length(WECC_lon)
            for year_i = 1:length(years)
                for week_i = 1:length(weeks)
            
                gridded_WECC_week_year_frac_devs_wrt_ann_mean(lat_i,lon_i,week_i,year_i,var_i) = ...
                (gridded_WECC_week_year(lat_i,lon_i,week_i,year_i,var_i) - nanmean(nanmean(gridded_WECC_week_year(lat_i,lon_i,:,:,var_i),3),4))...
                ./nanmean(nanmean(gridded_WECC_week_year(lat_i,lon_i,:,:,var_i),3),4);
            
                gridded_WECC_week_year_frac_devs_wrt_woy_mean(lat_i,lon_i,week_i,year_i,var_i) = ...
                (gridded_WECC_week_year(lat_i,lon_i,week_i,year_i,var_i) - nanmean(gridded_WECC_week_year(lat_i,lon_i,week_i,:,var_i),4))...
                ./nanmean(gridded_WECC_week_year(lat_i,lon_i,week_i,:,var_i),4);                            
                
                end
            end
        end
    end
end

%% save

save('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/David_ERA5_data/read_org_week_year_1950_2020.mat',...
     'years',...
     'weeks',...
     'WECC_daily_var_names',...
     'all_gridded_NH_data_varnames_desc',...
     'NH_lat_era5',...
     'NH_lon_era5',...
     'WECC_lat',...
     'WECC_lon',...
     'gridded_NH_week_year',...
     'gridded_NH_week_year_anoms',...
     'gridded_WECC_week_year',...
     'gridded_WECC_week_year_anoms',...
     'gridded_WECC_week_year_frac_devs_wrt_ann_mean',...
     'gridded_WECC_week_year_frac_devs_wrt_woy_mean')
