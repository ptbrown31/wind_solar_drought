close all
clear all

addpath(genpath('/Users/patrickbrown/Documents/MATLAB/'))
addpath(genpath('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/'))


%% load clim expl

% load('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/clim_expl_data/read_save_clim_expl_data.mat',...
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
%     'datetime_desc_time_daily')

%% load david data

load('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/David_ERA5_data/read_save_ERA5_data.mat',...
    'all_gridded_NH_data',...
    'gridded_WECC_all_var_daily',...
    'NH_lat_era5',...
    'NH_lon_era5',...
    'all_gridded_NH_data_varnames',...
    'WECC_daily_var_names',...
    'WECC_lat',...
    'WECC_lon',...
    'datetime_desc_time_daily');

%convert ssrd to W/m^2 (it is in J/m^2) (it was originally hour so divide
%by number of seconds in an hour)

gridded_WECC_all_var_daily(:,:,:,2) = gridded_WECC_all_var_daily(:,:,:,2)./3600;

%wind power appears to be in KW
%(https://www.nrel.gov/gis/assets/pdfs/windsmodel4pub1-1-9base200904enh.pdf)
gridded_WECC_all_var_daily(:,:,:,3) = gridded_WECC_all_var_daily(:,:,:,3).*1000;

%% make some derived variables

wind_speed_250 = sqrt(all_gridded_NH_data(:,:,:,5).^2 + all_gridded_NH_data(:,:,:,8).^2);
wind_speed_500 = sqrt(all_gridded_NH_data(:,:,:,6).^2 + all_gridded_NH_data(:,:,:,9).^2);
wind_speed_850 = sqrt(all_gridded_NH_data(:,:,:,7).^2 + all_gridded_NH_data(:,:,:,10).^2);
wind_speeds_gridded_NH_data = cat(4,wind_speed_250,wind_speed_500,wind_speed_850);

divergences = NaN(length(NH_lat_era5),length(NH_lon_era5),length(datetime_desc_time_daily),3);

for time_i = 1:length(datetime_desc_time_daily)

    divergences(:,:,time_i,1) = divergence(squeeze(all_gridded_NH_data(:,:,time_i,5)),squeeze(all_gridded_NH_data(:,:,time_i,8)));
    divergences(:,:,time_i,2) = divergence(squeeze(all_gridded_NH_data(:,:,time_i,6)),squeeze(all_gridded_NH_data(:,:,time_i,9)));
    divergences(:,:,time_i,3) = divergence(squeeze(all_gridded_NH_data(:,:,time_i,7)),squeeze(all_gridded_NH_data(:,:,time_i,10)));

end

vort_advection = NaN(length(NH_lat_era5),length(NH_lon_era5),length(datetime_desc_time_daily));

for time_i = 1:length(datetime_desc_time_daily)
    
    vort_advection(1:end-1,1:end-1,time_i) = -squeeze(all_gridded_NH_data(1:end-1,1:end-1,time_i,6)).*diff(squeeze(all_gridded_NH_data(:,1:end-1,time_i,11)),1,1) - ...
                                              squeeze(all_gridded_NH_data(1:end-1,1:end-1,time_i,8)).*diff(squeeze(all_gridded_NH_data(1:end-1,:,time_i,11)),1,2);

end

temp_advection = NaN(length(NH_lat_era5),length(NH_lon_era5),length(datetime_desc_time_daily));

for time_i = 1:length(datetime_desc_time_daily)
    
     temp_advection(1:end-1,1:end-1,time_i) = -squeeze(all_gridded_NH_data(1:end-1,1:end-1,time_i,6)).*diff(squeeze(all_gridded_NH_data(:,1:end-1,time_i,4)),1,1) - ...
                                               squeeze(all_gridded_NH_data(1:end-1,1:end-1,time_i,8)).*diff(squeeze(all_gridded_NH_data(1:end-1,:,time_i,4)),1,2);
 
end

all_gridded_NH_data_varnames_desc = {'250 hPa heights',...
                                     '500 hPa heights',...
                                     '850 hPa heights',...
                                     '850 hpa temperature',...
                                     '250 hPa U wind',...
                                     '500 hpa U wind',...
                                     '850 hpa U wind',...
                                     '250 hPa V wind',...
                                     '500 hpa V wind',...
                                     '850 hpa V wind',...
                                     '500 hpa vorticity',...
                                     '250 hPa wind speed',...
                                     '500 hpa wind speed',...
                                     '850 hpa wind speed',...
                                     '250 hpa divergence',...
                                     '500 hpa divergence',...
                                     '850 hpa divergence',...
                                     '500 hpa vorticity advection',...
                                     '850 hpa temperature advection'};
                                 
%% plot climatological stuff as a check
                                 
        all_gridded_NH_data = cat(4,all_gridded_NH_data,...
                                    wind_speeds_gridded_NH_data,...
                                    divergences,...
                                    vort_advection,...
                                    temp_advection);

        all_gridded_NH_data_clim = mean(all_gridded_NH_data,3);

        % min_color_ranges = [100 0 260 -20 0 0 0];
        % max_color_ranges = [275 400 300 20 25 20E6 1000000];

        row_num = 4;
        col_num = 3;

        coast = load('coastlines');

        NH_lat_era5 = double(NH_lat_era5);
        NH_lon_era5 = double(NH_lon_era5);

        FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
        set(gcf,'color',[1 1 1]);
        set(0, 'DefaultAxesFontSize',12);
        set(0,'defaultAxesFontName', 'helvetica')

        for var_i = 1:length(all_gridded_NH_data_varnames)

            h = subplot(row_num,col_num,var_i);
            hold on

            axesm('robinson',...
            'Frame', 'on',...
            'Grid', 'off',...
            'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
            'maplonlim',[min(NH_lon_era5) max(NH_lon_era5)])
            tightmap

            %if var_i == 4; colormap(redblue); end
        %      min_color_range = 0.0;
        %      max_color_range = 0.35;

            %contourfm(WECC_lat,WECC_lon,squeeze(gridded_WECC_all_var_daily_ann_clim(:,:,var_i)), 'LineStyle','none');
            pcolorm(NH_lat_era5,NH_lon_era5,squeeze(all_gridded_NH_data_clim(:,:,1,var_i)), 'LineStyle','none');        
            %caxis([min_color_ranges(var_i) max_color_ranges(var_i)])

            geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')

            t=contourcbar;
            %set(get(t,'ylabel'),'String', 'Fraction of days in season');

            title(all_gridded_NH_data_varnames_desc{var_i})

        end

        % better clim plots

        [NH_lat_era5_grid,NH_lon_era5_grid] = meshgrid(NH_lat_era5,NH_lon_era5);

        num_conts = 13;

        all_gridded_NH_data_clim = double(all_gridded_NH_data_clim);
        %wind_speeds_gridded_NH_data_clim = double(wind_speeds_gridded_NH_data_clim);

        height_vars = {'250 hPa heights, wind, isotachs',...
                       '500 hPa heights, wind, vorticity',...
                       '850 hPa heights, wind, temperature'};

        row_num = 3;
        col_num = 1;

        FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
        set(gcf,'color',[1 1 1]);
        set(0, 'DefaultAxesFontSize',12);
        set(0,'defaultAxesFontName', 'helvetica')

        for var_i = 1:length(height_vars)

            subplot(row_num,col_num,var_i);
            hold on

            axesm('robinson',...
            'Frame', 'on',...
            'Grid', 'off',...
            'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
            'maplonlim',[110 340])
           tightmap

            %if var_i == 4; colormap(redblue); end
        %      min_color_range = 0.0;
        %      max_color_range = 0.35;

            if var_i == 3

             contourfm(NH_lat_era5,NH_lon_era5,squeeze(all_gridded_NH_data_clim(:,:,1,4)),40,'Linestyle','none')

             caxis([min(min(squeeze(all_gridded_NH_data_clim(:,:,1,4)))) 1*max(max(squeeze(all_gridded_NH_data_clim(:,:,1,4))))])
             t=contourcbar;

            end

            if var_i == 2

             contourfm(NH_lat_era5,NH_lon_era5,squeeze(all_gridded_NH_data_clim(:,:,1,11)),80,'Linestyle','none')

             caxis([min(min(squeeze(all_gridded_NH_data_clim(:,:,1,11)))) 0.2*max(max(squeeze(all_gridded_NH_data_clim(:,:,1,11))))])
             t=contourcbar;

            end
            if var_i == 1

             contourfm(NH_lat_era5,NH_lon_era5,squeeze(all_gridded_NH_data_clim(:,:,1,12)),80,'Linestyle','none')

             caxis([min(min(squeeze(all_gridded_NH_data_clim(:,:,1,12)))) 1*max(max(squeeze(all_gridded_NH_data_clim(:,:,1,12))))])
             t=contourcbar;

            end

            contourm(NH_lat_era5,NH_lon_era5,squeeze(all_gridded_NH_data_clim(:,:,1,var_i)),num_conts,'LineColor','r')    
            quiverm(NH_lat_era5_grid,NH_lon_era5_grid,squeeze(all_gridded_NH_data_clim(:,:,1,var_i+7))',squeeze(all_gridded_NH_data_clim(:,:,1,var_i+4))')

            %pcolorm(NH_lat_era5,NH_lon_era5,squeeze(all_gridded_NH_data_clim(:,:,1,var_i)), 'LineStyle','none');        
            %caxis([min_color_ranges(var_i) max_color_ranges(var_i)])

            geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')

            %set(get(t,'ylabel'),'String', 'Fraction of days in season');

            title(height_vars{var_i})

        end

%% load gridded population

load('/Users/patrickbrown/Documents/MATLAB/seas_forecast_land_extremes/read_gridded_pop_GDP_crop.mat',...
    'gridded_pop_all_data',...
    'daily_SAT_land_mask_lat_lon',...
    'lat_BEST',...
    'lon_BEST',...
    'var_names_pop')

gridded_pop_all_data_2020 = squeeze(gridded_pop_all_data(:,:,5));

%gridded_pop_all_data_2020 = gridded_pop_all_data_2020';

%lon_BEST goes from -180 to 180, need to give WECC_lon in these terms
WECC_lon_adjust = WECC_lon-360;

%interpolate to ERA5 grid

[WECC_lat_grid,WECC_lon_grid] = meshgrid(WECC_lat,WECC_lon_adjust);

WECC_lat_grid = double(WECC_lat_grid);
WECC_lon_grid = double(WECC_lon_grid);
% WECC_lat = double(WECC_lat);
% WECC_lon = double(WECC_lon);

population_2020_regrid = interp2(lon_BEST,lat_BEST,gridded_pop_all_data_2020,WECC_lon_grid,WECC_lat_grid);

%% make date/time array

datetime_desc_time_daily_numeric = linspace(1979.0,2018+364/365,length(datetime_desc_time_daily));

%% make temperature departures and people departures

degree_day_thresh = 291.483;

temp_degree_departure = squeeze(gridded_WECC_all_var_daily(:,:,:,1)) - degree_day_thresh;
temp_degree_abs_departure = abs(temp_degree_departure);

people_temp_degree_abs_departure = NaN(size(gridded_WECC_all_var_daily,1),size(gridded_WECC_all_var_daily,2),size(gridded_WECC_all_var_daily,3));

for time_i = 1:size(people_temp_degree_abs_departure,3)
    
    people_temp_degree_abs_departure(:,:,time_i) = temp_degree_abs_departure(:,:,time_i).*population_2020_regrid';
    
end

gridded_WECC_temp_var_daily = cat(4,temp_degree_departure,temp_degree_abs_departure,people_temp_degree_abs_departure);

gridded_WECC_all_var_daily_reorder = cat(4,...
                                         gridded_WECC_all_var_daily(:,:,:,2:3),...
                                         gridded_WECC_all_var_daily(:,:,:,1),...
                                         gridded_WECC_temp_var_daily);

%% make the clim maps for WECC
% size(gridded_WECC_all_var_daily)
%           31          36       14610           3

gridded_WECC_all_var_daily_ann_clim = squeeze(mean(gridded_WECC_all_var_daily_reorder,3));

gridded_WECC_all_var_daily_ann_clim_plus_pop = cat(3,...
                                                   gridded_WECC_all_var_daily_ann_clim,...
                                                   population_2020_regrid');

WECC_daily_var_names_plus_pop = {'surface solar power (W/m^2)',...
                                 'wind power (W/m^2)',...
                                 'temperature (K)',...
                                 'temp departure from 291 (k)',...
                                 'temp abs departure from 291 (k)',...
                                 'people-degree-days',...
                                 'population (count)'};

min_color_ranges = [100 0 260 -20 0 0 0];
max_color_ranges = [275 400 300 20 25 20E6 1000000];

row_num = 3;
col_num = 2;

coast = load('coastlines');

lat_BEST = double(lat_BEST);
lon_BEST = double(lon_BEST);

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',12);
set(0,'defaultAxesFontName', 'helvetica')

for var_i = 1:length(WECC_daily_var_names_plus_pop)-1

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
    pcolorm(WECC_lat,WECC_lon,squeeze(gridded_WECC_all_var_daily_ann_clim_plus_pop(:,:,var_i)), 'LineStyle','none');        
    caxis([min_color_ranges(var_i) max_color_ranges(var_i)])
    
    if var_i == 4; colormap(h,redblue); end
    if var_i ~= 4; colormap(h,parula); end
    
    geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')

    t=contourcbar;
    %set(get(t,'ylabel'),'String', 'Fraction of days in season');
    
    title(WECC_daily_var_names_plus_pop{var_i})
    
end

%% make time series

daily_WECC_time_series = NaN(size(gridded_WECC_all_var_daily_reorder,3),size(gridded_WECC_all_var_daily_reorder,4));

for var_i = 1:size(gridded_WECC_all_var_daily_reorder,4)
    
    daily_WECC_time_series(:,var_i) = threed2oned_lat_lon_2(squeeze(gridded_WECC_all_var_daily_reorder(:,:,:,var_i)),WECC_lat,WECC_lon);
    
end

row_num = 3;
col_num = 2;

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',12);
set(0,'defaultAxesFontName', 'helvetica')

for var_i = 1:size(gridded_WECC_all_var_daily_reorder,4)

    subplot(row_num,col_num,var_i)
    hold on
    
    plot(datetime_desc_time_daily_numeric,daily_WECC_time_series(:,var_i),'-k','LineWidth',1)
    
    title(strcat('WECC Mean, ',WECC_daily_var_names_plus_pop{var_i}))
    
end

%% make more persistent measures

drought_timescale_days = 7;

% CSDI
    % Cold spell duration index (CSDI) 
    % is defined as annual or seasonal 
    % count of days with at least 6 
    % consecutive days when the daily 
    % minimum T fall below the 10th percentile 
    % in the calendar 5-day window for the 
    % base period 1979-2009.
    
    %Becasue of demand it makes sense to do everything on weekly (7 day)
    %schedule
    
daily_WECC_time_series_7_day_smooth = NaN(size(daily_WECC_time_series));

for var_i = 1:size(gridded_WECC_all_var_daily_reorder,4)
    
    daily_WECC_time_series_7_day_smooth(:,var_i) = smooth(squeeze(daily_WECC_time_series(:,var_i)),drought_timescale_days);
    
end

%% plot smooth time series

row_num = 3;
col_num = 2;

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',12);
set(0,'defaultAxesFontName', 'helvetica')

for var_i = 1:size(gridded_WECC_all_var_daily_reorder,4)

    subplot(row_num,col_num,var_i)
    hold on
    
    plot(datetime_desc_time_daily_numeric,daily_WECC_time_series_7_day_smooth(:,var_i),'-k','LineWidth',1)
    
    title(strcat('WECC Mean, 7-day mean',WECC_daily_var_names_plus_pop{var_i}))
    
end

%% line everything up in annual to make annual data

ann_desc_time_daily_numeric = 1979:2019;

daily_WECC_time_series_day_year = NaN(366,length(ann_desc_time_daily_numeric),size(daily_WECC_time_series,2));
daily_WECC_time_series_day_year_7_day_smooth = NaN(366,length(ann_desc_time_daily_numeric),size(daily_WECC_time_series,2));

datetime_desc_time_daily_numeric_day_year = datetime;

for var_i = 1:size(daily_WECC_time_series,2)
    for year_i = 1:length(ann_desc_time_daily_numeric)-1

        day_segment_for_this_year_inds = find(isbetween(datetime_desc_time_daily,datetime(ann_desc_time_daily_numeric(year_i),01,01),datetime(ann_desc_time_daily_numeric(year_i+1),01,01)) == 1);

        daily_WECC_time_series_day_year(1:length(day_segment_for_this_year_inds),year_i,var_i) = daily_WECC_time_series(day_segment_for_this_year_inds,var_i);
        daily_WECC_time_series_day_year_7_day_smooth(1:length(day_segment_for_this_year_inds),year_i,var_i) = daily_WECC_time_series_7_day_smooth(day_segment_for_this_year_inds,var_i);

        datetime_desc_time_daily_numeric_day_year(1:length(day_segment_for_this_year_inds),year_i) = datetime_desc_time_daily(day_segment_for_this_year_inds);

    end
end

%% plot trends in extremes for each thing individually

most_extreme_per_year = NaN(length(ann_desc_time_daily_numeric),size(daily_WECC_time_series,2));

for var_i = 1:size(daily_WECC_time_series,2)
    for year_i = 1:length(ann_desc_time_daily_numeric)
        if var_i <= 2
            
            most_extreme_per_year(year_i,var_i) = min(squeeze(daily_WECC_time_series_day_year_7_day_smooth(:,year_i,var_i)));
            
        end
        if var_i > 2
            
            most_extreme_per_year(year_i,var_i) = max(squeeze(daily_WECC_time_series_day_year_7_day_smooth(:,year_i,var_i)));
            
        end
    end
end

row_num = 3;
col_num = 2;

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',12);
set(0,'defaultAxesFontName', 'helvetica')

for var_i = 1:size(gridded_WECC_all_var_daily_reorder,4)

    subplot(row_num,col_num,var_i)
    hold on
    
    plot(ann_desc_time_daily_numeric,most_extreme_per_year(:,var_i),'-k','LineWidth',1)
    scatter(ann_desc_time_daily_numeric,most_extreme_per_year(:,var_i),25,'k','filled')
    lsline
    
    title(strcat('most extreme week that year, ',WECC_daily_var_names_plus_pop{var_i}))
    
end

%% plot concurrent extremes

%make these

%(fractional wind) + (fractional solar)
%(fractional wind) + (fractional people-degree-days)
%(fractional solar) + (fractional people-degree-days)

%(fractional wind) + (fractional solar) + (fractional people-degree-days)


frac_sum_wind_solar = daily_WECC_time_series_7_day_smooth(:,1)/mean(daily_WECC_time_series_7_day_smooth(:,1)) + daily_WECC_time_series_7_day_smooth(:,2)/mean(daily_WECC_time_series_7_day_smooth(:,2));
frac_sum_wind_PDD = daily_WECC_time_series_7_day_smooth(:,1)/mean(daily_WECC_time_series_7_day_smooth(:,1)) + 1./(daily_WECC_time_series_7_day_smooth(:,6)/mean(daily_WECC_time_series_7_day_smooth(:,6)));
frac_sum_solar_PDD = daily_WECC_time_series_7_day_smooth(:,2)/mean(daily_WECC_time_series_7_day_smooth(:,2)) + 1./(daily_WECC_time_series_7_day_smooth(:,6)/mean(daily_WECC_time_series_7_day_smooth(:,6)));

frac_sum_wind_solar_PDD = daily_WECC_time_series_7_day_smooth(:,1)/mean(daily_WECC_time_series_7_day_smooth(:,1)) + ...
                          daily_WECC_time_series_7_day_smooth(:,2)/mean(daily_WECC_time_series_7_day_smooth(:,2)) + ...
                          1./(daily_WECC_time_series_7_day_smooth(:,6)/mean(daily_WECC_time_series_7_day_smooth(:,6)));
                      
frac_sum_var_names = {'wind + solar','wind + PDD','solar + PDD','wind + solar + PDD'};

frac_sum_all = cat(2,frac_sum_wind_solar,frac_sum_wind_PDD,frac_sum_solar_PDD,frac_sum_wind_solar_PDD);

% do the lowest per year analysis

frac_sum_all_day_year = NaN(366,length(ann_desc_time_daily_numeric),length(frac_sum_var_names));

for var_i = 1:length(frac_sum_var_names)
    for year_i = 1:length(ann_desc_time_daily_numeric)-1

        day_segment_for_this_year_inds = find(isbetween(datetime_desc_time_daily,datetime(ann_desc_time_daily_numeric(year_i),01,01),datetime(ann_desc_time_daily_numeric(year_i+1),01,01)) == 1);

        frac_sum_all_day_year(1:length(day_segment_for_this_year_inds),year_i,var_i) = frac_sum_all(day_segment_for_this_year_inds,var_i);

    end
end

frac_sum_most_extreme_per_year = NaN(length(ann_desc_time_daily_numeric),size(daily_WECC_time_series,2));

for var_i = 1:length(frac_sum_var_names)
    for year_i = 1:length(ann_desc_time_daily_numeric)
            
            frac_sum_most_extreme_per_year(year_i,var_i) = min(squeeze(frac_sum_all_day_year(:,year_i,var_i)));

    end
end

row_num = 2;
col_num = 2;

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',12);
set(0,'defaultAxesFontName', 'helvetica')

for var_i = 1:length(frac_sum_var_names)

    subplot(row_num,col_num,var_i)
    hold on
    
    plot(datetime_desc_time_daily_numeric,frac_sum_all(:,var_i),'-k','LineWidth',1)
    scatter(ann_desc_time_daily_numeric,frac_sum_most_extreme_per_year(:,var_i),25,'r','filled')
    
    title(strcat('frac sum, ',frac_sum_var_names{var_i}))
    
    xlim([1979 2019])
    
    xlabel('year')
    ylabel('stress intensity, lower = more intense')
    
end

row_num = 2;
col_num = 2;

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',12);
set(0,'defaultAxesFontName', 'helvetica')

for var_i = 1:length(frac_sum_var_names)

    subplot(row_num,col_num,var_i)
    hold on
    
    plot(ann_desc_time_daily_numeric,frac_sum_most_extreme_per_year(:,var_i),'-k','LineWidth',1)
    scatter(ann_desc_time_daily_numeric,frac_sum_most_extreme_per_year(:,var_i),25,'r','filled')
    lsline
    
    title(strcat('most extreme week that year, ',frac_sum_var_names{var_i}))
    
    xlim([1979 2019])
    
    xlabel('year')
    ylabel('stress intensity, lower = more intense')
    
end

%% scatter the weeks

% FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
% set(gcf,'color',[1 1 1]);
% set(0, 'DefaultAxesFontSize',12);
% set(0,'defaultAxesFontName', 'helvetica')
% 
%     xlabel('solar for that week / annual mean solar')
%     ylabel('wind for that week / annual mean wind')
% 
%     xlim([0.2 2.2])
%     ylim([0.2 2.2])
% 
%     hold on
%     for week_i = 1:size(daily_WECC_time_series_7_day_smooth,1)
% 
%         %scatter(daily_WECC_time_series_7_day_smooth(1:7:end,1)/mean(daily_WECC_time_series_7_day_smooth(:,1)),daily_WECC_time_series_7_day_smooth(1:7:end,2)/mean(daily_WECC_time_series_7_day_smooth(:,2)),30,'k','filled')
%         plot(daily_WECC_time_series_7_day_smooth(1:week_i,1)/mean(daily_WECC_time_series_7_day_smooth(:,1)),daily_WECC_time_series_7_day_smooth(1:week_i,2)/mean(daily_WECC_time_series_7_day_smooth(:,2)),'-k')
% 
%         drawnow
%        % pause(0.1)
% 
% 
%     end

%% make mean seas cycles

seas_WECC_time_series_day_year = nanmean(daily_WECC_time_series_day_year,2);
seas_WECC_time_series_day_year_smooth = NaN(size(seas_WECC_time_series_day_year));

for var_i = 1:size(gridded_WECC_all_var_daily_reorder,4)
    
    seas_WECC_time_series_day_year_smooth(:,1,var_i) = smooth(squeeze(seas_WECC_time_series_day_year(:,1,var_i)),21,'lowess');
    
end

day_of_year = 1:366;

row_num = 3;
col_num = 2;

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',12);
set(0,'defaultAxesFontName', 'helvetica')

for var_i = 1:size(gridded_WECC_all_var_daily_reorder,4)

    subplot(row_num,col_num,var_i)
    hold on
        
    for year_i = 1:length(ann_desc_time_daily_numeric)
        
        plot(day_of_year,daily_WECC_time_series_day_year(:,year_i,var_i),'-k','LineWidth',1)
        
    end
    
    plot(day_of_year,seas_WECC_time_series_day_year_smooth(:,var_i),'-r','LineWidth',2)
    
    xlim([1 366])
    if var_i <= 2; ylim([0 500]); end
    
    xlabel('day of year')
    title(WECC_daily_var_names_plus_pop{var_i})
    
end

%% make daily clim and anomalies of all gridded data

%first must make day-ann arrays

%     'all_gridded_NH_data',...
%     'gridded_WECC_all_var_daily',...
%     'NH_lat_era5',...
%     'NH_lon_era5',...
%     'all_gridded_NH_data_varnames',...
%     'WECC_daily_var_names',...
%     'WECC_lat',...
%     'WECC_lon',...
%     'datetime_desc_time_daily');

gridded_WECC_all_var_daily_day_ann = NaN(length(WECC_lat),...
                                         length(WECC_lon),...
                                         366,...
                                         length(ann_desc_time_daily_numeric),...
                                         length(WECC_daily_var_names_plus_pop)-1);

for var_i = 1:length(WECC_daily_var_names_plus_pop)-1
    for year_i = 1:length(ann_desc_time_daily_numeric)-1

        day_segment_for_this_year_inds = find(isbetween(datetime_desc_time_daily,datetime(ann_desc_time_daily_numeric(year_i),01,01),datetime(ann_desc_time_daily_numeric(year_i+1),01,01)) == 1);

        %WECC
        gridded_WECC_all_var_daily_day_ann(:,:,1:length(day_segment_for_this_year_inds),year_i,var_i) = ...
            gridded_WECC_all_var_daily_reorder(:,:,day_segment_for_this_year_inds,var_i);


    end
end

all_gridded_NH_data_day_ann = NaN(length(NH_lat_era5),...
                                  length(NH_lon_era5),...
                                  366,...
                                  length(ann_desc_time_daily_numeric),...
                                  length(all_gridded_NH_data_varnames_desc));

for var_i = 1:length(all_gridded_NH_data_varnames_desc)
    for year_i = 1:length(ann_desc_time_daily_numeric)-1

        day_segment_for_this_year_inds = find(isbetween(datetime_desc_time_daily,datetime(ann_desc_time_daily_numeric(year_i),01,01),datetime(ann_desc_time_daily_numeric(year_i+1),01,01)) == 1);

        %all NH
        all_gridded_NH_data_day_ann(:,:,1:length(day_segment_for_this_year_inds),year_i,var_i) = ...
            all_gridded_NH_data(:,:,day_segment_for_this_year_inds,var_i);


    end
end

%make daily climatologies

all_gridded_NH_data_daily_clim = nanmean(all_gridded_NH_data_day_ann,4);
gridded_WECC_all_var_daily_clim = nanmean(gridded_WECC_all_var_daily_day_ann,4);

%% make gridded anoms

gridded_WECC_all_var_daily_anoms = NaN(length(WECC_lat),...
                                       length(WECC_lon),...
                                       length(datetime_desc_time_daily),...
                                       length(WECC_daily_var_names_plus_pop)-1);

for var_i = 1:length(WECC_daily_var_names_plus_pop)-1
    for year_i = 1:length(ann_desc_time_daily_numeric)-1

        day_segment_for_this_year_inds = find(isbetween(datetime_desc_time_daily,datetime(ann_desc_time_daily_numeric(year_i),01,01),datetime(ann_desc_time_daily_numeric(year_i+1),01,01)) == 1);
        
        gridded_WECC_all_var_daily_anoms(:,:,day_segment_for_this_year_inds,var_i) = ...
            gridded_WECC_all_var_daily_reorder(:,:,day_segment_for_this_year_inds,var_i) - ...
            gridded_WECC_all_var_daily_clim(:,:,1:length(day_segment_for_this_year_inds),1,var_i);
        
    end
end

all_gridded_NH_data_daily_anoms = NaN(length(NH_lat_era5),...
                                      length(NH_lon_era5),...
                                      length(datetime_desc_time_daily),...
                                      length(all_gridded_NH_data_varnames_desc));

for var_i = 1:length(all_gridded_NH_data_varnames_desc)
    for year_i = 1:length(ann_desc_time_daily_numeric)-1

        day_segment_for_this_year_inds = find(isbetween(datetime_desc_time_daily,datetime(ann_desc_time_daily_numeric(year_i),01,01),datetime(ann_desc_time_daily_numeric(year_i+1),01,01)) == 1);
        
        all_gridded_NH_data_daily_anoms(:,:,day_segment_for_this_year_inds,var_i) = ...
            all_gridded_NH_data(:,:,day_segment_for_this_year_inds,var_i) - ...
            all_gridded_NH_data_daily_clim(:,:,1:length(day_segment_for_this_year_inds),1,var_i);        
    end
end

%% smooth gridded 

% gridded_WECC_all_var_daily_anoms_smooth = NaN(size(gridded_WECC_all_var_daily_anoms));
% 
% for var_i = 1:length(WECC_daily_var_names_plus_pop)-1
%     var_i/(length(WECC_daily_var_names_plus_pop)-1)
%     for lat_i = 1:length(WECC_lat)
%         for lon_i = 1:length(WECC_lon)
%         
%             gridded_WECC_all_var_daily_anoms_smooth(lat_i,lon_i,:,var_i) = smooth(squeeze(gridded_WECC_all_var_daily_anoms(lat_i,lon_i,:,var_i)),drought_timescale_days);
%             
%         end
%     end
% end
% 
% all_gridded_NH_data_daily_anoms_smooth = NaN(size(all_gridded_NH_data_daily_anoms));
% 
% for var_i = 1:length(all_gridded_NH_data_varnames_desc)
%     var_i/length(all_gridded_NH_data_varnames_desc)
%     for lat_i = 1:length(NH_lat_era5)
%         for lon_i = 1:length(NH_lon_era5)
%         
%             all_gridded_NH_data_daily_anoms_smooth(lat_i,lon_i,:,var_i) = smooth(squeeze(all_gridded_NH_data_daily_anoms_smooth(lat_i,lon_i,:,var_i)),drought_timescale_days);
%             
%         end
%     end
% end


%% find the X largest 7-day extremes on the non-deseased data

% size(daily_WECC_time_series_7_day_smooth)
% ans =
%        14610           6

% [B,I] = sort(___) also returns a collection of
% index vectors for any of the previous syntaxes.
% I is the same size as A and describes the arrangement
% of the elements of A into B along the sorted dimension.
% For example, if A is a vector, then B = A(I).

daily_WECC_time_series_7_day_smooth_rank = NaN(size(daily_WECC_time_series_7_day_smooth));

for var_i = 1:size(gridded_WECC_all_var_daily_reorder,4)
    
    this_var_t_series = daily_WECC_time_series_7_day_smooth(:,var_i);
    
    [B,I] = sort(this_var_t_series); %default is ascend so 1 is lowest
    
    daily_WECC_time_series_7_day_smooth_rank(:,var_i) = I;
    
end

%% pick out the extremes

initial_cases_to_pick = 90;

extreme_inds = NaN(initial_cases_to_pick,2,size(gridded_WECC_all_var_daily_reorder,4));

extreme_inds_date_low = {};
extreme_inds_date_high = {};

for var_i = 1:size(gridded_WECC_all_var_daily_reorder,4)
    
    WECC_daily_var_names_plus_pop{var_i}
    
    date_inds_small = daily_WECC_time_series_7_day_smooth_rank(1:initial_cases_to_pick,var_i);
    date_inds_large = daily_WECC_time_series_7_day_smooth_rank(end-initial_cases_to_pick+1:end,var_i);
    
    %eliminate border indicies (no indicies 1-3 or 14608)
        date_inds_small(date_inds_small <= 3) = NaN;
        date_inds_small(date_inds_small >= 14608) = NaN;
        date_inds_large(date_inds_large <= 3) = NaN;
        date_inds_large(date_inds_large >= 14608) = NaN;
    
    %go through inds eliminating adjacent ones
        for ind_i = 1:length(date_inds_small)

            ind_now = date_inds_small(ind_i);

            if isnan(ind_now) == 0

                distance_from_ind = date_inds_small - ind_now;

                adjacent_inds = find(abs(distance_from_ind) <= 7 & abs(distance_from_ind) > 0);

                date_inds_small(adjacent_inds) = NaN;

            end
        end
        for ind_i = 1:length(date_inds_large)

            ind_now = date_inds_large(ind_i);

            if isnan(ind_now) == 0

                distance_from_ind = date_inds_large - ind_now;

                adjacent_inds = find(abs(distance_from_ind) <= 7 & abs(distance_from_ind) > 0);

                date_inds_large(adjacent_inds) = NaN;

            end
        end
        
    extreme_inds(:,1,var_i) = date_inds_small;
    extreme_inds(:,2,var_i) = date_inds_large;
    
    datetime_desc_time_daily(date_inds_small(isnan(date_inds_small) == 0))'
    datetime_desc_time_daily(date_inds_large(isnan(date_inds_large) == 0))'
    
    extreme_inds_date_low{var_i} = datetime_desc_time_daily(date_inds_small(isnan(date_inds_small) == 0));
    extreme_inds_date_high{var_i} = datetime_desc_time_daily(date_inds_large(isnan(date_inds_large) == 0));
    
end

%% make anom map composites conditioned on wind, sun, people-degree-days

% size(gridded_WECC_all_var_daily_anoms)
%           31          36       14610           6
% size(all_gridded_NH_data_daily_anoms)
%           21         144       14610          19

% extremes are low for wind and sun, high for people-degree days

extreme_vars = {'surface solar power (W/m^2)',...
                'wind power (W/m^2)',...
                'people-degree-days'};
            
            
        variable_to_condition_on = 2;
            % size(extreme_inds)
            %     40     2     6
        var_to_cond_on_comp_center_time_inds_w_NaN = extreme_inds(:,1,variable_to_condition_on);
        %var_to_cond_on_comp_center_time_inds_w_NaN = extreme_inds(:,2,variable_to_condition_on);
            


var_to_cond_on_comp_center_time_inds = var_to_cond_on_comp_center_time_inds_w_NaN(isnan(var_to_cond_on_comp_center_time_inds_w_NaN) == 0);


gridded_WECC_all_var_daily_anoms_cond_drought_or_flood = NaN(length(WECC_lat),...
                                                         length(WECC_lon),...
                                                         drought_timescale_days,...
                                                         length(var_to_cond_on_comp_center_time_inds),...
                                                         length(WECC_daily_var_names_plus_pop)-1);
                                   
all_gridded_NH_data_daily_anoms_cond_drought_or_flood = NaN(length(NH_lat_era5),...
                                                        length(NH_lon_era5),...
                                                        drought_timescale_days,...
                                                        length(var_to_cond_on_comp_center_time_inds),...
                                                        length(all_gridded_NH_data_varnames_desc));
                                                
all_gridded_NH_data_daily_cond_drought_or_flood = NaN(length(NH_lat_era5),...
                                                  length(NH_lon_era5),...
                                                  drought_timescale_days,...
                                                  length(var_to_cond_on_comp_center_time_inds),...
                                                  length(all_gridded_NH_data_varnames_desc));
                                              
daily_WECC_cond_drought_or_flood = NaN(drought_timescale_days,...
                                               length(var_to_cond_on_comp_center_time_inds),...
                                               length(WECC_daily_var_names_plus_pop)-1);

% daily_WECC_time_series
                                                
for event_i = 1:length(var_to_cond_on_comp_center_time_inds)
    
    gridded_WECC_all_var_daily_anoms_cond_drought_or_flood(:,:,:,event_i,:) = gridded_WECC_all_var_daily_anoms(:,:,var_to_cond_on_comp_center_time_inds(event_i)-3:var_to_cond_on_comp_center_time_inds(event_i)+3,:);
    all_gridded_NH_data_daily_anoms_cond_drought_or_flood(:,:,:,event_i,:) = all_gridded_NH_data_daily_anoms(:,:,var_to_cond_on_comp_center_time_inds(event_i)-3:var_to_cond_on_comp_center_time_inds(event_i)+3,:);
    all_gridded_NH_data_daily_cond_drought_or_flood(:,:,:,event_i,:) = all_gridded_NH_data(:,:,var_to_cond_on_comp_center_time_inds(event_i)-3:var_to_cond_on_comp_center_time_inds(event_i)+3,:);
    
    daily_WECC_cond_drought_or_flood(:,event_i,:) = daily_WECC_time_series(var_to_cond_on_comp_center_time_inds(event_i)-3:var_to_cond_on_comp_center_time_inds(event_i)+3,:);
    
end
    
WECC_anom_composite_cond_drought_or_flood = squeeze(nanmean(nanmean(gridded_WECC_all_var_daily_anoms_cond_drought_or_flood,3),4));
NH_anom_composite_cond_drought_or_flood = squeeze(nanmean(nanmean(all_gridded_NH_data_daily_anoms_cond_drought_or_flood,3),4));
NH_composite_cond_drought_or_flood = squeeze(nanmean(nanmean(all_gridded_NH_data_daily_cond_drought_or_flood,3),4));

NH_normal_comp = nanmean(all_gridded_NH_data_daily_cond_drought_or_flood,3);

daily_WECC_cond_drought_or_flood_event_mean = squeeze(nanmean(daily_WECC_cond_drought_or_flood,1));

%% count fraction of anoms that of the same sign for stippling

frac_thresh_bottom = 0.2;
frac_thresh_top = 1-frac_thresh_bottom;

%average over event days
WECC_anom_comp_for_stipple = nanmean(gridded_WECC_all_var_daily_anoms_cond_drought_or_flood,3);

WECC_anom_stipple_YN = zeros(size(WECC_anom_comp_for_stipple,1),...
                             size(WECC_anom_comp_for_stipple,2),...
                             size(WECC_anom_comp_for_stipple,5));
                       
for lat_i = 1:size(WECC_anom_comp_for_stipple,1)
    for lon_i = 1:size(WECC_anom_comp_for_stipple,2)
        for var_i = 1:size(WECC_anom_comp_for_stipple,5)
            
            anom_string = squeeze(WECC_anom_comp_for_stipple(lat_i,lon_i,1,:,var_i));
            num_pos_anoms = length(find(anom_string > 0));
            frac_pos_anoms = num_pos_anoms/length(anom_string);
            
            if frac_pos_anoms > frac_thresh_top || frac_pos_anoms < frac_thresh_bottom
                
                WECC_anom_stipple_YN(lat_i,lon_i,var_i) = 1;
                
            end
        end
    end
end

NH_anom_comp_for_stipple = nanmean(all_gridded_NH_data_daily_anoms_cond_drought_or_flood,3);

NH_anom_stipple_YN = zeros(size(NH_anom_comp_for_stipple,1),...
                           size(NH_anom_comp_for_stipple,2),...
                           size(NH_anom_comp_for_stipple,5));
                     
for lat_i = 1:size(NH_anom_comp_for_stipple,1)
    for lon_i = 1:size(NH_anom_comp_for_stipple,2)
        for var_i = 1:size(NH_anom_comp_for_stipple,5)
            
            anom_string = squeeze(NH_anom_comp_for_stipple(lat_i,lon_i,1,:,var_i));
            num_pos_anoms = length(find(anom_string > 0));
            frac_pos_anoms = num_pos_anoms/length(anom_string);
            
            if frac_pos_anoms > frac_thresh_top || frac_pos_anoms < frac_thresh_bottom
                
                NH_anom_stipple_YN(lat_i,lon_i,var_i) = 1;
                
            end
        end
    end
end

%% plot composites for WECC (conditioned on extremes in whatever variable is picked above)

row_num = 3;
col_num = 2;

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',12);
set(0,'defaultAxesFontName', 'helvetica')

min_color_ranges = [-45 -175 -7 -7 -3 -5E5];
max_color_ranges = [ 45  175  7  7 3 5E5];

colormap(redblue)

for var_i = 1:length(WECC_daily_var_names_plus_pop)-1

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

    %contourfm(WECC_lat,WECC_lon,squeeze(WECC_composite_cond_low_wind(:,:,var_i)), 'LineStyle','none');
    pcolorm(WECC_lat,WECC_lon,squeeze(WECC_anom_composite_cond_drought_or_flood(:,:,var_i)), 'LineStyle','none');        
    caxis([min_color_ranges(var_i) max_color_ranges(var_i)])
    
    %stipple
        mask = squeeze(WECC_anom_stipple_YN(:,:,var_i))<1;
        
        [WECC_lon_mesh,WECC_lat_mesh] = meshgrid(WECC_lon,WECC_lat);
        [x_mesh_trans,y_mesh_trans] = mfwdtran(WECC_lat_mesh,WECC_lon_mesh);
        
        stipple(x_mesh_trans,y_mesh_trans,mask,'density',50)
    
    
    map_mean = threed2oned_lat_lon_2(squeeze(WECC_anom_composite_cond_drought_or_flood(:,:,var_i)),WECC_lat,WECC_lon);
    
%     if var_i <= 4; colormap(h,redblue); end
%     if var_i > 4; colormap(h,parula); end
    
    geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')

    t=contourcbar;
    %set(get(t,'ylabel'),'String', 'Fraction of days in season');
    
    title(strcat(WECC_daily_var_names_plus_pop{var_i},', map mean anom =',num2str(round(map_mean,2))));
    
end

%% plot composites for NH (conditioned on extremes in whatever variable is picked above)

[NH_lat_era5_grid,NH_lon_era5_grid] = meshgrid(NH_lat_era5,NH_lon_era5);

num_conts = 13;

NH_composite_cond_drought_or_flood = double(NH_composite_cond_drought_or_flood);
%wind_speeds_gridded_NH_data_clim = double(wind_speeds_gridded_NH_data_clim);

height_vars = {'250 hPa heights, wind, isotachs',...
               '500 hPa heights, wind, vorticity',...
               '850 hPa heights, wind, temperature'};
           
colorbar_labes = {'wind speed (m/s)',...
                  'vorticity (s^-^1)',...
                  'temperature (k)'};
           
row_num = 3;
col_num = 2;

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',22);
set(0,'defaultAxesFontName', 'helvetica')

plot_pos_1 = 1:2:6;

for var_i = 1:length(height_vars)

    h = subplot(row_num,col_num,plot_pos_1(var_i));
    hold on

    axesm('robinson',...
    'Frame', 'on',...
    'Grid', 'off',...
    'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
    'maplonlim',[110 340])
   tightmap

    %if var_i == 4; colormap(redblue); end
%      min_color_range = 0.0;
%      max_color_range = 0.35;

    if var_i == 3
        
     contourfm(NH_lat_era5,NH_lon_era5,squeeze(NH_composite_cond_drought_or_flood(:,:,4)),40,'Linestyle','none')
     
     caxis([min(min(squeeze(NH_composite_cond_drought_or_flood(:,:,4)))) 1*max(max(squeeze(NH_composite_cond_drought_or_flood(:,:,4))))])
     t=contourcbar;
        
    end
    
    if var_i == 2
        
     contourfm(NH_lat_era5,NH_lon_era5,squeeze(NH_composite_cond_drought_or_flood(:,:,11)),80,'Linestyle','none')
     
     caxis([-0.5*max(max(abs(squeeze(NH_composite_cond_drought_or_flood(:,:,11))))) 0.5*max(max(abs(squeeze(NH_composite_cond_drought_or_flood(:,:,11)))))])
     t=contourcbar;
        
    end
    if var_i == 1
        
     contourfm(NH_lat_era5,NH_lon_era5,squeeze(NH_composite_cond_drought_or_flood(:,:,12)),80,'Linestyle','none')
     
     caxis([min(min(squeeze(NH_composite_cond_drought_or_flood(:,:,12)))) 1*max(max(squeeze(NH_composite_cond_drought_or_flood(:,:,12))))])
     t=contourcbar;
        
    end
    
    colormap(h,parula)

    contourm(NH_lat_era5,NH_lon_era5,squeeze(NH_composite_cond_drought_or_flood(:,:,var_i)),num_conts,'LineColor','r')    
    quiverm(NH_lat_era5_grid,NH_lon_era5_grid,squeeze(NH_composite_cond_drought_or_flood(:,:,var_i+7))',squeeze(NH_composite_cond_drought_or_flood(:,:,var_i+4))')

    %pcolorm(NH_lat_era5,NH_lon_era5,squeeze(NH_composite_cond_low_wind(:,:,1,var_i)), 'LineStyle','none');        
    %caxis([min_color_ranges(var_i) max_color_ranges(var_i)])
    
    geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')

    set(get(t,'ylabel'),'String', colorbar_labes{var_i});
    
    title(height_vars{var_i})
    
end

plot_pos_2 = 2:2:6;

height_vars_2 = {'250 hPa heights & height anomalies',...
                 '500 hPa heights & height anomalies',...
                 '850 hPa heights & height anomalies'};
             
colorbar_labes = {'height anomaly (meters)',...
                  'height anomaly (meters)',...
                  'height anomaly (meters)'};

for var_i = 1:length(height_vars_2)

    h = subplot(row_num,col_num,plot_pos_2(var_i));
    hold on
    
    axesm('robinson',...
    'Frame', 'on',...
    'Grid', 'off',...
    'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
    'maplonlim',[110 340])
   tightmap
        
    contourfm(NH_lat_era5,NH_lon_era5,squeeze(NH_anom_composite_cond_drought_or_flood(:,:,var_i)),40,'Linestyle','none')
     
     %caxis([min(min(squeeze(NH_anom_composite_cond_low_wind(:,:,var_i)))) 1*max(max(squeeze(NH_anom_composite_cond_low_wind(:,:,var_i))))])
     caxis([-1*max(max(abs(squeeze(NH_anom_composite_cond_drought_or_flood(:,:,var_i))))) 1*max(max(abs(squeeze(NH_anom_composite_cond_drought_or_flood(:,:,var_i)))))])

     t=contourcbar;
     
     colormap(h,redblue)
     
     %stipple
        mask = squeeze(NH_anom_stipple_YN(:,:,var_i))<1;

        [NH_lon_mesh,NH_lat_mesh] = meshgrid(NH_lon_era5,NH_lat_era5);
        [x_mesh_trans,y_mesh_trans] = mfwdtran(NH_lat_mesh,NH_lon_mesh);

        stipple(x_mesh_trans,y_mesh_trans,mask,'density',100)
     
    contourm(NH_lat_era5,NH_lon_era5,squeeze(NH_composite_cond_drought_or_flood(:,:,var_i)),num_conts,'LineColor','r')

    %pcolorm(NH_lat_era5,NH_lon_era5,squeeze(NH_composite_cond_low_wind(:,:,1,var_i)), 'LineStyle','none');        
    %caxis([min_color_ranges(var_i) max_color_ranges(var_i)])
    
    geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')

    set(get(t,'ylabel'),'String', colorbar_labes{var_i});
    
    title(height_vars_2{var_i})
    
end

%% plot the individual events ranked

% all_gridded_NH_data_varnames_desc'
%     {'250 hPa heights'              } 1
%     {'500 hPa heights'              } 2
%     {'850 hPa heights'              } 3
%     {'850 hpa temperature'          } 4
%     {'250 hPa U wind'               } 5
%     {'500 hpa U wind'               } 6
%     {'850 hpa U wind'               } 7
%     {'250 hPa V wind'               } 8
%     {'500 hpa V wind'               } 9
%     {'850 hpa V wind'               } 10
%     {'500 hpa vorticity'            } 11
%     {'250 hPa wind speed'           } 12
%     {'500 hpa wind speed'           } 13
%     {'850 hpa wind speed'           } 14
%     {'250 hpa divergence'           } 15
%     {'500 hpa divergence'           } 16
%     {'850 hpa divergence'           } 17
%     {'500 hpa vorticity advection'  } 18
%     {'850 hpa temperature advection'} 19

anom_var_to_plot = 2;

row_num = 5;
col_num = 3;

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',6);
set(0,'defaultAxesFontName', 'helvetica')

    var_mean_whole_time = squeeze(mean(seas_WECC_time_series_day_year,1));

% size(NH_anom_comp_for_stipple)
%     21   144     1    16    19

    date_set = extreme_inds_date_low{variable_to_condition_on};
    %date_set = extreme_inds_date_high{6};

for event_i = 1:row_num*col_num
    
    h=subplot(row_num,col_num,event_i);
    hold on
    
    axesm('robinson',...
    'Frame', 'on',...
    'Grid', 'off',...
    'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
    'maplonlim',[110 340])
     tightmap
        
    contourfm(NH_lat_era5,NH_lon_era5,squeeze(NH_anom_comp_for_stipple(:,:,1,event_i,anom_var_to_plot)),40,'Linestyle','none')
     
     %caxis([min(min(squeeze(NH_anom_composite_cond_low_wind(:,:,var_i)))) 1*max(max(squeeze(NH_anom_composite_cond_low_wind(:,:,var_i))))])
     caxis([-1*max(max(abs(squeeze(NH_anom_comp_for_stipple(:,:,1,event_i,anom_var_to_plot))))) 1*max(max(abs(squeeze(NH_anom_comp_for_stipple(:,:,1,event_i,anom_var_to_plot)))))])

     %t=contourcbar;
     
     colormap(h,redblue)
     
    contourm(NH_lat_era5,NH_lon_era5,squeeze(NH_normal_comp(:,:,1,event_i,anom_var_to_plot)),num_conts,'LineColor','r')

    %pcolorm(NH_lat_era5,NH_lon_era5,squeeze(NH_composite_cond_low_wind(:,:,1,var_i)), 'LineStyle','none');        
    %caxis([min_color_ranges(var_i) max_color_ranges(var_i)])
    
    geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')

    %set(get(t,'ylabel'),'String', colorbar_labes{var_i});
        
%     map_mean_perc_dev = 100.*...
%                         (daily_WECC_cond_drought_or_flood_event_mean(event_i,variable_to_condition_on) - var_mean_whole_time(variable_to_condition_on))./...
%                         var_mean_whole_time(variable_to_condition_on);
                    
    map_mean_perc_dev = 100.*...
                        (daily_WECC_cond_drought_or_flood_event_mean(event_i,6) - var_mean_whole_time(6))./...
                        var_mean_whole_time(6);
                    
    [y,m,d] = ymd(date_set(event_i));

    %title(strcat('WECC anom= ',num2str(map_mean_wind_solar_temp)))
    title(strcat(num2str(m),'/',num2str(d),'/',num2str(y),', WECC percent deviation= ',num2str(round(map_mean_perc_dev,1))))
    
    linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-k')
    
    %note that the plots dont appear to be in order of largest droguht
    %because the seasonal cycle is subtracted from the value displayed in
    %the the title (it is anomaly) while the largest droughts are
    %conditioned on absolute values
    
end

