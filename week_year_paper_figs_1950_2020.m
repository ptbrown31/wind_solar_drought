close all
clear all

addpath(genpath('/Users/patrickbrown/Documents/MATLAB/'))
addpath(genpath('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/'))


%% load SSTs

load('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/David_ERA5_data/read_save_ERA5_monthly_SST_data_1950_2020.mat',...
    'Global_SSTs_week_year_anom',...
    'Global_SSTs_week_year',...
    'Global_lat_era5',...
    'Global_lon_era5',...
    'desc_years',...
    'desc_weeks')

%can coursify SSTs if neccesary

%% load era-5 data

load('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/wind_solar_droughts/David_ERA5_data/read_org_week_year_1950_2020.mat',...
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

% 1  '250 hPa heights',...
% 2  '500 hPa heights',...
% 3  '700 hPa heights',...
% 4  '250 hPa U wind',...
% 5  '500 hpa U wind',...
% 6   '700 hpa U wind',...
% 7   '250 hPa V wind',...
% 8   '500 hpa V wind',...
% 9   '700 hpa V wind',...
% 10  '500 hpa vorticity',...
% 11  '700 hpa temperature',...
% 12  '700 Omega',...
% 13  'sea level pressure',...
% 14  '250 hPa wind speed',...
% 15  '500 hpa wind speed',...
% 16  '700 hpa wind speed',...
% 17  '250 hpa divergence',...
% 18  '500 hpa divergence',...
% 19  '700 hpa divergence',...
% 20  '500 hpa vorticity advection',...
% 21  '700 hpa temperature advection'};

    %peel off temperature deviations to keep them seperate
    
%     gridded_WECC_week_year_td = gridded_WECC_week_year(:,:,:,:,5);
%     gridded_WECC_week_year_anoms_td = gridded_WECC_week_year_anoms(:,:,:,:,5);
%     gridded_WECC_week_year_frac_devs_wrt_ann_mean_td = gridded_WECC_week_year_frac_devs_wrt_ann_mean(:,:,:,:,5);
%     gridded_WECC_week_year_frac_devs_wrt_woy_mean_td = gridded_WECC_week_year_frac_devs_wrt_woy_mean(:,:,:,:,5);
%     WECC_daily_var_names = WECC_daily_var_names{5};
%     
%     gridded_WECC_week_year = gridded_WECC_week_year(:,:,:,:,1:4);
%     gridded_WECC_week_year_anoms = gridded_WECC_week_year_anoms(:,:,:,:,1:4);
%     gridded_WECC_week_year_frac_devs_wrt_ann_mean = gridded_WECC_week_year_frac_devs_wrt_ann_mean(:,:,:,:,1:4);
%     gridded_WECC_week_year_frac_devs_wrt_woy_mean = gridded_WECC_week_year_frac_devs_wrt_woy_mean(:,:,:,:,1:4);
%     WECC_daily_var_names = WECC_daily_var_names{1:4};
        
%% option to NaN-out HDD and CDD areas without enough variation

gridded_WECC_week_year_frac_devs_wrt_ann_mean(gridded_WECC_week_year_frac_devs_wrt_ann_mean>10) = NaN;
gridded_WECC_week_year_frac_devs_wrt_woy_mean(gridded_WECC_week_year_frac_devs_wrt_woy_mean>10) = NaN;

%% define colors to use throughout

colormap_strings = {'wm',...   % wind raw
                    'kwm',...  % wind anomaly
                    'ky',...   % solar raw
                    'kwy',...  % solar anomaly
                    'wr',...   % CDD raw
                    'bwr',...  % CDD anomaly
                    'wb',...   % HDD raw
                    'rwb',...  % HDD anomaly
                    'bwr',...  % temp demperature anomalies
                    'cwg'};    % height anomalies
                
colormap_brewer_string = {'greens',...    % 1 wind raw
                          'PRGn',...      % 2 wind anomaly
                          '*YlOrRd',...   % 3 solar raw
                          'RdGy',...      % 4 solar anomaly
                          'wr',...        % 5 CDD raw
                          'bwr',...       % 6 CDD anomaly
                          'wb',...        % 7 HDD raw
                          'rwb',...       % 8 HDD anomaly
                          'bwr',...       % 9 temp demperature anomalies
                          '*BrBG',...     % 10 height anomalies
                          'PiYG',...      % 11 vorticity
                          '*PuOr'};       % 12 divergence
                
% [cmap]=buildcmap(colormap_strings{lin_count_i});
% colormap(h,cmap)
    
%% make WECC Clim Maps
% size(gridded_WECC_week_year)
%     31    36    52    40     5

var_names_2 = {'100m wind power', 'surface solar power','cooling degree-days','heating degree-days'};
var_units_2 = {'proportion of annual domain mean (%)', 'proportion of annual domain mean (%)','C*days','C*days'};

coast = load('coastlines');
WECC_lat = double(WECC_lat);
WECC_lon = double(WECC_lon);

summer_JJA_weeks = 22:35;
winter_DJF_weeks = horzcat(1:8,48:52);

        gridded_WECC_week_year_summer_JJA = squeeze(nanmean(nanmean(gridded_WECC_week_year(:,:,summer_JJA_weeks,:,:),4),3));
        gridded_WECC_week_year_winter_DJF = squeeze(nanmean(nanmean(gridded_WECC_week_year(:,:,winter_DJF_weeks,:,:),4),3));
        gridded_WECC_week_year_ann = squeeze(nanmean(nanmean(gridded_WECC_week_year,4),3));

        space_mean_wind = threed2oned_lat_lon_2(squeeze(gridded_WECC_week_year_ann(:,:,1)),WECC_lat,WECC_lon);
        space_mean_sun = threed2oned_lat_lon_2(squeeze(gridded_WECC_week_year_ann(:,:,2)),WECC_lat,WECC_lon);

        gridded_WECC_week_year_summer_JJA(:,:,1) = 100*gridded_WECC_week_year_summer_JJA(:,:,1)./space_mean_wind;
        gridded_WECC_week_year_winter_DJF(:,:,1) = 100*gridded_WECC_week_year_winter_DJF(:,:,1)./space_mean_wind;
        gridded_WECC_week_year_ann(:,:,1) = 100*gridded_WECC_week_year_ann(:,:,1)./space_mean_wind;
        
        gridded_WECC_week_year_summer_JJA(:,:,2) = 100*gridded_WECC_week_year_summer_JJA(:,:,2)./space_mean_sun;
        gridded_WECC_week_year_winter_DJF(:,:,2) = 100*gridded_WECC_week_year_winter_DJF(:,:,2)./space_mean_sun;
        gridded_WECC_week_year_ann(:,:,2) = 100*gridded_WECC_week_year_ann(:,:,2)./space_mean_sun;
        
plot_array = cat(4,gridded_WECC_week_year_winter_DJF,gridded_WECC_week_year_summer_JJA,gridded_WECC_week_year_ann);

plot_inds = [1 4 7 10 2 5 8 11 3 6 9 12];

    min_color_ranges = [0 0 0 0];
    %max_color_ranges = [4 288 30 30];
    max_color_ranges = [350 200 30 30];

    row_num = 4;
    col_num = 3;

    FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',10);
    set(0,'defaultAxesFontName', 'helvetica')
    
    plot_i = 1;

for seas_i = 1:3
    for var_i = 1:length(WECC_daily_var_names)-1
        
        if var_i == 1; cmap_ind = 1; end
        if var_i == 2; cmap_ind = 3; end
        if var_i == 3; cmap_ind = 5; end
        if var_i == 4; cmap_ind = 7; end
        
        h = subplot(row_num,col_num,plot_inds(plot_i));
        hold on
        
        [cmap]=buildcmap(colormap_strings{cmap_ind});

        axesm('robinson',...
        'Frame', 'on',...
        'Grid', 'off',...
        'maplatlim',[min(WECC_lat) max(WECC_lat)],...
        'maplonlim',[min(WECC_lon) max(WECC_lon)])
        tightmap

        levels = linspace(min_color_ranges(var_i),max_color_ranges(var_i),50);

        contourfm(WECC_lat,WECC_lon,squeeze(plot_array(:,:,var_i,seas_i)),levels, 'LineStyle','none');
        %pcolorm(WECC_lat,WECC_lon,squeeze(plot_array(:,:,var_i,seas_i)), 'LineStyle','none');        
        caxis([min_color_ranges(var_i) max_color_ranges(var_i)])

%         if var_i == 4; colormap(h,redblue); end

        geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
        
        if seas_i == 3
            
            t=colorbar;
        
            set(get(t,'ylabel'),'String', var_units_2{var_i});
            
        end
        
        map_mean = threed2oned_lat_lon_2(squeeze(plot_array(:,:,var_i,seas_i)),WECC_lat,WECC_lon);

       % title(strcat(var_names_2{var_i},', mean=',num2str(round(map_mean,0))))
       
        if var_i == 3 || var_i == 4
            title(strcat('mean = ',num2str(round(map_mean,1)),{' '},var_units_2{var_i}))
        end
%         if var_i == 2
%             title(strcat('mean = ',num2str(round(map_mean,1)),{' '},var_units_2{var_i}))
%         end        
        if var_i == 1
            colormap(h,brewermap([],colormap_brewer_string{cmap_ind}))
        end
        if var_i == 2 || var_i == 3 || var_i == 4
            colormap(h,cmap)
        end
        
        plot_i = plot_i + 1;

    end
end
 
%% make 1-D WECC data and 2d scatters

WECC_week_year_1d = NaN(length(weeks),length(years),length(WECC_daily_var_names));

for var_i = 1:length(WECC_daily_var_names)
    for year_i = 1:length(years)
        for week_i = 1:length(weeks)

            WECC_week_year_1d(week_i,year_i,var_i) = threed2oned_lat_lon_2(squeeze(gridded_WECC_week_year(:,:,week_i,year_i,var_i)),WECC_lat,WECC_lon);
            
        end
    end
end

%% calc WECC overall means

WECC_overall_means = nanmean(nanmean(WECC_week_year_1d,2),1);

%% normalize the 1D

WECC_week_year_1d_norm = NaN(size(WECC_week_year_1d));

for var_i = 1:length(WECC_daily_var_names)
    
    WECC_week_year_1d_norm(:,:,var_i) = WECC_week_year_1d(:,:,var_i)./(nanmean(nanmean(WECC_week_year_1d(:,:,var_i),2),1));
    
end

%% calculate overall variability & that explained by predictible seasonal cycle

TwoDArray = WECC_week_year_1d(:,:,1);
wind_raw_1d = reshape(TwoDArray',[1 size(TwoDArray,1)*size(TwoDArray,2)]);

wind_overall_variability = std(wind_raw_1d);

TwoDArray = WECC_week_year_1d(:,:,2);
solar_raw_1d = reshape(TwoDArray',[1 size(TwoDArray,1)*size(TwoDArray,2)]);

solar_overall_variability = std(solar_raw_1d);

solar_to_wind_variability_ratio = solar_overall_variability./wind_overall_variability

%% climatology

WECC_week_year_1d_norm_mean = squeeze(mean(WECC_week_year_1d_norm,2));
WECC_week_year_1d_norm_std = squeeze(std(WECC_week_year_1d_norm,[],2));

%% make 1D for scatters

    OneDArray_wind_norm = NaN(1,length(years)*length(weeks));
        week_total_i = 1;
        for year_i = 1:length(years)
            for week_i = 1:length(weeks)
                OneDArray_wind_norm(week_total_i) = WECC_week_year_1d_norm(week_i,year_i,1);
                week_total_i = week_total_i + 1;
            end
        end

    OneDArray_insolation_norm = NaN(1,length(years)*length(weeks));
        week_total_i = 1;
        for year_i = 1:length(years)
            for week_i = 1:length(weeks)
                OneDArray_insolation_norm(week_total_i) = WECC_week_year_1d_norm(week_i,year_i,2);
                week_total_i = week_total_i + 1;
            end
        end
    
    OneDArray_CDD_norm = NaN(1,length(years)*length(weeks));
        week_total_i = 1;
        for year_i = 1:length(years)
            for week_i = 1:length(weeks)
                OneDArray_CDD_norm(week_total_i) = WECC_week_year_1d_norm(week_i,year_i,3);
                week_total_i = week_total_i + 1;
            end
        end

    OneDArray_HDD_norm = NaN(1,length(years)*length(weeks));
        week_total_i = 1;
        for year_i = 1:length(years)
            for week_i = 1:length(weeks)
                OneDArray_HDD_norm(week_total_i) = WECC_week_year_1d_norm(week_i,year_i,4);
                week_total_i = week_total_i + 1;
            end
        end
        
solar_to_wind_variability_normed_ratio = std(OneDArray_insolation_norm)./std(OneDArray_wind_norm)

%% calc fraction of varience explained by seasonal cycle

%make continous repeating seasonal-cycle time series

repeating_wind_seasonal_cycle = NaN(size(OneDArray_wind_norm));
repeating_solar_seasonal_cycle = NaN(size(OneDArray_insolation_norm));

week_of_year_i = 1;

for week_i = 1:length(repeating_wind_seasonal_cycle)
    
    repeating_wind_seasonal_cycle(week_i) = WECC_week_year_1d_norm_mean(week_of_year_i,1);
    repeating_solar_seasonal_cycle(week_i) = WECC_week_year_1d_norm_mean(week_of_year_i,2);
    
    week_of_year_i = week_of_year_i + 1;
    
    if week_of_year_i == 53; week_of_year_i = 1; end
    
end

FigHandle = figure('Position', [100, 100, 1000, 1000]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',18);
set(0,'defaultAxesFontName', 'helvetica')
hold on

scatter(repeating_solar_seasonal_cycle,OneDArray_insolation_norm,50,'r','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
% h1 = lsline;
% s = h1.Color;
% h1.Color = [1 0 0];

scatter(repeating_wind_seasonal_cycle,OneDArray_wind_norm,50,'b','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
% h2 = lsline;
% s = h2.Color;
% h2.Color = [0 0 1];

%legend([h1 h2],'solar power','wind power')

lsline

title('varience explained by the seasonal cycle')
xlabel('typical value for that week-of-year')
ylabel('actual weekly value')

corr_mat = corrcoef(repeating_wind_seasonal_cycle,...
                    OneDArray_wind_norm,...
                    'rows',...
                    'pairwise');

r_sqrd_wind_week_with_mean_week = corr_mat(2,1).^2;
        
corr_mat = corrcoef(repeating_solar_seasonal_cycle,...
                    OneDArray_insolation_norm,...
                    'rows',...
                    'pairwise');

r_sqrd_sun_week_with_mean_week = corr_mat(2,1).^2;

%% do rankings

num_tot_weeks = (40*365.25)/7;

% = 2087 weeks
% 1st percentile = 20.87 weeks
% use 20 weeks for round numbers - composite over the 20 weeks

% most extreme in terms of each individually

% [B,I] = sort(___) also returns a collection of index vectors for any of the previous syntaxes. 
% I is the same size as A and describes the arrangement of the elements of A into B along the sorted dimension. 
% For example, if A is a vector, then B = A(I).

% F =
%     4     3     2     7     4    10

% [b,I] = sort(F)

% b =
%      2     3     4     4     7    10
% I =
%      3     2     1     5     4     6

[ranked_value_wind, ranked_index_wind] = sort(OneDArray_wind_norm);
[ranked_value_insolation, ranked_index_insolation] = sort(OneDArray_insolation_norm);

    % wind + sun
    OneDArray_wind_insolation_norm = OneDArray_wind_norm + OneDArray_insolation_norm;
    [ranked_value_wind_insolation, ranked_index_wind_insolation] = sort(OneDArray_wind_insolation_norm);

%% re-wrap rankings into week-time so that you can pull from the other arrays

week_year_2d_week_number_key = NaN(length(weeks),length(years));

week_total_i = 1;
for year_i = 1:length(years)
    for week_i = 1:length(weeks)
        week_year_2d_week_number_key(week_i,year_i,1) = week_total_i;
        week_total_i = week_total_i + 1;
    end
end
    
%% make conditionals

% round(0.01*length(OneDArray_wind_norm))
% ans =
%     37

% round(0.005*length(OneDArray_wind_norm))
% ans =
%     18

% round(0.02*length(OneDArray_wind_norm))
% ans =
%     74

num_weeks_for_drought = 37;

mean_wind_given_solar_drought = mean(OneDArray_wind_norm(ranked_index_insolation(1:num_weeks_for_drought)));
mean_wind_given_wind_drought = mean(OneDArray_wind_norm(ranked_index_wind(1:num_weeks_for_drought)));

mean_insolation_given_solar_drought = mean(OneDArray_insolation_norm(ranked_index_insolation(1:num_weeks_for_drought)));
mean_insolation_given_wind_drought = mean(OneDArray_insolation_norm(ranked_index_wind(1:num_weeks_for_drought)));

mean_wind_given_wind_solar_drought = mean(OneDArray_wind_norm(ranked_index_wind_insolation(1:num_weeks_for_drought)));
mean_insolation_given_wind_solar_drought = mean(OneDArray_insolation_norm(ranked_index_wind_insolation(1:num_weeks_for_drought)));

%% make and plot annual drought frequency over time

%make arrays of weeks where a week is a 1 if its a drought and 0 otherwise

label_wind_drought_weeks = zeros(length(ranked_index_insolation),1);
label_solar_drought_weeks = zeros(length(ranked_index_insolation),1);
label_wind_solar_drought_weeks = zeros(length(ranked_index_insolation),1);

label_wind_drought_weeks(ranked_index_wind(1:num_weeks_for_drought),1) = 1;
label_solar_drought_weeks(ranked_index_insolation(1:num_weeks_for_drought),1) = 1;
label_wind_solar_drought_weeks(ranked_index_wind_insolation(1:num_weeks_for_drought),1) = 1;

% wind
    num_ann_wind_droughts = NaN(length(years),1);

    week_start_now = 1;
    week_end_now = week_start_now + 52 -1;

    for year_i = 1:length(years)

        num_ann_wind_droughts(year_i) = sum(label_wind_drought_weeks(week_start_now:week_end_now));

        week_start_now = week_start_now + 52;
        week_end_now = week_end_now + 52;
        
    end
    
% sun
    num_ann_solar_droughts = NaN(length(years),1);

    week_start_now = 1;
    week_end_now = week_start_now + 52 -1;

    for year_i = 1:length(years)-1

        num_ann_solar_droughts(year_i) = sum(label_solar_drought_weeks(week_start_now:week_end_now));

        week_start_now = week_start_now + 52;
        week_end_now = week_end_now + 52;

    end
    
% wind+sun
    num_ann_wind_sun_droughts = NaN(length(years),1);

    week_start_now = 1;
    week_end_now = week_start_now + 52 -1;

    for year_i = 1:length(years)-1

        num_ann_wind_sun_droughts(year_i) = sum(label_wind_solar_drought_weeks(week_start_now:week_end_now));

        week_start_now = week_start_now + 52;
        week_end_now = week_end_now + 52;

    end
    
    
% plot

    FigHandle = figure('Position', [100, 100, 1500, 1000]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',16);
    set(0,'defaultAxesFontName', 'helvetica')

    subplot(3,1,1);
    hold on
    
        bar(years,num_ann_wind_droughts)

        %xlabel('years')
        ylabel('annual # of wind droughts')
        
    subplot(3,1,2);
    hold on
    
        bar(years,num_ann_solar_droughts)

        %xlabel('years')
        ylabel('annual # of solar droughts')
        
    subplot(3,1,3);
    hold on
    
        bar(years,num_ann_wind_sun_droughts)

        xlabel('years')
        ylabel('annual # of wind+solar droughts')     


%% plot scatters

xmin = 0;
ymin = 0;
xmax = 200;
ymax = 200;

FigHandle = figure('Position', [100, 100, 1500, 1000]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',20);
set(0,'defaultAxesFontName', 'helvetica')
hold on

    [cmap]=buildcmap('bwrr');
    %[cmap]=buildcmap('rwbb');
    
    % plot all weeks
    scatter(100*OneDArray_wind_norm,...
            100*OneDArray_insolation_norm,...
            100,...
            100*OneDArray_CDD_norm,...
            'filled','MarkerFaceAlpha',.5,'MarkerEdgeColor','k','MarkerEdgeAlpha',.2)
    
    % plot drought weeks
    
        %plot solar droughts
        scatter(100*OneDArray_wind_norm(ranked_index_insolation(1:num_weeks_for_drought)),...
                100*OneDArray_insolation_norm(ranked_index_insolation(1:num_weeks_for_drought)),...
                200,100*OneDArray_CDD_norm(ranked_index_insolation(1:num_weeks_for_drought)),...
                'filled','s','MarkerFaceAlpha',.5,'MarkerEdgeColor','k','MarkerEdgeAlpha',1)
        %plot solar droughts
        scatter(100*OneDArray_wind_norm(ranked_index_wind(1:num_weeks_for_drought)),...
                100*OneDArray_insolation_norm(ranked_index_wind(1:num_weeks_for_drought)),...
                200,100*OneDArray_CDD_norm(ranked_index_wind(1:num_weeks_for_drought)),...
                'filled','o','MarkerFaceAlpha',.5,'MarkerEdgeColor','k','MarkerEdgeAlpha',1)
            
        %plot wind+solar droughts
        scatter(100*OneDArray_wind_norm(ranked_index_wind_insolation(1:num_weeks_for_drought)),...
                100*OneDArray_insolation_norm(ranked_index_wind_insolation(1:num_weeks_for_drought)),...
                200,100*OneDArray_CDD_norm(ranked_index_wind_insolation(1:num_weeks_for_drought)),...
                'filled','d','MarkerFaceAlpha',.5,'MarkerEdgeColor','k','MarkerEdgeAlpha',1)
    
    % plot seasonal cycle
    plot(100*WECC_week_year_1d_norm_mean(:,1),100*WECC_week_year_1d_norm_mean(:,2),'-k','LineWidth',2)
    scatter(100*WECC_week_year_1d_norm_mean(:,1),100*WECC_week_year_1d_norm_mean(:,2),100,'k','filled','MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',1)
    
%     for week_i = 1:4:48
%         text(100*WECC_week_year_1d_norm_mean(week_i,1)+2,100*WECC_week_year_1d_norm_mean(week_i,2)+2,num2str(week_i),'Color','black','FontWeight','bold','FontSize',16)
%     end
    
    % plot condistions and wind and solar droughts
    scatter(100*mean_wind_given_wind_drought,100*mean_insolation_given_wind_drought,600,'k','filled','o','MarkerFaceAlpha',.7,'MarkerEdgeColor','k','MarkerEdgeAlpha',1)
    scatter(100*mean_wind_given_solar_drought,100*mean_insolation_given_solar_drought,600,'k','filled','s','MarkerFaceAlpha',.7,'MarkerEdgeColor','k','MarkerEdgeAlpha',1)
    scatter(100*mean_wind_given_wind_solar_drought,100*mean_insolation_given_wind_solar_drought,600,'k','filled','d','MarkerFaceAlpha',.7,'MarkerEdgeColor','k','MarkerEdgeAlpha',1)
    
    % plot conditional on CDD
    % scatter(100*mean_wind_given_high_CDD,100*mean_insolation_given_high_CDD,400,'r','filled','MarkerFaceAlpha',.7,'MarkerEdgeColor','k','MarkerEdgeAlpha',1)

    t=contourcbar;
    set(get(t,'ylabel'),'String', 'HDD proportion of long-term mean (%)');
    caxis([0 300])

    plot([100 100],[ymin ymax],'-k')
    plot([xmin xmax],[100 100],'-k')

    xlim([xmin xmax])
    ylim([xmin xmax])
    
    xlabel('proportion of long-term mean wind power (%)')
    ylabel('proportion of long-term mean solar power (%)')
    
    %title('All weeks 1979-2018')
    
    colormap(cmap)

% output count of weeks in each quandrant 

%% pull the lowest 37 weeks (1st percentile) for each variable individually

% % size(gridded_WECC_week_year)
% % ans =
% %     31    36    52    40     4
% 
% %WECC_week_year_1d_norm
% 
weeks_to_pull = num_weeks_for_drought;
% num_conts = 3;
% 
% coast = load('coastlines');
% %colormap(redblue)
% 
% % min_colors = [-80 -80];
% % max_colors = [80 80];
% 
% min_colors = [0 0];
% max_colors = [200 200];
% 
% for var_i = 1:2
%     
%     if var_i == 1
%         
%         [cmap]=buildcmap(colormap_strings{2});
%     
%     end
%     if var_i == 2
%         
%         [cmap]=buildcmap(colormap_strings{4});
%     
%     end
%     
%     FigHandle = figure('Position', [100, 100, 1300, 1000]); %[left bottom width height]
%     set(gcf,'color',[1 1 1]);
%     set(0, 'DefaultAxesFontSize',3);
%     set(0,'defaultAxesFontName', 'helvetica')
%     hold on
% 
%     for week_i = 1:weeks_to_pull
% 
%         if var_i == 1; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind(week_i)); end
%         if var_i == 2; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_insolation(week_i)); end
% 
%         h = subplot(10,4,week_i);
%         hold on
% 
%         axesm('robinson',...
%         'Frame', 'on',...
%         'Grid', 'off',...
%         'maplatlim',[min(WECC_lat) max(WECC_lat)],...
%         'maplonlim',[min(WECC_lon) max(WECC_lon)])
%         tightmap
%         
%         levels = linspace(min_colors(var_i),max_colors(var_i),50);
%                 
%         contourfm(WECC_lat,WECC_lon,100*squeeze(gridded_WECC_week_year_frac_devs_wrt_woy_mean(:,:,week_ind,year_ind,var_i))+100,levels,'LineStyle','none');
%         %pcolorm(WECC_lat,WECC_lon,squeeze(gridded_WECC_week_year_anoms(:,:,week_ind,year_ind,var_i)), 'LineStyle','none');
%         
%         if var_i == 1; contourm(WECC_lat,WECC_lon,gridded_WECC_week_year_ann(:,:,var_i),num_conts,'LineColor','w','LineWidth',1); end
%         if var_i == 2; contourm(WECC_lat,WECC_lon,gridded_WECC_week_year_ann(:,:,var_i),num_conts,'LineColor','k','LineWidth',1); end
% 
%         map_mean_anom =100*threed2oned_lat_lon_2(squeeze(gridded_WECC_week_year_frac_devs_wrt_woy_mean(:,:,week_ind,year_ind,var_i)),WECC_lat,WECC_lon)+100;
%         map_mean_raw = 100*WECC_week_year_1d_norm(week_ind,year_ind,var_i);
%         %map_mean_raw = 100*threed2oned_lat_lon_2(squeeze(gridded_WECC_week_year_frac_devs_wrt_ann_mean(:,:,week_ind,year_ind,var_i)),WECC_lat,WECC_lon);
%         
%              geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
%         
%         %t=colorbar;
%         %set(get(t,'ylabel'),'String', '%');
%         
%         caxis([min_colors(var_i) max_colors(var_i)])
% 
%         title(strcat(num2str(week_i),')',{' '},num2str(years(year_ind)),',',{' '},num2str(week_ind),{' '},'(',num2str(round(map_mean_raw,0)),'%,',num2str(round(map_mean_anom,0)),'%)'))
%         
%        
%         if var_i == 1
%             colormap(h,brewermap([],colormap_brewer_string{2}))
%         end
%         if var_i == 2
%             %colormap(h,brewermap([],colormap_brewer_string{4}))
%             colormap(h,cmap)
%         end       
%        
%        
%     end
% end

%% pull the lowest 37 weeks (1st percentile) for each variable individually

% % size(gridded_WECC_week_year)
% % ans =
% %     31    36    52    40     4
% 
% %WECC_week_year_1d_norm
% 
% weeks_to_pull = 37;
% num_conts = 3;
% 
% coast = load('coastlines');
% %colormap(redblue)
% 
% % min_colors = [-80 -80];
% % max_colors = [80 80];
% 
% min_colors = [0 0];
% max_colors = [200 200];
% 
% for var_i = 1:2
%     
%     if var_i == 1
%         
%         [cmap]=buildcmap(colormap_strings{2});
%     
%     end
%     if var_i == 2
%         
%         [cmap]=buildcmap(colormap_strings{4});
%     
%     end
%     
%     FigHandle = figure('Position', [100, 100, 1300, 1000]); %[left bottom width height]
%     set(gcf,'color',[1 1 1]);
%     set(0, 'DefaultAxesFontSize',6);
%     set(0,'defaultAxesFontName', 'helvetica')
%     hold on
% 
%     for week_i = 1:weeks_to_pull
% 
%         if var_i == 1; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind_insolation(week_i)); end
%         if var_i == 2; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind_insolation(week_i)); end
% 
%         h = subplot(10,4,week_i);
%         hold on
% 
%         axesm('robinson',...
%         'Frame', 'on',...
%         'Grid', 'off',...
%         'maplatlim',[min(WECC_lat) max(WECC_lat)],...
%         'maplonlim',[min(WECC_lon) max(WECC_lon)])
%         tightmap
%         
%         levels = linspace(min_colors(var_i),max_colors(var_i),50);
%                 
%         contourfm(WECC_lat,WECC_lon,100*squeeze(gridded_WECC_week_year_frac_devs_wrt_woy_mean(:,:,week_ind,year_ind,var_i))+100,levels,'LineStyle','none');
%         %pcolorm(WECC_lat,WECC_lon,squeeze(gridded_WECC_week_year_anoms(:,:,week_ind,year_ind,var_i)), 'LineStyle','none');
%         
%         if var_i == 1; contourm(WECC_lat,WECC_lon,gridded_WECC_week_year_ann(:,:,var_i),num_conts,'LineColor','w','LineWidth',1); end
%         if var_i == 2; contourm(WECC_lat,WECC_lon,gridded_WECC_week_year_ann(:,:,var_i),num_conts,'LineColor','k','LineWidth',1); end
% 
%         map_mean_anom =100*threed2oned_lat_lon_2(squeeze(gridded_WECC_week_year_frac_devs_wrt_woy_mean(:,:,week_ind,year_ind,var_i)),WECC_lat,WECC_lon)+100;
%         map_mean_raw = 100*WECC_week_year_1d_norm(week_ind,year_ind,var_i);
%         %map_mean_raw = 100*threed2oned_lat_lon_2(squeeze(gridded_WECC_week_year_frac_devs_wrt_ann_mean(:,:,week_ind,year_ind,var_i)),WECC_lat,WECC_lon);
%         
%              geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
%         
%         %t=colorbar;
%         %set(get(t,'ylabel'),'String', '%');
%         
%         caxis([min_colors(var_i) max_colors(var_i)])
% 
%         title(strcat(num2str(week_i),')',{' '},num2str(years(year_ind)),',',{' '},num2str(week_ind),{' '},'(',num2str(round(map_mean_raw,0)),'%,',num2str(round(map_mean_anom,0)),'%)'))
%         
%        
%         if var_i == 1
%             colormap(h,brewermap([],colormap_brewer_string{2}))
%         end
%         if var_i == 2
%             %colormap(h,brewermap([],colormap_brewer_string{4}))
%             colormap(h,cmap)
%         end       
%        
%        
%     end
% end

%% 500 mb all weeks

% 1  '250 hPa heights',...
% 2  '500 hPa heights',...
% 3  '700 hPa heights',...
% 4  '250 hPa U wind',...
% 5  '500 hpa U wind',...
% 6   '700 hpa U wind',...
% 7   '250 hPa V wind',...
% 8   '500 hpa V wind',...
% 9   '700 hpa V wind',...
% 10  '500 hpa vorticity',...
% 11  '700 hpa temperature',...
% 12  '700 Omega',...
% 13  'sea level pressure',...
% 14  '250 hPa wind speed',...
% 15  '500 hpa wind speed',...
% 16  '700 hpa wind speed',...
% 17  '250 hpa divergence',...
% 18  '500 hpa divergence',...
% 19  '700 hpa divergence',...
% 20  '500 hpa vorticity advection',...
% 21  '700 hpa temperature advection'};

weeks_to_pull = num_weeks_for_drought;
% 
% coast = load('coastlines');
% %colormap(redblue)
% 
% % min_colors = [-80 -80];
% % max_colors = [80 80];
% 
% [NH_lat_era5_grid,NH_lon_era5_grid] = meshgrid(NH_lat_era5,NH_lon_era5);
% 
% min_colors = [0 0];
% max_colors = [200 200];
% 
% for var_i = 1:3
%     
%     FigHandle = figure('Position', [100, 100, 1300, 1000]); %[left bottom width height]
%     set(gcf,'color',[1 1 1]);
%     set(0, 'DefaultAxesFontSize',6);
%     set(0,'defaultAxesFontName', 'helvetica')
%     hold on
% 
%     for week_i = 1:weeks_to_pull
% 
%           if var_i == 1; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind(week_i)); end
%           if var_i == 2; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_insolation(week_i)); end
%           if var_i == 3; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind_insolation(week_i)); end
% 
%         h = subplot(10,4,week_i);
%         hold on
% 
%          %500mb
% 
%              num_conts = 13;
% 
%              contour_var_to_plot = 2;
%              anom_var_to_plot = 2;
%              wind_var_to_plot_U = 5;
%              wind_var_to_plot_V = 8;
% 
%              axesm('robinson',...
%              'Frame', 'on',...
%              'Grid', 'off',...
%              'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
%              'maplonlim',[110 340])
%              tightmap
% 
%              level_min = -1300;
%              level_max = 1300;
%              num_levels = 50;
% 
%              levels = linspace(level_min,level_max,num_levels);
%              contourfm(NH_lat_era5,NH_lon_era5,squeeze(gridded_NH_week_year_anoms(:,:,week_ind,year_ind,anom_var_to_plot)),levels,'Linestyle','none')
% 
%              caxis([min(levels) max(levels)])
% 
%              contourm(NH_lat_era5,NH_lon_era5,squeeze(gridded_NH_week_year(:,:,week_ind,year_ind,contour_var_to_plot)),num_conts,'LineColor','k')
%              quiverm(NH_lat_era5_grid,...
%                      NH_lon_era5_grid,...
%                      squeeze(gridded_NH_week_year(:,:,week_ind,year_ind,wind_var_to_plot_V))',...
%                      squeeze(gridded_NH_week_year(:,:,week_ind,year_ind,wind_var_to_plot_U))',...
%                      'k')
% 
%              geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
%              linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-k','LineWidth',0.5)
% 
% %               t=colorbar;
% %               set(get(t,'ylabel'),'String', 'meters');
% 
% %                      [cmap]=buildcmap(colormap_strings{10});
% %                      colormap(h2,cmap)
%              colormap(h,brewermap([],colormap_brewer_string{10}))
% 
%              %title('500mb heights (contoured), wind (vectors), height anomalies (shaded)')  
% 
%        
%     end
% end

%% SST all weeks

weeks_to_pull = num_weeks_for_drought;
% num_conts = 3;
% 
% coast = load('coastlines');
% %colormap(redblue)
% 
% for var_i = 1:3
%     
%     FigHandle = figure('Position', [100, 100, 1300, 1000]); %[left bottom width height]
%     set(gcf,'color',[1 1 1]);
%     set(0, 'DefaultAxesFontSize',6);
%     set(0,'defaultAxesFontName', 'helvetica')
%     hold on
% 
%     for week_i = 1:weeks_to_pull
% 
%           if var_i == 1; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind(week_i)); end
%           if var_i == 2; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_insolation(week_i)); end
%           if var_i == 3; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind_insolation(week_i)); end
% 
%         h = subplot(10,4,week_i);
%         hold on
% 
%              num_conts = 13;
% 
%              axesm('robinson',...
%              'Frame', 'on',...
%              'Grid', 'off',...
%              'maplatlim',[-10 max(NH_lat_era5)],...
%              'maplonlim',[0 360])
%              tightmap
% 
%              level_min = -3;
%              level_max = 3;
%              num_levels = 50;
% 
%              levels = linspace(level_min,level_max,num_levels);
%              contourfm(Global_lat_era5,Global_lon_era5,squeeze(Global_SSTs_week_year_anom(:,:,week_ind,year_ind)),levels,'Linestyle','none')
% 
%              caxis([min(levels) max(levels)])
% 
%              geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
%              
%              [cmap]=buildcmap(colormap_strings{6});
%              colormap(h,cmap)
% 
%              %title('500mb heights (contoured), wind (vectors), height anomalies (shaded)')  
% 
%        
%     end
% end

%% plot events back on seasonal cycle and scatters

FigHandle = figure('Position', [100, 100, 1100, 800]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',20);
set(0,'defaultAxesFontName', 'helvetica')

var_markers = {'o','s','d'};
var_markers_colors = {'k','m','b'};
%var_markers_colors = {[.75, .0, .75],[0.75, 0.75, 0.75],[0.7, 0.3, 0.4]};
mkrsize = 120;
mkrtrans = 0.4;

h=subplot(2,1,1);
hold on

%set(h,'color',[0.95 0.95 0.95]);
set(h,'color',[1 1 1]);

    boundedline(weeks,100*WECC_week_year_1d_norm_mean(:,1),200*WECC_week_year_1d_norm_std(:,1),'-g','alpha','transparency',0.4)
    boundedline(weeks,100*WECC_week_year_1d_norm_mean(:,2),200*WECC_week_year_1d_norm_std(:,2),'-y','alpha','transparency',0.7)

    plot(weeks,100*ones(1,length(weeks)),'-k','LineWidth',3)
    plot([25 25],[0 100*5.5],'-k','LineWidth',3)
    
    h1 = plot(weeks,100*WECC_week_year_1d_norm_mean(:,1),'Color',[0, 0.5, 0],'linewidth',4);
    h2 = plot(weeks,100*WECC_week_year_1d_norm_mean(:,2),'Color',[0.9290, 0.6940, 0.1250],'linewidth',4);
    
%     plot(weeks,100*WECC_week_year_1d_norm_mean(:,2),'-k','linewidth',2);
%     plot(weeks,100*WECC_week_year_1d_norm_mean(:,2)+200*WECC_week_year_1d_norm_std(:,2),'-k','linewidth',.01);
%     plot(weeks,100*WECC_week_year_1d_norm_mean(:,2)-200*WECC_week_year_1d_norm_std(:,2),'-k','linewidth',.01);
    
    for var_i = 1:3
        for week_i = 1:weeks_to_pull

            if var_i == 1; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind(week_i)); end
            if var_i == 2; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_insolation(week_i)); end
            if var_i == 3; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind_insolation(week_i)); end

            scatter(week_ind,100*WECC_week_year_1d_norm(week_ind,year_ind,1),...
                mkrsize,var_markers_colors{var_i},var_markers{var_i},...
                'filled',...
                'MarkerFaceAlpha',mkrtrans,...
                'MarkerEdgeAlpha',mkrtrans,...
                'MarkerEdgeColor','k',...
                'MarkerEdgeAlpha',.1)
            
            scatter(week_ind,100*WECC_week_year_1d_norm(week_ind,year_ind,2),...
                mkrsize,var_markers_colors{var_i},var_markers{var_i},...
                'filled',...
                'MarkerFaceAlpha',mkrtrans,...
                'MarkerEdgeAlpha',mkrtrans,...
                'MarkerEdgeColor','k',...
                'MarkerEdgeAlpha',.1)
        end
    end
    
%     h3 = scatter(55,2,mkrsize,'k',var_markers{1},'filled','MarkerFaceAlpha',mkrtrans,'MarkerEdgeAlpha',mkrtrans);
%     h4 = scatter(55,2,mkrsize,'k',var_markers{2},'filled','MarkerFaceAlpha',mkrtrans,'MarkerEdgeAlpha',mkrtrans);
%     h5 = scatter(55,2,mkrsize,'k',var_markers{3},'filled','MarkerFaceAlpha',mkrtrans,'MarkerEdgeAlpha',mkrtrans);
%     legend([h1,h2,h3,h4,h5],'wind power','solar power','wind drought','solar drought','wind+solar drought')

     legend([h1,h2],'wind power','solar power')
        
    ylabel('proportion of long-term mean (%)')
    %xlabel('week of the year')
    title('power supply')
    
    ylim([0 200])
    xlim([1 52])
    
h2 = subplot(2,1,2);
hold on

%set(h2,'color',[0.95 0.95 0.95]);
set(h,'color',[1 1 1]);

    boundedline(weeks,100*WECC_week_year_1d_norm_mean(:,3),200*WECC_week_year_1d_norm_std(:,3),'-r','alpha','transparency',0.2)
    boundedline(weeks,100*WECC_week_year_1d_norm_mean(:,4),200*WECC_week_year_1d_norm_std(:,4),'-b','alpha','transparency',0.2)
        
    plot(weeks,100*ones(1,length(weeks)),'-k','LineWidth',3)
    plot([25 25],[0 100*5.5],'-k','LineWidth',3)
    
    h1 = plot(weeks,100*WECC_week_year_1d_norm_mean(:,3),'-r','linewidth',4);
    h2 = plot(weeks,100*WECC_week_year_1d_norm_mean(:,4),'-b','linewidth',4);
    
    for var_i = 1:3
        for week_i = 1:weeks_to_pull

            if var_i == 1; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind(week_i)); end
            if var_i == 2; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_insolation(week_i)); end
            if var_i == 3; [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind_insolation(week_i)); end

            scatter(week_ind,100*WECC_week_year_1d_norm(week_ind,year_ind,3),...
                mkrsize,var_markers_colors{var_i},var_markers{var_i},...
                'filled',...
                'MarkerFaceAlpha',mkrtrans,...
                'MarkerEdgeAlpha',mkrtrans,...
                'MarkerEdgeColor','k',...
                'MarkerEdgeAlpha',.1)
            
            scatter(week_ind,100*WECC_week_year_1d_norm(week_ind,year_ind,4),...
                mkrsize,var_markers_colors{var_i},var_markers{var_i},...
                'filled',...
                'MarkerFaceAlpha',mkrtrans,...
                'MarkerEdgeAlpha',mkrtrans,...
                'MarkerEdgeColor','k',...
                'MarkerEdgeAlpha',.1)

        end
    end
    
    h3 = scatter(55,2,mkrsize,var_markers_colors{1},var_markers{1},'filled','MarkerFaceAlpha',mkrtrans,'MarkerEdgeAlpha',mkrtrans,'MarkerEdgeColor','k','MarkerEdgeAlpha',.1);
    h4 = scatter(55,2,mkrsize,var_markers_colors{2},var_markers{2},'filled','MarkerFaceAlpha',mkrtrans,'MarkerEdgeAlpha',mkrtrans,'MarkerEdgeColor','k','MarkerEdgeAlpha',.1);
    h5 = scatter(55,2,mkrsize,var_markers_colors{3},var_markers{3},'filled','MarkerFaceAlpha',mkrtrans,'MarkerEdgeAlpha',mkrtrans,'MarkerEdgeColor','k','MarkerEdgeAlpha',.1);   
    legend([h1,h2,h3,h4,h5],'cooling degree-days','heating degree-days','wind drought','solar drought','wind+solar drought','location','NorthWest')
        
    ylabel('proportion of long-term mean (%)')
    xlabel('week of the year')
    title('power demand')
    
    ylim([0 100*5.5])
    xlim([1 52])

%% composite preproc (now add the wind+solar)

% wind comps

    extreme_array_gridded_WECC_week_year_anom_wrt_woy_wind = [];
    extreme_array_gridded_WECC_week_year_anom_wrt_ann_wind = [];
    
    extreme_array_gridded_WECC_week_year_wind = [];
    extreme_array_gridded_WECC_week_year_anoms_wind = [];
    
    extreme_array_mean_WECC_week_year_anom_wrt_ann_wind = [];

    extreme_array_gridded_NH_anoms_wind = [];
    extreme_array_gridded_Global_SSTs_anoms_wind = [];
    
    extreme_array_gridded_NH_wind = [];
    extreme_array_gridded_Global_SSTs_wind = [];

    for week_i = 1:weeks_to_pull

        [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind(week_i));

        extreme_array_gridded_WECC_week_year_anom_wrt_woy_wind = cat(4,extreme_array_gridded_WECC_week_year_anom_wrt_woy_wind,squeeze(gridded_WECC_week_year_frac_devs_wrt_woy_mean(:,:,week_ind,year_ind,:)));
        extreme_array_gridded_WECC_week_year_anom_wrt_ann_wind = cat(4,extreme_array_gridded_WECC_week_year_anom_wrt_ann_wind,squeeze(gridded_WECC_week_year_frac_devs_wrt_ann_mean(:,:,week_ind,year_ind,:)));
        
        extreme_array_gridded_WECC_week_year_wind = cat(4,extreme_array_gridded_WECC_week_year_wind,squeeze(gridded_WECC_week_year(:,:,week_ind,year_ind,:)));
        extreme_array_gridded_WECC_week_year_anoms_wind = cat(4,extreme_array_gridded_WECC_week_year_anoms_wind,squeeze(gridded_WECC_week_year_anoms(:,:,week_ind,year_ind,:)));
        
        extreme_array_mean_WECC_week_year_anom_wrt_ann_wind = cat(2,extreme_array_mean_WECC_week_year_anom_wrt_ann_wind,squeeze(WECC_week_year_1d_norm(week_ind,year_ind,:)));

        extreme_array_gridded_NH_anoms_wind = cat(4,extreme_array_gridded_NH_anoms_wind,squeeze(gridded_NH_week_year_anoms(:,:,week_ind,year_ind,:)));
        extreme_array_gridded_NH_wind = cat(4,extreme_array_gridded_NH_wind,squeeze(gridded_NH_week_year(:,:,week_ind,year_ind,:)));

        extreme_array_gridded_Global_SSTs_anoms_wind = cat(3,extreme_array_gridded_Global_SSTs_anoms_wind,squeeze(Global_SSTs_week_year_anom(:,:,week_ind,year_ind)));
        extreme_array_gridded_Global_SSTs_wind = cat(3,extreme_array_gridded_Global_SSTs_wind,squeeze(Global_SSTs_week_year(:,:,week_ind,year_ind)));
    end

% sun comps

    extreme_array_gridded_WECC_week_year_anom_wrt_woy_sun = [];
    extreme_array_gridded_WECC_week_year_anom_wrt_ann_sun = [];
    
    extreme_array_gridded_WECC_week_year_sun = [];
    extreme_array_gridded_WECC_week_year_anoms_sun = [];    
    
    extreme_array_mean_WECC_week_year_anom_wrt_ann_sun = [];

    extreme_array_gridded_NH_anoms_sun = [];
    extreme_array_gridded_Global_SSTs_anoms_sun = [];
    
    extreme_array_gridded_NH_sun = [];
    extreme_array_gridded_Global_SSTs_sun = [];


    for week_i = 1:weeks_to_pull

        [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_insolation(week_i));

        extreme_array_gridded_WECC_week_year_anom_wrt_woy_sun = cat(4,extreme_array_gridded_WECC_week_year_anom_wrt_woy_sun,squeeze(gridded_WECC_week_year_frac_devs_wrt_woy_mean(:,:,week_ind,year_ind,:)));
        extreme_array_gridded_WECC_week_year_anom_wrt_ann_sun = cat(4,extreme_array_gridded_WECC_week_year_anom_wrt_ann_sun,squeeze(gridded_WECC_week_year_frac_devs_wrt_ann_mean(:,:,week_ind,year_ind,:)));
        
        extreme_array_gridded_WECC_week_year_sun = cat(4,extreme_array_gridded_WECC_week_year_sun,squeeze(gridded_WECC_week_year(:,:,week_ind,year_ind,:)));
        extreme_array_gridded_WECC_week_year_anoms_sun = cat(4,extreme_array_gridded_WECC_week_year_anoms_sun,squeeze(gridded_WECC_week_year_anoms(:,:,week_ind,year_ind,:)));        
        
        extreme_array_mean_WECC_week_year_anom_wrt_ann_sun = cat(2,extreme_array_mean_WECC_week_year_anom_wrt_ann_sun,squeeze(WECC_week_year_1d_norm(week_ind,year_ind,:)));

        extreme_array_gridded_NH_anoms_sun = cat(4,extreme_array_gridded_NH_anoms_sun,squeeze(gridded_NH_week_year_anoms(:,:,week_ind,year_ind,:)));
        extreme_array_gridded_NH_sun = cat(4,extreme_array_gridded_NH_sun,squeeze(gridded_NH_week_year(:,:,week_ind,year_ind,:)));

        extreme_array_gridded_Global_SSTs_anoms_sun = cat(3,extreme_array_gridded_Global_SSTs_anoms_sun,squeeze(Global_SSTs_week_year_anom(:,:,week_ind,year_ind)));
        extreme_array_gridded_Global_SSTs_sun = cat(3,extreme_array_gridded_Global_SSTs_sun,squeeze(Global_SSTs_week_year(:,:,week_ind,year_ind)));
    end
    
% sun+wind comps

    extreme_array_gridded_WECC_week_year_anom_wrt_woy_wind_sun = [];
    extreme_array_gridded_WECC_week_year_anom_wrt_ann_wind_sun = [];
    
    extreme_array_gridded_WECC_week_year_wind_sun = [];
    extreme_array_gridded_WECC_week_year_anoms_wind_sun = [];
    
    extreme_array_mean_WECC_week_year_anom_wrt_ann_wind_sun = [];

    extreme_array_gridded_NH_anoms_wind_sun = [];
    extreme_array_gridded_Global_SSTs_anoms_wind_sun = [];
    
    extreme_array_gridded_NH_wind_sun = [];
    extreme_array_gridded_Global_SSTs_wind_sun = [];

    for week_i = 1:weeks_to_pull

        [week_ind,year_ind] = find(week_year_2d_week_number_key == ranked_index_wind_insolation(week_i));

        extreme_array_gridded_WECC_week_year_anom_wrt_woy_wind_sun = cat(4,extreme_array_gridded_WECC_week_year_anom_wrt_woy_wind_sun,squeeze(gridded_WECC_week_year_frac_devs_wrt_woy_mean(:,:,week_ind,year_ind,:)));
        extreme_array_gridded_WECC_week_year_anom_wrt_ann_wind_sun = cat(4,extreme_array_gridded_WECC_week_year_anom_wrt_ann_wind_sun,squeeze(gridded_WECC_week_year_frac_devs_wrt_ann_mean(:,:,week_ind,year_ind,:)));
        
        extreme_array_gridded_WECC_week_year_wind_sun = cat(4,extreme_array_gridded_WECC_week_year_wind_sun,squeeze(gridded_WECC_week_year(:,:,week_ind,year_ind,:)));
        extreme_array_gridded_WECC_week_year_anoms_wind_sun = cat(4,extreme_array_gridded_WECC_week_year_anoms_wind_sun,squeeze(gridded_WECC_week_year_anoms(:,:,week_ind,year_ind,:)));         
        
        extreme_array_mean_WECC_week_year_anom_wrt_ann_wind_sun = cat(2,extreme_array_mean_WECC_week_year_anom_wrt_ann_wind_sun,squeeze(WECC_week_year_1d_norm(week_ind,year_ind,:)));

        extreme_array_gridded_NH_anoms_wind_sun = cat(4,extreme_array_gridded_NH_anoms_wind_sun,squeeze(gridded_NH_week_year_anoms(:,:,week_ind,year_ind,:)));
        extreme_array_gridded_NH_wind_sun = cat(4,extreme_array_gridded_NH_wind_sun,squeeze(gridded_NH_week_year(:,:,week_ind,year_ind,:)));

        extreme_array_gridded_Global_SSTs_anoms_wind_sun = cat(3,extreme_array_gridded_Global_SSTs_anoms_wind_sun,squeeze(Global_SSTs_week_year_anom(:,:,week_ind,year_ind)));
        extreme_array_gridded_Global_SSTs_wind_sun = cat(3,extreme_array_gridded_Global_SSTs_wind_sun,squeeze(Global_SSTs_week_year(:,:,week_ind,year_ind)));
    end   

%% main figure

num_conts_2 = 3;

for var2condon_i = 1:3
    
    if var2condon_i == 1
        %anom
            extreme_array_gridded_WECC_week_year_raw_anoms_var2co = extreme_array_gridded_WECC_week_year_anoms_wind;
            extreme_array_gridded_WECC_week_year_anoms_var2co = extreme_array_gridded_WECC_week_year_anom_wrt_woy_wind;
            extreme_array_gridded_NH_anoms_var2co = extreme_array_gridded_NH_anoms_wind;
            extreme_array_gridded_Global_SSTs_anoms_var2co = extreme_array_gridded_Global_SSTs_anoms_wind;
            
        %raw
            extreme_array_gridded_WECC_week_year_raw_var2co = extreme_array_gridded_WECC_week_year_wind;
            extreme_array_gridded_WECC_week_year_var2co = extreme_array_gridded_WECC_week_year_anom_wrt_ann_wind;
            extreme_array_mean_WECC_week_year_var2co = extreme_array_mean_WECC_week_year_anom_wrt_ann_wind;
            extreme_array_gridded_NH_var2co = extreme_array_gridded_NH_wind;
            extreme_array_gridded_Global_SSTs_var2co = extreme_array_gridded_Global_SSTs_wind;
    end
    if var2condon_i == 2
        %anom
            extreme_array_gridded_WECC_week_year_raw_anoms_var2co = extreme_array_gridded_WECC_week_year_anoms_sun;
            extreme_array_gridded_WECC_week_year_anoms_var2co = extreme_array_gridded_WECC_week_year_anom_wrt_woy_sun;
            extreme_array_gridded_NH_anoms_var2co = extreme_array_gridded_NH_anoms_sun;
            extreme_array_gridded_Global_SSTs_anoms_var2co = extreme_array_gridded_Global_SSTs_anoms_sun;
        %raw
            extreme_array_gridded_WECC_week_year_raw_var2co = extreme_array_gridded_WECC_week_year_sun;
            extreme_array_gridded_WECC_week_year_var2co = extreme_array_gridded_WECC_week_year_anom_wrt_ann_sun;
            extreme_array_mean_WECC_week_year_var2co = extreme_array_mean_WECC_week_year_anom_wrt_ann_sun;
            extreme_array_gridded_NH_var2co = extreme_array_gridded_NH_sun;
            extreme_array_gridded_Global_SSTs_var2co = extreme_array_gridded_Global_SSTs_sun;
    end
    if var2condon_i == 3
        %anom
            extreme_array_gridded_WECC_week_year_raw_anoms_var2co = extreme_array_gridded_WECC_week_year_anoms_wind_sun;
            extreme_array_gridded_WECC_week_year_anoms_var2co = extreme_array_gridded_WECC_week_year_anom_wrt_woy_wind_sun;
            extreme_array_gridded_NH_anoms_var2co = extreme_array_gridded_NH_anoms_wind_sun;
            extreme_array_gridded_Global_SSTs_anoms_var2co = extreme_array_gridded_Global_SSTs_anoms_wind_sun;
        %raw
            extreme_array_gridded_WECC_week_year_raw_var2co = extreme_array_gridded_WECC_week_year_wind_sun;
            extreme_array_gridded_WECC_week_year_var2co = extreme_array_gridded_WECC_week_year_anom_wrt_ann_wind_sun;
            extreme_array_mean_WECC_week_year_var2co = extreme_array_mean_WECC_week_year_anom_wrt_ann_wind_sun;
            extreme_array_gridded_NH_var2co = extreme_array_gridded_NH_wind_sun;
            extreme_array_gridded_Global_SSTs_var2co = extreme_array_gridded_Global_SSTs_wind_sun;
    end
    
        frac_for_thresh = 0.7;

        %WECC

            stipple_lat_lons_WECC_var2co = NaN(length(WECC_lat)*length(WECC_lon),2,length(WECC_daily_var_names));

            for var_i = 1:length(WECC_daily_var_names)

                lin_loc = 1;

                for lat_i = 1:length(WECC_lat)
                    for lon_i = 1:length(WECC_lon)

                        if var_i <= 2
                            anoms_across_extreme_weeks = squeeze(extreme_array_gridded_WECC_week_year_anoms_var2co(lat_i,lon_i,var_i,:));
                        end
                        if var_i >2
                            anoms_across_extreme_weeks = squeeze(extreme_array_gridded_WECC_week_year_raw_anoms_var2co(lat_i,lon_i,var_i,:));
                        end

                        num_pos = length(find(anoms_across_extreme_weeks > 0.1));

                        if sum(isnan(anoms_across_extreme_weeks)) ~= length(anoms_across_extreme_weeks)
                            if isempty(num_pos) == 0
                                if num_pos/length(anoms_across_extreme_weeks) >= frac_for_thresh || num_pos/length(anoms_across_extreme_weeks) <= (1-frac_for_thresh)

                                    stipple_lat_lons_WECC_var2co(lin_loc,1,var_i) = WECC_lat(lat_i);
                                    stipple_lat_lons_WECC_var2co(lin_loc,2,var_i) = WECC_lon(lon_i);

                                end
                            end
                        end

                        lin_loc = lin_loc + 1;

                    end
                end
            end

        %NH

            stipple_lat_lons_NH_var2co = NaN(length(NH_lat_era5)*length(NH_lon_era5),2,length(all_gridded_NH_data_varnames_desc));

            for var_i = 1:length(all_gridded_NH_data_varnames_desc)

                lin_loc = 1;

                for lat_i = 1:length(NH_lat_era5)
                    for lon_i = 1:length(NH_lon_era5)
                        
                        if var_i == 11 || var_i == 15 || var_i == 16 || var_i == 17 || var_i == 18 || var_i == 19   %certain variables you want to stipple on abs pos/neg not anoms
                            anoms_across_extreme_weeks = squeeze(extreme_array_gridded_NH_var2co(lat_i,lon_i,var_i,:));
                        else
                            anoms_across_extreme_weeks = squeeze(extreme_array_gridded_NH_anoms_var2co(lat_i,lon_i,var_i,:));
                        end

                        num_pos = length(find(anoms_across_extreme_weeks > 0));

                        if sum(isnan(anoms_across_extreme_weeks)) ~= length(anoms_across_extreme_weeks)
                            if isempty(num_pos) == 0
                                if num_pos/length(anoms_across_extreme_weeks) >= frac_for_thresh || num_pos/length(anoms_across_extreme_weeks) <= (1-frac_for_thresh)

                                    stipple_lat_lons_NH_var2co(lin_loc,1,var_i) = NH_lat_era5(lat_i);
                                    stipple_lat_lons_NH_var2co(lin_loc,2,var_i) = NH_lon_era5(lon_i);

                                end
                            end
                        end

                        lin_loc = lin_loc + 1;

                    end
                end
            end

        %Global SSTs

            stipple_lat_lons_Global_SSTs_var2co = NaN(length(Global_lat_era5)*length(Global_lon_era5),2);

                lin_loc = 1;

                for lat_i = 1:length(Global_lat_era5)
                    for lon_i = 1:length(Global_lon_era5)

                        anoms_across_extreme_weeks = squeeze(extreme_array_gridded_Global_SSTs_anoms_var2co(lat_i,lon_i,:));

                        num_pos = length(find(anoms_across_extreme_weeks > 0));

                        if sum(isnan(anoms_across_extreme_weeks)) ~= length(anoms_across_extreme_weeks)
                            if isempty(num_pos) == 0
                                if num_pos/length(anoms_across_extreme_weeks) >= frac_for_thresh || num_pos/length(anoms_across_extreme_weeks) <= (1-frac_for_thresh)

                                    stipple_lat_lons_Global_SSTs_var2co(lin_loc,1) = Global_lat_era5(lat_i);
                                    stipple_lat_lons_Global_SSTs_var2co(lin_loc,2) = Global_lon_era5(lon_i);

                                end
                            end
                        end

                        lin_loc = lin_loc + 1;

                    end
                end

        extreme_array_gridded_WECC_week_year_anoms_comp_var2co = mean(extreme_array_gridded_WECC_week_year_anoms_var2co,4);
        extreme_array_gridded_WECC_week_year_raw_anoms_comp_var2co = mean(extreme_array_gridded_WECC_week_year_raw_anoms_var2co,4);
        extreme_array_gridded_WECC_week_year_raw_comp_var2co = mean(extreme_array_gridded_WECC_week_year_raw_var2co,4);
        extreme_array_gridded_WECC_week_year_comp_var2co = mean(extreme_array_gridded_WECC_week_year_var2co,4);
        extreme_array_mean_WECC_week_year_comp_var2co = mean(extreme_array_mean_WECC_week_year_var2co,2);

        extreme_array_gridded_NH_anoms_comp_var2co = mean(extreme_array_gridded_NH_anoms_var2co,4);
        extreme_array_gridded_NH_comp_var2co = mean(extreme_array_gridded_NH_var2co,4);

        extreme_array_gridded_Global_SSTs_anoms_comp_var2co = mean(extreme_array_gridded_Global_SSTs_anoms_var2co,3);
        extreme_array_gridded_Global_SSTs_comp_var2co = mean(extreme_array_gridded_Global_SSTs_var2co,3);

        % plot composites - wind

        % plot WECC composites

%                 min_colors = [-170 -50 -5 -5];
%                 max_colors = [170 50 5 5];

%                  min_colors = [-80 -80 -80 -80];
%                  max_colors = [80 80 80 80];

%                   min_colors = [-100 -100 -100 -100];
%                   max_colors = [100 100 100 100];

                  min_colors = [0 0 0 0 -7];
                  max_colors = [200 200 200 200 7];

                FigHandle = figure('Position', [100, 100, 1000, 1000]); %[left bottom width height]
                set(gcf,'color',[1 1 1]);
                set(0, 'DefaultAxesFontSize',11);
                set(0,'defaultAxesFontName', 'helvetica')
                hold on

                for var_i = 1:2
                    
                        if var_i == 1
                            [cmap]=buildcmap(colormap_strings{2});
                        end
                        if var_i == 2
                            [cmap]=buildcmap(colormap_strings{4});
                        end

                            h = subplot(3,1,var_i);
                            hold on

                            axesm('robinson',...
                            'Frame', 'on',...
                            'Grid', 'off',...
                            'maplatlim',[min(WECC_lat) max(WECC_lat)],...
                            'maplonlim',[min(WECC_lon) max(WECC_lon)])
                            tightmap

                            levels = linspace(min_colors(var_i),max_colors(var_i),30);

                            %pcolorm(WECC_lat,WECC_lon,squeeze(extreme_array_gridded_WECC_week_year_anoms_comp_var2co(:,:,var_i)), 'LineStyle','none');
                            contourfm(WECC_lat,WECC_lon,100*squeeze(extreme_array_gridded_WECC_week_year_anoms_comp_var2co(:,:,var_i))+100,levels,'LineStyle','none');
                            
                            contourm(WECC_lat,WECC_lon,gridded_WECC_week_year_ann(:,:,var_i),num_conts_2,'LineColor','k','LineWidth',1)

                            scatterm(stipple_lat_lons_WECC_var2co(:,1,var_i),stipple_lat_lons_WECC_var2co(:,2,var_i),6,'k','filled')

                            map_mean_anom = 100*threed2oned_lat_lon_2(squeeze(extreme_array_gridded_WECC_week_year_anoms_comp_var2co(:,:,var_i)),WECC_lat,WECC_lon)+100;
                            map_mean = 100*extreme_array_mean_WECC_week_year_comp_var2co(var_i);

                            caxis([min_colors(var_i) max_colors(var_i)])

                            geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')

                            t=colorbar;
                            set(get(t,'ylabel'),'String', '%');
                            
                            if var_i == 1
                                colormap(h,brewermap([],colormap_brewer_string{2}))
                            end                            
                            if var_i == 2
                                colormap(h,cmap);
                            end


                            title(strcat(var_names_2{var_i},{' '},'(',num2str(round(map_mean,0)),'%,',num2str(round(map_mean_anom,0)),'%)'))

                end
                
                %plot deviations
                
                var_i = 5;
                
                    [cmap]=buildcmap(colormap_strings{6});
                
                            h = subplot(3,1,3);
                            hold on

                            axesm('robinson',...
                            'Frame', 'on',...
                            'Grid', 'off',...
                            'maplatlim',[min(WECC_lat) max(WECC_lat)],...
                            'maplonlim',[min(WECC_lon) max(WECC_lon)])
                            tightmap
                            
                            contour_map = extreme_array_gridded_WECC_week_year_raw_comp_var2co(:,:,var_i);
                            
                                contour_map_pos = zeros(size(contour_map));
                                pos_inds = find(contour_map >= 0);
                                contour_map_pos(pos_inds) = contour_map(pos_inds);

                                contour_map_neg = zeros(size(contour_map));
                                neg_inds = find(contour_map < 0);
                                contour_map_neg(neg_inds) = contour_map(neg_inds);
                            
                            levels = linspace(min_colors(var_i),max_colors(var_i),30);

                            contourfm(WECC_lat,WECC_lon,extreme_array_gridded_WECC_week_year_raw_anoms_comp_var2co(:,:,var_i),levels,'LineStyle','none');
                            
                            contourm(WECC_lat,WECC_lon,contour_map_pos,3,'LineColor','r','LineWidth',1)
                            contourm(WECC_lat,WECC_lon,contour_map_neg,3,'LineColor','b','LineWidth',1)

                            scatterm(stipple_lat_lons_WECC_var2co(:,1,var_i),stipple_lat_lons_WECC_var2co(:,2,var_i),6,'k','filled')

                            map_mean_anom = threed2oned_lat_lon_2(squeeze(extreme_array_gridded_WECC_week_year_raw_anoms_comp_var2co(:,:,var_i)),WECC_lat,WECC_lon);
                            map_mean = threed2oned_lat_lon_2(squeeze(extreme_array_gridded_WECC_week_year_raw_comp_var2co(:,:,var_i)),WECC_lat,WECC_lon);

                            caxis([min_colors(var_i) max_colors(var_i)])

                            geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')

                            t=colorbar;
                            set(get(t,'ylabel'),'String', 'temperature anomaly');
                            
                            colormap(h,cmap);

                            title(strcat('temp departure from 18C & anomaly (',num2str(round(map_mean,0)),',',num2str(round(map_mean_anom,0)),')'))
                
                
        % plot NH composites

        [NH_lat_era5_grid,NH_lon_era5_grid] = meshgrid(NH_lat_era5,NH_lon_era5);

            % 1  '250 hPa heights',...
            % 2  '500 hPa heights',...
            % 3  '700 hPa heights',...
            % 4  '250 hPa U wind',...
            % 5  '500 hpa U wind',...
            % 6   '700 hpa U wind',...
            % 7   '250 hPa V wind',...
            % 8   '500 hpa V wind',...
            % 9   '700 hpa V wind',...
            % 10  '500 hpa vorticity',...
            % 11  '700 hpa temperature',...
            % 12  '700 Omega',...
            % 13  'sea level pressure',...
            % 14  '250 hPa wind speed',...
            % 15  '500 hpa wind speed',...
            % 16  '700 hpa wind speed',...
            % 17  '250 hpa divergence',...
            % 18  '500 hpa divergence',...
            % 19  '700 hpa divergence',...
            % 20  '500 hpa vorticity advection',...
            % 21  '700 hpa temperature advection'};

                 FigHandle = figure('Position', [100, 100, 1000, 1000]); %[left bottom width height]
                 set(gcf,'color',[1 1 1]);
                 set(0, 'DefaultAxesFontSize',14);
                 set(0,'defaultAxesFontName', 'helvetica')
                 hold on

                 %250mb

                     h1 = subplot(4,1,1);
                     hold on

                     num_conts = 13;

                     contour_var_to_plot = 1;
                     anom_var_to_plot = 14;
                     wind_var_to_plot_U = 4;
                     wind_var_to_plot_V = 7;

                     axesm('robinson',...
                     'Frame', 'on',...
                     'Grid', 'off',...
                     'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
                     'maplonlim',[110 340])
                     tightmap

                     levels = linspace(0,70,15);
                     contourfm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot)),levels,'Linestyle','none')
                     %contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot)),5,'LineWidth',2,'LineColor','w')

                        %stipple
                        %scatterm(stipple_lat_lons_NH_var2co(:,1,anom_var_to_plot),stipple_lat_lons_NH_var2co(:,2,anom_var_to_plot),1,'k','filled')

                     caxis([min(levels) max(levels)])

                     contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,contour_var_to_plot)),num_conts,'LineColor','k')
                     quiverm(NH_lat_era5_grid,...
                             NH_lon_era5_grid,...
                             squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_V))',...
                             squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_U))',...
                             'k')

                     geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
                     linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-k','LineWidth',3)

                      t=colorbar;
                      set(get(t,'ylabel'),'String', 'meters per second');
                      
                     colormap(h1,brewermap([],colormap_brewer_string{1}))
                      
%                      [cmap]=buildcmap(colormap_strings{1});
%                      colormap(h1,cmap)

                     title('250mb heights (contoured), wind (vectors), wind speed (shaded)')

                 %500mb

                     h2 = subplot(4,1,2);
                     hold on

                     num_conts = 13;

                     contour_var_to_plot = 2;
                     anom_var_to_plot = 2;
                     wind_var_to_plot_U = 5;
                     wind_var_to_plot_V = 8;

                     axesm('robinson',...
                     'Frame', 'on',...
                     'Grid', 'off',...
                     'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
                     'maplonlim',[110 340])
                     tightmap

                     level_min = -1300;
                     level_max = 1300;
                     num_levels = 50;

                     levels = linspace(level_min,level_max,num_levels);
                     contourfm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_anoms_comp_var2co(:,:,anom_var_to_plot)),levels,'Linestyle','none')

                     %stipple
                     scatterm(stipple_lat_lons_NH_var2co(:,1,anom_var_to_plot),stipple_lat_lons_NH_var2co(:,2,anom_var_to_plot),10,'k','filled')

                     caxis([min(levels) max(levels)])

                     contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,contour_var_to_plot)),num_conts,'LineColor','k')
                     quiverm(NH_lat_era5_grid,...
                             NH_lon_era5_grid,...
                             squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_V))',...
                             squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_U))',...
                             'k')

                     geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
                     linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-k','LineWidth',3)

                      t=colorbar;
                      set(get(t,'ylabel'),'String', 'meters');
                      
%                      [cmap]=buildcmap(colormap_strings{10});
%                      colormap(h2,cmap)
                     colormap(h2,brewermap([],colormap_brewer_string{10}))

                     title('500mb heights (contoured), wind (vectors), height anomalies (shaded)')

                 %850mb - primary

%                      h3 = subplot(4,1,3);
%                      hold on
% 
%                      num_conts = 13;
% 
%                      contour_var_to_plot = 4;
%                      anom_var_to_plot = 10; %4 is 850mb
%                      wind_var_to_plot_U = 7;
%                      wind_var_to_plot_V = 10;
% 
%                      axesm('robinson',...
%                      'Frame', 'on',...
%                      'Grid', 'off',...
%                      'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
%                      'maplonlim',[110 340])
%                      tightmap
% 
%                      %levels (temp)
%                          level_min = -5;
%                          level_max = 5;
% 
%                      num_levels = 50;
% 
%                      levels = linspace(level_min,level_max,num_levels);
%                      contourfm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_anoms_comp_var2co(:,:,anom_var_to_plot)),levels,'Linestyle','none')
% 
%                      %stipple
%                      scatterm(stipple_lat_lons_NH_var2co(:,1,anom_var_to_plot),stipple_lat_lons_NH_var2co(:,2,anom_var_to_plot),10,'k','filled')
% 
%                      caxis([min(levels) max(levels)])
% 
%                      colormap(h,redblue)
% 
%                      contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,contour_var_to_plot)),num_conts,'LineColor','k')
%                      quiverm(NH_lat_era5_grid,...
%                              NH_lon_era5_grid,...
%                              squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_V))',...
%                              squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_U))',...
%                              'k')
% 
%                      geoshow(coast.lat, coast.long, 'Color', 'black')
%                      linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-k','LineWidth',3)
% 
%                       t=colorbar;
%                       set(get(t,'ylabel'),'String', 'degrees C');
% 
%                      [cmap]=buildcmap(colormap_strings{6});
%                      colormap(h3,cmap)
%                       
%                      %title('850mb winds and temperature anomalies')
%                      title('sea level pressure (contoured), 850mb wind (vectors), 850mb temperature anomalies (shaded)')
                     
                     %SLP
                     
                     h3 = subplot(4,1,3);
                     hold on

                     num_conts = 10;

                     contour_var_to_plot = 13;
                     anom_var_to_plot = 13; 
                     wind_var_to_plot_U = 6;
                     wind_var_to_plot_V = 9;

                     axesm('robinson',...
                     'Frame', 'on',...
                     'Grid', 'off',...
                     'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
                     'maplonlim',[110 340])
                     tightmap

                     %levels (temp)
%                          level_min = -1000;
%                          level_max = 1000;
                         level_min = -1000;
                         level_max = 1000;

                     num_levels = 50;

                     levels = linspace(level_min,level_max,num_levels);
                     %contourfm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_anoms_comp_var2co(:,:,anom_var_to_plot)),levels,'Linestyle','none')
                     contourfm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_anoms_comp_var2co(:,:,contour_var_to_plot)),levels,'Linestyle','none')

                     %stipple
                     scatterm(stipple_lat_lons_NH_var2co(:,1,anom_var_to_plot),stipple_lat_lons_NH_var2co(:,2,anom_var_to_plot),10,'k','filled')

                     caxis([min(levels) max(levels)])

                     colormap(h,redblue)

                     contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_anoms_comp_var2co(:,:,contour_var_to_plot)),num_conts,'LineColor','k')
%                      quiverm(NH_lat_era5_grid,...
%                              NH_lon_era5_grid,...
%                              squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_V))',...
%                              squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_U))',...
%                              'k')

                     geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
                     linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-k','LineWidth',3)

                      t=colorbar;
                      set(get(t,'ylabel'),'String', 'mb');

%                      [cmap]=buildcmap(colormap_strings{6});
%                      colormap(h3,cmap)
                     colormap(h3,brewermap([],colormap_brewer_string{10}))
                      
                     %title('850mb winds and temperature anomalies')
                     title('sea level pressure anomalies')

                  %SSTs

                     h4 = subplot(4,1,4);
                     hold on

                     axesm('robinson',...
                     'Frame', 'on',...
                     'Grid', 'off',...
                     'maplatlim',[-10 max(NH_lat_era5)],...
                     'maplonlim',[0 360])
                     tightmap

                     level_min = -2;
                     level_max = 2;
                     num_levels = 50;

                        levels = linspace(level_min,level_max,num_levels);

                        %pcolorm(Global_lat_era5,Global_lon_era5,extreme_array_gridded_Global_SSTs_anoms_comp_var2co, 'LineStyle','none');        
                        contourfm(Global_lat_era5,Global_lon_era5,extreme_array_gridded_Global_SSTs_anoms_comp_var2co,levels,'Linestyle','none')

                        %stipple
                        scatterm(stipple_lat_lons_Global_SSTs_var2co(:,1),stipple_lat_lons_Global_SSTs_var2co(:,2),10,'k','filled')

                        map_mean_anom = threed2oned_lat_lon_2(extreme_array_gridded_Global_SSTs_anoms_comp_var2co,Global_lat_era5,Global_lon_era5);
                        map_mean_raw = threed2oned_lat_lon_2(extreme_array_gridded_Global_SSTs_comp_var2co,Global_lat_era5,Global_lon_era5);

                       caxis([min(levels) max(levels)])

                        geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')

                        t=colorbar;
                        set(get(t,'ylabel'),'String', 'degree C');

                     %linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-g','LineWidth',4)
                     %linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-k','LineWidth',3)

                     [cmap]=buildcmap(colormap_strings{6});
                     colormap(h4,cmap)
                     
                        title('sea surface temperature anomalies (shaded)')
                        
        % plot NH composites - Synoptic Dynamics
        
        second_per_day = 24*60*60;

                % 1  '250 hPa heights',...
                % 2  '500 hPa heights',...
                % 3  '700 hPa heights',...
                % 4  '250 hPa U wind',...
                % 5  '500 hpa U wind',...
                % 6   '700 hpa U wind',...
                % 7   '250 hPa V wind',...
                % 8   '500 hpa V wind',...
                % 9   '700 hpa V wind',...
                % 10  '500 hpa vorticity',...
                % 11  '700 hpa temperature',...
                % 12  '700 Omega',...
                % 13  'sea level pressure',...
                % 14  '250 hPa wind speed',...
                % 15  '500 hpa wind speed',...
                % 16  '700 hpa wind speed',...
                % 17  '250 hpa divergence',...
                % 18  '500 hpa divergence',...
                % 19  '700 hpa divergence',...
                % 20  '500 hpa vorticity advection',...
                % 21  '700 hpa temperature advection'};

                 FigHandle = figure('Position', [100, 100, 1000, 1000]); %[left bottom width height]
                 set(gcf,'color',[1 1 1]);
                 set(0, 'DefaultAxesFontSize',14);
                 set(0,'defaultAxesFontName', 'helvetica')
                 hold on
                     
                 % 500 divergence

                     h1 = subplot(4,1,1);
                     hold on

                     num_conts = 13;

                     contour_var_to_plot = 2;
                     anom_var_to_plot = 18;
                     wind_var_to_plot_U = 5;
                     wind_var_to_plot_V = 8;

                     axesm('robinson',...
                     'Frame', 'on',...
                     'Grid', 'off',...
                     'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
                     'maplonlim',[110 340])
                     tightmap
                     
                     %max_val = 0.9*max(max(abs(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot))));

                     %levels = linspace(-max_val,max_val,30);
                     levels = linspace(-1E-5,1E-5,30);
                     contourfm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot)),levels,'Linestyle','none')
                     %contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot)),5,'LineWidth',2,'LineColor','w')

                     %stipple
                     scatterm(stipple_lat_lons_NH_var2co(:,1,anom_var_to_plot),stipple_lat_lons_NH_var2co(:,2,anom_var_to_plot),10,'k','filled')

                     caxis([min(levels) max(levels)])

                     contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,contour_var_to_plot)),num_conts,'LineColor','k')
                     quiverm(NH_lat_era5_grid,...
                             NH_lon_era5_grid,...
                             squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_V))',...
                             squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_U))',...
                             'k')

                     geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
                     linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-k','LineWidth',3)

                      t=colorbar;
                      set(get(t,'ylabel'),'String', 'per second');
                      
%                      [cmap]=buildcmap(colormap_strings{1});
%                      colormap(h1,cmap)
%                      [cmap]=buildcmap(colormap_strings{10});
%                      colormap(h3,cmap)
                     
                     colormap(h1,brewermap([],colormap_brewer_string{12}))

                     title('500mb divergence (shaded)')
                     
                 % 500 vorticity

                     h2 = subplot(4,1,2);
                     hold on

                     num_conts = 13;

                     contour_var_to_plot = 2;
                     anom_var_to_plot = 10;
                     wind_var_to_plot_U = 5;
                     wind_var_to_plot_V = 8;

                     axesm('robinson',...
                     'Frame', 'on',...
                     'Grid', 'off',...
                     'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
                     'maplonlim',[110 340])
                     tightmap
                     
                     %max_val = 0.9*max(max(abs(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot))));

                     %levels = linspace(-max_val,max_val,30);
                     levels = linspace(-3E-5,3E-5,30);
                     contourfm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot)),levels,'Linestyle','none')
                     %contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot)),5,'LineWidth',2,'LineColor','w')

                     %stipple
                     scatterm(stipple_lat_lons_NH_var2co(:,1,anom_var_to_plot),stipple_lat_lons_NH_var2co(:,2,anom_var_to_plot),10,'k','filled')

                     caxis([min(levels) max(levels)])

                     contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,contour_var_to_plot)),num_conts,'LineColor','k')
                     quiverm(NH_lat_era5_grid,...
                             NH_lon_era5_grid,...
                             squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_V))',...
                             squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_U))',...
                             'k')

                     geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
                     linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-k','LineWidth',3)

                      t=colorbar;
                      set(get(t,'ylabel'),'String', 'per second');
                      
%                      [cmap]=buildcmap(colormap_strings{1});
%                      colormap(h1,cmap)
%                      colormap(h1,parula)

%                      [cmap]=buildcmap(colormap_strings{10});
%                      colormap(h1,cmap)

                     colormap(h2,brewermap([],colormap_brewer_string{11}))
                     
                     title('500mb relative vorticity (shaded)')
                     
                 % 700 temp

                     h4 = subplot(4,1,3);
                     hold on

                     num_conts = 13;

                     contour_var_to_plot = 3;
                     anom_var_to_plot = 11;
                     wind_var_to_plot_U = 6;
                     wind_var_to_plot_V = 9;

                     axesm('robinson',...
                     'Frame', 'on',...
                     'Grid', 'off',...
                     'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
                     'maplonlim',[110 340])
                     tightmap
                     
%                      max_val = 0.9*max(max(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot)));
%                      min_val = 1.1*max(max(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot)));

                     levels = linspace(245,300,20);
                     contourfm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot)),levels,'Linestyle','none')
                     %contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot)),5,'LineWidth',2,'LineColor','w')

                        %stipple
                        %scatterm(stipple_lat_lons_NH_var2co(:,1,anom_var_to_plot),stipple_lat_lons_NH_var2co(:,2,anom_var_to_plot),1,'k','filled')

                     caxis([min(levels) max(levels)])

                     contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,contour_var_to_plot)),num_conts,'LineColor','k')
                     quiverm(NH_lat_era5_grid,...
                             NH_lon_era5_grid,...
                             squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_V))',...
                             squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_U))',...
                             'k')

                     geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
                     linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-k','LineWidth',3)

                      t=colorbar;
                      set(get(t,'ylabel'),'String', 'degrees K');
                      
%                      [cmap]=buildcmap(colormap_strings{1});
%                      colormap(h1,cmap)
%                      colormap(h4,parula)
                       colormap(h4,jet)

                     title('700mb temperature (shaded)')
                     
                 % 700 mb vertical velocity
                     
                     h5 = subplot(4,1,4);
                     hold on
                     
                     contour_var_to_plot = 3;
                     anom_var_to_plot = 12;
                     wind_var_to_plot_U = 6;
                     wind_var_to_plot_V = 9;

                     vv_plot_field = squeeze(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot));
                     
                      K = (1/9)*ones(2);
                      vv_plot_field_smooth = conv2(vv_plot_field,K,'same');
                     
                     num_conts = 13;

                     axesm('robinson',...
                     'Frame', 'on',...
                     'Grid', 'off',...
                     'maplatlim',[min(NH_lat_era5) max(NH_lat_era5)],...
                     'maplonlim',[110 340])
                     tightmap
                     
                    % max_val = 0.9*max(max(abs(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot))));

                     levels = linspace(-0.2,0.2,30);
                     contourfm(NH_lat_era5,NH_lon_era5,vv_plot_field_smooth,levels,'Linestyle','none')
                     %contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,anom_var_to_plot)),5,'LineWidth',2,'LineColor','w')

                        %stipple
                        %scatterm(stipple_lat_lons_NH_var2co(:,1,anom_var_to_plot),stipple_lat_lons_NH_var2co(:,2,anom_var_to_plot),1,'k','filled')

                     caxis([min(levels) max(levels)])

                     contourm(NH_lat_era5,NH_lon_era5,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,contour_var_to_plot)),num_conts,'LineColor','k')
                     quiverm(NH_lat_era5_grid,NH_lon_era5_grid,squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_V))',squeeze(extreme_array_gridded_NH_comp_var2co(:,:,wind_var_to_plot_U))')

                     geoshow(coast.coastlat, coast.coastlon, 'Color', 'black')
                     linem([min(WECC_lat); min(WECC_lat); max(WECC_lat); max(WECC_lat); min(WECC_lat)],[min(WECC_lon); max(WECC_lon); max(WECC_lon); min(WECC_lon); min(WECC_lon)],'-k','LineWidth',3)

                      t=colorbar;
                      set(get(t,'ylabel'),'String', 'hPa/s');
                      
%                      [cmap]=buildcmap(colormap_strings{1});
%                      colormap(h1,cmap)
%                      colormap(h5,parula)
                     [cmap]=buildcmap(colormap_strings{6});
                     colormap(h5,cmap)

                    %colormap(h5,brewermap([],'BrBG'))
                     
                     title('Omega (shaded)')
                     
end
