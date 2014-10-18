function hazard=climada_meteorite_hazard_set
% climada
% NAME:
%   climada_meteorite_hazard_set
% PURPOSE:
%   generate the meteorite hazard event set
%
%   first, create a centroids file prior to calling this code:
%   ...data/system/Meteorite_centroids.xls
%
%   for 1'000'000 events:
%       generation of 1000000 events took 0.107830 sec
%       latitude/longitude conversion took 0.620342 sec
%       diameter conversion took 1.090989 sec
%       CalculationUnit interpolation took 13996.156000 sec (almost 4h)
%       37272 (100.00%) CalculationUnits hit (by 1390444 events)
%       matrix density 0.0037%
%       storing as C:\Data\climada\climada_data\hazards\MEXX_B_Probabilistic_1e6.mat ...
%       for the Test ptf: loss calculation took 10.275969 seconds, based
%       upon 5000000 events @ 279 locations
% CALLING SEQUENCE:
%   hazard=climada_meteorite_hazard_set
% EXAMPLE:
%   hazard=climada_meteorite_hazard_set
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   hazard: a climada hazard set structure, stored automatically also 
%       in ...data/hazards... struct fields are:
%       lon(centroid_i): the longitude of each centroid
%       lat(centroid_i): the latitude of each centroid
%       centroid_ID(centroid_i): a unique ID for each centroid
%       peril_ID: just an ID identifying the peril, e.g. 'TC' for
%       tropical cyclone or 'ET' for extratropical cyclone
%       comment: a free comment, normally containing the time the hazard
%           event set has been generated
%       orig_years: the original years the set is based upon
%       orig_event_count: the original events
%       event_count: the total number of events in the set, including all
%           probabilistic ones, hence event_count>=orig_event_count
%       orig_event_flag(event_i): a flag for each event, whether it's an original
%           (1) or probabilistic (0) one - here all=0
%       event_ID: a unique ID for each event
%       date: the creation date of the set
%       arr(event_i,centroid_i),sparse: the hazard intensity of event_i at
%           centroid_i
%       frequency(event_i): the frequency of each event
%       matrix_density: the density of the sparse array hazard.intensity
%       filename: the filename of the hazard event set (if passed as a
%           struct, this is often useful)
% MODIFICATION HISTORY: 
% David N. Bresch, david.bresch@gmail.com, 20130317, based on work back in 20080801
%-

hazard=[]; % init output

global climada_global
if ~climada_init_vars(1),return;end % init/import global variables

% poor man's version to check arguments
%%if ~exist('param1','var'),param1=[];end

% PARAMETERS
%
hazard.peril_ID='ME'; % ME for MEteorite
hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep 'ME_hazard.mat'];
centroids_file='Meteorite_centroids_small.xls'; % next line only to add path
centroids_file=[climada_global.modules_dir filesep 'meteorite' filesep 'data' filesep 'system' filesep centroids_file]; % add path
%
% whether we show the check plots (default=1)
% if =1, only key plots, if =2 all plots, if=0 no plots
% no plots, if more than 100000 events
show_check_plots=2; % default=1
%
matrix_density=0.0001; % decimal, hence 0.01 means 1% density of hazard.intensity
%
% some hazard parameters
hazard.event_count=100000; % default 1000000
hazard.event_years=hazard.event_count/10; % to produce at least a 1000yr loss
hazard.comment=sprintf('ME meteorite hazard event set, (c) David N. Bresch, %s',datestr(now));
%
rand_state=0; % see help rand

if hazard.event_count>100000
    show_check_plots=0;
    fprintf('WARNING: plots suppressed due to too many events (>100000)\n');
end

% init hazard
hazard.event_ID=1:hazard.event_count;
hazard.orig_event_flag=hazard.event_ID*0;
hazard.frequency=ones(1,hazard.event_count)/hazard.event_years;
hazard.reference_year=2000; % hard-wired, unimportant
hazard.date = datestr(now);

% generate impact location and diameter
tic;
rand('state',rand_state); % init
meteor_impact.raw=rand(hazard.event_count,3);
fprintf('generation of %i events took %f sec\n',hazard.event_count,toc)

% convert random numbers to intervals
% -----------------------------------

if show_check_plots>1,figure('Name','Conversions','Color',[1 1 1]);end

tic
% longitude (easy)
meteor_impact.lon=meteor_impact.raw(:,1)*360-180; % -180..180 degrees
if show_check_plots>1
    subplot(3,2,1);
    hist(meteor_impact.lon,11);title('Longitude histogram');
    xlabel('Longitude');ylabel('number events');
end % show_check_plots

% latitude (cos phi is latitude element)
x0=0:pi/2000:pi/2;y0=cos(x0).^(1/2);
x2=x0;x1=x0-x0(end);
y2=y0;y1=2-y0(end:-1:1);
x=[x1 x2];
y=[y1 y2];
[x,unique_pos]=unique(x);y=y(unique_pos);
x=(x-min(x));x=x/max(x); % normalize to 1
y=(y-min(y));y=y/max(y); % normalize to 1
y=-(y*180-90); % convert to latitude
meteor_impact.lat=interp1(x,y,meteor_impact.raw(:,2));
fprintf('latitude/longitude conversion took %f sec\n',toc)
if show_check_plots>1
    subplot(3,2,2);
    plot(y,x); % to check the latitude conversion
    title('Latitude conversion (CDF)');
    axis tight;xlabel('latitude');ylabel('random input 0..1');
    subplot(3,2,3);
    hist(meteor_impact.lat,11);title('Latitude histogram');
    xlabel('Latitude');ylabel('number events');
    subplot(3,2,4);
    plot(meteor_impact.lon,meteor_impact.lat,'.r','MarkerSize',1);
    title('Impact center map');
    fprintf('map borders ...\n');
    hold on;climada_plot_world_borders;
end % show_check_plots

% diameter
tic
x=0.0001:0.0001:1;
y=-(log2(x).^2);
x=(x-min(x));x=x/max(x); % normalize to 1
y=(y-min(y));y=y/max(y); % normalize to 1
x=x*1000; % max 1000 km
meteor_impact.diam_km=interp1(y,x,meteor_impact.raw(:,3));
fprintf('diameter conversion took %f sec\n',toc)
if show_check_plots>1
    subplot(3,2,5);
    plot(x,y);title('Diameter conversion (CDF)');
    axis tight;xlabel('diameter [km]');ylabel('random input 0..1');
    subplot(3,2,6);
    hist(meteor_impact.diam_km,100);title('Diameter histogram');
    xlabel('diameter [km]');ylabel('number events');
end % show_check_plots

if show_check_plots>0
    % and the final figure
    MeteoriteMap_figure=figure('Name','Meteorite map','Color',[1 1 1]);
    circle_diam=[];
    circle_format=[];
    marker_size=1;
    marker_format='.b';
    climada_circle_plot(meteor_impact.diam_km,meteor_impact.lon,meteor_impact.lat,'Meteorite impact map',circle_diam,circle_format,marker_size,marker_format);
    hold on
end % show_check_plots

% load CalculationUnits, interpolate to
if exist(centroids_file,'file')
    centroids=climada_centroids_read(centroids_file);

    %     if show_check_plots>0
    %         % add centroids
    %         hold on
    %         plot(centroids.Longitude,centroids.Latitude,'.g','MarkerSize',1);
    %     end % show_check_plots

    hazard.centroid_ID=centroids.centroid_ID';
    hazard.centroids_file=centroids_file;
    hazard.lat=centroids.Latitude';
    hazard.lon=centroids.Longitude';
    hazard.severity=meteor_impact.diam_km'; % the original diameter
    % next line is a bit redundant (a target for memory optimization if needed)
    hazard.meteor_impact=meteor_impact; % original data for storage (as long as space allows)
    n_centroids=length(hazard.centroid_ID);
    hazard.intensity=spalloc(hazard.event_count,n_centroids,fix(n_centroids*hazard.event_count*matrix_density));

    % interpolate hazard to CalculationUnit
    msgstr=sprintf('Interpolation to %i centroids...',n_centroids);
    h = waitbar(0,msgstr);
    fprintf('%s\n',msgstr);
    cent_hit=0;cent_hit_abs=0;
    cent_hit_pos=zeros(1,n_centroids);
    t0 = clock;
    for cent_i=1:n_centroids
        waitbar(cent_i/n_centroids,h);
        
        if mod(cent_i,1000)==0 % show time est. every 1000th CU
            est_time_rem_min=(n_centroids-cent_i)*etime(clock,t0)/cent_i/60;
            fprintf('estimated time remaining %3.2f min\n',est_time_rem_min)
            msgstr=sprintf('time remaining %3.2f min',est_time_rem_min);
            waitbar(cent_i/n_centroids,h,msgstr);

        else
            waitbar(cent_i/n_centroids,h);
        end
                       
        dist_km=climada_distance_km(hazard.lon(cent_i),hazard.lat(cent_i),meteor_impact.lon,meteor_impact.lat);
        pos=find((dist_km-meteor_impact.diam_km)<0); % CU within impact distance
        if length(pos)>0
            hazard.intensity(pos,cent_i)=dist_km(pos);
            cent_hit=cent_hit+1;
            cent_hit_abs=cent_hit_abs+length(pos);
            cent_hit_pos(cent_i)=1;
            %fprintf('CU %i (%i): %i impacts\n',hazard.centroid_ID(cent_i),cent_i,length(pos));
        end
    end % cent_i
    fprintf('centroids interpolation took %f sec\n',etime(clock,t0))
    close(h)
    hazard.comment2=sprintf('%i (%2.2f%%) centroids hit (by %i events)',cent_hit,cent_hit/n_centroids*100,cent_hit_abs);
    fprintf('%s\n',hazard.comment2);

    if show_check_plots>0 && sum(cent_hit_pos)>0% add hit centroids to plot
        figure(MeteoriteMap_figure)
        cent_hit_pos=logical(cent_hit_pos); % to use as index
        plot(hazard.lon(cent_hit_pos),hazard.lat(cent_hit_pos),'xg','MarkerSize',3);
    end % show_check_plots

    hazard.matrix_density=nnz(hazard.intensity)/numel(hazard.intensity);
    fprintf('matrix density %2.4f%%\n',hazard.matrix_density*100);

    if show_check_plots>0
        figure('Name','Meteorite hazard array','Color',[1 1 1])
        spy(hazard.intensity); hold on
        title('Meteorite hazard array');
        xlabel('events');ylabel('centroids');
    end % show_check_plots

    % save hazard event set
    hazard.filename=hazard_set_file;
    fprintf('storing as %s ...\n',hazard_set_file);
    save(hazard_set_file,'hazard','-v6'); % v6: un-zipped for faster access (load)
else
    fprintf('WARNING: centroids file %s does not exist\n',centroids_file);
end

return