clear all; clc; close all;
global snd_FN
% attempt to create profile to test RRTM4PBL set up;

% 1. verify again whether or not the monthly sounding fits the RICO profile
% I feel like it doesn't fit at all... I think that is okay though because
% it is not my purpose to have accurate forcing with this case, it is more
% like to enable the code to have ability to run with setup for the PBL
% itself.


workdir = '/Users/xchen/Documents/SAM/RICO/RRTM4PBL_setup';
ERA5_forcingFN= 'era5_pressure_levels_monthly_means_rico_region_Jan2005.nc';
absFN = [workdir filesep ERA5_forcingFN];
ERA5 = read_netCDF_into_matlab_structure(absFN);
p_mat = repmat(ERA5.level*100,1,length(ERA5.longitude), length(ERA5.latitude), 2);
ERA5.theta = compute_potential_temp( ERA5.t, permute(double(p_mat),[2,3,1,4]));

grav = 9.8; % m/s^2;
% obtain an area mean:
lonmask = ERA5.longitude==-62;
latmask = ERA5.latitude==17.5;
ERA5_parms = fieldnames(ERA5);
for i = 1:length(ERA5_parms)
    FN = ERA5_parms{i};
    if ~strcmp(FN,'longitude') & ~strcmp(FN,'latitude') & ~strcmp(FN, 'level') &  ~strcmp(FN, 'time') & ~strcmp(FN,'time_num')
        tmp = ERA5.(FN)(lonmask,latmask,:,:);
        ERA5_areamean.(FN) = squeeze(mean(tmp,[1,2],'omitnan'));
    end
end

% convert temperature to potential temperature:

 snd=read_SAM_snd_file('./RICO/snd');
% % so now, I will try making use of both profile:
 dz_RICO  = snd.z(2) - snd.z(1);
% % build another vertical axis with the same vertical increment as RICO
% % sounding;
% % interpolate the ERA5 sounding to the high resolution vertical axis.
% Z_common = min(snd.z):dz_RICO:40E3;

path2grdfile = './SAM_input/grd_streched_1.05.txt';
Z_common = readmatrix(path2grdfile)';

%Z_common = read_SAM_grd_file(path2grdfile);   % needs to be developed.

% interpolation:
for i = 1:length(ERA5_parms)
    FN = ERA5_parms{i};
    if ~strcmp(FN,'longitude') & ~strcmp(FN,'latitude') & ~strcmp(FN, 'level') &  ~strcmp(FN, 'time') & ~strcmp(FN,'time_num')
        ERA5_intp.(FN) = interp1(ERA5_areamean.z(:,1)./grav, ERA5_areamean.(FN)(:,1), Z_common','linear');
    elseif strcmp(FN,'level')
        ERA5_intp.(FN) = interp1( ERA5_areamean.z(:,1)./grav, (double(ERA5.level)),Z_common', 'linear');
    end
end

clear RICO_intp;
snd_parms = fieldnames(snd);
for i = 1:length(snd_parms)
    FN = snd_parms{i};
    if ~strcmp(FN, 'day') & ~strcmp(FN,'pres0')
        for iz = 1:2
            RICO_intp.(FN)(:,iz) = interp1(snd.z(:,iz), snd.(FN)(:,iz),Z_common, 'linear','extrap');
        end
    end
end


hfig_era = figure(1);clf;
set(hfig_era, 'name','ERA5 monthly profiles');
% subplot(2,3,1);
% plot(ERA5_intp.theta(:,1), Z_common./1E3,'.-');
% hold on
% plot(RICO_intp.theta(:,1), Z_common./1E3, '.-r');
% plot(theta, Z_common./1E3, '.-');

for i = 1:5
    if i==1
        PN = 'level';
        
    else
        
        PN  = snd_parms{3+i};
    end
    subplot(2,3,i)
    plot(ERA5_intp.(PN)(:,1), Z_common,'.-');
    hold on;
    title(PN)
    ylabel('geopotential (m)');
    
    %legend('day0','day100');
end


% plot sounding data:

% hfig_snd = figure(2); clf;
% set(hfig_snd, 'name','RICO snd');
for i = 1:5
    PN  = snd_parms{3+i};
    subplot(2,3,i)
    plot(snd.(PN), snd.z,'.-');
    title(PN)
    ylabel('geopotential (m)');
    legend('day0','day100');
end
    
    

%% now merge the two profiles:
clear snd_ext

wgt_ERA = min( 1, max(0, (Z_common-2000)/(6500-2000)) );   % 4500 to 6500

%wgt_ERA = reshape(wgt_ERA,size(T_RICO));
for s = 1:size(snd.z,2)
    for i = 1:4
        PN = snd_parms{4+i};
        snd_ext.(PN)(:,s) = RICO_intp.(PN)(:,s) + wgt_ERA'.*(ERA5_intp.(PN) - RICO_intp.(PN)(:,s));
    end
    
    % now merge pressure separatly: taking directly from the ERA5 since
    % RICO pressure is -999 (NaN).
    snd_ext.p(:,s) = ERA5_intp.level;
end

nd = datenum('2004-12-01') - datenum('2004-01-01');
snd_ext.day = [1+nd ,  1+nd+31];
snd_ext.pres0 = max(snd_ext.p,[],1);
snd_ext.z = repmat(Z_common',1,2);




% display the merged profile:
figure(1);
hold on;
for i = 1:5
    PN  = snd_parms{3+i};
    subplot(2,3,i)
    plot(snd_ext.(PN)(:,1), snd_ext.z,'.-r');
    title(PN)
    ylabel('geopotential (m)');
    %legend('day0','day100');
end
   
% looking good here.



% write the extended data into a snd file again;
svdir = 'SAM_input';
if ~exist(svdir,'dir')
    mkdir(svdir)
end
nlevs = size(snd_ext.z,1);
headerstr = 'z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]';
%snd_FN = 'snd_RRTM4PBL';
snd_FN = 'snd_RRTM_40E3_streched';

fid = fopen([svdir filesep snd_FN], 'w');
for s = 1:size(snd_ext.z,2)
    % first line
    if s==1
    fprintf(fid, '%s\n', headerstr);
    end
    % second line
    fprintf(fid, '%4.1f, %i, %7.2f ', snd_ext.day(s), nlevs, 1015.40);  %snd_ext.pres0(s)
    fprintf(fid, '%s, %s, %s\n', 'day', 'levels','pres0');
    % start to put the records in here:
    snd_matrix = [snd_ext.z(:,s), snd_ext.p(:,s)*0-999, snd_ext.theta(:,s), ...
                  snd_ext.q(:,s)*1000, snd_ext.u(:,s), snd_ext.v(:,s)];
    fprintf(fid, '%6.1f  %7.2f  %6.2f  %5.2f  %6.2f  %6.2f\n', snd_matrix');
end
fclose(fid);

type([svdir filesep snd_FN])

% save the merged data just in case:
save([svdir filesep 'ERA5_merged_RICO_sounding_40E3_streched_grid.mat'], 'snd_ext','ERA5_intp');


% then apply the matlab code shared by peter to generate the netCDF file.
