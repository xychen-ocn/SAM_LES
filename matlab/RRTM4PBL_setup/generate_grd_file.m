% This script will set up the vertical grid level for SAM. 
% in the "grd" file, I will list all the grid levels out. 
% the number of grd levels will be equal to nz_gl defined in the
% SRC/domain.f90 
% note, the grd level specified here is the mid-levels (in meters) for the
% Arakawa-C grid, not the interfacial level. 

clear all; close all;
optname = 'streched';
dz = 40;           % units: m (same for the rest of the height related params)
zmid0= 20;
PBL_lev = 4E3;
domain_top = 6E3;

%% method 1: set different resolution in different layers:
if strcmp(optname,'discrete')
    % require input of number of discrete vertical resolution
    vres=[dz, 100, 200];
    %z_mid_PBL = zmid0:dz:PBL_lev;
    a=10E3;
    layers_top = [PBL_lev,a, domain_top];
    z0_layers=[zmid0, PBL_lev, a];  % meant to be input to the function.
    
    nlayers = length(vres);
    for i=1:nlayers       
        layers(i).zmid = z0_layers(i):vres(i):layers_top(i);
    end
    %z_mid_upper = z_mid_tropo(end):vres(3):domain_top;
    
    z_mid_all = unique([layers.zmid]);
    grid_parm = [num2str(len(vres)) 'vertRes'];
    
    %% method 2: stretching vertical resolution above certain level:
elseif strcmp(optname,'streched')
    % require input of the streching coeff
    strcoeff=1.05;
    z_mid_PBL = zmid0:dz:PBL_lev;
    i = 2;
    z_mid_upper=[];
    z_mid_upper(1) = z_mid_PBL(end);
    
    while i>0
        z_mid_upper(i) = z_mid_upper(i-1) + dz*strcoeff^(i-1);
        dz_streched(i-1) = dz*strcoeff^(i-1);
        if z_mid_upper(i)>=domain_top
            break
        else
            i=i+1;
        end
    end
    
    z_mid_all = unique([z_mid_PBL, z_mid_upper]);
    grid_parm = num2str(strcoeff*100);

end

%% make a plot to confirm:
x = repmat([0;1],1,length(z_mid_all));
y = repmat(z_mid_all,2, 1);

figure()
subplot(2,1,1)
plot(x,y,'-k');
hold on;
set(gca,'yscale','linear');

subplot(2,1,2)
hold on;
plot(z_mid_PBL, ones(size(z_mid_PBL))*dz,'-b');
plot(z_mid_upper(2:end), dz_streched)
xlabel('height (m)')
ylabel('dz (m)');


%% save data
svdir = ['/Users/xchen/Documents/GitHub/SAM_LES/input'];
if ~exist(svdir,'dir')
    mkdir(svdir)
end
grd_FN = ['grd_dz' num2str(dz) 'm' '_' num2str(domain_top/1E3) 'km' '_' optname '_' grid_parm];
fid = fopen([svdir filesep grd_FN ],'w');
fprintf(fid,'%6.1f\n',z_mid_all);
fclose(fid)

type([svdir filesep grd_FN])