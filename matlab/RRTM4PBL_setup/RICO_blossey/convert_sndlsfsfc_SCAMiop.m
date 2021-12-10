clear all;
thermo_constants;

% get sounding from snd file
snd_FN = 'snd_RRTM_40E3_streched';
fid = fopen(snd_FN,'r');
fgetl(fid)

n=0;
while ~feof(fid)
  n=n+1
  in = fgetl(fid)
  in(in==',') = ' ';
  in_values = sscanf(in,'%g',3);
  if isempty(in_values); break; end
  % NOTE: snd is a structure which holds the various fields
  %         in the snd files (day, pres0, z, p, etc.)
  snd.day(n) = in_values(1);
  snd.pres0(n) = in_values(3);
  a = fscanf(fid,'%g',[6 in_values(2)])';  % 6 fields
  snd.z(:,n) = a(:,1);
  snd.p(:,n) = a(:,2);
  snd.theta(:,n) = a(:,3); % NOTE: theta, not T.
  snd.q(:,n) = a(:,4); % NOTE: q in g/kg.
  snd.u(:,n) = a(:,5);
  snd.v(:,n) = a(:,6);
  fgetl(fid);
end
snd.q = 1e-3*snd.q; % scale q into kg/kg
 
fclose(fid);

% get forcings from lsf file
fid = fopen('lsf','r');
fgetl(fid);

n=0;
while ~feof(fid)
  n=n+1;
  in = fgetl(fid);
  in(in==',') = ' ';
  in_values = sscanf(in,'%g',3);
  if isempty(in_values); break; end
  lsf.day(n) = in_values(1);
  lsf.pres0(n) = in_values(3);
  a = fscanf(fid,'%g',[7 in_values(2)])';
  lsf.z(:,n) = a(:,1);
  lsf.p(:,n) = a(:,2);
  lsf.dt(:,n) = a(:,3); % NOTE: K/s
  lsf.dq(:,n) = a(:,4); % NOTE: kg/kg/s
  lsf.ug(:,n) = a(:,5); % GEOSTROPHIC WIND
  lsf.vg(:,n) = a(:,6); % GEOSTROPHIC WIND
  lsf.wg(:,n) = a(:,7); % LARGE-SCALE VERTICAL MOTION
  fgetl(fid);
end
 
fclose(fid);

% get forcings from sfc file
fid = fopen('sfc','r');
fgetl(fid);

n=0;
a = fscanf(fid,'%g',[5 inf])';
sfc.day = a(:,1);
sfc.sst = a(:,2);
sfc.shf = a(:,3);
sfc.lhf = a(:,4);
sfc.tau = a(:,5);

fclose(fid);

%% Harmonize time dimensions and delete duplicate times
calday = [MakeRowVector(sfc.day) MakeRowVector(snd.day) MakeRowVector(lsf.day)];
calday = sort(calday,'Ascend');
keep = [1];
for n = 2:length(calday)
  if abs(diff(calday(n-1:n)))>0.5/86400
    keep = [keep n];
  end
end
calday = calday(keep);
Ntime = length(calday);
disp(calday')

if min([snd.z(:); lsf.z(:)])>-1
  zsounding = true;
  zz = [snd.z(:); lsf.z(:)];
  zz = sort(zz,'Ascend');
  keep = [1];
  for n = 2:length(zz)
    if abs(diff(zz(n-1:n)))>0.5
      keep = [keep n];
    end
  end
  zz = zz(keep);
  disp(zz);

  Nlev = length(zz);
  zi(1) = 0;
  zi(2:Nlev) = 0.5*(zz(1:Nlev-1)+zz(2:Nlev));
  zi(Nlev+1) = 1.5*zz(Nlev)-0.5*zz(Nlev-1);

  thetav = interp1(snd.z(:,1),snd.theta(:,1).*(1+0.61*snd.q(:,1)),zz, ...
                   'linear','extrap');

  % code taken from pressz.f90 in SAM
  presr(1) = (snd.pres0(1)/1000).^(Rd/Cp);
  presi(1) = snd.pres0(1);
  for k = 2:Nlev+1
    presr(k) = presr(k-1) - g/Cp/thetav(k-1)*(zi(k)-zi(k-1));
    presi(k) = 1000*presr(k)^(Cp/Rd);
    pp(k-1) = exp( log(presi(k-1)) + log(presi(k)/presi(k-1)) ...
                   *(zz(k-1)-zi(k-1))/(zi(k)-zi(k-1)) );
  end

elseif min([snd.p(:); lsf.p(:)])>-1
  zsounding = false;
  
  pp = [psnd(:); plsf(:)];
  pp = sort(pp,'Descend');
  keep = [1];
  for n = 2:length(pp)
    if abs(diff(pp(n-1:n)))>0.01
      keep = [keep n];
    end
  end
  pp = pp(keep);

else
  disp('snd and lsf files do not use same vertical coordinate');
  error(['This will work in SAM but I haven''t coded this up yet for ' ...
         'the conversion of a SCAM netcdf IOP file.'])
end
Nlev = length(pp);

% Re-order pressure coordinate, so that pp(1) is the top of the
% sounding, and convert to Pa.
pp = fliplr(reshape(100*pp,[1 Nlev]));

%%%  REGRID FIELDS ONTO NEW COMBINED TIME/PRESSURE GRID

% First, sfc quantities
wh = {'sst','shf','lhf','tau'};
for m = 1:length(wh)
  out.(wh{m}) = interp1(sfc.day,sfc.(wh{m}),calday,'linear','extrap');
end

% Next, snd quantities
wh = {'theta','q','u','v'};
for m = 1:length(wh)
  tmp = snd.(wh{m});
  % First in height
  clear tmp2
  for k = 1:length(snd.day)
    if zsounding
      tmp2(:,k) = interp1(snd.z(:,k),tmp(:,k),zz);
      tmp2(:,k) = tmp2(end:-1:1,k);
    else
      tmp2(:,k) = interp1(100*snd.p(:,k),tmp(:,k),pp);
    end
  end
  out.(wh{m}) = interp1(snd.day,tmp2',calday,'linear','extrap')';
end
out.pres0 = interp1(snd.day,snd.pres0,calday,'linear','extrap');

% Last, lsf quantities
wh = {'dt','dq','ug','vg','wg'};
for m = 1:length(wh)
  tmp = lsf.(wh{m});
  % First in height
  clear tmp2
  for k = 1:length(lsf.day)
    if zsounding
      tmp2(:,k) = interp1(lsf.z(:,k),tmp(:,k),zz);
      tmp2(:,k) = tmp2(end:-1:1,k);
    else
      tmp2(:,k) = interp1(100*lsf.p(:,k),tmp(:,k),pp);
    end
  end
  out.(wh{m}) = interp1(lsf.day,tmp2',calday,'linear','extrap')';
end

%%%%% OPEN NETCDF FILE FOR FORCINGS %%%%%%%%%%
tmp = strsplit(snd_FN,'_');
caseID = strjoin({tmp{2}, tmp{3}},'_');
nc = ['RICO_forcing_' datestr(now,'dd-mmm-yyyy') '_' caseID '_test.nc'];
comment = ['Forcings for RICO based on idealized case distributed ' ...
           'with SAM and also including a sounding from ERA5 patched ' ...
           'on top of the idealized sounding.  Prepared by X.Chen'];


iyear = 2004;
create_scam_netcdf_file_new(nc,comment,pp,298,18,iyear,calday,0) 

%%%%% TIMESERIES OF SURFACE/TOA FIELDS (FLUXES, SURFACE TEMP, ETC.) %%%%%%%%%%
%% Note that all variables have dimensions {'time','lat','lon'}
%%   and are single precision.
Variables = { ...
    {'Ps','Pa','Surface Pressure','surface_air_pressure'}, ...
    {'Ptend','Pa/s','Surface Pressure Tendency','tendency_of_surface_air_pressure'}, ...
    {'Tg','K','Surface Temperature (SST if over water)','surface_temperature'}, ...
    {'shflx','W/m2','Surface Sensible Heat Flux','surface_upward_sensible_heat_flux'}, ...
    {'lhflx','W/m2','Surface Latent Heat Flux','surface_upward_latent_heat_flux'}, ...
    };
for n = 1:length(Variables)
  disp(Variables{n}{1})
  nccreate(nc,Variables{n}{1}, ...
           'Dimensions',{'lon','lat','time'}, ...
           'Datatype','single')
  ncwriteatt(nc,Variables{n}{1}, ...
             'units',Variables{n}{2})
  ncwriteatt(nc,Variables{n}{1}, ...
             'long_name',Variables{n}{3})
  ncwriteatt(nc,Variables{n}{1}, ...
             'standard_name',Variables{n}{4})
% $$$   ncwrite(nc,Variables{n}{1},...
% $$$           reshape(Variables{n}{5},[Ntime 1 1]));
end

ncwrite(nc,'Ps',reshape(100*out.pres0,[1 1 Ntime]));
ncwrite(nc,'Ptend',reshape(deriv_nonuniform(86400*calday,100*out.pres0),[1 ...
                    1 Ntime]));
ncwrite(nc,'Tg',reshape(out.sst,[1 1 Ntime]));
ncwrite(nc,'shflx',reshape(out.shf,[1 1 Ntime]));
ncwrite(nc,'lhflx',reshape(out.lhf,[1 1 Ntime]));

%%%%% TIMESERIES OF VERTICALLY-VARYING FIELDS (T,q,etc.) %%%%%%%%%%
%% Note that all variables have dimensions {'time','lev','lat','lon'}
%%   and are single precision.
Variables = { ...
    {'u','m/s','Zonal Wind','eastward_wind'}, ...
    {'v','m/s','Meridional Wind','northward_wind'}, ...
    {'ug','m/s','Geostrophic Zonal Wind','geostrophic_eastward_wind'}, ...
    {'vg','m/s','Geostrophic Meridional Wind','geostrophic_northward_wind'}, ...
    {'omega','Pa/s','Vertical Pressure Velocity','lagrangian_tendency_of_air_pressure'}, ...
    {'T','K','Absolute Temperature','air_temperature'}, ...
    {'q','kg/kg','Water Vapor Mass Mixing Ratio',''}, ...
    {'divT','K/s','Large-scale Horizontal Temperature Advection',''}, ...
    {'divq','K/s','Large-scale Horizontal Advection of Water Vapor Mass Mixing Ratio',''}, ...
    };
for n = 1:length(Variables)
  disp(Variables{n}{1})
  nccreate(nc,Variables{n}{1}, ...
           'Dimensions',{'lon','lat','lev','time'}, ...
           'Datatype','single')
  ncwriteatt(nc,Variables{n}{1}, ...
             'units',Variables{n}{2})
  ncwriteatt(nc,Variables{n}{1}, ...
             'long_name',Variables{n}{3})
  ncwriteatt(nc,Variables{n}{1}, ...
             'standard_name',Variables{n}{4})
% $$$   ncwrite(nc,Variables{n}{1},...
% $$$           reshape(Variables{n}{5},[Ntime 1 1]));
end

ncwrite(nc,'u',reshape(out.u,[1 1 Nlev Ntime]))
ncwrite(nc,'v',reshape(out.v,[1 1 Nlev Ntime]))
ncwrite(nc,'ug',reshape(out.ug,[1 1 Nlev Ntime]))
ncwrite(nc,'vg',reshape(out.vg,[1 1 Nlev Ntime]))

pres = pp'*ones(1,Ntime);
out.T = out.theta.*(pres/1e5).^(Rd/Cp);
out.Tv = out.T.*(1 + 0.61*out.q);
out.rho = pres./(Rd*out.Tv);
out.omega = -g*out.rho.*out.wg;

ncwrite(nc,'omega',reshape(out.omega,[1 1 Nlev Ntime]))
ncwrite(nc,'T',reshape(out.T,[1 1 Nlev Ntime]))
ncwrite(nc,'q',reshape(out.q,[1 1 Nlev Ntime]))
ncwrite(nc,'divT',reshape(out.dt,[1 1 Nlev Ntime]))
ncwrite(nc,'divq',reshape(out.dq,[1 1 Nlev Ntime]))

