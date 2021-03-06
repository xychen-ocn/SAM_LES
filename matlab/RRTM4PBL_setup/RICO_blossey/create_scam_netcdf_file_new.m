function [out] = create_scam_netcdf_file_new(ncfile,comment,...
                                         lev,lon,lat,iyear,calday, ...
                                         phis)
                                         

% FUNCTION [OUT] = CREATE_SCAM_NETCDF_FILE_NEW(NCFILE,COMMENT,...
%                                         LEV,LON,LAT,IYEAR,CALDAY, ...
%                                         PHIS)
% 
% creates a netcdf file suitable for using as an IOP forcing
% dataset for SCAM, the single column version of the Community
% Atmosphere Model (CAM).  Uses Matlab netcdf interface (circa 2019).
%
% Written by Peter Blossey, 200X using SNCTOOL netcdf interface.
% Revised 2019-11 to bring Matlab netcdf interface up to date.

if exist(ncfile)
  eval(['!rm -i ' ncfile])
end

mySchema.Name   = '/';
mySchema.Format = 'classic';

Attributes = { {'author','Peter Blossey, pblossey@u.washington.edu'}, ...
               {'institution', 'University of Washington'}, ...
               {'Conventions', 'CF-1.3'}, ...
               {'date', date}, ...
               {'comment', comment} };
for n = 1:length(Attributes)
  mySchema.Attributes(n).Name = Attributes{n}{1};
  mySchema.Attributes(n).Value = Attributes{n}{2};
end

Dimensions = {{'lat',1,false}, ...
              {'lon',1,false}, ...
              {'lev',length(lev),false}, ...
              {'time',length(calday),false}}; 
for n = 1:length(Dimensions)
  mySchema.Dimensions(n).Name = Dimensions{n}{1};
  mySchema.Dimensions(n).Length = Dimensions{n}{2};
  mySchema.Dimensions(n).Unlimited= Dimensions{n}{3};
end

ncwriteschema(ncfile,mySchema);

% define time variables
time = round(86400*(calday - floor(calday(1))));
whos calday iyear

dv = datevec(calday+datenum(iyear,1,1)); % create date vector.

% Set nbdate
yy = mod(dv(1,1),100);
mm = dv(1,2);
dd = dv(1,3);
nbdate = 1e4*yy + 1e2*dv(1,2) + dv(1,3);

%%%%% ADD DIMENSION VARIABLE (LATITUDE) %%%%%%%%%%
Variables = {...
    {'lat','double',{'lat'},'degrees_north','Latitiude','latitude',lat}, ...
    {'lon','double',{'lon'},'degrees_east','Lonitiude','longitude',lon}, ...
    {'lev','double',{'lev'},'Pa','Pressure level','air_pressure',lev}, ...
    {'time','int32',{'time'},'s','Time in seconds after 00Z on nbdate','time',time}, ...
    {'calday','double',{'time'},'d', ...
     sprintf('Time in days after 00Z on 31 Dec %d',iyear-1),'',calday}, ...
    {'year','int32',{'time'},'year','Year','',dv(:,1)}, ...
    {'month','int32',{'time'},'month','Month','',dv(:,2)}, ...
    {'day','int32',{'time'},'day','Day','',dv(:,3)}, ...
    {'hour','double',{'time'},'hour','Hour','',dv(:,4)}, ...
    {'nbdate','int32',[],'yymmdd','Base date (Note that only two digit year is permitted)','',nbdate}, ...
    {'bdate','int32',[],'yymmdd','Base date (Note that only two digit year is permitted)','',nbdate}, ...
    {'phis','double',{'lat','lon'},'m2/s2','Surface geopotential', ...
     'surface_geopotential',phis} };

for n = 1:length(Variables)
  vname = Variables{n}{1};
  disp(vname)
  nccreate(ncfile,vname, ...
           'Dimensions',Variables{n}{3}, ...
           'Datatype',Variables{n}{2})
  ncwriteatt(ncfile,vname, ...
               'units',Variables{n}{4})
  ncwriteatt(ncfile,vname, ...
             'long_name',Variables{n}{5})
  ncwriteatt(ncfile,vname, ...
             'standard_name',Variables{n}{6})
  ncwrite(ncfile,vname,Variables{n}{7})
end

