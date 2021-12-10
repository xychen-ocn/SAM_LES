function snd=read_SAM_snd_file(absFN)
% function: read SAM sounding file
fid = fopen(absFN,'r');
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

end