%% Water savings with 30% coverage and max area = 30 km2

vap_all = ncread('cru_ts4.05.2001.2010.vap.dat.nc','vap');
tmx_all = ncread('cru_ts4.05.2001.2010.tmx.dat.nc','tmx');
tmn_all = ncread('cru_ts4.05.2001.2010.tmn.dat.nc','tmn');
tmp_all = ncread('cru_ts4.05.2001.2010.tmp.dat.nc','tmp');
cld_all = ncread('cru_ts4.05.2001.2010.cld.dat.nc','cld');

vap_all = cat(3, vap_all, ncread('cru_ts4.05.2011.2020.vap.dat.nc','vap'));
tmx_all = cat(3, tmx_all, ncread('cru_ts4.05.2011.2020.tmx.dat.nc','tmx'));
tmn_all = cat(3, tmn_all, ncread('cru_ts4.05.2011.2020.tmn.dat.nc','tmn'));
tmp_all = cat(3, tmp_all, ncread('cru_ts4.05.2011.2020.tmp.dat.nc','tmp'));
cld_all = cat(3, cld_all, ncread('cru_ts4.05.2011.2020.cld.dat.nc','cld'));

lon_vap = ncread('cru_ts4.05.2011.2020.vap.dat.nc','lon');
lat_vap = ncread('cru_ts4.05.2011.2020.vap.dat.nc','lat');

yr = 2001; % 2001~2020

vap = squeeze(vap_all(:,:,(yr-2001)*12+1:(yr-2001)*12+12));
tmx = squeeze(tmx_all(:,:,(yr-2001)*12+1:(yr-2001)*12+12));
tmn = squeeze(tmn_all(:,:,(yr-2001)*12+1:(yr-2001)*12+12));
tmp = squeeze(tmp_all(:,:,(yr-2001)*12+1:(yr-2001)*12+12));

Rs_Name = ['CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed4.1_Subset_',num2str(yr),'01-',num2str(yr),'12.nc'];
sw = ncread(Rs_Name, 'ini_sfc_sw_down_all_mon');
lon_sw = ncread(Rs_Name, 'lon');
lat_sw = ncread(Rs_Name, 'lat');

uv_Name = ['uv10_monthly_',num2str(yr),'.nc'];
u10 = ncread(uv_Name, 'u10');
v10 = ncread(uv_Name, 'v10');
lon_u = ncread(uv_Name, 'longitude');
lat_u = ncread(uv_Name, 'latitude');

Epen = nan(2,12);

lon = X_centroid;
lat = Y_centroid;
x = PanelArea/ReservoirArea;

for t = 1:12

    % data from CRU
    J = find(lon_vap >= lon, 1, 'first');
    I = find(lat_vap >= lat, 1, 'first');
    ea = vap(J, I, t)/10;
    Tmax = tmx(J, I, t);
    Tmin = tmn(J, I, t);
    Tm = tmp(J, I, t);
    if Tm < 0
        Epen(:,t) = 0;
        continue
    end

    % data from CERES
    J = find(lon_sw-180 >= lon, 1, 'first');
    I = find(lat_sw >= lat, 1, 'first');
    Rs = sw(J, I, t)*3600*24*10^-6;

    % data from ERA5
    if lon >= 0
        J = find(lon_u >= lon, 1, 'first');
    else
        J = find(lon_u >= min(lon+360,lon_u(end)), 1, 'first');
    end
    I = find(lat_u <= lat, 1, 'first');
    U = u10(J, I, t);
    V = v10(J, I, t);
    WS = sqrt(U^2+V^2);
    u2 = WS*4.87/(log(67.8*10-5.42));

    if isnan(Tm + Rs + u2)
        continue
    end

    % other parameters
    e_Tmax = 0.6108*exp(17.27*Tmax/(Tmax+237.3));
    e_Tmin = 0.6108*exp(17.27*Tmin/(Tmin+237.3));
    es = (e_Tmax+e_Tmin)/2;
    D = es-ea;
    D(D<0) = 0;

    DELTA = 4098*(0.6108*exp(17.27*Tm/(Tm+237.3)))/(Tm+237.3)^2;
    gamma = 0.665*10^-3*101.325;
    lamda = 2.501-0.002361*Tm;
    fu = 1+0.536*u2;

    Gsc = 0.082; J = t*31-15; phi = pi/180*lat;
    dr = 1+0.033*cos(2*pi/365*J);
    delta = 0.409*sin(2*pi/365*J-1.39);
    X = 1-(tan(phi)*tan(delta))^2;
    if X <= 0
        X = 0.00001;
    end
    omegas = pi/2-atan(-tan(phi)*tan(delta)/X^0.5);
    Ra = 24*60/pi * Gsc * dr * (omegas*sin(phi)*sin(delta) + cos(phi)*cos(delta)*sin(omegas));
    Ra = max(Ra,0);
    Rso = (0.25+0.5)*Ra; % as = 0.25; bs = 0.5
    f = 1.35*min(Rs/Rso,1)-0.35; % ac = 1.35; bc = -0.35

    alpha = 0.08; sigma = 4.903*10^(-9);
    Rns_free = (1-alpha)*Rs;
    Rnl_free = -f*(0.34-0.14*sqrt(ea))*sigma*(Tm+273.2)^4;
    Rn_free = Rns_free + Rnl_free;
    Rns_fcover = 0;
    Rnl_fcover = -0.1*(0.34-0.14*sqrt(ea))*sigma*(Tm+273.2)^4;
    Rn_fcover = Rns_fcover+Rnl_fcover;

    Epen(1,t) = DELTA/(DELTA+gamma)*Rn_free/lamda +...
        gamma/(DELTA+gamma)*6.43*fu*D/lamda; % free water
    Epen(2,t) = DELTA/(DELTA+gamma)*(Rn_free*(1-x)+Rn_fcover*x)/lamda +...
        gamma/(DELTA+gamma)*6.43*fu*D/lamda; % after covering

end


