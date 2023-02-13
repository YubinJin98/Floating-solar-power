function [Results, ACPower, TMYData] = FSPCalculate(res, T, WS, CERES, CountryNumber)
%Calculate floating solar power
%   need data of reservoir, weather condition from ERA5, solar radiation
%   from CERES and CountryNumber


%% Prepare reservoir dataset
for i = 1:length(res)
    res(i).ID = i-1;
end
ID = double([res.ID])';
X_centroid = double([res.X_centroid])';
Y_centroid = double([res.Y_centroid])';
ResArea = double([res.Area])'; % square km
X_centroid2 = X_centroid;
X_centroid2(X_centroid<0) = X_centroid2(X_centroid<0) + 360;


%% Prepare solar radiation data
DNI_raw = CERES.DNI_raw;
DHI_raw = CERES.DHI_raw;
lon = CERES.lon;
lat = CERES.lat;
time = CERES.time;
time_char = datestr(double(time)+datenum('03/01/2000'),'dd-mm-yyyy HH:MM:SS');
Time = struct('year',[],'month',[],'day',[],'hour',[],'minute',[],'second',[],'UTCOffset',[]);
for t = 1:length(time)
    Time.year(t,1) = str2double(time_char(t,7:10));
    Time.month(t,1) = str2double(time_char(t,4:5));
    Time.day(t,1) = str2double(time_char(t,1:2));
    Time.hour(t,1) = str2double(time_char(t,12:13));
    Time.minute(t,1) = str2double(time_char(t,15:16));
    Time.second(t,1) = str2double(time_char(t,18:19));
    Time.UTCOffset(t,1) = 0;
end

DNI = nan(length(ID),length(time)); DHI = DNI;
for n = 1:length(ID)
    I = find(lat >= Y_centroid(n), 1, 'first');
    J = find(lon >= min(X_centroid2(n),lon(end)), 1, 'first');
    DNI(n,:) = squeeze(DNI_raw(J,I,:));
    DHI(n,:) = squeeze(DHI_raw(J,I,:));
end

Location = struct('latitude',[],'longitude',[]);
Location = repmat(Location,[length(ID),1]);
SunEl = nan(length(ID),length(time));
for n = 1:length(ID)
    Location(n,1).latitude = Y_centroid(n);
    Location(n,1).longitude = X_centroid2(n);
    [~, SunEl(n,:), ~, ~] = pvl_ephemeris(Time, Location(n,1));
end

GHI = DHI + DNI.*cosd(90-SunEl);
DNI(SunEl<0) = 0;
DHI(SunEl<0) = 0;
GHI(SunEl<0) = 0;

TMYData = struct('DateNumber',[],'SiteTimeZone',[],'DNI',[],'DHI',[],'GHI',[]...
    ,'SiteLatitude',[],'SiteLongitude',[],'Pressure',[],'DryBulb',[],'Wspd',[]);
TMYData = repmat(TMYData, [length(ID),1]);
for n = 1:length(ID)
    TMYData(n,1).DateNumber = double(time)+datenum('03/01/2000');
    TMYData(n,1).SiteTimeZone = 0;
    TMYData(n,1).DNI = DNI(n,:)';
    TMYData(n,1).DHI = DHI(n,:)';
    TMYData(n,1).GHI = GHI(n,:)';
    TMYData(n,1).SiteLatitude = Y_centroid(n);
    TMYData(n,1).SiteLongitude = X_centroid(n);
    TMYData(n,1).Pressure = 1013.25;
    TMYData(n,1).DryBulb = T(n,:)';
    TMYData(n,1).Wspd = WS(n,:)';
end
idx = isnan(sum(T,2));
TMYData(idx) = [];
ID(idx) = [];
X_centroid(idx) = [];
Y_centroid(idx) = [];
ResArea(idx) = [];


%% Calculate floating solar power
ModuleParameters = pvl_sapmmoduledb(523,'SandiaModuleDatabase_20120925.xlsx');

InverterDatabase = load('SandiaInverterDatabaseSAM2014.1.14.mat');
Inverter = InverterDatabase.SNLInverterDB(1461);
clear InverterDatabase

ACPower = nan(length(time),length(ID));
Array_Tilt = nan(length(ID),1); Array_Ms = nan(length(ID),1); Array_Mp = nan(length(ID),1);
for n = 1:length(ID)
    [ACPower(:,n), Array, ~] = PV_LIB_Model(TMYData(n,1), ModuleParameters, Inverter);
    Array_Tilt(n,1) = Array.Tilt;
    Array_Ms(n,1) = Array.Ms;
    Array_Mp(n,1) = Array.Mp;
end
ACPower(ACPower < 0) = 0;
dEnergy = sum(ACPower,1)'./(Array_Ms.*Array_Mp.*ModuleParameters.area.*(cosd(Array_Tilt)+1.2*sind(Array_Tilt)))/10^3; % kWh


%% Organize results
Results = struct('CountryNumber',[],'Res_ID',[],'X_centroid',[],'Y_centroid',[],...
    'ResArea',[],'dEnergy',[],'Time',[]);
Results.CountryNumber = CountryNumber;
Results.Res_ID = ID;
Results.X_centroid = X_centroid;
Results.Y_centroid = Y_centroid;
Results.ResArea = ResArea;
Results.dEnergy = dEnergy;
Results.Time = Time.year(1);


end

