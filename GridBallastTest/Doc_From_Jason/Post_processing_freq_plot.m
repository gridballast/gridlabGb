%% Post Processing for Dynamics Model
clear all
clc
close all
format short g

%Files
%Generator
fHandle=fopen('Gen_1_Speed.csv','rt');
RawDataGen1=textscan(fHandle,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',9,'Delimiter',',');
fclose(fHandle);
RawDataGen1{2}=RawDataGen1{2}/(2*pi);

%GFA devices
fHandle=fopen('load_l671_freq.csv','rt');
RawDataL671=textscan(fHandle,'%s %s %f %f %f %f %f %f %f','HeaderLines',9,'Delimiter',',');
fclose(fHandle);

fHandle=fopen('triplex_load_freq.csv','rt');
RawDataTriplexLoad=textscan(fHandle,'%s %s %f %f %f %f %f','HeaderLines',9,'Delimiter',',');
fclose(fHandle);

%Frequency values
fHandle=fopen('node_n630_freq.csv','rt');
RawDataN630=textscan(fHandle,'%s %f','HeaderLines',9,'Delimiter',',');
fclose(fHandle);

fHandle=fopen('swing_node_freq.csv','rt');
RawDataSwing=textscan(fHandle,'%s %f','HeaderLines',9,'Delimiter',',');
fclose(fHandle);

fHandle=fopen('triplex_second_node_freq.csv','rt');
RawDataSecondTrip=textscan(fHandle,'%s %f','HeaderLines',9,'Delimiter',',');
fclose(fHandle);

fHandle=fopen('triplex_meter_freq.csv','rt');
RawDataTripM=textscan(fHandle,'%s %f','HeaderLines',9,'Delimiter',',');
fclose(fHandle);

fHandle=fopen('triplex_main_node_freq.csv','rt');
RawDataOtherTrip=textscan(fHandle,'%s %f','HeaderLines',9,'Delimiter',',');
fclose(fHandle);

%Extract date information
tempData=RawDataGen1{1};
Data1Lengths=cellfun(@length,tempData);
Data1NS=Data1Lengths>23;
DataSpeedTimes=zeros(length(tempData),1);

if (any(Data1NS))
    DataSpeedTimes(Data1NS)=datenum(tempData(Data1NS),'yyyy-mm-dd HH:MM:SS.FFF');
end

if (any(~Data1NS))
    DataSpeedTimes(~Data1NS)=datenum(tempData(~Data1NS),'yyyy-mm-dd HH:MM:SS');
end

DataSpeedTimes=(DataSpeedTimes-DataSpeedTimes(1,1))*(24*3600);

%Create voltages
Voltages671=[complex(RawDataL671{4},RawDataL671{5}) complex(RawDataL671{6},RawDataL671{7}) complex(RawDataL671{8},RawDataL671{9})]/2401.7771;
VoltagesTrip=[complex(RawDataTriplexLoad{4},RawDataTriplexLoad{5}) complex(RawDataTriplexLoad{6},RawDataTriplexLoad{7})]/120;

%% Plots

%Frequency plots
figure('Name','Frequencies','NumberTitle','off');
plot(DataSpeedTimes,[RawDataGen1{2} RawDataL671{3} RawDataTriplexLoad{3} RawDataN630{2} RawDataSwing{2} RawDataSecondTrip{2} RawDataTripM{2} RawDataOtherTrip{2}],'.');
xlabel('Time, seconds');
ylabel('Frequency, Hz');
legend('Gen 1 rotor','Load 671','Triplex Load','Node 630','Swing node','S. Triplex','Triplex Meter','Other Triplex');

%Voltage plots
figure('Name','Voltages - 671','NumberTitle','off');
plot(DataSpeedTimes,abs(Voltages671),'.');
xlabel('Time, seconds');
ylabel('Voltage, puV');

figure('Name','Voltages - Triplex','NumberTitle','off');
plot(DataSpeedTimes,abs(VoltagesTrip),'.');
xlabel('Time, seconds');
ylabel('Voltage, puV');


%% GFA Status
GFAStatusLoadA=strcmp(RawDataL671{2},'TRUE');
GFAStatusLoadB=strcmp(RawDataTriplexLoad{2},'TRUE');

%Plot that
figure('Name','GFA Status','NumberTitle','off');
Axvals=plotyy(DataSpeedTimes,[RawDataL671{3} RawDataTriplexLoad{3}],DataSpeedTimes,[GFAStatusLoadA GFAStatusLoadB]);
ylim(Axvals(2),[-0.5 1.5]);
xlabel('Time, seconds');
ylabel(Axvals(1),'Frequency, Hz');
ylabel(Axvals(2),'On/Off');
legend('L671 Freq','Trip Freq','L671 Status','Trip Status');
