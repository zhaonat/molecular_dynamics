%% Diffusion Coefficient Analysis
cd('D:\Nathan\Documents\StanfordYearOne\AdvancedNumericalMethods\Homework4\Problem2')
load('AutocorrData50Temps')

LinFunc = @(x, m, b) m*x+b
tf = 30; t0 =0; h = (tf-t0)/3000;
 Time = 0:h:tf-h;
ValidAVG = AVG(:,1:48);
figure;
plot(ValidAVG)
title('Autocorrelation plots for Multiple Temperatures')
xlabel('Time')
ylabel('Autocorrelation Values')

%% Get Temperature Data
temps = [];
NumParticles = 64;
for T =  0.05:0.02:2.0
    vxVec = T*randn(NumParticles,1);
    vyVec = T*randn(NumParticles,1);
    vxinit = vxVec; vyinit = vyVec;

    vxMeanSq = mean(vxinit.^2);
    vyMeanSq = mean(vyinit.^2);
    vxrms = sqrt(vxMeanSq); vyrms = sqrt(vxMeanSq);
    vrms = sqrt(vxrms^2 + vyrms^2);
    T = vrms^2/3
    temps = [temps, T];
end

Tused = temps(1:48);
%% Extract Slope for Diffusion Coefficients
figure
diffCoeffs = [];
for i = 1:48
   data = AVG(:, i);
   [r,m,b] = regression(Time.',data);
   y = mean(m)*Time.'+b;
   plot(Time, y)
   hold on;
   diffCoeffs = [diffCoeffs, mean(m)];
   
end

figure;
plot(Tused, diffCoeffs, '.', 'markersize', 15)