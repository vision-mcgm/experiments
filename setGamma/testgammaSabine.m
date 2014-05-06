%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Karl Gegenfurtner
% Am besten mittlere EInstellungen. Es sollte im ganz hellen bereich nicht zu 
% Sättigung kommen, und im dunklen Bereich nicht ewig lange bei null rumschwirren.
% Du misst die Gamma-Kurven getrennt für jeden Phosphor, und dann zur Kontrolle 
% nochmal die Summe. Das sollte additiv sein, wenn Du die Hintergrundhelligkeit 
% (bei 0) von jeder Messung abziehst. Die Kurevn für die einzelnen Phosphore sind 
% meist nahezu gleich, es schadet aber nicht das getrennt zu machen. Die mittlere 
% Leuchtdichte (Output 128 mit Gamma-Korrektur, circa 192 ohne) sollte bei circa 
% 30 cd/m2 liegen, auf alle Fälle über 20. Wenn's höher ist (50 oder 60) dann 
% macht es auch nichts, aber dann musst du auf Additivität bei hohen INtensitäten 
% achten.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thorsten Hansen
% Optimal ist eine Helligkeit von 30cd/m^2 bei einem Grauwert von 
% RGB = (190, 190, 190).
% 190 ist der Wert, auf den 128 nach der Gamma-Korrektur abgebildet wird, also
% das mittlere Grau. Das gilt aber nur für CRT monitore bei einem angenommenen 
% Gamma von 2.4: 255*(0.5^(1/2.4))=191.03
% Wenn du ein LCD verwendest, machst du es erstmal genauso, auch wenn es 
% wahrscheinlich falsch ist. ;-) Nee, die LCDs sind meist so eingestellt, das 
% ein Bild liefern, wie es die Leute von ihrem CRT gewöhnt sind.
% Du zeigt also ein Grau von RGB = (190, 190, 190) und fummelst solange an
% dem Kontrast und der Helligkeit herum, bis dein Photometer 30cd/m^2 anzeigt.
% Wir hatten neulich
%     Contrast 80%
%     Brightness 60%
% Schlecht sind Werte von 100%.
% Als nächstes solltest du die Farbtemperatur des Monitors so einstellen, dass 
% das mittlere Grau RGB = (190, 190, 190) einen CIE x,y Wert von 0.33, 0.33 
% liefert. Das kannst du ja mit deiner Greytag Wunderwaffe feststellen.
% Oder du stellst einfach eine Farbtemperatur von 6400 K ein, die bei unserem 
% Iiyama die beste Annäherung an 0.33 0.33 lieferte.
% Für jeden Kanal einzeln an jeweils 22 Punkten:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information für Elo-Touch (IRIS)
%%%%%%%%%%%%%%%%%%%%%%%%%%% Setting b0c76 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 127 -> 30.52 cd/m2; 128 -> 30.52 cd/m2
% gammaValue = 2.4155
% isoluminance (in original lookup table): 
%    green: 0,142,0 (10.50 cd/m2)
%    gray:  123,123,123 (10.50 cd/m2)
%    red:   0,0,255 (10.50 cd/m2)
%    physikalisch: 10.6-10.8 cd/m2
%%%%%%%%%%%%%%%%%%%%%%%%%%% Setting b0c45 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 127 -> 8.44 cd/m2; 128 -> 8.55 cd/m2
% gammaValue = 2.6139
% isoluminance (in original lookup table): 
%   green:  0,146,0 (2.95 cd/m2))
%   gray:   130,130,130 (2.95 cd/m2)
%   red:    0,0,255 (2.95 cd/m2)

%%%%%%%%%%%%%%%%%%%%%%%%%% leily-NEC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gammaValue = 1.3833
% red-green isoluminance (in original setting): 139 (22.94 cd/m2))
% red-gray isoluminance (in original setting): 84 (22.94 cd/m2)

%%%%%%%%%%%%%%%%%%%%%%%%%% gamma-Mitsubishi_IRIS.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gammaValue = 1.9874
% red-green isoluminance (in original setting): 147 (25.26 cd/m2))
% red-gray isoluminance (in original setting): 120 (25.26 cd/m2)

%%%%%%%%%%%%%%%%%%%%%%%%%% gamma-Mitsubishi_3.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gammaValue = 1.8554
% red-green isoluminance (in original setting): 138 (26.20 cd/m2))
% red-gray isoluminance (in original setting): 113 (26.20 cd/m2)
% max: 118 cd/m2

%%%%%%%%%%%%%%%%%%%%%%%%%% gamma-SonyG520.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gammaValue = 2.6362
% red-green isoluminance (in original setting): 161 (20.72 cd/m2))
% red-gray isoluminance (in original setting): 141 (20.72 cd/m2)
%%%%%%%%%%%%%%%%%%%%%%%%%% Dell-Zebris.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gammaValue = 2.0392
% red-green isoluminance (in original setting): 153 (29.95 cd/m2))
% red-gray isoluminance (in original setting): 126 (29.95 cd/m2)
%%%%%%%%%%%%%%%%%%%%%%%%%% NEC EyeLInk1000.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gammaValue =   2.0471
% red(255)-green(x) isoluminance (in original setting): x=140 (24.07 cd/m2))
% red(255)-gray(x) isoluminance (in original setting): x=117 (24.07 cd/m2)
%%%%%%%%%%%%%%%%%%%%%%%%%% elot-touch-flat.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gammaValue = 2.3363
% red(255)-green(x) isoluminance (in original setting): x=153 (46.81 cd/m2))
% red(255)-gray(x) isoluminance (in original setting): x=131 (46.81 cd/m2)
%%%%%%%%%%%%%%%%%%%%%%%%%% elot-touch-classic.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gammaValue =    2.0108
% red(255)-green(x) isoluminance (in original setting): x=128 (19.97 cd/m2))
% red(255)-gray(x) isoluminance (in original setting): x=109 (19.97 cd/m2)
clear;
file_name   = 'SONY_Monitor_ColorCal_Readings';
data        = textread([file_name '.txt']);
volt        = 1;
red         = 2;
green       = 3;
blue        = 4;
gray        = 5;
maxIndex    = 255;  %  +1 = total num of indices
indexNormal = (0:0.001:1)';

opts = fitoptions;
for i=red:gray,
    fresult = fit(data(2:end,volt), data(2:end,i),'power1',opts)
    predict(:,i) = fresult.a*(indexNormal*max(data(:,volt))).^fresult.b;
end

% plot raw data and adjusted curve
subplot(2,2,1);
plot(maxIndex*(0:0.001:1)', predict(:,gray));
hold on;
scatter(data(:,volt), data(:,gray));
title('Raw Data and Power Function');

% create a gamma-corrected table
% Inverse function is calculated numerically
% independent of fitting curves
orgTable = [maxIndex*indexNormal  predict(:,:)];
correctedTable = [];
for i = 0: 255
    for j=red:gray,
        eqLum = i * max(data(:,j))/maxIndex;
        numInd = max(find(orgTable(:,j+1) <= eqLum));
        numTable(j) = round(orgTable(numInd,1));
    end
    correctedTable = [correctedTable ;  i numTable(red:gray)];
end
%%

% confirm linearity
subplot(2,2,2)
linearTable = [];
for i=0:255
    for j = red:gray,
        numInd = max(find(orgTable(:,1) <= correctedTable(i+1, j)));
        numTable(j) = orgTable(numInd, j+1);
    end
    linearTable = [linearTable;  i  numTable(red:gray) ];
end


%%
sumTable = linearTable(:,red) + linearTable(:,green) + linearTable(:,blue);

plot(linearTable(:,1), linearTable(:,red)');
hold on;
plot(linearTable(:,1), linearTable(:,green)');
hold on;
plot(linearTable(:,1), linearTable(:,blue)');
hold on;
plot(linearTable(:,1), linearTable(:,gray)');
hold on;
plot(linearTable(:,1), sumTable, 'r');
title('adjusted LUT-Table');
xlabel('Voltage (0-255)');
ylabel('Luminance cd/m2');

dlmwrite([file_name '.r'], correctedTable(:,red));
dlmwrite([file_name '.g'], correctedTable(:,green));
dlmwrite([file_name '.b'], correctedTable(:,blue));
dlmwrite([file_name '.gray'], correctedTable(:,gray));

%correctedTable
%gammaValue = fresult.b

% plot rgb and gray
subplot(2,2,3);
plot(data(:,volt),data(:,red),'r*');
xlabel('Voltage (0-255)');
ylabel('Luminance cd/m2');
title('red, green, blue, gray');
hold on;
plot(data(:,volt),data(:,green),'g*');
plot(data(:,volt),data(:,blue),'b*');
plot(data(:,volt),data(:,gray),'k*');

% confirm additivity
subplot(2,2,4);
plot(data(:,volt),data(:,red)+data(:,green)+data(:,blue),'k');
hold on
scatter(data(:,volt),data(:,gray),'k+');
title('Gray vs Sum of rgb');
xlabel('Voltage (0-255)');
ylabel('Luminance cd/m2');

%%
cyan = 6;
yellow = 7;
magenta = 8;

predict(:,cyan) = predict(:,green) + predict(:,blue);
predict(:,yellow) = predict(:,red) + predict(:,green);
predict(:,magenta) = predict(:,red) + predict(:,blue);


% equiluminant to maximal blue
eqLum = [];
for i = red:magenta,
%     if i == blue,
%         continue;
%     end
    numInd = max(find(predict(:,i) <= data(end,blue)));
    eqLum(i) = round(orgTable(numInd,1));
    %colorVec = {'dummy', 'For red:', 'For green:', 'For blue:', 'For gray:','For cyan:', 'For yellow:', 'For magenta:'};
    %curCol = colorVec(i);
    %sprintf(curCol);
    fprintf('blue(255)-color(x) isoluminance (in original setting): x=%d (%.2f cd/m2)\n', eqLum(i), data(end,blue));
end

dlmwrite('eqLum.log', eqLum);

%%

% output gamma-table
fid = fopen(['../' file_name],'w');
fprintf(fid,'%d %d %d %d %d\n',correctedTable');
fclose(fid);
 