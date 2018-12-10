%% -----------------------------------------------------------------------
%
% Title       : test_average.m
% Author      : Alexander Kapitanov	
% Company     : AO Insys
% E-mail      : sallador@bk.ru 
% Version     : 1.0	 
%
%-------------------------------------------------------------------------
%
% Description : 
%    Moving average filters
%
%-------------------------------------------------------------------------
%
% Version     : 1.0 
% Date        : 2018.12.09
%
%-------------------------------------------------------------------------

% Preparing to work
close all;
clear all;

set(0, 'DefaultAxesFontSize', 14, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 14, 'DefaultTextFontName', 'Times New Roman'); 

% Settings
Num = 2^8;      % Number of points
SNR = 11-2;      % Signal noise ratio

F1 = 5;         % Signal freq1
F2 = 4;         % Signal freq2


% Create input data vector
for i = 1:Num
  Dat(i,1) = (cos(F1 * (i-1) * 2 * pi / Num));
  
  Dat(i,1) = cos(F1 * (i-1) * 2 * pi / Num) + 4*i/Num;

  %  Dat(i,1) = (cos(F1 * (i-1) * 2 * pi / Num)) + ((1/4) * cos(F2 * (i-1) * 2 * pi / Num));
%  Dat(i,1) = i-1;

%  if (i < Num/2)
%    Dat(i,1) = cos(F1 * (i-1) * 2 * pi / Num) + 4*i/Num;
%  else
%    Dat(i,1) = cos(F1 * (i-1) * 2 * pi / Num) + (4*Num/2 - 4*i/Num);
%  endif

%  if (i < Num/2)
%    Dat(i,1) = 4*i/Num;
%  else
%    Dat(i,1) = 4 - 4*i/Num;
%  endif

endfor

% Adding noise to real signal 
DatAwgn = awgn(Dat, SNR, 0, 1);     

% Add peaks
DatAwgn(1*Num/4,1) = 8;
DatAwgn(2*Num/4,1) = 6;
DatAwgn(3*Num/4,1) = 8;


% Execute moving average filters
function fil_out = fn_filter(Data, Ndat, Nord)
  fil_dat = zeros(Ndat, 1);
  for i = 1:Ndat
    for j = 0:Nord-1
      if ((i-1) < j)
        fil_dat(i,1) = fil_dat(i,1);
      else 
        fil_dat(i,1) = fil_dat(i,1) + Data(i-j,1);
      endif
    endfor
  endfor
  fil_out = fil_dat / Nord;
endfunction

% Calculate some MAF filters
Fl2 = fn_filter(DatAwgn, Num, 4);
Fl3 = fn_filter(DatAwgn, Num, 8);
Fl4 = fn_filter(DatAwgn, Num, Num/16);
Fl5 = fn_filter(DatAwgn, Num, Num/8);

FlN = fn_filter(DatAwgn, Num, 16);

% Find sorting data from min to max
DatSort = sort(DatAwgn);

MedCalc = 2 * median(DatAwgn) / Num;
AvgCalc = 2 * mean(DatAwgn) / Num;


for i = 1:Num
  DatAvg(i,1) = (i-1) * AvgCalc;
endfor


% Plot average data
figure('Name', 'Moving average filter', 'NumberTitle','off');
plot(DatAwgn, '-', 'LineWidth', 1, 'Color',[0.8 0 0]); hold on;
plot(FlN, '-*', 'LineWidth', 1, 'Color', [0 0 1]); hold on; 
plot(DatAvg, '-', 'LineWidth', 2, 'Color', [0.4 0.4 0.8]); hold on; 
grid on; axis tight;
legend({'Input Signal', 'Filter N = 16', 'Average Line'}, 'Location', 'northwest');
title(['Moving average filter']) 


% Plot some filters
figure('Name', 'Moving average filters N order', 'NumberTitle','off');
subplot(2,2,1)
plot(DatAwgn, '-', 'LineWidth', 1, 'Color',[1 0 0]); hold on;
plot(DatAvg, '-', 'LineWidth', 2, 'Color', [0.4 0.4 0.8]); hold on; 
plot(Fl2, '*', 'LineWidth', 1, 'Color', [0 0 1]); grid on; axis tight;
title(['Order = 4'])

subplot(2,2,2)
plot(DatAwgn, '-', 'LineWidth', 1, 'Color',[1 0 0]); hold on;
plot(DatAvg, '-', 'LineWidth', 2, 'Color', [0.4 0.4 0.8]); hold on; 
plot(Fl3, '*', 'LineWidth', 1, 'Color', [0 0.2 1]); grid on; axis tight;
title(['Order = 8'])

subplot(2,2,3)
plot(DatAwgn, '-', 'LineWidth', 1, 'Color',[1 0 0]); hold on;
plot(DatAvg, '-', 'LineWidth', 2, 'Color', [0.4 0.4 0.8]); hold on; 
plot(Fl4, '*', 'LineWidth', 1, 'Color', [0 0.4 1]); grid on; axis tight;
title(['Order = N / 8'])

subplot(2,2,4)
plot(DatAwgn, '-', 'LineWidth', 1, 'Color',[1 0 0]); hold on;
plot(DatAvg, '-', 'LineWidth', 2, 'Color', [0.4 0.4 0.8]); hold on; 
plot(Fl5, '*', 'LineWidth', 1, 'Color', [0 0.6 1]); grid on; axis tight;
title(['Order = N / 4'])
