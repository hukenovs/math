%% -----------------------------------------------------------------------
%
% Title       : test_parallel_fir.m
% Author      : Alexander Kapitanov	
% Company     : AO "Insys"
% E-mail      : sallador@bk.ru 
% Version     : 1.0	 
%
% ------------------------------------------------------------------------
%
% Description : 
%    Create Parallel Finite Impulse Response (FIR) filters for High-Freq signals. 
%    It is like polyphase FIR but is without decimation after last adder!
%
%    For example: IF the input freq > DSP freq you can use interleave mode on the inputs.
%    So you have M = round( Freq[Input] / Freq[DSP] ) points of input signal,
%    You can implement parallel scheme as is. But this method uses AREA resources
%    Improve IR and Input Signal by simple ADD/SUB operation you can optimize it!
%
%    See "FFA FIR filters" for more information.
%      * FFA - fast filter architecture.
%
% Example: L = 2 and NFIR = 128 taps:
%
%  > Common method:
%    Y0 = H0 * X0 + (Z^-2) * H1 * X1
%    Y1 = H0 * X1 + H1 * X0
%    
%    Resources: 4*128 + 2 = 514 DSPs.
%
%  > FFA method (Less Area resources):
%    Y0 = H0 * X0 + (Z^-2) * H1 * X1
%    Y1 = (H0+H1)*(X0+X1) - H0*X0 - H1*X1
%
%    Resources: 3*128 + 5 = 389 DSPs.
%
%  where:
%    > Y0 and Y1 - two parts of output signal (Interleave-2 mode)
%    > X0 and X1 - two parts of input signal  (Even / Odd)
%    > H0 and H1 - two parts of FIR filters coefficients
%  
%  As you can see 1st method takes 4 FIRs and 2 ADD operations. Total DSPs = 514.
%  But the 2nd method takes only 3 FIRs and 5 ADD operations.  Total DSPs = 389.
%  
%  You can save ~24.3% DSPs resources of FPGA for L = 2 and NFIR = 128 taps!!
% 
%  | NFIR | Saved DSPs % |
%  |  64  | 23.643       |
%  | 128  | 24.319       |
%  | 256  | 24.659       |
%  | 512  | 24.829       |
%
% ------------------------------------------------------------------------
%
% Version     : 1.0 
% Date        : 2017.11.29 
%
%% ----------------------------------------------------------------------- 
%
%	GNU GENERAL PUBLIC LICENSE
% Version 3, 29 June 2007
%
%	Copyright (c) 2018 Kapitanov Alexander
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
% APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT 
% HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY 
% OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, 
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM 
% IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF 
%  ALL NECESSARY SERVICING, REPAIR OR CORRECTION. 
%
%% -----------------------------------------------------------------------		   

% Preparing to work
close all;
clear all;

set(0, 'DefaultAxesFontSize', 14, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 14, 'DefaultTextFontName', 'Times New Roman'); 

% Settings
NFFT = 2^9;            % Sampling Frequency
t = 0:1/NFFT:1-1/NFFT;  % Time vector #1
tt = 1:NFFT;            % Time vector #2

Asig = 2^14-1;
Fsig = 2;
Fm = 1*(NFFT/4);
B = Fm / NFFT;
Ffm = 1;

% For testing FORWARD and INVERSE FFT: FWT
STAGE = log2(NFFT);

N = 128; % Filter order
COE_WIDTH = 16; % Real width for FIR coeffs (Implementation)

L = 2;
%% -------------------------------------------------------------------------- %%
% ---------------- 0: CREATE INPUT DATA FOR CPP/RTL -------------------------- % 
%% -------------------------------------------------------------------------- %%

Rep = 200; % Number of clones for input signal

for i = 1:NFFT
    
  Dat(i,1) = Asig * cos((Fsig*i + B*i*i/2) * 2*pi/NFFT) * sin(i * Ffm * pi / NFFT);
    
  if (mod(i+Rep-1, Rep+1) == Rep)
    Dat(i,1) = Asig;
  else
    Dat(i,1) = 0;
  end   
    
end

% Adding noise to real signal 
%SNR = -50;
SNR = -35;

DatX = awgn(Dat, SNR, 0, 1);        

%Mre = max(abs(DatRe));
%Mim = max(abs(DatIm));
%Mdt = max(Mre, Mim);

%Din(:,1) = round(((2^15 - 1)/Mdt)*DatRe);
%Din(:,2) = round(((2^15 - 1)/Mdt)*DatIm);
%Dcm(:,1) = Din(:,1) + j*Din(:,2);

%% -------------------------------------------------------------------------- %%
% ---------------- 1:  INPUT DATA (TEST MATH SIGNAL) ------------------------- % 
%% -------------------------------------------------------------------------- %%

figure(1) % Plot loaded data in Time Domain
  subplot(2,1,1);
  plot(tt(1:NFFT), DatX, '-', 'LineWidth', 2, 'Color',[1 0  0]);
  grid on;
  hold on;
  axis tight;      
  title(['Test Data in Time Domain']);   

%% -------------------------------------------------------------------------- %%
% ---------------- 2:  FIR FILTER IMPLEMENTATION ----------------------------- % 
%% -------------------------------------------------------------------------- %%

Fcut = 20; % First freq (passband)
Fs = 100; % Sampling freq

BETA = 3; % Beta (Kaiser)
WIND = kaiser(N, BETA); % KAISER WINDOW IS USED!

t = 1:N;
f =  [Fcut]/(Fs/2);

% Filter type: 'low', 'high', 'stop', 'pass';
Hc = fir1(round(N)-1, f, 'low', WIND);
MHc = max(abs(Hc));
Hc = Hc / MHc;
Hc_n = (Hc/max(Hc))*(2^(COE_WIDTH-1)-1);
Hc_r = round(Hc_n);
Hf = 20 * log10(abs(fftshift(fft(Hc_n, 10000))));
Sp_err = 20 * log10(abs(fftshift(fft(Hc_r, 10000))));

ff = -0.5:1/10000:0.5-1/10000;

% Plot FIR in T/F Domains

%% -------------------------------------------------------------------------- %%
% ---------------- 3:  FILTER INPUT SIGNAL L = 3 ----------------------------- % 
%% -------------------------------------------------------------------------- %%

if (L == 2) 
  % 1) Model: 
  Yc_s = filter(Hc, 1, DatX);
  % 2) L=2 (4FIR implemented) 
  Xc = DatX;
  
  % Impulse repsonse:
  Hc0 = Hc(1:2:end); % each 2i
  Hc1 = Hc(2:2:end); % each 2i+1, where i=0:N/2-1

  % Input data:
  Xc0 = Xc(1:2:end);
  Xc1 = Xc(2:2:end);

  Ac0 = filter(Hc0, 1, Xc0);
  Ac1 = filter(Hc1, 1, Xc1);
  
  AAc1(1,1) = 0;
  AAc1(2:NFFT/2,1) = Ac1(1:end-1,1);

  Bc0 = filter(Hc0, 1, Xc1);
  Bc1 = filter(Hc1, 1, Xc0);

  Yc0 = Ac0 + AAc1;
  Yc1 = Bc0 + Bc1;

  for i = 1:NFFT/2
    for j = 1:2
      if (j == 1)  
        Yc(2*i-1,1) = Yc0(i,1);
      else
        Yc(2*i-0,1) = Yc1(i,1); 
      end
    end   
  end

  % L = 2 (improved)
  Xc01 = filter((Hc0+Hc1), 1, (Xc0+Xc1));
  Zc1 = Xc01 - Ac0 - Ac1;

  for i = 1:NFFT/2
    for j = 1:2
      if (j == 1)  
        Zc(2*i-1,1) = Yc0(i,1);
      else
        Zc(2*i-0,1) = Zc1(i,1); 
      end
    end   
  end

  % L = 2 (Proposed: symmetric)
  Qc01 = filter((Hc0(1:length(Hc0)/2)+Hc1(1:length(Hc1)/2)), 1, (Xc0+Xc1));
  Qc10 = filter((Hc0(1:length(Hc0)/2)-Hc1(1:length(Hc1)/2)), 1, (Xc0-Xc1));
  Qc11 = Ac1 + AAc1;
  
  Qc0 = Qc01 + Qc10 - Qc11;
  Qc1 = Qc01 - Qc10;

  for i = 1:NFFT/2
    for j = 1:2
      if (j == 1)  
        Qc(2*i-1,1) = Qc0(i,1);
      else
        Qc(2*i-0,1) = Qc1(i,1); 
      end
    end   
  end  

  Yc_diff = Yc-Yc_s;
  Zc_diff = Zc-Yc_s;
  Qc_diff = Qc-Yc_s;
  
end

%% -------------------------------------------------------------------------- %%
% ---------------- 3:  FILTER INPUT SIGNAL L = 3 ----------------------------- % 
%% -------------------------------------------------------------------------- %%
if (L == 3)
  % 1) Model: 
  Yc_s = filter(Hc, 1, DatX(1:length(DatX)-rem(length(DatX), 3)));
  % 2) L=2 (4FIR implemented) 
  Xc = DatX(1:length(DatX)-rem(length(DatX), 3));
  % Impulse repsonse:
  Hc0 = Hc(1:3:end); % each 3i
  Hc1 = Hc(2:3:end); % each 3i+1,
  Hc2 = Hc(3:3:end); % each 3i+2,

  % Input data:
  Xc0 = Xc(1:3:end);
  Xc1 = Xc(2:3:end);
  Xc2 = Xc(3:3:end);

  Ac0 = filter(Hc0, 1, Xc0);
  Ac1 = filter(Hc0, 1, Xc1);
  Ac2 = filter(Hc0, 1, Xc2);

  Bc0 = filter(Hc1, 1, Xc0);
  Bc1 = filter(Hc1, 1, Xc1);
  Bc2 = filter(Hc1, 1, Xc2);  
  
  Cc0 = filter(Hc2, 1, Xc0);
  Cc1 = filter(Hc2, 1, Xc1);
  Cc2 = filter(Hc2, 1, Xc2);
  
%  AAc1(1,1) = 0;
%  AAc1(2:round(NFFT/3),1) = Ac1(1:round(NFFT/3)-1,1);

  Yc0 = Ac0 + Bc2 + Cc1;
  Yc1 = Ac1 + Bc0 + Cc2;
  Yc2 = Ac2 + Bc1 + Cc0;

  for i = 1:floor(NFFT/3)
    for j = 1:3
      if (j == 1)  
        Yc(3*i-2,1) = Yc0(i,1);
      elseif (j == 2)
        Yc(3*i-1,1) = Yc1(i,1);
      else
        Yc(3*i-0,1) = Yc2(i,1);        
      end
    end   
  end

  % L = 3 (improved)
  
  Xc01  = filter((Hc0+Hc1), 1, (Xc0+Xc1));
  Xc12  = filter((Hc1+Hc2), 1, (Xc1+Xc2));
  Xc012 = filter((Hc0+Hc1+Hc2), 1, (Xc0+Xc1+Xc2));
  
  Tco = Ac0 - Cc2;
  Tc1 = Xc12 - Bc1;
  Tc2 = Xc01 - Bc1;
  
  Zc0 = Tco + Tc1;
  Zc1 = Tc2 - Tco;
  Zc2 = Xc012 - Tc2 - Tc1;

  for i = 1:floor(NFFT/3)
    for j = 1:3
      if (j == 1)  
        Zc(3*i-2,1) = Zc0(i,1);
      elseif (j == 2)
        Zc(3*i-1,1) = Zc1(i,1);
      else
        Zc(3*i-0,1) = Zc2(i,1);        
      end
    end   
  end
  
  Yc_diff = Yc - Yc_s(1:length(Zc));
  Zc_diff = Zc - Yc_s(1:length(Zc));  
 
end


figure(1) % Plot loaded data in Time Domain  
  subplot(3,1,1);
  plot(Yc_s, '-', 'LineWidth', 2, 'Color',[0 0  1])
  grid on
  hold on
  axis tight      
  title(['Test Data in Time Domain'])

figure(1)
  subplot(3,1,2);
  plot(Hc, '*-', 'LineWidth', 1, 'Color',[1 0 0]);
  axis tight;
  title(['Filter IR, Order = ',num2str(N)]);
  xlabel ('Time');
  ylabel ('Magnitude');
  grid on;

  subplot(3,1,3);
  plot(Hf-max(Hf), '-', 'LineWidth', 1, 'Color',[1 0 0]);
  hold on
  
  plot(ff, Sp_err-max(Sp_err), '-', 'LineWidth', 1, 'Color',[0 0 1]);
  grid on;
  axis([0 ff(10000) -120 0]) 
  title('FIR Spectrum');
  xlabel ('Freq ( x rad / samples)');
  ylabel ('Magnitude (dB)');  
  legend ('Double dtype', ['CoeWidth: ',num2str(COE_WIDTH)], 'location', 'east'); 


figure(2) % Plot loaded data in Time Domain  
  subplot(5,1,1);
  plot(Yc_s, '-', 'LineWidth', 1, 'Color',[1 0 0])
  grid on
  hold on
  axis tight      
  title(['Yc (algorithm as is)'])   
  
  subplot(5,1,2);
  plot(Yc, '-', 'LineWidth', 1, 'Color',[0 0 1])
  grid on
  hold on
  axis tight      
  title(['Yc (parallel)'])  
  
  subplot(5,1,3);
  plot(Yc_diff, '-', 'LineWidth', 2, 'Color',[0 1 0])
  grid on
  hold on
  axis tight      
  title(['Yc difference'])  
  
  subplot(5,1,4);
  plot(Zc_diff, '-', 'LineWidth', 2, 'Color',[0 1 0])
  grid on
  hold on
  axis tight      
  title(['Yc difference (II method)'])    
  
  subplot(5,1,5);
  plot(Qc_diff, '-', 'LineWidth', 2, 'Color',[0 1 0])
  grid on
  hold on
  axis tight      
  title(['Yc difference (Proposed)'])    
  
  
%disp("Y res   : ="), disp(Yc_s(145:149,1));
%disp("Yc1     : ="), disp(Yc1(71:75,1));
%disp("Ac0     : ="), disp(Ac0(71:75,1));
%disp("Ac1     : ="), disp(Ac1(71:75,1));
%disp("Ac0+Ac1 : ="), disp(Ac0(73,1)+Ac1(72,1));
%disp("Ac0+Ac1 : ="), disp(Ac0(74,1)+Ac1(73,1));
%disp("Ac0+Ac1 : ="), disp(Ac0(75,1)+Ac1(74,1));

clear Asig;
clear B;
clear BETA;
clear Hf;
clear Hc_r;
clear MHc;
clear SNR;
clear Sp_err;
clear Hc_n;
clear f;
clear i;
clear j;
clear N;
clear Dat;
clear Fcut;
clear Ffm;
clear Fm;
clear Fs;
clear Fsig;
clear WIND;
clear t;
clear ff;
clear STAGE;
clear SP_err;
clear Rep;
clear COE_WIDTH;