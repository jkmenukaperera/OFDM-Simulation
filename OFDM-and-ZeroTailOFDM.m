clear all
close all

clc
nbits = 208000;
modlevel = 2 ;
nbitpersym  = 52;   % number of bits per qam OFDM symbol (same as the number of subcarriers for 16-qam)
nsym        = 10^4; % number of symbols
len_fft     = 64;   % fft size
sub_car     = 52;   % number of data subcarriers
EbNo        = 0:2:15;

EsNo= EbNo+10*log10(52/64)+ 10*log10(64/80) +10*log10(4);%code rate +___+modulation rate%

snr=EsNo - 10*log10((64/80));


M = modem.qammod('M',16); % modulation object

% Generating data

t_data=randint(nbitpersym*nsym*4,1,2);


qamdata=bi2de(reshape(t_data,4,520000).','left-msb');

maping = bin2gray(qamdata,'qam',16);

% modulating data

mod_data =1/sqrt(10)* modulate(M,maping);

% serial to parallel conversion

par_data = reshape(mod_data,nbitpersym,nsym).';


% pilot insertion

pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;

% fourier transform time doamain data

IFFT_data =ifft(fftshift(pilot_ins_data.')).';
a=max(max(abs(IFFT_data)));
IFFT_data=IFFT_data./a; % normalization

% addition cyclic prefix

 cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';
 cylic_add_data_1 = [zeros(10000,16) IFFT_data].';
%   size(cylic_add_data)
%   size(cylic_add_data_1)
% parallel to serial coversion

ser_data = reshape(cylic_add_data,80*nsym,1);
% 
zero1 = reshape(cylic_add_data_1,80*nsym,1);


% passing thru channel

no_of_error=[];
ratio=[];

for ii=1:length(snr)
  
chan_awgn = awgn(ser_data,snr(ii),'measured'); % awgn addition

ser_to_para = reshape(chan_awgn,80,nsym).'; % serial to parallel coversion

cyclic_pre_rem = ser_to_para(:,[17:80]);   %cyclic prefix removal

FFT_recdata =a*fftshift(fft(cyclic_pre_rem.')).';    % freq domain transform

rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); %pilot removal

ser_data_1 =sqrt(10)* reshape(rem_pilot.',nbitpersym*nsym,1);  % serial coversion

z=modem.qamdemod('M',16);

demod_Data = demodulate(z,ser_data_1);  %demodulatin the data

demaping = gray2bin(demod_Data,'qam',16);

data1 = de2bi(demaping,'left-msb');

data2 = reshape(data1.',nbitpersym*nsym*4,1);

[no_of_error(ii),ratio(ii)]=biterr(t_data , data2) ; % error rate calculation

end


for iii=1:length(snr)
  
chan_awgn_1 = awgn(zero1,snr(iii),'measured'); % awgn addition

ser_to_para_1 = reshape(chan_awgn_1,80,nsym).'; % serial to parallel coversion

cyclic_pre_rem_1 = ser_to_para_1(:,[17:80]);   %cyclic prefix removal

FFT_recdata_1 =a*fftshift(fft(cyclic_pre_rem_1.')).';    % freq domain transform

rem_pilot_1 = FFT_recdata_1 (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); %pilot removal

ser_data_1_1 =sqrt(10)* reshape(rem_pilot_1.',nbitpersym*nsym,1);  % serial coversion

z=modem.qamdemod('M',16);

demod_Data_1 = demodulate(z,ser_data_1_1);

demaping_1 = gray2bin(demod_Data_1,'qam',16);

data1_1 = de2bi(demaping_1,'left-msb');

data2_1 = reshape(data1_1.',nbitpersym*nsym*4,1);

[no_of_error_1(iii),ratio_1(iii)]=biterr(t_data , data2_1) ; % error rate calculation

end

% plotting the result
semilogy(EbNo,ratio,'--*r','linewidth',2);
hold on;
semilogy(EbNo,ratio_1,'--*g','linewidth',2);
theoryBer = (1/4)*3/2*erfc(sqrt(4*0.1*(10.^(EbNo/10))));
semilogy(EbNo,theoryBer ,'--b','linewidth',2);
axis([0 15 10^-5 1])
legend('Simulated OFDM','Simulated ZERO-TAIL-OFDM','Theoritical')
grid on
xlabel('EbNo');
ylabel('BER')
title('Bit error probability curve for qam');
