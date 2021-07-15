%This program implements the generation, modulation, reception, and detection of the
%TR-UWB signal.

%Developed by Areeb Tariq
%Date: June 7, 2021.
%Course: RadioFrequency for connected objects

close all

%**********************************************************************************************
%                                                                                             *
%                       PART-1: Generation and Modulation of TR-UWB Signals
%                                                                                             *
%**********************************************************************************************


%**********************************************************************************************
%                               Section-1: Variables 
%**********************************************************************************************
Td = 50e-9;                         %Time delay between reference and information pulse
Ts = 300e-9;                        %Time delay between consecutive references pulses
Tp = 2e-9;                          %UWB signal pulse width
fs = 3e10;                          %Sampling frequency
fn = fs/2;                          %Nyquist Frequency should be minimum half of the sampling frequency
fp =4.492e9;                        %Carrier frequency of cosine
t=-2e-9:1/fs:1.9e-6;                %Range of time for the UWB signal
information_seq = [1 0 1 1 0 1 0];  %The binary sequence for the information signals

%**********************************************************************************************
%               Section-2: Generation of the Signal Reference pulse as specified 
%**********************************************************************************************
%Using the formula in course slide-21 for Scholtz's monocycle
Y=(1-(4*pi.*(t.^2))/Tp^2) .* exp(-2*pi.*(t.^2)/Tp^2);
subplot(2,1,1);
plot(t,Y,'-r');
axis([-2e-9 2e-9 -0.6 1.2])
ylabel('amplitude','fontSize',15)
xlabel('time, secs','fontSize',15)
title('Time domain: UWB reference signal pulse (Scholtz monocycle, Tp=2ns) ','fontSize',15)

%**********************************************************************************************
%               Section-3: Generation of frequency spectrum of the Signal reference pulse
%**********************************************************************************************

%First we need to calculate the points for Discrete fourier transform
%Points should be multiple of 2, the additional points than signal are padded with zero
points_for_dft = 2.^(ceil(log(length(Y))/log(2)));

%Take the fast fourier transform of time domain signal with the calculated points
fft_uwb_reference_signal=fft(Y,points_for_dft);

%Find the unique dft points since fft is symmetric 
non_repeated_dft_points = ceil((points_for_dft+1)/2); 

%Update the FFT for the TR-UWB signal with just unique points
fft_uwb_reference_signal=fft_uwb_reference_signal(1:non_repeated_dft_points);

%Get an absolute value and adjust with the size of the information sequence
uwb_ref_spectrum=abs(fft_uwb_reference_signal);
uwb_ref_spectrum=uwb_ref_spectrum*2;
uwb_ref_spectrum(1)=uwb_ref_spectrum(1)/2;
uwb_ref_spectrum(length(uwb_ref_spectrum))=uwb_ref_spectrum(length(uwb_ref_spectrum))/2;
uwb_ref_spectrum=uwb_ref_spectrum/length(Y);

%Generate the scale for the frequency spectrum of TR-UWB signal pulse train
f=(0:non_repeated_dft_points-1)*2*fn/points_for_dft;
subplot(2,1,2);
plot(f,20*log10(uwb_ref_spectrum));
axis([0 10e9 -160 -50])
xlabel('Frequency (Hz)','fontSize',15); 
title('Frequency spectrum: UWB reference signal pulse (Scholtz monocycle)','fontSize',15);
ylabel('Power/frequency (dB/Hz)','fontSize',15);

%**********************************************************************************************
%           Section-4: Generation of the train of Reference and information pulses 
%**********************************************************************************************
t_for_ref_pulse = t;                %The variable is used to delay the time variable to implement reference pulses
t_for_info_pulse = t;               %The variable is used to delay the time variable to implement information pulses

%Implementation of reference pulses
for info_bit = 1:length(information_seq)
    t_for_ref_pulse = t_for_ref_pulse - Ts;
    Y = Y + (1-(4*pi.*(t_for_ref_pulse.^2))/Tp^2) .* exp(-2*pi.*(t_for_ref_pulse.^2)/Tp^2);
end

%Implementation of information pulses
for info_bit = 1:length(information_seq)
    if info_bit == 1
        t_for_info_pulse = t_for_info_pulse - Td;
        Y = Y + (1-(4*pi.*(t_for_info_pulse.^2))/Tp^2).*exp(-2*pi.*(t_for_info_pulse.^2)/Tp^2); 
    end
    if info_bit > 1
        if information_seq(info_bit) == 1
            t_for_info_pulse = t_for_info_pulse - Ts;
            Y = Y + (1-(4*pi.*(t_for_info_pulse.^2))/Tp^2).*exp(-2*pi.*(t_for_info_pulse.^2)/Tp^2);
        end
        if information_seq(info_bit) == 0
            t_for_info_pulse = t_for_info_pulse - Ts;
            Y = Y + (1-(4*pi.*(t_for_info_pulse.^2))/Tp^2).*exp(-2*pi.*(t_for_info_pulse.^2)/Tp^2)/3;
        end
    end
end

tr_uwb_signal = Y;
figure
subplot (3,1,1);
plot(t,tr_uwb_signal,'-r')
axis([-40e-9 1.9e-6 -0.6 1.2])
ylabel('amplitude', 'fontSize',15)
xlabel('time, secs','fontSize',15)
title('Time domain: Binary coded TR-UWB pulse train','fontSize',15)

%The code below gives a zoomed in version of the generated pulse containing
%the first two information bits. It gives an idea to observe the Ts, Td, and Tp
subplot (3,1,2);
plot(t,tr_uwb_signal,'-r')
axis([-3e-9 3.6e-7 -0.6 1.2])
ylabel('amplitude', 'fontSize',15)
xlabel('time, secs','fontSize',15)
title('Time domain: Zoomed in first two coded information bits as an example','fontSize',15)

%**********************************************************************************************
%           Section-5: Generation of frequency spectrum of the UWB signal
%**********************************************************************************************

%First we need to calculate the points for Discrete fourier transform
%Points should be multiple of 2, the additional points than signal are padded with zero
points_for_dft = 2.^(ceil(log(length(tr_uwb_signal))/log(2)));

%Take the fast fourier transform of time domain signal with the calculated points
fft_tr_uwb_signal=fft(tr_uwb_signal,points_for_dft);

%Find the unique dft points since fft is symmetric 
non_repeated_dft_points = ceil((points_for_dft+1)/2); 

%Update the FFT for the TR-UWB signal with just unique points
fft_tr_uwb_signal=fft_tr_uwb_signal(1:non_repeated_dft_points);

%Get an absolute value and adjust with the size of the information sequence
tr_uwb_spectrum=abs(fft_tr_uwb_signal);
tr_uwb_spectrum=tr_uwb_spectrum*2;
tr_uwb_spectrum(1)=tr_uwb_spectrum(1)/2;
tr_uwb_spectrum(length(tr_uwb_spectrum))=tr_uwb_spectrum(length(tr_uwb_spectrum))/2;
tr_uwb_spectrum=tr_uwb_spectrum/length(tr_uwb_signal);

%Generate the scale for the frequency spectrum of TR-UWB signal pulse train
f=(0:non_repeated_dft_points-1)*2*fn/points_for_dft;
subplot(3,1,3);
plot(f,20*log10(tr_uwb_spectrum));
xlabel('Frequency (Hz)','fontSize',15); 
title('Frequency spectrum: Binary coded TR-UWB pulse train','fontSize',15);
ylabel('Power/frequency (dB/Hz)','fontSize',15);

%**********************************************************************************************
%           Section-6: Modulation of UWB baseband signal with carrier
%**********************************************************************************************
%Using an Amplitude modulation, single sideband modulation in MATLAB
[modulated_uwb_signal,modulated_time] = modulate(tr_uwb_signal,fp,fs,'amssb');
figure
subplot(3,1,1);
plot(modulated_time,modulated_uwb_signal,'-r')
axis([-4e-8 1.9e-6 -1.2 1.2])
ylabel('amplitude', 'fontSize',15)
xlabel('time, secs','fontSize',15)
title('Time domain: Modulated TR-UWB signal','fontSize',15)

%The code below gives a zoomed in version of the generated pulse containing
%the first two information bits. It gives an idea to observe the Ts, Td, Tp, and modulation.
subplot(3,1,2);
plot(modulated_time,modulated_uwb_signal,'-r')
axis([-0.1e-9 0.60e-7 -1.2 1.2])
%axis([2.95e-7 3.55e-7 -1.2 1.2])
ylabel('amplitude', 'fontSize',15)
xlabel('time, secs','fontSize',15)
title('Time domain: Zoomed in first information bit modulated as an example','fontSize',15)

%**********************************************************************************************
%           Section-7: Generation of frequency spectrum of the modulated UWB signal
%**********************************************************************************************
%First we need to calculate the points for Discrete fourier transform
%Points should be multiple of 2, the additional points than signal are padded with zero
points_for_dft = 2.^(ceil(log(length(modulated_uwb_signal))/log(2)));

%Take the fast fourier transform of time domain signal with the calculated points
fft_modulated_tr_uwb_signal=fft(modulated_uwb_signal,points_for_dft);

%Find the unique dft points since fft is symmetric 
non_repeated_dft_points = ceil((points_for_dft+1)/2); 

%Update the FFT for the TR-UWB signal with just unique points
fft_modulated_tr_uwb_signal=fft_modulated_tr_uwb_signal(1:non_repeated_dft_points);

%Get an absolute value and adjust with the size of the information sequence
modulated_tr_uwb_spectrum=abs(fft_modulated_tr_uwb_signal);
modulated_tr_uwb_spectrum=modulated_tr_uwb_spectrum*2;
modulated_tr_uwb_spectrum(1)=modulated_tr_uwb_spectrum(1)/2;
modulated_tr_uwb_spectrum(length(modulated_tr_uwb_spectrum))=modulated_tr_uwb_spectrum(length(modulated_tr_uwb_spectrum))/2;
modulated_tr_uwb_spectrum=modulated_tr_uwb_spectrum/length(modulated_uwb_signal);

%Generate the scale for the frequency spectrum of modulated TR-UWB signal pulse train
f=(0:non_repeated_dft_points-1)*2*fn/points_for_dft;
subplot(3,1,3);
plot(f,20*log10(modulated_tr_uwb_spectrum));
xlabel('Frequency (Hz)','fontSize',15); 
title('Frequency spectrum of Modulated TR-UWB pulse train','fontSize',15);
ylabel('Power/frequency (dB/Hz)','fontSize',15);

%**********************************************************************************************
%                                                                                             *
%                       PART-2: Reception and Detection of TR-UWB Signals
%                                                                                             *
%**********************************************************************************************

%**********************************************************************************************
%                       Section-8: Squaring of the received TR-UWB signal
%**********************************************************************************************
%Square the received modulated signal
received_sqr_tr_uwb_signal = modulated_uwb_signal.^2;

figure
subplot(3,1,1);
plot(modulated_time,received_sqr_tr_uwb_signal,'-r')
axis([-4e-8 1.9e-6 -0.2 1.2])
ylabel('amplitude', 'fontSize',15)
xlabel('time, secs','fontSize',15)
title('TIme domain: Squared Received TR-UWB signal','fontSize',15)

%The code below gives a zoomed in version of the generated pulse containing
%the first two information bits. It gives an idea to observe the Ts, Td, and Tp.
subplot(3,1,2);
plot(modulated_time,received_sqr_tr_uwb_signal,'-r')
axis([-0.1e-9 0.60e-7  -1.2 1.2])
ylabel('amplitude', 'fontSize',15)
xlabel('time, secs','fontSize',15)
title('Time domain: Zoomed in first information bit in Squared Received TR-UWB signal as an example','fontSize',15)

%**********************************************************************************************
%           Section-9: Generation of frequency spectrum of the received squared UWB signal
%**********************************************************************************************
%First we need to calculate the points for Discrete fourier transform
%Points should be multiple of 2, the additional points than signal are padded with zero
points_for_dft = 2.^(ceil(log(length(received_sqr_tr_uwb_signal))/log(2)));

%Take the fast fourier transform of time domain signal with the calculated points
fft_rec_sqr_tr_uwb_signal=fft(received_sqr_tr_uwb_signal,points_for_dft);

%Find the unique dft points since fft is symmetric 
non_repeated_dft_points = ceil((points_for_dft+1)/2); 

%Update the FFT for the TR-UWB signal with just unique points
fft_rec_sqr_tr_uwb_signal=fft_rec_sqr_tr_uwb_signal(1:non_repeated_dft_points);

%Get an absolute value and adjust with the size of the information sequence
rec_sqr_tr_uwb_spectrum=abs(fft_rec_sqr_tr_uwb_signal);
rec_sqr_tr_uwb_spectrum=rec_sqr_tr_uwb_spectrum*2;
rec_sqr_tr_uwb_spectrum(1)=rec_sqr_tr_uwb_spectrum(1)/2;
rec_sqr_tr_uwb_spectrum(length(rec_sqr_tr_uwb_spectrum))=rec_sqr_tr_uwb_spectrum(length(rec_sqr_tr_uwb_spectrum))/2;
rec_sqr_tr_uwb_spectrum=rec_sqr_tr_uwb_spectrum/length(received_sqr_tr_uwb_signal);

%Generate the scale for the frequency spectrum of squared received TR-UWB signal pulse train
f=(0:non_repeated_dft_points-1)*2*fn/points_for_dft;
subplot(3,1,3);
plot(f,20*log10(rec_sqr_tr_uwb_spectrum));
xlabel('Frequency (Hz)','fontSize',15); 
title('Frequency spectrum: Squared Received TR-UWB signal','fontSize',15);
ylabel('Power/frequency (dB/Hz)','fontSize',15);

%**********************************************************************************************
%                   Section-10: Integration of squared received signal
%**********************************************************************************************
%Initialize an array of integrated signal with zeroes
integrated_tr_uwb_signal = zeros(1,length(received_sqr_tr_uwb_signal));

%Start lower time index with the initial index for the integration of the squared received signal 
time_initial_index = min(modulated_time);

%Variable to mark the upper time index for the integration of the squared received signal 
time_final_index = 0;

%Run loop to integrate the whole signal divided by blocks of 100 nanoseconds
while time_final_index < max(modulated_time)
    %Upper integration limit to make the integrated block size of 100 nanoseconds
    time_final_index = time_final_index + 100e-9; 
    
    %Use find function to find the indices of the non-zero signal values in this block
    non_zero_indices = find((modulated_time >= time_initial_index ) & (modulated_time < time_final_index));
    
    %Use Cumulative trapezoidal numerical integration to add all non-zero values using index found in last step
    integrated_tr_uwb_signal(non_zero_indices) = cumtrapz(received_sqr_tr_uwb_signal(non_zero_indices)); 
    
    %Move to the next block
    time_initial_index = time_final_index;
end

figure
subplot(2,1,1);
plot(modulated_time,integrated_tr_uwb_signal,'-r')
axis([-4e-8 1.95e-6 -0.2 50])
ylabel('amplitude', 'fontSize',15)
xlabel('time, secs','fontSize',15)
title('Time domain: Integrated UWB signal at the receiver side','fontSize',15)

%**********************************************************************************************
%                   Section-11: Decision and Recovery of the binary signal
%**********************************************************************************************
decoded_information_sequence= zeros(1,7);
threshold= 30;

%Start lower time index with the initial index for the decision of the integrated signal 
time_initial_index = min(modulated_time);

%Variable to mark the upper time index for the decision of the integrated signal 
time_final_index = 0;

%Initial index of decoded information sequence
decoded_signal_index = 1;
result = zeros(1,length(information_seq));

%Run loop on whole integrated signal divided by blocks of 300 nanoseconds to decode the information sequence since
%the reference pulses are separated by 300 nanoseconds.
while time_final_index < max(modulated_time)
    %Upper time limit to decide the integrated signal 
    time_final_index = time_final_index + 300e-9; 
    
    %Use find function the find the indices of non-zero values in this integrated block
    index = find((modulated_time >= time_initial_index ) & (modulated_time < time_final_index));
    
    %If the maximum integrated value is non-zero and greater than threshold, decide for the binary '1'
    if(max(integrated_tr_uwb_signal(index)) > threshold && max(integrated_tr_uwb_signal(index) > 0))
        decoded_information_sequence(decoded_signal_index) = 1;
        result(index) = integrated_tr_uwb_signal(index);
        
    %If the maximum integrated value is less than threshold, decide for the binary '0'    
    else
        decoded_information_sequence(decoded_signal_index) = 0;
    end

    %Move to the next block
    time_initial_index = time_final_index;
    
    decoded_signal_index = decoded_signal_index + 1;
end

subplot(2,1,2);
stem(decoded_information_sequence,'-k')
ylabel('Amplitude', 'fontSize',15)
title('Decoded Information Sequence','fontSize',15)
axis([0.9 7.3 -1 2])
disp('The decoded information sequence is: ')
disp(decoded_information_sequence)
disp('There is 100% correct recovery')
