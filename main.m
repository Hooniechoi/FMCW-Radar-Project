clear 
close all
clc

%% Data load 
% radar data load 
D = dir();

% 파일 번호 설정  
file_num = 1;

filename = D(file_num).name;
file = [D(file_num).folder,'/', filename];

% data 길이 
Nscan = 1200;

% 샘플링 수
Sample_cnt = 2^13 ; 

cnt = 1;
for j = 1:Nscan
    SetName1 = ['/SCAN_', num2str(j,'%05d'), '/Sim_TimeData'];
    try
        tmp = reshape(h5read(file,SetName1), Sample_cnt, 8, []);
        for ch = 1:8
            timeData{ch}(:, : , cnt) = reshape(tmp(:, ch, :), Sample_cnt, []);
        end
        cnt = cnt + 1;
    catch
        disp('Error')
        break;
    end
end

num_scan = Nscan;

% channel 선택 
selected_channel = 4;
radar_data = double(reshape(timeData{selected_channel}(:,1,:), Sample_cnt, [])); 

%% Parameter
BW=3*10^9;
Tc=256*(1/(2*10^6));
Fs=10;
Fc=61*10^9;
c=3*10^8;
scan_duration=50*10^(-3); % scan 간 간격

NFFT=Sample_cnt;
freqindex = (0:NFFT/2-1)*(Fs/NFFT); 
rang_res = c/(2*BW);
time = [1:num_scan]*scan_duration;

% 1차 Bandpass Filtering
low_cutoff_1st = 0.1;       
high_cutoff_1st = 0.6;        
[b1, a1] = butter(4, [low_cutoff_1st high_cutoff_1st] / (Fs / 2), 'stop');
filtered_signal_1st = filtfilt(b1, a1, radar_data);

% 2차 Bandpass Filtering
low_cutoff_2nd = 0.9;
high_cutoff_2nd = 2.5;
[b2, a2] = butter(4, [low_cutoff_2nd high_cutoff_2nd] / (Fs / 2), 'bandpass');
filtered_signal = filtfilt(b2, a2, filtered_signal_1st);

% 필터링된 데이터로 FFT 수행
w = hamming(Sample_cnt);
signal_fft_once = zeros(size(filtered_signal(1:end/2,:)));
selected_range = 2; 

for scan_idx = 1:size(filtered_signal,2) 
    radar_data_EC = (filtered_signal(:,scan_idx) - mean(filtered_signal(:,scan_idx))) / std(filtered_signal(:,scan_idx));
    fft_inst_peak = fft((radar_data_EC - mean(radar_data_EC)).*w,NFFT);
    signal_fft_once(:,scan_idx) = fft_inst_peak(1:NFFT/2);
    magnitude(:,scan_idx) = abs(signal_fft_once(:,scan_idx)); 
    inst_phase_map(:,scan_idx) = angle(signal_fft_once(:,scan_idx)); 
end 
phase_map = unwrap(inst_phase_map,[],2);

%Phase map을 통한 거리 환산
sample_phase_index = 1:size(phase_map, 1);
selected_phase = phase_map(selected_range, :);
range_offset = (c / (4 * pi * Fc)) .* selected_phase;

freq_vector = find(freqindex >= low_cutoff_2nd & freqindex <= high_cutoff_2nd);

[~, max_idx] = max(magnitude(freq_vector, selected_range));
HR_frequency = freqindex(freq_vector(max_idx));
HR = HR_frequency * 60;


% range_offset 배열의 크기
N = length(range_offset);
range_offset_filtered = range_offset;
threshold = 0.009; 

for i = 2:N-1  
    diff_prev = abs(range_offset(i) - range_offset(i-1));
    diff_next = abs(range_offset(i) - range_offset(i+1));
    
    if diff_prev > threshold || diff_next > threshold
        range_offset_filtered(i) = mean([range_offset(i-2), range_offset(i+2)]);
    end
end

figure;
plot(time, range_offset_filtered, 'b-', 'DisplayName', 'Range Offset');
xlabel('시간 [s]','fontsize',15);
ylabel('거리 보정값','fontsize',15);
title('심박수 피크 검출 그래프');
set(gca,'fontsize',10) 

% 3차 Bandpass Filtering
low_cutoff_3rd = 0.9;
high_cutoff_3rd = 2.5;
[b3, a3] = butter(4, [low_cutoff_3rd high_cutoff_3rd] / (Fs / 2), 'bandpass');
range_offset_filtered_signal = filtfilt(b3, a3, range_offset_filtered);

% range_offset을 FFT
N_2nd = length(range_offset_filtered_signal); 
range_offset_fft = fft(range_offset_filtered_signal);

f = (0:N_2nd-1)*(Fs/N_2nd);  

magnitude_fft = abs(range_offset_fft);

half_idx = floor(N_2nd/2);

%% 주파수 스펙트럼 시각화
figure;
plot(f(2:half_idx), magnitude_fft(2:half_idx));
xlabel('Frequency (Hz)','fontsize',15);
ylabel('Magnitude','fontsize',15);
title('주파수 스펙트럼');

%% 가장 큰 주파수 성분의 위치 찾기
[max_magnitude, max_idx] = max(magnitude_fft(2:half_idx));  % 최대값과 그 인덱스 찾기
HR_real_frequency = f(max_idx + 1);  % 주파수 벡터에서 해당 주파수 찾기 (인덱스 보정 필요)

HR_real = HR_real_frequency * 60;
disp(['HR: ', num2str(HR_real), ' BPM']);

%% Filtered Range Offset and Peak Detection
figure;
plot(time, range_offset_filtered_signal, 'b-', 'DisplayName', 'Range Offset');
hold on;
[fil_peaks, fil_locs] = findpeaks(range_offset_filtered_signal, 'MinPeakDistance', Fs* 90 / HR);
plot(time(fil_locs), fil_peaks, 'ro', 'DisplayName', 'Detected Peaks');
xlabel('시간 [s]','fontsize',15); 
ylabel('거리 보정값 [m]','fontsize',15); 
title('심장 박동 파형');
set(gca,'fontsize',10) 
xlim([1 10]);
grid on;


num_detected_filt_peaks = numel(fil_peaks);
%disp(['심박 피크 수 : ', num2str(num_detected_filt_peaks)]);

if numel(fil_locs) > 1
    % 피크 간 간격 계산
    peak_intervals_filt = diff(time(fil_locs)); 
    
    % 피크 간 간격 시각화
    figure;
    plot(1:length(peak_intervals_filt), peak_intervals_filt, 'b-', 'DisplayName', '심박 간격');
    xlabel('시간 인덱스','fontsize',15);
    ylabel('심박 간격 [s]', 'fontsize',15);
    title('심박 간격 변동성');
    legend('show');
    ylim([0.3 1.8]);
    grid on;

else
    disp('피크 검출 안됨');
end

% 피크 간 간격의 표준 편차 계산
std_peak_intervals_int = std(peak_intervals_filt);
std_peak_intervals = std_peak_intervals_int*100;
disp(['SDNN : ', num2str(std_peak_intervals), 'ms']);

% RMSSD 계산
rmssd_int = sqrt(mean(diff(peak_intervals_filt).^2));
rmssd = rmssd_int*100;
disp(['RMSSD : ', num2str(rmssd), 'ms']);

HR_min = 50;   % HR 최솟값 (스트레스가 거의 없는 상태)
HR_max = 180;  % HR 최댓값 (매우 높은 스트레스 상태)

RMSSD_min = 10;  % RMSSD 최솟값 (매우 높은 스트레스 상태)
RMSSD_max = 50; % RMSSD 최댓값 (스트레스가 거의 없는 상태)

std_peak_min = 10;  % 피크 간 간격 표준 편차 최솟값 (매우 높은 스트레스 상태)
std_peak_max = 100;   % 피크 간 간격 표준 편차 최댓값 (스트레스가 거의 없는 상태)

%% HR 정규화
HR_norm = (HR_real - HR_min) / (HR_max - HR_min);
HR_norm = min(max(HR_norm, 0), 1);  % 0~1로 한정

%% RMSSD 정규화 (역수로 처리)
RMSSD_norm = (RMSSD_max - rmssd) / (RMSSD_max - RMSSD_min);
RMSSD_norm = min(max(RMSSD_norm, 0), 1);  % 0~1로 한정

%% 피크 간 표준 편차 정규화 (역수로 처리)
std_peak_norm = (std_peak_max - std_peak_intervals) / (std_peak_max - std_peak_min);
std_peak_norm = min(max(std_peak_norm, 0), 1);  % 0~1로 한정

%% 각 지표에 대한 가중치 설정
w_HR = 0.2;         % HR 가중치
w_RMSSD = 0.4;      % RMSSD 가중치
w_std_peak = 0.4;   % 피크 간 표준 편차 가중치

%% 스트레스 점수 계산 (0~100)
stress_score = 100 * (w_HR * HR_norm + w_RMSSD * RMSSD_norm + w_std_peak * std_peak_norm);

%% 결과 출력
disp(['스트레스 점수: ', num2str(stress_score), ' / 100']);

magfor = magnitude_fft * 10^6;
ranfor = range_offset_filtered_signal * 10^4;

