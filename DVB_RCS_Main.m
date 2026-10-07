% 按照DVB-RCS标准，实现返回链路信号处理流程：
% 支持突发类型 TRF ACQ SYNC CSC
%支持调制方式 pi/2 BPSK QPSK 8PSK 16QAM(TC-LM) CPM(CC-CPM)
%TC-LM模式：突发成帧->能量扩散->CRC->Turbo->调制->突发构造->基带成型
%CC-CPM模式：突发成帧->能量扩散->CRC->RS+卷积->交织->CPM调制->突发构造->基带成型
clear;clc;close all;
%##############################仿真参数配置##########################
burst_type = 'TRF';
%编码方案 'concatenated'，'turbo'
coding_scheme = 'concatenated';
%滚降因子（0.2,0.25,0.35）
rolloff = 0.2;
%过采样率（每个符号采样点数）
sps = 4;
%调制阶数（TC-LM模式） 1=pi/2 BPSK 2=QPSK 3=8PSK 4=16QAM
mod_order = 2;

%载波频率
fc = 14e5; %1.4M波段

%采样率
Rs = 2e6;%符号率
fs = sps * Rs;%TC-LM 采样率（sps 点/符号）
waveform_cfg.waveform_id = 1;
%获取DVB-RCS参数
params = get_rcs_params(burst_type,coding_scheme,waveform_cfg);
input_len = params.input_len;

rng(42);%固定种子保证结果可重复性
burst_frame = randi([0,1],input_len,1);
%能量扩散
scrambled_data = rcs_energy_dispersal(burst_frame);

%CRC校验
if strcmp(params.crc_type,'crc32')
    [data_with_crc,crc_value] = rcs_crc32_encode(scrambled_data);
else
    [data_with_crc,crc_value] = rcs_crc16_encode(scrambled_data);
end

%信道编码
coded_data = rcs_channel_coding(data_with_crc,params);

%比特交织
if strcmp(params.mod_type,'cpm')
    interleaved_data = rcs_bit_interleaving(coded_data,params);
else
    interleaved_data = coded_data;%turbo编码不交织
end
%调制
if strcmp(params.mod_type,'cpm')
    [alpha_burst, burst_info] = rcs_burst_construction_cpm(interleaved_data, params);
    %CC-CPM 突发符号序列：UW|数据1|尾符号|UW(中导)|数据2|尾符号（§7.3.6.2）
    [symbols, cpm_info] = rcs_cpm_modulate(alpha_burst, params);
    burst_signal = symbols;
else
    symbols = rcs_modulation(interleaved_data, params, false);
    [burst_signal,burst_info] = rcs_burst_construction(symbols,params);
end



%基带成型
I_data = real(burst_signal);
Q_data = imag(burst_signal);
if strcmp(params.mod_type, 'cpm')
    %CC-CPM：脉冲成形（q_AV/g_AV）已在 rcs_cpm_modulation 内完成，
    %        burst_signal 就是 cpm_samples_per_sym 点/符号的复基带，
    %        这里不能再插值/滤波（否则变成 32*4=128 点/符号）
    I_shaped = I_data;
    Q_shaped = Q_data;
    h_rrc = 1;
    fs_use = params.cpm_samples_per_sym * Rs;%CPM 实际采样率
else

    [I_shaped,Q_shaped,h_rrc] = baseband_shaping(I_data, Q_data, rolloff, sps);
    fs_use = fs;
end

%正交调制
rf_signal = quadrature_modulation(I_shaped, Q_shaped, fc, fs_use);
fprintf('RF 信号长度：%d 采样点\n',length(rf_signal));
fprintf('RF 信号功率：%.4f \n\n',mean(rf_signal.^2));

%************************画图************************************
%星座图
figure;
plot(real(symbols), imag(symbols), '*');
title('星座图');

%滤波器
t_filt = (-floor(length(h_rrc)/2):floor(length(h_rrc)/2))/sps;
figure;
plot(t_filt, h_rrc, 'b-', 'LineWidth', 1.5);
title('滤波器-脉冲响应');


[H, f] = freqz(h_rrc, 1, 1024, 'whole');
H = fftshift(H);
f_norm = f / (2*pi) *sps;
figure;
plot(f_norm, 20*log10(abs(H)/max(abs(H))), 'r-', 'LineWidth', 1.5);
title('滤波器-频率响应');

%基带信号
len = length(I_shaped);
figure;
plot(I_shaped(1:len), 'b-', 'LineWidth', 0.8);
title('I路基带信号');

figure;
plot(Q_shaped(1:len), 'b-', 'LineWidth', 0.8);
title('Q路基带信号');

%射频信号
figure;
n_show_rf = min(2000, length(rf_signal));
t_rf = (0:n_show_rf -1) / fs_use*1e6;%微秒
plot(t_rf, rf_signal(1:n_show_rf),'-b');
title('射频时域信号波形');

[Pxx_rf, f_rf] = pwelch(rf_signal, hanning(512), 256, 1024, fs_use, 'centered');
figure;
plot((f_rf + fc)/1e6, 10*log10(Pxx_rf), 'r-');
title('射频信号功率谱密度');






