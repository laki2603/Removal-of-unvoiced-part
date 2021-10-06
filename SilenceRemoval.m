close all; clear all;

% read sound 
[data, fs] = audioread('intro.mp3');
plot(data); title('given data');
% normalize data
data = data / abs(max(data));
figure;
plot(data); title('normalised data');
% do framing
f_d = 0.025;

f_size = round(f_d * fs);
n = length(data)/f_size;
n = floor(n);
temp = 0;
for i = 1 : n
    frames(i, :) = data(temp +1 : temp + f_size);
    temp = temp + f_size
end




%% finding ZCR of all frames
[r, c] = size(frames);

for i = 1 : r
    
    x = frames(i, :);
    
    ZCR(i) = 0;
    for k = 1:length(x) - 1
    
        if ((x(k) < 0) && (x(k + 1) > 0 ))
            ZCR(i) = ZCR(i) + 1;
   
        elseif ((x(k) > 0) && (x(k + 1) < 0))
            ZCR(i) = ZCR(i) + 1;
        end
    end

    
end

% calculating rate
ZCRr = ZCR/length(x);
ZCRr = ZCRr/max(ZCRr); %within range of 1
f_size = round(f_d * fs);
zcr_wave = 0;
for j = 1 : length(ZCRr)
    l = length(zcr_wave);
    zcr_wave(l : l + f_size) = ZCRr(j);
end
figure;plot(zcr_wave); title('zcr wave');

% plot the ZCR with Signal
t = [0 : 1/fs : length(data)/fs]; % time in sec
t = t(1:end - 1);
t1 = [0 : 1/fs : length(zcr_wave)/fs];
t1 = t1(1:end - 1);
figure;
plot(t,data); hold on;
plot(t1,zcr_wave,'r','LineWidth',2); title('zcr with signal');
legend('Speach signal', 'ZCR');

% Silence Removal
id = find(ZCRr <= 0.14); 
fr_ws = frames(id,:); % frames without silence

% reconstruct signal
data_r = reshape(fr_ws',1,[]);

figure;
plot(data);hold on;
plot(data_r,'g'); title('speech without silence using zcr');
legend('Original signal', 'Signal after silence removal');

% By pre-emphasis
for i1 = 1 : r
   temp = frames(i1,:);
    sum1(i1) = 0;
    for j1 = 2 : length(temp)
        
        sum1(i1) = sum1(i1) + abs(temp(j1) - temp(j1 - 1));        
    end
    
    enr(i1) = sum1(i1)/sum(temp.^2); % pre-emphasized energy ratio
end

enr = enr./max(enr);

f_size = round(f_d * fs);
enr_wave = 0;
for j1 = 1 : length(enr)
    l1 = length(enr_wave);
    enr_wave(l1 : l1 + f_size) = enr(j1);
end


% plot the enr with Signal
T = [0 : 1/fs : length(data)/fs]; % time in sec
T = T(1:end - 1);
T1 = [0 : 1/fs : length(enr_wave)/fs];
T1 = T1(1:end - 1);
figure;
plot(T,data'); hold on;
plot(T1,enr_wave,'r','LineWidth',2); title('enr wave with signal')
legend('Speech Signal','Pre-Emph. Energy Ratio');

% Silence Removal
id1 = find(enr <= 0.1);
fr_ws1 = frames(id1,:); % frames without silence

% reconstruct signal
data_r1 = reshape(fr_ws1',1,[]);

figure;
plot(data);hold on;
plot(data_r1,'g'); title('speech without silence using pre-emphasis');
legend('Original signal', 'Signal after silence removal');