%%
%% load the data as p, q, t and etc
unit_m3_to_ml = 1e6;
unit_Pa_to_mmHg = 1/133.32;
data = dlmread('./55arteryNetwork_age19.csv',';',1,0);
t = data(:,1);
p = data(:,2)*unit_Pa_to_mmHg;
q = data(:,3)*unit_m3_to_ml;
T = t(2)-t(1);
Fs = 1/T;
L = length(q);

%%
Rp = 1.023;
C = 1.36;

% 1.c
f= 0:0.1:20;
w = 2*pi*f;
Z = ZWK(w,Rp,C);

figure
plot(f,log(abs(Z)))
ylabel('|Z|')
xlabel('Freq [Hz]')

figure
plot(f,angle(Z)*180/pi)
xlabel('Freq [Hz]')
ylabel('\angle Z')




Qm = fft(q);
f = Fs*(0:L-1)/L;
w = 2*pi*f;
Z = ZWK(w,Rp,C);
pwk=real_ifft(transpose(Z).*Qm, w, t);

% plot
figure()
plot(t, pwk);
hold all;
plot(t,p);
xlabel('Time')
ylabel('Pressure')
%legend('Pwk','Pdata')



%%
Rp = 1.023;
C = 1.06;

pwk1 = @(C) real_ifft(transpose(ZWK(w,Rp,C)).*Qm, w, t);
res = @(C) sum((pwk1(C) - p).^2);

C_opt = fminsearch(res,C);

% plot
plot(t,p,'ko');
hold all
plot(t,pwk1(C));
plot(t,pwk1(C_opt));

xlabel('Time')
ylabel('Pressure')
legend('Pwk','Pdata')


%% 
Pm = fft(p);
Zm = Pm/Qm;
Zdc = Zm(1);
pwk2 = @(C) real_ifft(transpose(ZWK(w,Zdc,C)).*Qm, w, t);
res = @(C) sum((pwk2(C) - p).^2);

C_opt = fminsearch(res,C);

% plot
figure
plot(t,p,'ko');
hold all
plot(t,pwk2(C));
figure
plot(t,pwk2(C_opt),t,pwk1(C_opt));

xlabel('Time')
ylabel('Pressure')
legend('Pwk','Pdata')
