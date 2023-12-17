clear all; close all; clc;
%kod radi za proivoljan broj M, D, f, rezoluciju ugla, ali samo za
%omnidirekcionu karakteristku posto se to meni trazilo, a vec je dovoljno
%sporo radio pa nisam htela da ga prosirujem

%ulazni parametri
bbbb=144;
gggg=2020;
M=3+mod(bbbb,8);
D=2-mod(bbbb,10)/10;
U=1+mod(bbbb,4);
f=gggg*3-bbbb*8; 
dteta=1; %menjati rezoluciju
teta=[0:dteta:359];

%racunanje pozicija zvucnika na y osi
d=D/(M-1);
for br=1:M
    Tx(1,br)=(-1)*D/2+d*(br-1); % y osa
    Tx(2,br)=0; % x osa
end

r=10;       
fs=48000;
c=340;
jedinicni_impuls=zeros(1,fs); 
%jedinicni_impuls(1)=1;

%pozicije merenja 
for br=1:length(teta)
    Rx(1,br)=r*sin(teta(br)/180*pi); % y osa
    Rx(2,br)=r*cos(teta(br)/180*pi); % x osa
end

%prikaz pozicije izvora i prijemnih tačaka
figure(1), plot(Tx(2,:),Tx(1,:),'xr', LineWidth=2), hold on,
plot(Rx(1,:),Rx(2,:),'og', LineWidth=1);  
xlim([-r-1 r+1]);ylim([-r-1 r+1]);axis equal; grid on;
title('prikaz pozicije izvora i prijemnih tačaka');
legend('pozicije izvora','pozicije prijemnih tačaka');

%rastojanje svakog zvucnika od mesta merenja i kasnjenje
for br1=1:M
    for br2=1:length(teta)
        r(br1, br2)=sqrt(Rx(2,br2)^2 +(Rx(1,br2)-Tx(1,br1))^2);
        t(br1, br2)=ceil((r(br1, br2)./c).*fs);
    end
end

%odabir karakteristike usmerenosti
switch U % ne menjati U jer sam napravila samo omnidirekcionu
    case 1
        u_f=1;
    case 2
        u_f=cos(teta_u); %nije u funkciji
    case 3
        u_f=1/2*(1+cos(teta_u)); %nije u funkciji
    case 4
        u_f=1/4*(1+3*cos(teta_u)); %nije u funkciji
end

%fft
for br1=1:length(teta)
    %za svaki Rx
    x(br1,:)=jedinicni_impuls(:);
    for br2=1:M
        x(br1,t(br2, br1))=(sqrt(413/4*pi)/r(br2, br1)).*u_f;
    end
    X(br1,:)=fft(x(br1,:));
    usmerenost_za_f(br1)=X(br1,f);
end

%prikaz usmerenosti za frekvenciju 4908 Hz
figure(2), 
mmpolar(teta/180*pi,20*log10(usmerenost_za_f/max(abs(usmerenost_za_f))),'TTickDelta',15','RLimit',[0 -15],'TLimit',[-pi pi]); 
title('Prikaz usmerenosti za frekvenciju 4908 Hz'); 

%Oktavni spektar
br_oktava=7;
for br=1:br_oktava
    fc(br)=125*2^(br-1);
end
fg=round(fc/sqrt(2)); fg(br_oktava+1)=round(fc(br_oktava)*sqrt(2));

for br=1:7
    for br1=1:length(teta)
        y_r(br1,br)=sum(abs(X(br1,fg(br):fg(br+1))).^2);
    end
    figure,
    mmpolar(teta/180*pi,20*log10(y_r(:,br)/max(abs(y_r(:,br)))),'TTickDelta',15','RLimit',[0 -15],'TLimit',[-pi pi]); 
end