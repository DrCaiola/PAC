%%Parallel Phase-Amplitude Coupling a la Tort
%Inputs:    x - input signal (as column or row vector)
%           freq - 1) [a1 a2;b1 b2] makes 1-D MI up to 2 diagrams
%           (depending on TOL). a1 and b1 are lowcutoff frequencies and a2
%           and b2 are hicutoff frequencies. The order doesn't matter as
%           all combinations are checked.
%                  2) n x 2 matrix allows for all phase amplitude pairs
%                  combinations to be compared.
%                  3) [] makes 2-D MI diagram
%           fs - sampling rate
%           Num_bin - number of bins (default is 18)
%Optional: 'average', <cut> - average data over <cut> seconds (will turn
%               off 'savemore')
%Optional: 'plotless' - does not display plot
%Optional: 'filtorder' - change filtorder
%Optional: 'savemore' - saves all filtered signals and hilbert signals
%               currently only on 2-D code
%Optional: 'meanless' - turns off mean correcting feature (does not work on
%               pop_eegfiltnew)
%Optional: 'debug' - allows interactive filter view at individual PA pairs
%Optional: 'trueplot' - plots overlapping patches per bin size
%
%Optional: 'filter', <filter name> - change filter from pop_eegfiltnew
%   Filter options: 'old'       - uses eegfilt.
%                   'firls'     - uses firls (requires 'filtorder').
%                   'fir1'      - uses fir1 with hamming windowing
%                                   (requires 'filtorder').
%                   'kaiser'    - uses fir1 with kaiser windowing.
%                   'cheby'     - uses cheby1.
%                   'ellip'     - uses ellip.
%                   'ellipp'    - uses ellip with order calculation.
%                   'butter'    - uses butter filter.
%                   'variable'  - uses larger Amplitude bins to account for
%                                   increasing phase and butter for filter.
%                   'wavelet'    - Morlet wavelet with optional paramters.
%
%Example 1: Take signal x (sampled at 1000 Hz) and compute the 1-D PAC for
%   phase 20 and amplitude 100 with a buffer of 5 Hz on both sides.
%
%   MI=PAC_par(x,[15,25;95,105],1000);
%   Output: MI vector with two entries (20 phase and 100 amp & 100 phase
%   and 20 amp) and up to two 1-D MI graphs (depending if the MI is over
%   TOL).
%
%Example 2: Take signal x (sampled at 1000 Hz) and compute the 1-D PAC for
%   phase 20 and amplitude 100 and also for phase 40 and amplitude 100 with
%   5 Hz buffer.
%
%   MI=PAC_par(x,[15,25;95,105;35,45],1000);
%   Output: MI vector with 6 entries ([20,100], [20,40], [100,20],
%   [100,40], [40,20], [40,100]) and up to 6 1-D MI graphs (depending if
%   the MI is over TOL).
%
%Example 3: Take signal x (sampled at 1000 Hz) and compute the 2-D PAC.
%
%   MI=PAC_par(x,[],1000);
%   Output: MI and 2-D PAC with default parameters.
%
%Example 4: Take signal x (sampled at 1000 Hz) and compute the 2-D PAC with
%   Phase from 4-50 in intervals of 2Hz and Amplitude from 45-200 with
%   intervals of 5Hz.
%
%   MI=PAC_par(x,[],1000,[],'param',[4 50 2 45 100 5]);
%   Output: MI and 2-D PAC with given parameters.
%
%Example 5: Take 100 second signal, x (sampled at 1000 Hz), and compute the
%   MI for every 10 seconds and average them together into one 2-D plot and
%   one average MI.
%
%   [MI,MI_all]=PAC_par(x,[],1000,[],'average',10);
%   Output: Average MI and vector of 10 individual MIs in MI_all. Plot
%   output is a 2-D PAC.
%
%Example 6: Take signal x (sampled at 1000 Hz) and compute the 2-D PAC with
%   Phase from 2-50 in intervals of 2Hz and AMplitude from 5-400 with
%   intervales of 5Hz with a filter order of 4 with no output of plots.
%
%   MI=PAC_par(x,[],1000,[],'param',[2 50 2 5 400 5], 'filter','variable','filtorder',4,'plotless');
%   Output: MI calculated with variable bin filter.
%
%  Code based on the MI procedures by Tort et al.
%
%
%LAST UPDATED: 10.24.2017 by Mike Caiola
%   changelog: 09.07.16 - Began code branch to add parallel option
%              09.13.16 - Bug fixes
%              11.10.16 - Added Variable filter and documentation
%              11.21.16 - Consolidated 1D into 2D code
%              12.06.16 - Added 'debug' option to allow interactive looks
%                           at the psd of the completed PSDs
%              12.07.16 - Added legacy option to run older versions of code
%              12.26.16 - Allowed 'Average' to truncate signal
%              03.02.17 - Began to remove uneeded documentation
%              09.26.17 - Added MI bar graph into 'debug' option
%              10.06.17 - Bug Fix
%              10.10.17 - Added 'trueplot' option
%              10.19.17 - Added wavelet filter
%              10.20.17 - Added wavlet options, removed TOL and hybrid
%              10.24.17 - Added wavelet3d option, quite slow still
%              10.26.17 - Disabled wavelet3d and added superior PAC3D
function [MI,MI_all]=PAC_par(x,freq,fs,Num_bin,varargin)
%% Setup
%Rearrange signal
if size(x,2)==1
    x=x';
end
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary')
%Optional parameters
opt_all=0;
opt_avg=0;
opt_diag=0;
opt_plot=1;
opt_param=0;
opt_circle=0;
opt_meanc=1;
opt_filter=0;
opt_old=0;
opt_firls=0;
opt_fir1=0;
opt_kaiser=0;
opt_butter=0;
opt_cheby=0;
opt_ellip=0;
opt_ellipp=0;
opt_wavelet=0;
opt_save=0;
opt_var=0;
opt_debug=0;
opt_runtime=0;
opt_trueplot=0;
opt_w3D=0;
woptions={};
dfp=[]; df=[];
cut=length(x);
filtorder=[];
if nargin<4 || isempty(Num_bin)
    Num_bin=18;
end
numvarargs= length(varargin);
if numvarargs~=0; k=1;
    if isempty(varargin{1})
        k=k+1;
    end
    while k<=numvarargs
        switch varargin{k}
            case 'plotless'
                opt_plot=0;
            case 'param'
                opt_param=1;
                param=varargin{k+1};
                k=k+1;
            case 'filtorder'
                filtorder=varargin{k+1};
                k=k+1;
            case 'circle'
                opt_circle=1;
            case 'meanless'
                opt_meanc=0;
            case 'runtime'
                opt_runtime=1;
                tstart=tic;
            case 'filter'
                opt_filter=1;
                switch varargin{k+1}
                    case 'old'
                        opt_old=1;
                        filtername=varargin{k+1};
                    case 'firls'
                        opt_firls=1;
                        filtername=varargin{k+1};
                    case 'fir1'
                        opt_fir1=1;
                        filtername=varargin{k+1};
                    case 'kaiser'
                        opt_kaiser=1;
                        filtername=varargin{k+1};
                    case 'butter'
                        opt_butter=1;
                        filtername=varargin{k+1};
                    case 'cheby'
                        opt_cheby=1;
                        filtername=varargin{k+1};
                    case 'ellip'
                        opt_ellip=1;
                        filtername=varargin{k+1};
                    case 'ellipp'
                        opt_ellipp=1;
                        filtername='ellip plus';
                    case 'variable'
                        opt_var=1;
                        filtername='Butter w/ variable bins';
                    case 'wavelet'
                        opt_wavelet=1;
                        filtername=varargin{k+1};
                end
                k=k+1;
            case 'savemore'
                opt_save=1;
            case 'wavelettime'
                woptions=[woptions,'time',varargin{k+1}];
                k=k+1;
            case 'waveletcycles'
                woptions=[woptions,'cycles',varargin{k+1}];
                k=k+1;
            case 'waveletlog'
                woptions=[woptions,'log'];
            case 'debug'
                opt_debug=1;
            case 'trueplot'
                opt_trueplot=1;
            case 'PAC3D'
                opt_w3D=1;
                m_sec=30;
                dt=1;
                plotnum=0;
                if k~=numvarargs && isnumeric(varargin{k+1})
                    w=varargin{k+1};
                    if length(w)==1
                        plotnum=w;
                    elseif length(w)==2
                        m_sec=w(1);
                        dt=w(2);
                    else
                        m_sec=w(1);
                        dt=w(2);
                        plotnum=w(3);
                    end
                    k=k+1;
                end
            case 'average'
                opt_avg=1;
                cut=varargin{k+1}*fs;
                k=k+1;
                if mod(length(x),cut)
                    x=x(1:end-mod(length(x),cut));
                    warning('cut does not divide equally: signal truncated')
                end
                if cut==0 || cut>length(x)
                    cut=length(x);
                    fprintf('Removed cut \n')
                end
                opt_save=0;
        end
        k=k+1;
    end
end
if isempty(freq)
    opt_all=1;
elseif isequal(freq(2,1),freq(2,2)) && ~opt_var
    opt_var=1;
    opt_filter=1;
    filtername='Butter w/ variable bins';
    warning('Using Variable Amplitude Bins')
else
    opt_save=0;
end
if ~opt_filter
    filtername='new';
end
if (opt_firls || opt_fir1) && isempty(filtorder)
    filtorder=1000;
    warning('Default filtorder of 1000 used');
end
if (opt_butter || opt_ellip || opt_cheby || opt_ellipp || opt_var) && isempty(filtorder)
    filtorder=2;
    warning('Default filtorder of 2 used');
end
if ~(opt_firls || opt_fir1 || opt_butter || opt_ellip || opt_cheby || opt_ellipp || opt_old) && opt_meanc==0
    warning('Mean may still be corrected')
end
if opt_circle && opt_var
    opt_circle=0;
    warning('Circle option turned off')
end
%Subtract Mean
if opt_meanc
    x=x-mean(x);
end
[pxx,f0]=pwelch(x,[],[],[],fs);

x_original=x;   %in case of multiple passes

bin_length=2*pi/Num_bin;

%% Code
for cuti=1:length(x_original)/cut    %Large for-loop for multiple run-throughs
    x=x_original(cut*(cuti-1)+1:cut*cuti);
    if opt_all
        if opt_param
            Phase_min=param(1);
            Phase_max=param(2);
            Phase_step=param(3);
            Amp_min=param(4);
            Amp_max=param(5);
            Amp_step=param(6);
        else
            Phase_min=4; %Default 2-D Parameters
            Phase_max=52;
            Phase_step=2;
            Amp_min=15;
            Amp_max=400;
            Amp_step=5;
            param=[Phase_min Phase_max Phase_step Amp_min Amp_max Amp_step];
        end
    else
        Phase_min=(freq(1,1)+freq(1,2))/2;
        Phase_max=Phase_min;
        Phase_step=freq(1,2)-freq(1,1);
        Amp_min=(freq(2,1)+freq(2,2))/2;
        Amp_max=Amp_min;
        if opt_var
            Amp_step=2*Phase_min;
            freq(2,:)=[Amp_min-Phase_min Amp_min+Phase_min];
        else
            Amp_step=freq(2,2)-freq(2,1);
        end
    end
    pb_s=Phase_min:Phase_step:Phase_max;
    lpb_s=length(pb_s);
    Phase_bin=zeros(lpb_s,2);
    theta=zeros(lpb_s,length(x));
    ab_s=Amp_min:Amp_step:Amp_max;
    lab_s=length(ab_s);
    Amp_bin=zeros(lab_s,2);
    %Ampp=zeros(lab_s,length(x));
    Amp=zeros(lab_s,length(x));
    dbg_bin=zeros(lab_s,2);
    dbg_binp=zeros(lpb_s,2);
    for i=1:lpb_s
        Phase_bin(i,:)=[pb_s(i)-Phase_step pb_s(i)+Phase_step];
    end
    for i=1:lab_s
        if ~opt_var
            Amp_bin(i,:)=[ab_s(i)-Amp_step ab_s(i)+Amp_step];
        end
    end
    if cuti==1
        MI=zeros((Amp_max-Amp_min)/Amp_step+1,(Phase_max-Phase_min)/Phase_step+1,length(x_original)/cut);
    end
    if opt_wavelet
        theta=Wavletfilter(x,fs,pb_s,'phase',woptions{:});
    else
        EEG=struct('data',x,'srate',fs,'pnts',length(x));       %Structure needed for pop_eegfiltnew
        parfor ip=1:lpb_s
            if opt_old
                [y,bp]=eegfilt(x,fs,Phase_bin(ip,1),Phase_bin(ip,1)+2*Phase_step,0,filtorder);
            elseif opt_firls
                bp=firls(filtorder,[0 (1-.15)*Phase_bin(ip,1)/(fs/2) Phase_bin(ip,1)/(fs/2) (Phase_bin(ip,1)+2*Phase_step)/(fs/2) (1+.15)*(Phase_bin(ip,1)+2*Phase_step)/(fs/2) 1],[0 0 1 1 0 0]);
                y=filtfilt(bp,1,x);
            elseif opt_fir1
                bp=fir1(filtorder,[Phase_bin(ip,1) (Phase_bin(ip,1)+2*Phase_step)]/(fs/2));
                y=filtfilt(bp,1,x);
            elseif opt_butter
                [bp,ap]=butter(filtorder,[max([Phase_bin(ip,1),.1]) (Phase_bin(ip,1)+2*Phase_step)]/(fs/2));
                y=filtfilt(bp,ap,x);
            elseif opt_ellip
                [bp,ap]=ellip(filtorder,.05,50,[Phase_bin(ip,1) (Phase_bin(ip,1)+2*Phase_step)]/(fs/2));
                y=filtfilt(bp,ap,x);
            elseif opt_ellipp
                [n,Wp] = ellipord([Phase_bin(ip,1) (Phase_bin(ip,1)+2*Phase_step)]/(fs/2),[Phase_bin(ip,1)-(min(max(Phase_bin(ip,1)/4, 2),Phase_bin(ip,1))-.00001),(Phase_bin(ip,1)+2*Phase_step)+(min(max(Phase_bin(ip,1)/4, 2),Phase_bin(ip,1))-.00001)]/(fs/2),.05,50);
                [bp,ap] = ellip(n,.05,50,Wp);
                y=filter(bp,ap,x);
            elseif opt_cheby
                [bp,ap]=cheby1(filtorder,.05,[Phase_bin(ip,1) (Phase_bin(ip,1)+2*Phase_step)]/(fs/2));
                y=filtfilt(bp,ap,x);
            elseif opt_kaiser
                [n,Wn,beta,ftype] = kaiserord([Phase_bin(ip,1)-(min(max(Phase_bin(ip,1)/4, 2),Phase_bin(ip,1))),Phase_bin(ip,1), (Phase_bin(ip,1)+2*Phase_step),(Phase_bin(ip,1)+2*Phase_step)+(min(max(Phase_bin(ip,1)/4, 2),Phase_bin(ip,1)))],[0 1 0],[.01 .05 .01],fs);
                n = n + rem(n,2);
                bp = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
                y=kfilt(x,bp);
            elseif opt_var
                [bp,ap]=butter(filtorder,[max([Phase_bin(ip,1),.1]) (Phase_bin(ip,1)+2*Phase_step)]/(fs/2));
                y=filtfilt(bp,ap,x);
            else
                [y,~,bp]=pop_eegfiltnew(EEG,Phase_bin(ip,1), (Phase_bin(ip,1)+2*Phase_step),filtorder);
            end
            if opt_debug
                dbg_binp(ip,:)=[Phase_bin(ip,1) (Phase_bin(ip,1)+2*Phase_step)];
            end
            theta(ip,:)=angle(hilbert(y));
        end
    end
    %close(h);
    if ~opt_var
        if opt_wavelet
            Amp=Wavletfilter(x,fs,ab_s,'amplitude',woptions{:});
        else
            parfor jp=1:lab_s
                if opt_old
                    [y,ba]=eegfilt(x,fs,(Amp_bin(jp,1)), (Amp_bin(jp,1)+2*Amp_step),0,filtorder);
                elseif opt_firls
                    ba=firls(filtorder,[0 (1-.15)*(Amp_bin(jp,1))/(fs/2) (Amp_bin(jp,1))/(fs/2) (Amp_bin(jp,1)+2*Amp_step)/(fs/2) (1+.15)*(Amp_bin(jp,1)+2*Amp_step)/(fs/2) 1],[0 0 1 1 0 0]);
                    y=filtfilt(ba,1,x);
                elseif opt_fir1
                    ba=fir1(filtorder,[(Amp_bin(jp,1)) (Amp_bin(jp,1)+2*Amp_step)]/(fs/2));
                    y=filtfilt(ba,1,x);
                elseif opt_butter
                    [ba,aa]=butter(filtorder,[(Amp_bin(jp,1)) (Amp_bin(jp,1)+2*Amp_step)]/(fs/2));
                    y=filtfilt(ba,aa,x);
                elseif opt_ellip
                    [ba,aa]=ellip(filtorder,.05,50,[(Amp_bin(jp,1)) (Amp_bin(jp,1)+2*Amp_step)]/(fs/2));
                    y=filtfilt(ba,aa,x);
                elseif opt_ellipp
                    [n,Wp] = ellipord([(Amp_bin(jp,1)) (Amp_bin(jp,1)+2*Amp_step)]/(fs/2),[(Amp_bin(jp,1))-(min(max((Amp_bin(jp,1))/4, 2),(Amp_bin(jp,1)))-.00001),(Amp_bin(jp,1)+2*Amp_step)+(min(max((Amp_bin(jp,1))/4, 2),(Amp_bin(jp,1)))-.00001)]/(fs/2),.05,50);
                    [ba,aa] = ellip(n,.05,50,Wp);
                    y=filter(ba,aa,x);
                elseif opt_cheby
                    [ba,aa]=cheby1(filtorder,.05,[(Amp_bin(jp,1)) (Amp_bin(jp,1)+2*Amp_step)]/(fs/2));
                    y=filtfilt(ba,aa,x);
                elseif opt_kaiser
                    [n,Wn,beta,ftype] = kaiserord([(Amp_bin(jp,1))-(min(max((Amp_bin(jp,1))/4, 2),(Amp_bin(jp,1)))),(Amp_bin(jp,1)), (Amp_bin(jp,1)+2*Amp_step), (Amp_bin(jp,1)+2*Amp_step)+(min(max((Amp_bin(jp,1))/4, 2),(Amp_bin(jp,1))))],[0 1 0],[.01 .05 .01],fs);
                    n = n + rem(n,2);
                    ba = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
                    y=kfilt(x,ba);
                else
                    [y,~,ba]=pop_eegfiltnew(EEG,(Amp_bin(jp,1)), (Amp_bin(jp,1)+2*Amp_step),filtorder);
                end
                if opt_debug
                    if ~(opt_butter || opt_ellip || opt_ellipp || opt_cheby)
                        aa=1;
                    end
                    dbg_coeffa{jp}=aa;
                    dbg_coeffb{jp}=ba;
                    dbg_bin(jp,:)=[max([Amp_bin(jp,1),.1]) (Amp_bin(jp,1)+2*Amp_step)];
                end
                Amp(jp,:)=abs(hilbert(y));
            end
        end
    end
    if ~opt_all
        MeanA=zeros(1,Num_bin);
    end
    if opt_var
        dbg_bin=zeros(lab_s,lpb_s,2);
        if opt_w3D
            tAmp=zeros(lab_s,lpb_s,length(x));
        end
    end
    parfor i=1:lpb_s
        for j=1:lab_s
            if opt_circle
                MI(j,i,cuti)=abs(mean(Amp(j,:).*exp(1i*theta(i,:))));
            elseif opt_var
                [b,a]=butter(filtorder,[max([.1,ab_s(j)-pb_s(i)]) ab_s(j)+pb_s(i)]/(fs/2));
                y=filtfilt(b,a,x);
                if opt_debug
                    dbg_coeffbp{j,i}=b;
                    dbg_coeffap{j,i}=a;
                end
                dbg_bin(j,i,:)=[max([.1,ab_s(j)-pb_s(i)]) ab_s(j)+pb_s(i)];
                Ampp=abs(hilbert(y));
                if opt_w3D
                    tAmp(j,i,:)=Ampp;
                end
                MeanAp=zeros(1,Num_bin);
                if ~isnan(Ampp)
                    for kp=1:Num_bin
                        in=find(theta(i,:)>=(kp-1)*bin_length-pi & theta(i,:)<kp*bin_length-pi);
                        if ~isempty(in)
                            MeanAp(kp)=mean(Ampp(in));
                        else
                            MeanAp(kp)=0;
                        end
                    end
                end
                MI(j,i,cuti)=sum(MeanAp/sum(MeanAp).*(log(MeanAp/sum(MeanAp))/log(Num_bin)))+1;
                if isnan(MI(j,i,cuti))
                    MI_inf=MeanAp/sum(MeanAp).*(log(MeanAp/sum(MeanAp))/log(Num_bin)+1);
                    MI_inf(isnan(MI_inf))=0;
                    MI(j,i,cuti)=sum(MI_inf);
                end
                if ~opt_all
                    MeanA(i,:)=MeanAp;
                end
            else
                MeanAp=zeros(1,Num_bin);
                if ~isnan(Amp(j,:))
                    for kp=1:Num_bin
                        in=find(theta(i,:)>=(kp-1)*bin_length-pi & theta(i,:)<kp*bin_length-pi);
                        if ~isempty(in)
                            MeanAp(kp)=mean(Amp(j,in));
                        else
                            MeanAp(kp)=0;
                        end
                    end
                end
                MI(j,i,cuti)=sum(MeanAp/sum(MeanAp).*(log(MeanAp/sum(MeanAp))/log(Num_bin)))+1;
                if isnan(MI(j,i,cuti))
                    MI_inf=MeanAp/sum(MeanAp).*(log(MeanAp/sum(MeanAp))/log(Num_bin)+1);
                    MI_inf(isnan(MI_inf))=0;
                    MI(j,i,cuti)=sum(MI_inf);
                end
                if ~opt_all
                    MeanA(i,:)=MeanAp;
                end
            end
        end
    end
end

MI_all=MI;      %Average Multi run-throughs
MI=mean(MI,3);
if opt_plot     %Plotting procedure
    if opt_all  %2-D
        if opt_trueplot
            P=zeros(4,lpb_s);
            for i=1:lpb_s
                P(:,i)=[Phase_bin(i,1); Phase_bin(i,2); Phase_bin(i,2); Phase_bin(i,1)];
            end
            MI_tp=zeros(size(MI,2),1);
            for i=1:size(MI,2)
                MI_tp(1+(i-1)*size(MI,1):i*size(MI,1))=MI(:,i);
            end
            figure; hold all;
            k=0;
            if opt_var
                for i=1:lpb_s
                    for j=1:lab_s
                        A1=[dbg_bin(j,i,1);dbg_bin(j,i,1);dbg_bin(j,i,2);dbg_bin(j,i,2)];
                        k=k+1;
                        patch(P(:,i),A1,MI_tp(k),'FaceAlpha',.2,'EdgeColor','none')
                    end
                end
                A=[dbg_bin(1,1,1),dbg_bin(end,end,end)];
            else
                A=zeros(4,lab_s);
                for i=1:lab_s
                    A(:,i)=[Amp_bin(i,1); Amp_bin(i,1); Amp_bin(i,2); Amp_bin(i,2)];
                end
                for i=1:lpb_s
                    for j=1:lab_s
                        k=k+1;
                        patch(P(:,i),A(:,j),MI_tp(k),'FaceAlpha',.2,'EdgeColor','none')
                    end
                end
            end
            if opt_filter
                title(['True Plot (Filter: ', filtername,')'])
            else
                title('True Plot')
            end
            axis([P(1,1),P(end,end),A(1,1),A(end,end)])
        else
            figure; hold all; contourf(Phase_min:Phase_step:Phase_max,Amp_min:Amp_step:Amp_max,MI,30,'lines','none')
            if opt_filter
                title(['(Filter: ' filtername,')'])
            end
        end
        xlabel('Phase')
        ylabel('Amplitude')
        colorbar;
        plot(ceil(Amp_min/2):Phase_step:Phase_max,2*(ceil(Amp_min/2):Phase_step:Phase_max),'w--')
        if opt_debug
            button=1;
            while button<3
                title({'DEBUG MODE';'Click to investigate point'})
                [xin,yin] = ginput(1);
                if 2*mod(xin,Phase_step)<Phase_step %Round xin
                    dbg_x=xin-mod(xin,Phase_step);
                else
                    dbg_x=xin+(Phase_step-mod(xin,Phase_step));
                end
                if 2*mod(yin,Amp_step)<Amp_step %Round yin
                    dbg_y=yin-mod(yin,Amp_step);
                else
                    dbg_y=yin+(Amp_step-mod(yin,Amp_step));
                end
                dbg_xind=find(pb_s==dbg_x);
                dbg_yind=find(ab_s==dbg_y);
                if opt_var
                    filterview(x,dbg_coeffbp{dbg_yind,dbg_xind},dbg_coeffap{dbg_yind,dbg_xind},pxx,f0,fs,dbg_bin(dbg_yind,dbg_xind,:),dbg_binp(dbg_xind,:),filtername,filtorder,Num_bin)
                    title({['Amplitude',' : ',num2str(dbg_bin(dbg_yind,dbg_xind,1)),' - ',num2str(dbg_bin(dbg_yind,dbg_xind,2)), ' Hz : ',filtername];['@ Phase = ', num2str(dbg_x)]});
                else
                    filterview(x,dbg_coeffb{dbg_yind},dbg_coeffa{dbg_yind},pxx,f0,fs,dbg_bin(dbg_yind,:),dbg_binp(dbg_xind,:),filtername,filtorder,Num_bin)
                end
                button2 = uicontrol('Style','pushbutton','String',...
                    'Close','Callback',@callbackfun);
                uiwait(gcf);
                title({'DEBUG MODE';'Left click to continue. Right click to close.'})
                [~,~,button] = ginput(1);
            end
        end
    else        %1-D
        bin_lengthD=360/Num_bin;
        figure;
        if opt_circle
            hold on
            plot(Amp(1,:).*exp(1i*theta(1,:)));
            quiver(0,0,real(MI),imag(MI))
            xlabel('Real Part')
            ylabel('Imaginary Part')
            title({[num2str(freq(1,1)), ' - ' num2str(freq(1,2)) ,' Hz Phase vs ' ,num2str(freq(2,1)), ' - ' num2str(freq(2,2)),' Hz Amplitude'];
                ['Mean Vector Length = ', num2str(abs(MI))]})
            axis equal
        else
            bar(0:bin_lengthD:720-bin_lengthD,[MeanA(1,:) MeanA(1,:)]/sum(MeanA(1,:)))
            xlim([-10 720])
            if opt_filter
                title({[num2str(freq(1,1)), ' - ' num2str(freq(1,2)) ,' Hz Phase vs ' ,num2str(freq(2,1)), ' - ' num2str(freq(2,2)),' Hz Amplitude'];
                    ['MI = ', num2str(MI)];
                    ['(Filter: ' filtername,')']})
            else
                title({[num2str(freq(1,1)), ' - ' num2str(freq(1,2)) ,' Hz Phase vs ' ,num2str(freq(2,1)), ' - ' num2str(freq(2,2)),' Hz Amplitude'];
                    ['MI = ', num2str(MI)]})
            end
            xlabel('Phase')
            ylabel('Amplitude')
        end
    end
end
if opt_w3D
    if opt_var
        Amp=tAmp;
    end
    PAC3D(x,theta,Amp,fs,m_sec,dt,param,Num_bin,plotnum);
end
if opt_save
    MIs.MI=MI;
    MI=MIs;
end
if opt_runtime
    toc(tstart)
end
end

%% SUBROUTINES
function callbackfun(src,event)
closereq
uiresume(gcbf)
end
function filterview(x,coeffb,coeffa,pxx,f0,fs,bina,binp,filtername,filtorder,Num_bin)
% filterview -  Displays original psd, filtered psd and filter
%y is amp signal
%z is phase signal
switch filtername
    case 'new'
        EEG=struct('data',x,'srate',fs,'pnts',length(x));
        [y,~,~]=pop_eegfiltnew(EEG,bina(1),bina(2),filtorder);
        [z,~,~]=pop_eegfiltnew(EEG,binp(1),binp(2),filtorder);
    case 'old'
        [y,~]=eegfilt(x,fs,bina,0,filtorder);
        [z,~]=eegfilt(x,fs,binp,0,filtorder);
    case 'kaiser'
        [n,Wn,beta,ftype] = kaiserord([binp(1)-(min(max((binp(1))/4, 2),(binp(1)))),(binp(1)), (binp(2)), (binp(2))+(min(max((binp(2))/4, 2),(binp(1))))],[0 1 0],[.01 .05 .01],fs);
        n = n + rem(n,2);
        y=kfilt(x,coeffb);
        z=kfilt(x,fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale'));
    case 'ellip'
        [b,a]=ellip(filtorder,.05,50,[(binp(1)) (binp(2))]/(fs/2));
        y=filter(coeffb,coeffa,x);
        z=filter(b,a,x);
    case 'ellip plus'
        [n,Wp] = ellipord([(binp(1)) (binp(2))]/(fs/2),[(binp(1))-(min(max((binp(1))/4, 2),(binp(1)))-.00001),(binp(2))+(min(max((binp(1))/4, 2),binp(1))-.00001)]/(fs/2),.05,50);
        [b,a] = ellip(n,.05,50,Wp);
        y=filter(coeffb,coeffa,x);
        z=filter(b,a,x);
    case 'firls'
        b=firls(filtorder,[0 (1-.15)*(binp(1))/(fs/2) (binp(1))/(fs/2) (binp(2))/(fs/2) (1+.15)*(binp(2))/(fs/2) 1],[0 0 1 1 0 0]);
        y=filtfilt(coeffb,coeffa,x);
        z=filtfilt(b,1,x);
    case 'firl'
        b=fir1(filtorder,[(binp(1)) (binp(2))]/(fs/2));
        y=filtfilt(coeffb,coeffa,x);
        z=filtfilt(b,1,x);
    case 'cheby'
        [b,a]=cheby1(filtorder,.05,[(binp(1)) (binp(2))]/(fs/2));
        y=filtfilt(coeffb,coeffa,x);
        z=filtfilt(b,a,x);
    case 'Butter w/ variable bins'
        [b,a]=butter(filtorder,[binp(1) (binp(2))]/(fs/2));
        y=filtfilt(coeffb,coeffa,x);
        z=filtfilt(b,a,x);
end
[pxy,f]=pwelch(y,[],[],[],fs);
[h,f1]=freqz(coeffb,coeffa,length(f),fs);
figure;
subplot(1,2,1)
plot(f0,pxx/max(pxx),'--');
hold on
plot(f,pxy/max(pxx));
plot(f1,abs(h));
xlim([0 max(200,bina(2)+50)]);
title(['Amplitude',' : ',num2str(bina(1)),' - ',num2str(bina(2)), ' Hz : ',filtername]);
subplot(1,2,2)
theta=angle(hilbert(z));
Amp=abs(hilbert(y));
MeanAp=zeros(1,Num_bin);
bin_length=2*pi/Num_bin;
for kp=1:Num_bin
    in=find(theta>=(kp-1)*bin_length-pi & theta<kp*bin_length-pi);
    if ~isempty(in)
        MeanAp(kp)=mean(Amp(in));
    else
        MeanAp(kp)=0;
    end
end
MI=sum(MeanAp/sum(MeanAp).*(log(MeanAp/sum(MeanAp))/log(Num_bin)))+1;
bin_lengthD=360/Num_bin;
bar(0:bin_lengthD:720-bin_lengthD,[MeanAp(1,:) MeanAp(1,:)]/sum(MeanAp(1,:)))
xlim([-10 720])
title({[num2str(binp(1)), ' - ' num2str(binp(2)) ,' Hz Phase vs ' ,num2str(bina(1)), ' - ' num2str(bina(2)),' Hz Amplitude'];
    ['MI = ', num2str(MI)];
    ['(Filter: ' filtername,')']})
xlabel('Phase')
ylabel('Amplitude')
end
function Y=Wavletfilter(x,fs,bins,varargin)
opt_phase=0;
opt_amp=0;
opt_log=0;
wtime=4; %8 for high accuracy
cycles=[3,10];
numvarargs=length(varargin);k=1;
if numvarargs
    while k<=numvarargs
        switch varargin{k}
            case 'phase'
                opt_phase=1;
            case 'amplitude'
                opt_amp=1;
            case 'time'
                wtime=varargin{k+1};
                k=k+1;
            case 'log'
                opt_log=1;
            case 'cycles'
                cycles=varargin{k+1};
                k=k+1;
        end
        k=k+1;
    end
end

num_frex=length(bins);

% define wavelet parameters
time=-wtime:1/fs:wtime;
if length(cycles)>1
    if opt_log
        s=logspace(log10(cycles(1)),log10(cycles(2)),num_frex)./(2*pi*bins);
    else
        s=linspace(cycles(1),cycles(2),num_frex)./(2*pi*bins);
    end
else
    s=cycles.*ones(1,num_frex)./(2*pi*bins);
end

% definte convolution parameters
n_wavelet=length(time);
n_data=length(x);
n_total=n_wavelet+n_data-1;
n_conv_pow2=pow2(nextpow2(n_total));
half_of_wavelet_size=(n_wavelet-1)/2;

% get FFT of data
eegfft = fft(x,n_conv_pow2);

Z=zeros(num_frex,n_data);
% loop through frequencies and compute synchronization
for fi=1:num_frex
    wavelet = fft(sqrt(1/(s(fi)*sqrt(pi)))*exp(2*1i*pi*bins(fi).*time).*exp(-time.^2./(2*(s(fi)^2))),n_conv_pow2);
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_total);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    Z(fi,:)=eegconv;
end
if opt_phase
    Y=angle(Z);
elseif opt_amp
    Y=abs(Z);
else
    Y=Z;
end
end
function PAC3D(x,theta,Amp,fs,m_sec,dt,param,Num_bin,plotnum)
if nargin<9
    plotnum=0;
end
if nargin<8
    Num_bin=18;
end
bin_length=2*pi/Num_bin;
frames=floor((length(x)-m_sec*fs-1)/(dt*fs));
MI3=zeros(size(Amp,1),size(theta,1),frames);
sA=size(Amp,1);
sP=size(theta,1);
if size(sA,3)>1
    for k=1:frames
        btemp=zeros(sA,sP);
        kamp=Amp(:,:,(1+(k-1)*(dt*fs)):(1+(k-1)*(dt*fs)+m_sec*fs));
        ktheta=theta(:,(1+(k-1)*(dt*fs)):(1+(k-1)*(dt*fs)+m_sec*fs));
        parfor i=1:sP
            atemp=zeros(sA,1);
            for j=1:sA
                MeanAp=zeros(1,Num_bin);
                if ~isnan(kamp(j,i,:))
                    for kp=1:Num_bin
                        in=find(ktheta(i,:)>=(kp-1)*bin_length-pi & ktheta(i,:)<kp*bin_length-pi);
                        if ~isempty(in)
                            MeanAp(kp)=mean(kamp(j,i,in));
                        else
                            MeanAp(kp)=0;
                        end
                    end
                end
                atemp(j)=sum(MeanAp/sum(MeanAp).*(log(MeanAp/sum(MeanAp))/log(Num_bin)))+1;
                if isnan(atemp(j))
                    MI_inf=MeanAp/sum(MeanAp).*(log(MeanAp/sum(MeanAp))/log(Num_bin)+1);
                    MI_inf(isnan(MI_inf))=0;
                    atemp(j)=sum(MI_inf);
                end
            end
            btemp(:,i)=atemp;
        end
        MI3(:,:,k)=btemp;
    end
else
    for k=1:frames
        btemp=zeros(sA,sP);
        kamp=Amp(:,(1+(k-1)*(dt*fs)):(1+(k-1)*(dt*fs)+m_sec*fs));
        ktheta=theta(:,(1+(k-1)*(dt*fs)):(1+(k-1)*(dt*fs)+m_sec*fs));
        parfor i=1:sP
            atemp=zeros(sA,1);
            for j=1:sA
                MeanAp=zeros(1,Num_bin);
                if ~isnan(kamp(j,:))
                    for kp=1:Num_bin
                        in=find(ktheta(i,:)>=(kp-1)*bin_length-pi & ktheta(i,:)<kp*bin_length-pi);
                        if ~isempty(in)
                            MeanAp(kp)=mean(kamp(j,in));
                        else
                            MeanAp(kp)=0;
                        end
                    end
                end
                atemp(j)=sum(MeanAp/sum(MeanAp).*(log(MeanAp/sum(MeanAp))/log(Num_bin)))+1;
                if isnan(atemp(j))
                    MI_inf=MeanAp/sum(MeanAp).*(log(MeanAp/sum(MeanAp))/log(Num_bin)+1);
                    MI_inf(isnan(MI_inf))=0;
                    atemp(j)=sum(MI_inf);
                end
            end
            btemp(:,i)=atemp;
        end
        MI3(:,:,k)=btemp;
    end
end
time_array=m_sec:dt:m_sec+(frames-1)*dt;
if plotnum==0 || plotnum==1
    figure;
    contourslice(param(1):param(3):param(2),time_array,param(4):param(6):param(5),permute(MI3,[3,2,1]),[],time_array,[]);
    alpha(.9);
    xlabel('Phase'); ylabel('Time (s)'); zlabel('Amplitude');
    view(34,39);
    caxis([0 max(max(max(MI3)))])
end
if plotnum==0 || plotnum==2
    figure;
    isosurface(param(1):param(3):param(2),time_array,param(4):param(6):param(5),permute(MI3,[3,2,1]));
    alpha(.9);
    xlabel('Phase'); ylabel('Time (s)'); zlabel('Amplitude');
    view(34,39);
    caxis([0 max(max(max(MI3)))])
end
end
function [EEG, com, b] = pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder, revfilt, usefft, plotfreqz)
% pop_eegfiltnew() - Filter data using Hamming windowed sinc FIR filter
%
% Usage:
%   >> [EEG, com, b] = pop_eegfiltnew(EEG); % pop-up window mode
%   >> [EEG, com, b] = pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder,
%                                     revfilt, usefft, plotfreqz);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   locutoff  - lower edge of the frequency pass band (Hz)
%               {[]/0 -> lowpass}
%   hicutoff  - higher edge of the frequency pass band (Hz)
%               {[]/0 -> highpass}
%
% Optional inputs:
%   filtorder - filter order (filter length - 1). Mandatory even
%   revfilt   - [0|1] invert filter (from bandpass to notch filter)
%               {default 0 (bandpass)}
%   usefft    - ignored (backward compatibility only)
%   plotfreqz - [0|1] plot filter's frequency and phase response
%               {default 0}
%
% Outputs:
%   EEG       - filtered EEGLAB EEG structure
%   com       - history string
%   b         - filter coefficients
%
% Note:
%   pop_eegfiltnew is intended as a replacement for the deprecated
%   pop_eegfilt function. Required filter order/transition band width is
%   estimated with the following heuristic in default mode: transition band
%   width is 25% of the lower passband edge, but not lower than 2 Hz, where
%   possible (for bandpass, highpass, and bandstop) and distance from
%   passband edge to critical frequency (DC, Nyquist) otherwise. Window
%   type is hardcoded to Hamming. Migration to windowed sinc FIR filters
%   (pop_firws) is recommended. pop_firws allows user defined window type
%   and estimation of filter order by user defined transition band width.
%
% Author: Andreas Widmann, University of Leipzig, 2012
%
% See also:
%   firfilt, firws, windows

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2008 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
com = '';

if nargin < 1
    help pop_eegfiltnew;
    return
end
if isempty(EEG.data)
    error('Cannot filter empty dataset.');
end

% GUI
if nargin < 2
    
    geometry = {[3, 1], [3, 1], [3, 1], 1, 1, 1};
    geomvert = [1 1 1 2 1 1];
    
    uilist = {{'style', 'text', 'string', 'Lower edge of the frequency pass band (Hz)'} ...
        {'style', 'edit', 'string', ''} ...
        {'style', 'text', 'string', 'Higher edge of the frequency pass band (Hz)'} ...
        {'style', 'edit', 'string', ''} ...
        {'style', 'text', 'string', 'FIR Filter order (Mandatory even. Default is automatic*)'} ...
        {'style', 'edit', 'string', ''} ...
        {'style', 'text', 'string', {'*See help text for a description of the default filter order heuristic.', 'Manual definition is recommended.'}} ...
        {'style', 'checkbox', 'string', 'Notch filter the data instead of pass band', 'value', 0} ...
        {'style', 'checkbox', 'string', 'Plot frequency response', 'value', 1}};
    
    result = inputgui('geometry', geometry, 'geomvert', geomvert, 'uilist', uilist, 'title', 'Filter the data -- pop_eegfiltnew()', 'helpcom', 'pophelp(''pop_eegfiltnew'')');
    
    if isempty(result), return; end
    
    locutoff = str2num(result{1});
    hicutoff = str2num(result{2});
    filtorder = str2num(result{3});
    revfilt = result{4};
    plotfreqz = result{5};
    usefft = [];
    
else
    
    if nargin < 3
        hicutoff = [];
    end
    if nargin < 4
        filtorder = [];
    end
    if nargin < 5 || isempty(revfilt)
        revfilt = 0;
    end
    if nargin < 6
        usefft = [];
    elseif usefft == 1
        error('FFT filtering not supported. Argument is provided for backward compatibility in command line mode only.')
    end
    if nargin < 7
        plotfreqz = 0;
    end
    
end

% Constants
TRANSWIDTHRATIO = 0.25;
fNyquist = EEG.srate / 2;

% Check arguments
if locutoff == 0, locutoff = []; end
if hicutoff == 0, hicutoff = []; end
if isempty(hicutoff) % Convert highpass to inverted lowpass
    hicutoff = locutoff;
    locutoff = [];
    revfilt = ~revfilt;
end
edgeArray = sort([locutoff hicutoff]);

if isempty(edgeArray)
    error('Not enough input arguments.');
end
if any(edgeArray < 0 | edgeArray >= fNyquist)
    error('Cutoff frequency out of range');
end

if ~isempty(filtorder) && (filtorder < 2 || mod(filtorder, 2) ~= 0)
    error('Filter order must be a real, even, positive integer.')
end

% Max stop-band width
maxTBWArray = edgeArray; % Band-/highpass
if revfilt == 0 % Band-/lowpass
    maxTBWArray(end) = fNyquist - edgeArray(end);
elseif length(edgeArray) == 2 % Bandstop
    maxTBWArray = diff(edgeArray) / 2;
end
maxDf = min(maxTBWArray);

% Transition band width and filter order
if isempty(filtorder)
    
    % Default filter order heuristic
    if revfilt == 1 % Highpass and bandstop
        df = min([max([maxDf * TRANSWIDTHRATIO 2]) maxDf]);
    else % Lowpass and bandpass
        df = min([max([edgeArray(1) * TRANSWIDTHRATIO 2]) maxDf]);
    end
    
    filtorder = 3.3 / (df / EEG.srate); % Hamming window
    filtorder = ceil(filtorder / 2) * 2; % Filter order must be even.
    %filtorder = 1650;      %For constant order
else
    
    df = 3.3 / filtorder * EEG.srate; % Hamming window
    filtorderMin = ceil(3.3 ./ ((maxDf * 2) / EEG.srate) / 2) * 2;
    filtorderOpt = ceil(3.3 ./ (maxDf / EEG.srate) / 2) * 2;
    if filtorder < filtorderMin
        error('Filter order too low. Minimum required filter order is %d. For better results a minimum filter order of %d is recommended.', filtorderMin, filtorderOpt)
    elseif filtorder < filtorderOpt
        warning('firfilt:filterOrderLow', 'Transition band is wider than maximum stop-band width. For better results a minimum filter order of %d is recommended. Reported might deviate from effective -6dB cutoff frequency.', filtorderOpt)
    end
    
end

filterTypeArray = {'lowpass', 'bandpass'; 'highpass', 'bandstop (notch)'};
% fprintf('pop_eegfiltnew() - performing %d point %s filtering.\n', filtorder + 1, filterTypeArray{revfilt + 1, length(edgeArray)})
% fprintf('pop_eegfiltnew() - transition band width: %.4g Hz\n', df)
% fprintf('pop_eegfiltnew() - passband edge(s): %s Hz\n', mat2str(edgeArray))

% Passband edge to cutoff (transition band center; -6 dB)
dfArray = {df, [-df, df]; -df, [df, -df]};
cutoffArray = edgeArray + dfArray{revfilt + 1, length(edgeArray)} / 2;
% fprintf('pop_eegfiltnew() - cutoff frequency(ies) (-6 dB): %s Hz\n', mat2str(cutoffArray))

% Window
winArray = windows('hamming', filtorder + 1);

% Filter coefficients
if revfilt == 1
    filterTypeArray = {'high', 'stop'};
    b = firws(filtorder, cutoffArray / fNyquist, filterTypeArray{length(cutoffArray)}, winArray);
else
    b = firws(filtorder, cutoffArray / fNyquist, winArray);
end

% Plot frequency response
if plotfreqz
    freqz(b, 1, 8192, EEG.srate);
end

% Filter
% disp('pop_eegfiltnew() - filtering the data');
EEG = firfilt(EEG, b);

% History string
com = sprintf('%s = pop_eegfiltnew(%s, %s, %s, %s, %s, %s, %s);', inputname(1), inputname(1), mat2str(locutoff), mat2str(hicutoff), mat2str(filtorder), mat2str(revfilt), mat2str(usefft), mat2str(plotfreqz));


% firfilt() - Pad data with DC constant, filter data with FIR filter,
%             and shift data by the filter's group delay
%
% Usage:
%   >> EEG = firfilt(EEG, b, nFrames);
%
% Inputs:
%   EEG           - EEGLAB EEG structure
%   b             - vector of filter coefficients
%
% Optional inputs:
%   nFrames       - number of frames to filter per block {default 1000}
%
% Outputs:
%   EEG           - EEGLAB EEG structure
%
% Note:
%   Higher values for nFrames increase speed and working memory
%   requirements.
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   filter, findboundaries

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    function y = firfilt(EEG, b, nFrames)
        
        if nargin < 2
            error('Not enough input arguments.');
        end
        if nargin < 3 || isempty(nFrames)
            nFrames = 1000;
        end
        
        % Filter's group delay
        if mod(length(b), 2) ~= 1
            error('Filter order is not even.');
        end
        groupDelay = (length(b) - 1) / 2;
        
        % Find data discontinuities and reshape epoched data
        % if EEG.trials > 1 % Epoched data
        %     EEG.data = reshape(EEG.data, [EEG.nbchan EEG.pnts * EEG.trials]);
        %     dcArray = 1 : EEG.pnts : EEG.pnts * (EEG.trials + 1);
        % else % Continuous data
        %     dcArray = [findboundaries(EEG.event) EEG.pnts + 1];
        % end
        
        dcArray = [1 EEG.pnts + 1];
        % Initialize progress indicator
        nSteps = 20;
        step = 0;
        % fprintf(1, 'firfilt(): |');
        % strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
        %tic
        
        for iDc = 1:(length(dcArray) - 1)
            
            % Pad beginning of data with DC constant and get initial conditions
            ziDataDur = min(groupDelay, dcArray(iDc + 1) - dcArray(iDc));
            [temp, zi] = filter(b, 1, double([EEG.data(:, ones(1, groupDelay) * dcArray(iDc)) ...
                EEG.data(:, dcArray(iDc):(dcArray(iDc) + ziDataDur - 1))]), [], 2);
            
            blockArray = [(dcArray(iDc) + groupDelay):nFrames:(dcArray(iDc + 1) - 1) dcArray(iDc + 1)];
            for iBlock = 1:(length(blockArray) - 1)
                
                % Filter the data
                [EEG.data(:, (blockArray(iBlock) - groupDelay):(blockArray(iBlock + 1) - groupDelay - 1)), zi] = ...
                    filter(b, 1, double(EEG.data(:, blockArray(iBlock):(blockArray(iBlock + 1) - 1))), zi, 2);
                
                % Update progress indicator
                %             [step, strLength] = mywaitbar((blockArray(iBlock + 1) - groupDelay - 1), size(EEG.data, 2), step, nSteps, strLength);
            end
            
            % Pad end of data with DC constant
            temp = filter(b, 1, double(EEG.data(:, ones(1, groupDelay) * (dcArray(iDc + 1) - 1))), zi, 2);
            EEG.data(:, (dcArray(iDc + 1) - ziDataDur):(dcArray(iDc + 1) - 1)) = ...
                temp(:, (end - ziDataDur + 1):end);
            
            % Update progress indicator
            %         [step, strLength] = mywaitbar((dcArray(iDc + 1) - 1), size(EEG.data, 2), step, nSteps, strLength);
            
        end
        
        % Reshape epoched data
        % if EEG.trials > 1
        %     EEG.data = reshape(EEG.data, [EEG.nbchan EEG.pnts EEG.trials]);
        % end
        
        % Deinitialize progress indicator
        % fprintf(1, '\n')
        y=EEG.data;
    end

% function [step, strLength] = mywaitbar(compl, total, step, nSteps, strLength)
%
% progStrArray = '/-\|';
% tmp = floor(compl / total * nSteps);
% if tmp > step
%     fprintf(1, [repmat('\b', 1, strLength) '%s'], repmat('=', 1, tmp - step))
%     step = tmp;
%     ete = ceil(toc / step * (nSteps - step));
%     strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '%s %3d%%, ETE %02d:%02d'], progStrArray(mod(step - 1, 4) + 1), floor(step * 100 / nSteps), floor(ete / 60), mod(ete, 60));
% end
%
% end

%firws() - Designs windowed sinc type I linear phase FIR filter
%
% Usage:
%   >> b = firws(m, f);
%   >> b = firws(m, f, w);
%   >> b = firws(m, f, t);
%   >> b = firws(m, f, t, w);
%
% Inputs:
%   m - filter order (mandatory even)
%   f - vector or scalar of cutoff frequency/ies (-6 dB;
%       pi rad / sample)
%
% Optional inputs:
%   w - vector of length m + 1 defining window {default blackman}
%   t - 'high' for highpass, 'stop' for bandstop filter {default low-/
%       bandpass}
%
% Output:
%   b - filter coefficients
%
% Example:
%   fs = 500; cutoff = 0.5; tbw = 1;
%   m  = pop_firwsord('hamming', fs, tbw);
%   b  = firws(m, cutoff / (fs / 2), 'high', windows('hamming', m + 1));
%
% References:
%   Smith, S. W. (1999). The scientist and engineer's guide to digital
%   signal processing (2nd ed.). San Diego, CA: California Technical
%   Publishing.
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   pop_firws, pop_firwsord, pop_kaiserbeta, windows

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    function [b a] = firws(m, f, t, w)
        
        a = 1;
        
        if nargin < 2
            error('Not enough input arguments');
        end
        if length(m) > 1 || ~isnumeric(m) || ~isreal(m) || mod(m, 2) ~= 0 || m < 2
            error('Filter order must be a real, even, positive integer.');
        end
        f = f / 2;
        if any(f <= 0) || any(f >= 0.5)
            error('Frequencies must fall in range between 0 and 1.');
        end
        if nargin < 3 || isempty(t)
            t = '';
        end
        if nargin < 4 || isempty(w)
            if ~isempty(t) && ~ischar(t)
                w = t;
                t = '';
            else
                w = windows('blackman', (m + 1));
            end
        end
        w = w(:)'; % Make window row vector
        
        b = fkernel(m, f(1), w);
        
        if length(f) == 1 && strcmpi(t, 'high')
            b = fspecinv(b);
        end
        
        if length(f) == 2
            b = b + fspecinv(fkernel(m, f(2), w));
            if isempty(t) || ~strcmpi(t, 'stop')
                b = fspecinv(b);
            end
        end
        
        % Compute filter kernel
        function b = fkernel(m, f, w)
            m = -m / 2 : m / 2;
            b(m == 0) = 2 * pi * f; % No division by zero
            b(m ~= 0) = sin(2 * pi * f * m(m ~= 0)) ./ m(m ~= 0); % Sinc
            b = b .* w; % Window
            b = b / sum(b); % Normalization to unity gain at DC
        end
        
        % Spectral inversion
        function b = fspecinv(b)
            b = -b;
            b(1, (length(b) - 1) / 2 + 1) = b(1, (length(b) - 1) / 2 + 1) + 1;
        end
    end
% windows() - Returns handle to window function or window
%
% Usage:
%   >> h = windows(t);
%   >> h = windows(t, m);
%   >> h = windows(t, m, a);
%
% Inputs:
%   t - char array 'rectangular', 'bartlett', 'hann', 'hamming',
%       'blackman', 'blackmanharris', or 'kaiser'
%
% Optional inputs:
%   m - scalar window length
%   a - scalar or vector with window parameter(s)
%
% Output:
%   h - function handle or column vector window
%
% Author: Andreas Widmann, University of Leipzig, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    function h = windows(t, m, a)
        
        if nargin < 1
            error('Not enough input arguments.');
        end
        h = str2func(t);
        switch nargin
            case 2
                h = h(m);
            case 3
                h = h(m, a);
        end
    end

    function w = rectangular(m)
        w = ones(m, 1);
    end

    function w = bartlett(m)
        w = 1 - abs(-1:2 / (m - 1):1)';
    end

% von Hann
    function w = hann(m);
        w = hamming(m, 0.5);
    end

% Hamming
    function w = hamming(m, a)
        if nargin < 2 || isempty(a)
            a = 25 / 46;
        end
        m = [0:1 / (m - 1):1]';
        w = a - (1 - a) * cos(2 * pi * m);
    end

% Blackman
    function w = blackman(m, a)
        if nargin < 2 || isempty(a)
            a = [0.42 0.5 0.08 0];
        end
        m = [0:1 / (m - 1):1]';
        w = a(1) - a(2) * cos (2 * pi * m) + a(3) * cos(4 * pi * m) - a(4) * cos(6 * pi * m);
    end

% Blackman-Harris
    function w = blackmanharris(m)
        w = blackman(m, [0.35875 0.48829 0.14128 0.01168]);
    end

% Kaiser
    function w = kaiser(m, a)
        if nargin < 2 || isempty(a)
            a = 0.5;
        end
        m = [-1:2 / (m - 1):1]';
        w = besseli(0, a * sqrt(1 - m.^2)) / besseli(0, a);
    end

% Tukey
    function w = tukey(m, a)
        if nargin < 2 || isempty(a)
            a = 0.5;
        end
        if a <= 0
            w = ones(m, 1);
        elseif a >= 1
            w = hann(m);
        else
            a = (m - 1) / 2 * a;
            tapArray = (0:a)' / a;
            w = [0.5 - 0.5 * cos(pi * tapArray); ...
                ones(m - 2 * length(tapArray), 1); ...
                0.5 - 0.5 * cos(pi * tapArray(end:-1:1))];
        end
    end
end