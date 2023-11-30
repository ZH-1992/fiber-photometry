clear;close;clc;
cdir=dir('*.mat');
matname={cdir.name};
seginital=strfind(matname,'a000');
fn=0;
dur=0.5;
step_i=50;
for n=1:length(seginital)
    if ~isempty(seginital{n})
        fn=fn+1;
        filest(fn)=n;
    end
end
fileet=[filest(2:end)-1 length(matname)];    %为什么filest(2:end)后面要减去1？
clear cdir fn n seginital;

for m=1:length(filest)
    sn=2;     %record the output position     %为什么sn的初始值设为2？因为抬头一行用来写标题。
    num=0;    %count the number of the sweep
    thres=20; %this value is the dv/dt threshold
    for n=filest(m):fileet(m)
        num=num+1;
        load(matname{n});
        varia=[];v=[];dv=[];
        varia=who('*Ch*');
        v=eval([varia{1} '.values']);
        dt=eval([varia{1} '.interval']);
        clear(varia{:});
        dv=(diff(v(1:end-1))+diff(v(2:end)))/(2000*dt);
        v(1)=[];v(end)=[];
        row=[];bj=[];lbj=[];rbj=[];
        row=find(v>0);                            %%限定action potential的peak在0 mV以上。
        if isempty(row)
            continue;
        end
        row=[];
        row=find(dv>=thres);
        if isempty(row)
            continue;
        end
        bj=find(diff(row)>10);      %%??????????????????????? As sampling frequency is 100 kHz. 防止一个AP找到多个顶点
        lbj=[row(1);row(bj+1)];
        rbj=[row(bj);row(end)];
        
        
        for k=1:length(lbj)
            row=[];
            row=find(v(lbj(k):rbj(k))>0);         %%限定action potential的peak在0 mV以上。
            if isempty(row)
                lbj(k)=0;
            end
        end
        rbj(find(lbj==0))=[];
        lbj(find(lbj==0))=[];
        %find the max dv/dt and its position
        maxdv=zeros(length(lbj),1);
        maxdvpos=zeros(length(lbj),1);
        for k=1:length(lbj)
            [maxdv(k),maxdvpos(k)]=max(dv(lbj(k):rbj(k)));
        end
        maxdvpos=lbj+maxdvpos-1;
        
        %if the distance between the neighbor-spike less than 0.001 second
        %we will delete the smaller one
        erow=[];
        erow=find(diff(lbj)<fix(0.001/dt));
        for k=1:length(erow)
            if maxdv(erow(k))<maxdv(erow(k)+1)
                ;
            else
                erow(k)=erow(k)+1;
            end
        end
        lbj(erow)=[];                  
        rbj(erow)=[];
        maxdv(erow)=[];
        maxdvpos(erow)=[];
        
        %find the min dv/dt and it's position
        dlbj=[];drbj=[];
        if length(rbj)==1
            dlbj=rbj;
            drbj=rbj+fix(0.001/dt);
        else
            dlbj=rbj;
            if rbj(end)+lbj(end)-rbj(end-1)<size(dv,1)
                drbj=[lbj(2:end);rbj(end)+lbj(end)-rbj(end-1)];
            else
                drbj=[lbj(2:end);size(dv,1)];
            end
        end
        mindv=zeros(length(dlbj),1);
        mindvpos=zeros(length(dlbj),1);
        for k=1:length(dlbj)
            [mindv(k),mindvpos(k)]=min(dv(dlbj(k):drbj(k)));
        end
        mindvpos=dlbj+mindvpos-1;
        
        %find the value and position of every peak
        threshold=zeros(length(lbj),1);
        threshold=v(lbj);
        peak=zeros(length(maxdvpos),1);
        peakpos=zeros(length(maxdvpos),1);
        for k=1:length(maxdvpos)
            [peak(k),peakpos(k)]=max(v(maxdvpos(k):mindvpos(k)));
        end
        peakpos=maxdvpos+peakpos-1;
        
        %calculate the amplitude
        amp=[];
        amp=peak-threshold;
        %calculate the ahp
        ahp=zeros(length(peakpos),1);
        for k=1:length(peakpos)-1
            %%ahp(k)=min(v(peakpos(k):peakpos(k)+fix(0.002/dt))); %%%20110329 限定AHP的寻找范围是peakpos后的2秒， 这和PC的情况不符合， 所以进行修改，改为两个spike之间的最低点。
            ahp(k)=min(v(peakpos(k):peakpos(k+1)));
        end
        ahp=ahp-threshold(1:end);
        ahp(end)=0;
        %calculate the half-width
        halfheight=[];
        halfheight=threshold+amp/2;
        halflbj=zeros(length(halfheight),1);
        halfrbj=zeros(length(halfheight),1);
        for k=1:length(halfheight)
            lnrow=[];lnbj=[];lnlbj=[];lnrbj=[];
            lnrow=find(v>=halfheight(k));
            lnbj=find(diff(lnrow)>1);
            lnlbj=[lnrow(1);lnrow(lnbj+1)];
            lnrbj=[lnrow(lnbj);lnrow(end)];
            for k1=1:length(lnlbj)
                if lnlbj(k1)<peakpos(k)&lnrbj(k1)>peakpos(k)
                    halflbj(k)=lnlbj(k1);
                    halfrbj(k)=lnrbj(k1);
                    break;
                end
            end
        end
        halfwide=[];
        halfwide=(halfrbj-halflbj)*dt*1000;
        %find the end of every spike
        npheight=[];
        npheight=threshold+0.1*amp;
        nprbj=zeros(length(npheight),1);
        for k=1:length(npheight)
            lnrow=[];lnbj=[];lnlbj=[];lnrbj=[];
            lnrow=find(v>=npheight(k));
            lnbj=find(diff(lnrow)>1);
            lnlbj=[lnrow(1);lnrow(lnbj+1)];
            lnrbj=[lnrow(lnbj);lnrow(end)];
            for k1=1:length(lnlbj)
                if lnlbj(k1)<peakpos(k)&lnrbj(k1)>peakpos(k)
                    nprbj(k)=lnrbj(k1);
                    break;
                end
            end
        end
        %calculate the curve area between the threshold position and the end of the spike
        area=zeros(length(nprbj),1);   
        for k=1:length(nprbj)
            area(k)=polyarea((lbj(k):nprbj(k))'*dt*1000,v(lbj(k):nprbj(k)));
        end
        %calculate inter-spike intervals
        IS_interval=zeros(length(peakpos),1); 
        if length(peakpos)>1
          for j=1:length(peakpos)-1;
          IS_interval(j)=(peakpos(j+1)-peakpos(j))*1000*dt;
          end
          IS_interval(length(peakpos))=0;
        end
        %output the result
        res=[];
        res=[ones(length(area),1)*num peakpos*1000*dt IS_interval threshold amp peak halfwide ...
            maxdv mindv area ahp];
        res(find(peak<0),:)=[];                    %%限定action potential的peak在0 mV以上。
        res_interval=[ones(size(peakpos,1)-1,1)*num IS_interval(1:end-1,1)];
        % mean inter-spike interval for each sweep
        spike_num=length(peakpos);
        step_duration=dur; % the unit was second.  Current injection has a duration of 500 ms.
        spike_freq=spike_num/step_duration;
        mean_interval=[num mean(IS_interval(1:end-1,1)) spike_num spike_freq (num-1)*step_i+0];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        warning off matlab:xlswrite:addsheet;
        xlswrite('res.xls',res,matname{n}(1:end-8),['a' num2str(sn)]);
        if ~isempty(res_interval)
        xlswrite('res_ISI.xls',res_interval,matname{n}(1:end-8),['a' num2str(sn)]);
        end
        xlswrite(strcat('mean_ISI_',matname{n}(1:end-8)),mean_interval,matname{n}(1:end-8),['a' num2str(num)]);
        sn=sn+size(res,1);
        
    end
    xlswrite('res.xls',{'sweepnum','peakpos','ISI', 'threshold','amp','peak',...
        'halfwide','导数最大值','导数最小值','面积','ahp'},matname{n}(1:end-8));
    xlswrite('res_ISI.xls',{'sweepnum','inter-spike interval'},matname{n}(1:end-8));
    xlswrite(strcat('mean_ISI_',matname{n}(1:end-8)),{'sweepnum','mean ISI for each sweep','spike_num','spike_freq (Hz)','Injected current (pA)'},matname{n}(1:end-8));
        % mean value of threshold,amp,peak,halfwidth,面积,ahp
    NUMERIC=xlsread('res',matname{n}(1:end-8));
    mean_values=[mean(NUMERIC(:,4:end-1)) mean(NUMERIC(find(NUMERIC(:,11)<0),11)) mean(NUMERIC(find(NUMERIC(:,3)>0),3))];  % mean(NUMERIC(find(NUMERIC(:,3)>0),3))用于计算ISI， mean(NUMERIC(find(NUMERIC(:,11)<0),11))用于计算AHP。
    xlswrite(strcat('Mean_',matname{n}(1:end-8)),mean_values,matname{n}(1:end-8),'a2:i2');
    xlswrite(strcat('Mean_',matname{n}(1:end-8)),{'threshold','amp','peak',...
        'halfwidth','导数最大值','导数最小值','面积','ahp','mean_ISI','Current threhold','Cm','Rm','Ra','Tau','RMP','1st spike Vthresh','sustained spike Vthresh'},matname{n}(1:end-8));
    
    %%%%%%  只对每一个sweep的第一个spike的参数进行统计平均
    spike_1st=xlsread('res',matname{n}(1:end-8));
    spike_1st_num(1)=1;
    for i=1:size(spike_1st,1)
        if spike_1st(1,1)+i<=spike_1st(end,1)
        spike_1st_num(i+1)=min(find(spike_1st==spike_1st(1,1)+i));
        else
        break
        end
    end
    spike_1st=spike_1st(spike_1st_num,:);
    mean_1st=[mean(spike_1st(:,4:end-1))];
    xlswrite(strcat('Spike_1st_',matname{n}(1:end-8)),spike_1st,matname{n}(1:end-8),'a2');
    xlswrite(strcat('Spike_1st_',matname{n}(1:end-8)),mean_1st,matname{n}(1:end-8),'d30');
    xlswrite(strcat('Spike_1st_',matname{n}(1:end-8)),{'sweep num','peakpos',...
        'ISI','threshold','amplitude (mV)','peak','half width','导数最大值','导数最小值','area','ahp','Current threhold','Cm','Rm','Tau','RMP'},matname{n}(1:end-8));
    
    xlswrite(strcat('iN_1st_',matname{n}(1:end-8)),mean_1st,matname{n}(1:end-8),'a2');
    xlswrite(strcat('iN_1st_',matname{n}(1:end-8)),{'threshold','amplitude (mV)','peak','half width','导数最大值','导数最小值','area','Current threhold','Cm','Rm','Tau','RMP'},matname{n}(1:end-8));
    
end
