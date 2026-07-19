function specific_rates_lin_only(is, nis, ind, gal, ver, tmax_window, Mmax_window, zform_window, tdyn, leg)

% gal is the simulation index: V**
% gal==0 means stack all simulations together

% ver=1 --> M_in and M_out from Rc with Vr>0
% ver=2 --> M_in and M_out from average with Vr>0
% ver=3 --> M_out from Rc with Vr>Vesc and M_in from mass conservation
% ver=4 --> M_out from average with Vr>Vesc and M_in from mass conservation
% ver=5 --> M_out from Frederic's method and M_in from mass conservation

dirname = './specific_rates/';
mkdir(dirname);

Tdown = tmax_window(1);
Tup   = tmax_window(2);

Mdown = Mmax_window(1);
Mup   = Mmax_window(2);

zdown = zform_window(1);
zup   = zform_window(2);

if(gal~=0)
    dirname = strcat(dirname,'V',num2str(gal,'%02i'),'/');
    tit = strcat('$V',num2str(gal,'%02i'),'$');
    mkdir(dirname);
else
    dirname = strcat(dirname,'all/');
    tit = 'all';
    mkdir(dirname);
end

dirname = strcat(dirname,num2str(Mdown,'%3.1f'),'_M50_',num2str(Mup,'%3.1f'),'/');
dirname = strcat(dirname,num2str(zdown,'%3.1f'),'_zform_',num2str(zup,'%3.1f'),'/');
mkdir(dirname);

if(ver==1) % M_in and M_out from Rc with Vr>0
    prop2 = 31; % M_in
    prop3 = 34; % M_out
    dirname = strcat(dirname,'Rc_Vr_0/');
elseif(ver==2) % M_in and M_out from average with Vr>0
    prop2 = 51; % M_in
    prop3 = 52; % M_out
    dirname = strcat(dirname,'average_Vr_0/');
elseif(ver==3) % M_out from Rc with Vr>Vesc and M_in from mass conservation
    prop2 = 75; % M_in
    prop3 = 37; % M_out
    dirname = strcat(dirname,'Rc_Vr_Vesc_gas_cons/');
elseif(ver==4) % M_out from Rc with Vr>0 V>Vesc and M_in from mass conservation
    prop2 = 76; % M_in
    prop3 = 41; % M_out
    dirname = strcat(dirname,'Rc_Vr_0_V_Vesc_gas_cons/');
elseif(ver==5) % M_out from average with Vr>Vesc and M_in from mass conservation
    prop2 = 77; % M_in
    prop3 = 53; % M_out
    dirname = strcat(dirname,'avg_Vr_Vesc_gas_cons/');
elseif(ver==6) % M_out from average with Vr>0 V>Vesc and M_in from mass conservation
    prop2 = 78; % M_in
    prop3 = 54; % M_out
    dirname = strcat(dirname,'avg_Vr_0_V_Vesc_gas_cons/');
elseif(ver==7) % M_out from Frederic's method and M_in from mass conservation
    prop2 = 79; % M_in
    prop3 = 40; % M_out
    dirname = strcat(dirname,'Frederic_gas_cons/');    
elseif(ver==8) % M_out from Rc with Vr>Vesc and M_in from mass conservation
    prop2 = 86; % M_in
    prop3 = 37; % M_out
    dirname = strcat(dirname,'Rc_Vr_Vesc_bar_cons/');
elseif(ver==9) % M_out from Rc with Vr>0 V>Vesc and M_in from mass conservation
    prop2 = 87; % M_in
    prop3 = 41; % M_out
    dirname = strcat(dirname,'Rc_Vr_0_V_Vesc_bar_cons/');
elseif(ver==10) % M_out from average with Vr>Vesc and M_in from mass conservation
    prop2 = 88; % M_in
    prop3 = 53; % M_out
    dirname = strcat(dirname,'avg_Vr_Vesc_bar_cons/');
elseif(ver==11) % M_out from average with Vr>0 V>Vesc and M_in from mass conservation
    prop2 = 89; % M_in
    prop3 = 54; % M_out
    dirname = strcat(dirname,'avg_Vr_0_V_Vesc_bar_cons/');
elseif(ver==12) % M_out from Frederic's method and M_in from mass conservation
    prop2 = 90; % M_in
    prop3 = 40; % M_out
    dirname = strcat(dirname,'Frederic_bar_cons/');
end
mkdir(dirname);

ylab = '${\rm specific\:\:rates\:\:[Gyr^{-1}]}$';
ylab2 = '${\rm specific\:\:\:rates\:\:\:[Gyr^{-1}\:]}$';
yl = 0;
yu = 40;
dy = 5;
ytick = [0.01 0.1 1 10 40];
ystr = {num2str(0.01,'%4.2f'),num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f'),num2str(40,'%2.0f')};

xlab = '$t\:{\rm [Myr]}$';
xlab2 = '$t\:\:{\rm [Myr]}$';
xl = 0;
if(gal==19)
    xu = 200;
    dx = 20;
    xtick = [1 10 100 200];
    xstr = {num2str(1,'%1.0f'),num2str(10,'%2.0f'),num2str(100,'%3.0f'),num2str(200,'%3.0f')};
    leg_pos = [0.6, 0.14, 0.16, 0.12];
    leg_pos2 = [0.58, 0.14, 0.16, 0.12];
else
    xu = 1000;
    dx = 100;
    xtick = [1 10 100 1000];
    xstr = {num2str(1,'%1.0f'),num2str(10,'%2.0f'),num2str(100,'%3.0f'),num2str(1000,'%4.0f')};
    leg_pos = [0.6, 0.74, 0.16, 0.12]
    leg_pos2 = [0.58, 0.77, 0.16, 0.12];
end
tmed = [2*tdyn 6*tdyn];
dt = tdyn/4;
tH = 0.5*tdyn;

tmax = 4000;
tmin = 0;
tbin = (tmin+dt/2):dt:(tmax-dt/2);
Nt = length(tbin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(gal~=0)
    m = find(ind==gal);
    vec = (nis(m)+1):(nis(m+1)-1);
    vec2 = vec+1;
    
    b = find(is(vec,62)>=Tdown & is(vec,62)<Tup & is(vec,1)~=is(vec2,1));
    end_index = vec(b);
    if( is(nis(m+1),62)>=Tdown & is(nis(m+1),62)<Tup )
        end_index = [end_index, nis(m+1)];
    end
    
    b = find(is(vec2,62)>=Tdown & is(vec2,62)<Tup & is(vec2,1)~=is(vec,1));
    start_index = vec2(b);
    if( is(nis(m)+1,62)>=Tdown & is(nis(m)+1,62)<Tup )
        start_index = [nis(m)+1, start_index];
    end
else    
    end_index = [];
    start_index = [];
    for m=1:length(ind)
        vec = (nis(m)+1):(nis(m+1)-1);
        vec2 = vec+1;
        
        b = find(is(vec,62)>=Tdown & is(vec,62)<Tup & is(vec,1)~=is(vec2,1));
        end_index = [end_index,vec(b)];
        if( is(nis(m+1),62)>=Tdown & is(nis(m+1),62)<Tup )
            end_index = [end_index, nis(m+1)];
        end        
        
        if( is(nis(m)+1,62)>=Tdown & is(nis(m)+1,62)<Tup )
            start_index = [nis(m)+1, start_index];
        end
        b = find(is(vec2,62)>=Tdown & is(vec2,62)<Tup & is(vec2,1)~=is(vec,1));
        start_index = [start_index,vec2(b)];
    end
end
b = find(is(start_index,2)>=zdown & is(start_index,2)<zup);
start_index = start_index(b);
end_index = end_index(b);
Mavg = 0.*end_index;
for i=1:length(start_index)
    tvec = start_index(i):end_index(i);
    b = find(abs(is(tvec,48)-50)<=10);
    if(~isempty(b))
        Mavg(i) = mean(log10(is(tvec(b),6)));
    end
%     b = find(abs(is(tvec,48)-50)<=20);
%     if(~isempty(b))
%         Mavg(i) = max(Mavg(i), mean(log10(is(tvec(b),6))));
%     end
%     b = find(abs(is(tvec,48)-50)<=30);
%     if(~isempty(b))
%         Mavg(i) = max(Mavg(i), mean(log10(is(tvec(b),6))));
%     end
end
b = find(Mavg>=Mdown & Mavg<Mup);
start_index = start_index(b);
end_index = end_index(b);
Mavg = Mavg(b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( length(start_index) ~= length(end_index) )
    [0,0,length(start_index), length(end_index)]
    return
elseif( min(end_index-start_index)<=0 )
    [0,0,min(end_index-start_index)]
    return
else
    
    Nclump = length(start_index);
    Nindex = max(end_index-start_index) + 1;
    ydat1  = 1e50.*ones(Nclump,Nindex);
    ydat2  = ydat1;
    ydat3  = ydat1;
    ydat4  = ydat1;
    clump_time = ydat1;
    stack1 = 1e50.*ones(3,Nt+1);
    stack2 = stack1;
    stack3 = stack1;
    stack4 = stack1;
    min(is(end_index(1:Nclump),2))
    [Nclump, Nindex]
    
    % Get individual clump data and smooth with Gaussian
    for i=1:Nclump
        j = start_index(i);
        k = end_index(i);
        for n=j:k
            ydat1(i,n-j+1) = 1e9 .* is(n,15)./is(n,6); % sSFR Gyr^{-1}
            ydat2(i,n-j+1) = 1e9 .* is(n,prop2)./is(n,6); % sMg_in Gyr^{-1}
            ydat3(i,n-j+1) = 1e9 .* is(n,prop3)./is(n,6); % sMg_out Gyr^{-1}
            ydat4(i,n-j+1) = 1e9.*(is(n,45)-is(n,44))./is(n,6); % sM*_out_net Gyr^{-1}            
            clump_time(i,n-j+1) = is(n,48);
        end
        temp1 = ydat1(i,1:(k-j+1));
        temp2 = ydat2(i,1:(k-j+1));
        temp3 = ydat3(i,1:(k-j+1));
        temp4 = ydat4(i,1:(k-j+1));
        xdat = clump_time(i,1:(k-j+1));
        for n=1:(k-j+1)
            b1 = find(abs(xdat-xdat(n))<=4.*tH);  % Gaussian with HWHM=tH out to 4*HWHM
            temp1(n) = sum( ydat1(i,b1).*0.5.^(((xdat(b1)-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(n))./tH).^2) );
            temp2(n) = sum( ydat2(i,b1).*0.5.^(((xdat(b1)-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(n))./tH).^2) );
            temp3(n) = sum( ydat3(i,b1).*0.5.^(((xdat(b1)-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(n))./tH).^2) );
            temp4(n) = sum( ydat4(i,b1).*0.5.^(((xdat(b1)-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(n))./tH).^2) );
        end
        ydat1(i,1:(k-j+1)) = temp1;
        ydat2(i,1:(k-j+1)) = temp2;
        ydat3(i,1:(k-j+1)) = temp3;
        ydat4(i,1:(k-j+1)) = max( 1e-6, temp4 );
    end
    clear temp1 temp2 temp3 temp4 xdat
    
    % Stack smoothed clump data
    for i=1:Nt
        b = find(clump_time>(tbin(i)-dt/2) & clump_time<=(tbin(i)+dt/2));
        temp1 = sort(ydat1(b));
        temp2 = sort(ydat2(b));
        temp3 = sort(ydat3(b));
        temp4 = sort(ydat4(b));
        Ntemp = length(b);
        if(Ntemp>=3)
            Nlow  = floor(Ntemp/6) + 1;
            Nhigh = floor(5*Ntemp/6);
            
            stack1(1,i) = median(temp1);
            stack1(2,i) = median(temp1) - temp1(Nlow);
            stack1(3,i) = temp1(Nhigh) - median(temp1);
            
            stack2(1,i) = median(temp2);
            stack2(2,i) = median(temp2) - temp2(Nlow);
            stack2(3,i) = temp2(Nhigh) - median(temp2);
            
            stack3(1,i) = median(temp3);
            stack3(2,i) = median(temp3) - temp3(Nlow);
            stack3(3,i) = temp3(Nhigh) - median(temp3);
            
            stack4(1,i) = median(temp4);
            stack4(2,i) = median(temp4) - temp4(Nlow);
            stack4(3,i) = temp4(Nhigh) - median(temp4);
        end
    end
    clear temp1 temp2 temp3 temp4 Ntemp Nlow Nhigh
    
    b = find( stack1(1,1:Nt)~=1e50 );
    temp = stack1(1,b);
    xdat = tbin(b);
    for n=1:length(b)
        b1 = find(abs(xdat-xdat(n))<=4*tH);
        temp(n) = sum( stack1(1,b(b1)).*0.5.^(((xdat(b1)-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(n))./tH).^2) );
    end
    stack1(1,b) = temp;
    
    b = find( stack2(1,1:Nt)~=1e50 );
    temp = stack2(1,b);
    xdat = tbin(b);
    for n=1:length(b)
        b1 = find(abs(xdat-xdat(n))<=4*tH);
        temp(n) = sum( stack2(1,b(b1)).*0.5.^(((xdat(b1)-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(n))./tH).^2) );
    end
    stack2(1,b) = temp;
    
    b = find( stack3(1,1:Nt)~=1e50 );
    temp = stack3(1,b);
    xdat = tbin(b);
    for n=1:length(b)
        b1 = find(abs(xdat-xdat(n))<=4*tH);
        temp(n) = sum( stack3(1,b(b1)).*0.5.^(((xdat(b1)-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(n))./tH).^2) );
    end
    stack3(1,b) = temp;
        
    b = find( stack4(1,1:Nt)~=1e50 );
    temp = stack4(1,b);
    xdat = tbin(b);
    for n=1:length(b)
        b1 = find(abs(xdat-xdat(n))<=4*tH);
        temp(n) = sum( stack4(1,b(b1)).*0.5.^(((xdat(b1)-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(n))./tH).^2) );
    end
    stack4(1,b) = temp;
    clear temp xdat b b1
    
    b = find(tbin>=tmed(1) & tbin<=tmed(2));
    tmed2 = (tmed(1)+tmed(2))/2;
    deltmed = (tmed(2)-tmed(1))/10;
    b1 = find(abs(tbin-tmed2)==min(abs(tbin-tmed2)));
    b2 = find(abs(tbin-tmed2-3*deltmed)==min(abs(tbin-tmed2-3*deltmed)));
    b3 = find(abs(tbin-tmed2-2*deltmed)==min(abs(tbin-tmed2-2*deltmed)));
    b4 = find(abs(tbin-tmed2-deltmed)==min(abs(tbin-tmed2-deltmed)));
    b1 = b1(1);
    b2 = b2(1);
    b3 = b3(1);
    b4 = b4(1);
    
    stack1(1,Nt+1) = stack1(1,b1);
    stack1(2,Nt+1) = stack1(1,b1) - median( stack1(2,b) );
    stack1(3,Nt+1) = stack1(1,b1) + median( stack1(3,b) );
    
    stack2(1,Nt+1) = stack2(1,b2);
    stack2(2,Nt+1) = stack2(1,b2) - median( stack2(2,b) );
    stack2(3,Nt+1) = stack2(1,b2) + median( stack2(3,b) );
    
    stack3(1,Nt+1) = stack3(1,b3);
    stack3(2,Nt+1) = stack3(1,b3) - median( stack3(2,b) );
    stack3(3,Nt+1) = stack3(1,b3) + median( stack3(3,b) );
    
    stack4(1,Nt+1) = stack4(1,b4);
    stack4(2,Nt+1) = stack4(1,b4) - max( median( stack4(2,b) ), 0.1*stack4(1,b4) );
    stack4(3,Nt+1) = stack4(1,b4) + max( median( stack4(3,b) ), 0.1*stack4(1,b4) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename1 = strcat(dirname,'LIN__',num2str(Tdown,'%04i'),'_Tmax_',num2str(Tup,'%04i'));
    filename2 = strcat(dirname,'LOG__',num2str(Tdown,'%04i'),'_Tmax_',num2str(Tup,'%04i'));
    filename1e = strcat(filename1,'.eps');
    filename1j = strcat(filename1,'.jpg');
    filename2e = strcat(filename2,'.eps');
    filename2j = strcat(filename2,'.jpg');    
    dirname = strcat(dirname,'INDIVIDUAL_CLUMPS_',num2str(Tdown,'%04i'),'_Tmax_',num2str(Tup,'%04i'),'/');
    mkdir(dirname)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:Nclump
        j = start_index(i);
        k = end_index(i);
        
        filename = strcat(dirname,num2str(i,'%02i'),'_LIN.jpg');
        figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],'Visible','off');
        set(figure1,'WindowStyle','docked');
        axes1 = axes('PareNt',figure1,'YTick',yl:dy:yu,...
            'YMinorTick','on',...
            'XTick',xl:dx:xu,...
            'XMinorTick','on',...
            'TickLength',[0.02 0.04],...
            'PlotBoxAspectRatio',[1 1 1],...
            'LineWidth',1.5,...
            'FoNtSize',14,...
            'FontName','Arial');
        xlim(axes1,[xl xu]);
        ylim(axes1,[yl yu]);
        box(axes1,'on');
        hold(axes1,'all');
        xlabel(xlab2,'Interpreter','latex','FoNtSize',16,...
            'FoNtName','Times New Roman','units','normalized',...
            'position',[0.5 -0.06 0]);
        ylabel(ylab2,'Interpreter','latex','FoNtSize',16,...
            'FoNtName','Times New Roman','units','normalized',...
            'position',[-0.10 0.45 0]);
        set(gcf,'renderer','painters')        
        l(1) = plot(clump_time(i,1:(k-j+1)),ydat1(i,1:(k-j+1)),'marker','none',...
            'linestyle','-','linewidth',3,'color','k',...
            'DisplayName','${\rm sSFR}$');        
        l(2) = plot(clump_time(i,1:(k-j+1)),ydat4(i,1:(k-j+1)),'marker','none',...
            'linestyle',':','linewidth',3,'color','r',...
            'DisplayName','${\rm stellar\:\:mass\:\:loss}$');        
        l(3) = plot(clump_time(i,1:(k-j+1)),ydat3(i,1:(k-j+1)),'marker','none',...
            'linestyle','-.','linewidth',3,'color','b',...
            'DisplayName','${\rm gas\:\:outflow}$');        
        l(4) = plot(clump_time(i,1:(k-j+1)),ydat2(i,1:(k-j+1)),'marker','none',...
            'linestyle','--','linewidth',3,'color','g',...
            'DisplayName','${\rm gas\:\:inflow}$');
        ti = get(gca,'TightInset');
        set(gca,'Position',[ti(1)+0.02 ti(2)+0.08 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
        set(gca,'units','ceNtimeters')
        pos = get(gca,'Position');
        set(gcf, 'PaperUnits','ceNtimeters');
        set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)-1 pos(4)+ti(2)+ti(4)+2]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)+2]);
%         legend1 = legend(gca,l(1:4));
%         set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
%             'Location','Best','FoNtSize',14,...
%             'Interpreter','latex');
%         legend boxoff   
        print(gcf,'-djpeg',filename);
        close all
        fclose all
        
        filename = strcat(dirname,num2str(i,'%02i'),'_LOG.jpg');
        figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],'Visible','off');
        set(figure1,'WindowStyle','docked');
        axes1 = axes('PareNt',figure1,'YTick',ytick,...
            'YScale','log',...
            'YtickLabel',ystr,...
            'XTick',xl:dx:xu,...
            'XMinorTick','on',...
            'TickLength',[0.02 0.04],...
            'PlotBoxAspectRatio',[1 1 1],...
            'LineWidth',1.5,...
            'FoNtSize',14,...
            'FoNtName','Arial');
        xlim(axes1,[xl xu]);
        ylim(axes1,[ytick(1) ytick(length(ytick))]);
        box(axes1,'on');
        hold(axes1,'all');
        xlabel(xlab2,'Interpreter','latex','FoNtSize',16,...
            'FoNtName','Times New Roman','units','normalized',...
            'position',[0.5 -0.06 0]);
        ylabel(ylab2,'Interpreter','latex','FoNtSize',16,...
            'FoNtName','Times New Roman','units','normalized',...
            'position',[-0.10 0.45 0]);
        set(gcf,'renderer','paiNters')        
        l(1) = plot(clump_time(i,1:(k-j+1)),ydat1(i,1:(k-j+1)),'marker','none',...
            'linestyle','-','linewidth',3,'color','k',...
            'DisplayName','${\rm sSFR}$');        
        l(2) = plot(clump_time(i,1:(k-j+1)),ydat4(i,1:(k-j+1)),'marker','none',...
            'linestyle',':','linewidth',3,'color','r',...
            'DisplayName','${\rm stellar\:\:mass\:\:loss}$');        
        l(3) = plot(clump_time(i,1:(k-j+1)),ydat3(i,1:(k-j+1)),'marker','none',...
            'linestyle','-.','linewidth',3,'color','b',...
            'DisplayName','${\rm gas\:\:outflow}$');        
        l(4) = plot(clump_time(i,1:(k-j+1)),ydat2(i,1:(k-j+1)),'marker','none',...
            'linestyle','--','linewidth',3,'color','g',...
            'DisplayName','${\rm gas\:\:inflow}$');
        ti = get(gca,'TightInset');
        set(gca,'Position',[ti(1)+0.02 ti(2)+0.08 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
        set(gca,'units','ceNtimeters')
        pos = get(gca,'Position');
        set(gcf, 'PaperUnits','ceNtimeters');
        set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)-1 pos(4)+ti(2)+ti(4)+2]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)+2]);
%         legend1 = legend(gca,l(1:4));
%         set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
%             'Location','Best','FoNtSize',14,...
%             'Interpreter','latex');
%         legend boxoff
        print(gcf,'-djpeg',filename);
        close all
        fclose all
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
    set(figure1,'WindowStyle','docked');
    axes1 = axes('PareNt',figure1,'YTick',yl:dy:yu,...
        'YMinorTick','on',...
        'XTick',xl:dx:xu,...
        'XMinorTick','on',...
        'TickLength',[0.02 0.04],...
        'PlotBoxAspectRatio',[1 1 1],...
        'LineWidth',1.5,...
        'FoNtSize',16,...
        'FoNtName','Arial');    
    xlim(axes1,[xl xu]);
    ylim(axes1,[yl yu]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');    
    xlabel(xlab,'INterpreter','latex','FoNtSize',24,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[0.5 -0.06 0]);    
    ylabel(ylab,'INterpreter','latex','FoNtSize',24,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);    
    title(tit,'INterpreter','latex','FoNtSize',24,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[0.5 1.02 0]);    
    set(gcf,'renderer','paiNters')
    
    b = find(stack1(1,1:Nt) ~= 1e50);
    l(1) = plot(tbin(b),stack1(1,b),'marker','none',...
        'linestyle','-','linewidth',3,'color','k',...
        'DisplayName','${\rm sSFR}$');
    
    b = find(stack4(1,1:Nt) ~= 1e50);
    l(2) = plot(tbin(b),stack4(1,b),'marker','none',...
        'linestyle',':','linewidth',3,'color','r',...
        'DisplayName','${\rm stellar\:\:mass\:\:loss}$');
    
    b = find(stack3(1,1:Nt) ~= 1e50);
    l(3) = plot(tbin(b),stack3(1,b),'marker','none',...
        'linestyle','-.','linewidth',3,'color','b',...
        'DisplayName','${\rm gas\:\:outflow}$');
    
    b = find(stack2(1,1:Nt) ~= 1e50);
    l(4) = plot(tbin(b),stack2(1,b),'marker','none',...
        'linestyle','--','linewidth',3,'color','g',...
        'DisplayName','${\rm gas\:\:inflow}$');    
    
    y = linspace(stack1(2,Nt+1),stack1(3,Nt+1),10);
    x = tmed2.*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','k');
    x = linspace(tmed2-0.4*deltmed,tmed2+0.4*deltmed,10);
    y = stack1(3,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','k');
    y = stack1(2,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','k');
    
    y = linspace(stack4(2,Nt+1),stack4(3,Nt+1),10);
    x = (tmed2+deltmed).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','r');
    x = linspace(tmed2+0.6*deltmed,tmed2+1.4*deltmed,10);
    y = stack4(3,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','r');
    y = stack4(2,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','r');
    
    y = linspace(stack3(2,Nt+1),stack3(3,Nt+1),10);
    x = (tmed2+2*deltmed).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','b');
    x = linspace(tmed2+1.6*deltmed,tmed2+2.4*deltmed,10);
    y = stack3(3,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','b');
    y = stack3(2,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','b');
    
    y = linspace(stack2(2,Nt+1),stack2(3,Nt+1),10);
    x = (tmed2+3*deltmed).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','g');
    x = linspace(tmed2+2.6*deltmed,tmed2+3.4*deltmed,10);
    y = stack2(3,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','g');
    y = stack2(2,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','g');
    if(leg==1)
        legend1 = legend(gca,l(1:4));
        set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
            'Position',leg_pos,'FoNtSize',14,...
            'Interpreter','latex');
        legend boxoff
    end
    print(gcf,'-depsc',filename1e);
    
    xlabel(xlab2,'INterpreter','latex','FoNtSize',20,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[0.5 -0.06 0]);    
    ylabel(ylab2,'INterpreter','latex','FoNtSize',20,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);    
    title(tit,'Interpreter','latex','FoNtSize',20,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[0.5 1.02 0]);
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1)+0.05 ti(2)+0.08 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','ceNtimeters')
    pos = get(gca,'Position');
    set(gcf, 'PaperUnits','ceNtimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    
    if(leg==1)
        set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
            'Position',leg_pos2,'FoNtSize',14,...
            'INterpreter','latex');
        legend boxoff
    end
    print(gcf,'-djpeg',filename1j);
    close all
    fclose all
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
    set(figure1,'WindowStyle','docked');
    axes1 = axes('PareNt',figure1,'YTick',ytick,...
        'YScale','log',...
        'YtickLabel',ystr,...
        'XTick',xl:dx:xu,...
        'XMinorTick','on',...
        'TickLength',[0.02 0.04],...
        'PlotBoxAspectRatio',[1 1 1],...
        'LineWidth',1.5,...
        'FoNtSize',16,...
        'FoNtName','Arial');    
    xlim(axes1,[xl xu]);
    ylim(axes1,[ytick(1) ytick(length(ytick))]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');    
    xlabel(xlab,'INterpreter','latex','FoNtSize',24,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[0.5 -0.06 0]);    
    ylabel(ylab,'INterpreter','latex','FoNtSize',24,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);    
    title(tit,'INterpreter','latex','FoNtSize',24,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[0.45 1.02 0]);    
    set(gcf,'renderer','paiNters')
    
    b = find(stack1(1,1:Nt) ~= 1e50);
    l(1) = plot(tbin(b),stack1(1,b),'marker','none',...
        'linestyle','-','linewidth',3,'color','k',...
        'DisplayName','${\rm sSFR}$');
    
    b = find(stack4(1,1:Nt) ~= 1e50);
    l(2) = plot(tbin(b),stack4(1,b),'marker','none',...
        'linestyle',':','linewidth',3,'color','r',...
        'DisplayName','${\rm stellar\:\:mass\:\:loss}$');
    
    b = find(stack3(1,1:Nt) ~= 1e50);
    l(3) = plot(tbin(b),stack3(1,b),'marker','none',...
        'linestyle','-.','linewidth',3,'color','b',...
        'DisplayName','${\rm gas\:\:outflow}$');
    
    b = find(stack2(1,1:Nt) ~= 1e50);
    l(4) = plot(tbin(b),stack2(1,b),'marker','none',...
        'linestyle','--','linewidth',3,'color','g',...
        'DisplayName','${\rm gas\:\:inflow}$');
    
    y = linspace(max(0.01,stack1(2,Nt+1)),stack1(3,Nt+1),10);
    x = tmed2.*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','k');
    x = linspace(tmed2-0.4*deltmed,tmed2+0.4*deltmed,10);
    y = stack1(3,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','k');
    y = max(0.01,stack1(2,Nt+1)).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','k');
    
    y = linspace(max(0.01,stack4(2,Nt+1)),stack4(3,Nt+1),10);
    x = (tmed2+deltmed).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','r');
    x = linspace(tmed2+0.6*deltmed,tmed2+1.4*deltmed,10);
    y = stack4(3,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','r');
    y = max(0.01,stack4(2,Nt+1)).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','r');
    
    y = linspace(max(0.01,stack3(2,Nt+1)),stack3(3,Nt+1),10);
    x = (tmed2+2*deltmed).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','b');
    x = linspace(tmed2+1.6*deltmed,tmed2+2.4*deltmed,10);
    y = stack3(3,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','b');
    y = max(0.01,stack3(2,Nt+1)).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','b');
    
    y = linspace(max(0.01,stack2(2,Nt+1)),stack2(3,Nt+1),10);
    x = (tmed2+3*deltmed).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','g');
    x = linspace(tmed2+2.6*deltmed,tmed2+3.4*deltmed,10);
    y = stack2(3,Nt+1).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','g');
    y = max(0.01,stack2(2,Nt+1)).*ones(size(y));    
    plot(x,y,'marker','none','linestyle','-','linewidth',3,'color','g');
    
    if(leg==1)
        legend1 = legend(gca,l(1:4));
        set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
            'Position',leg_pos,'FoNtSize',14,...
            'INterpreter','latex');
        legend boxoff
    end
    priNt(gcf,'-depsc',filename2e);
    
    xlabel(xlab2,'INterpreter','latex','FoNtSize',20,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[0.5 -0.06 0]);    
    ylabel(ylab2,'INterpreter','latex','FoNtSize',20,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);    
    title(tit,'INterpreter','latex','FoNtSize',20,...
        'FoNtName','Times New Roman','units','normalized',...
        'position',[0.5 1.02 0]);    
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1)+0.05 ti(2)+0.08 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','ceNtimeters')
    pos = get(gca,'Position');
    set(gcf, 'PaperUnits','ceNtimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    
    if(leg==1)
        set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
            'Position',leg_pos2,'FoNtSize',14,...
            'INterpreter','latex');
        legend boxoff
    end
    priNt(gcf,'-djpeg',filename2j);
    close all
    fclose all
end