function specific_rates(is, nis, ind, gal, ver, tmax_window, Mmax_window, zform_window)

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
ylab_log = '${\rm log(specific\:\:rates\:\:[Gyr^{-1}])}$';
ylab_log2 = '${\rm log(specific\:\:\:rates\:\:\:[Gyr^{-1}\:])}$';
log_yu = 2;
log_yl = -2;
log_dy = 0.4;

xlab = '$t\:{\rm [Myr]}$';
xlab2 = '$t\:\:{\rm [Myr]}$';
xlab_log = '${\rm log(}t\:{\rm [Myr])}$';
xlab_log2 = '${\rm log(}\:\:t\:\:{\rm [Myr])}$';
xl = 0;
log_xl = 0.6;
if(gal==19)
    yu = 20;
    yl = 0;
    dy = 2;
    xu = 200;
    dx = 20;
    dt = 5; %Myr
    log_xu = 2.4;
    log_dx = 0.2;
else
    yu = 8;
    yl = 0;
    dy = 1;
    xu = 1000;
    dx = 100;
    dt = 10; %Myr
    log_xu = 3;
    log_dx = 0.2;
end

tmax = 4000;
tmin = 0;
tbin = (tmin+dt/2):dt:(tmax-dt/2);
Nt = length(tbin);

dlogt = 0.1;
logtmax = 4;
logtmin = -2;
logtbin = (logtmin+dlogt/2):dlogt:(logtmax-dlogt/2);
Nlogt = length(logtbin);

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
    if(isempty(b))
        b = find(abs(is(tvec,48)-50)<=20);
        if(isempty(b))
            b = find(abs(is(tvec,48)-50)<=30);
            if(isempty(b))
                b = find(abs(is(tvec,48)-50)<=40);
                if(isempty(b))
                    b = find(abs(is(tvec,48)-50)<=50);
                end
            end
        end
    end
    Mavg(i) = mean(log10(is(tvec(b),6)));
end
b = find(Mavg>=Mdown & Mavg<Mup);
start_index = start_index(b);
end_index = end_index(b);
Mavg = Mavg(b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[length(start_index), length(end_index)]
if( length(start_index) ~= length(end_index) )
    [0,0,length(start_index), length(end_index)]
    return
elseif( min(end_index-start_index)<=0 )
    [0,0,min(end_index-start_index)]
    return
else
    
    Nclump = length(start_index);
    ydat1  = 1e50.*ones(Nclump,max(end_index-start_index));
    ydat2  = ydat1;
    ydat3  = ydat1;
    ydat4  = ydat1;
    clump_time = ydat1;
    stack1 = 1e50.*ones(3,Nt);
    stack2 = stack1;
    stack3 = stack1;
    stack4 = stack1;
    
    log_ydat1  = 1e50.*ones(Nclump,max(end_index-start_index));
    log_ydat2  = log_ydat1;
    log_ydat3  = log_ydat1;
    log_ydat4  = log_ydat1;
    log_clump_time = log_ydat1;
    log_stack1 = 1e50.*ones(3,Nlogt);
    log_stack2 = log_stack1;
    log_stack3 = log_stack1;
    log_stack4 = log_stack1;
    min(is(end_index(1:Nclump),2))
    
    for i=1:Nclump
        j = start_index(i);
        k = end_index(i);
        for n=j:k
            ydat1(i,n) = max( 1e-2, 1e9 .* is(n,15)./is(n,6) ); % sSFR Gyr^{-1}
            ydat2(i,n) = max( 1e-2, 1e9 .* is(n,prop2)./is(n,6) ); % sMg_in Gyr^{-1}
            ydat3(i,n) = max( 1e-2, 1e9 .* is(n,prop3)./is(n,6) ); % sMg_out Gyr^{-1}
            ydat4(i,n) = max( 1e-2, 1e9.*(is(n,45)-is(n,44))./is(n,6) ); % sM*_out_net Gyr^{-1}
            
            log_ydat1(i,n) = log10(ydat1(i,n));
            log_ydat2(i,n) = log10(ydat2(i,n));
            log_ydat3(i,n) = log10(ydat3(i,n));
            log_ydat4(i,n) = log10(ydat4(i,n));
            
            clump_time(i,n) = is(n,48);
            log_clump_time(i,n) = log10(clump_time(i,n));
        end
    end
    for i=1:Nt
        b = find(clump_time>(tbin(i)-dt/2) & clump_time<=(tbin(i)+dt/2));
        temp_vec1 = sort(ydat1(b));
        temp_vec2 = sort(ydat2(b));
        temp_vec3 = sort(ydat3(b));
        temp_vec4 = sort(ydat4(b));
        Ntemp_vec = length(b);
        if(Ntemp_vec>=3)
            Nlow  = floor(Ntemp_vec/6) + 1;
            Nhigh = floor(5*Ntemp_vec/6);
            
            stack1(1,i) = median(temp_vec1);
            stack1(2,i) = median(temp_vec1) - temp_vec1(Nlow);
            stack1(3,i) = temp_vec1(Nhigh) - median(temp_vec1);
            
            stack2(1,i) = median(temp_vec2);
            stack2(2,i) = median(temp_vec2) - temp_vec2(Nlow);
            stack2(3,i) = temp_vec2(Nhigh) - median(temp_vec2);
            
            stack3(1,i) = median(temp_vec3);
            stack3(2,i) = median(temp_vec3) - temp_vec3(Nlow);
            stack3(3,i) = temp_vec3(Nhigh) - median(temp_vec3);
            
            stack4(1,i) = median(temp_vec4);
            stack4(2,i) = median(temp_vec4) - temp_vec4(Nlow);
            stack4(3,i) = temp_vec4(Nhigh) - median(temp_vec4);
        end
    end
    for i=1:Nlogt
        b = find(log_clump_time>(logtbin(i)-dlogt/2) & log_clump_time<=(logtbin(i)+dlogt/2));
        temp_vec1 = sort(log_ydat1(b));
        temp_vec2 = sort(log_ydat2(b));
        temp_vec3 = sort(log_ydat3(b));
        temp_vec4 = sort(log_ydat4(b));
        Ntemp_vec = length(b);
        if(Ntemp_vec>=3)
            Nlow  = floor(Ntemp_vec/6) + 1;
            Nhigh = floor(5*Ntemp_vec/6);
            
            log_stack1(1,i) = median(temp_vec1);
            log_stack1(2,i) = median(temp_vec1) - temp_vec1(Nlow);
            log_stack1(3,i) = temp_vec1(Nhigh) - median(temp_vec1);
            
            log_stack2(1,i) = median(temp_vec2);
            log_stack2(2,i) = median(temp_vec2) - temp_vec2(Nlow);
            log_stack2(3,i) = temp_vec2(Nhigh) - median(temp_vec2);
            
            log_stack3(1,i) = median(temp_vec3);
            log_stack3(2,i) = median(temp_vec3) - temp_vec3(Nlow);
            log_stack3(3,i) = temp_vec3(Nhigh) - median(temp_vec3);
            
            log_stack4(1,i) = median(temp_vec4);
            log_stack4(2,i) = median(temp_vec4) - temp_vec4(Nlow);
            log_stack4(3,i) = temp_vec4(Nhigh) - median(temp_vec4);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename1 = strcat(dirname,'LINEAR__',num2str(Tdown,'%04i'),'_Tmax_',num2str(Tup,'%04i'));
    filename2 = strcat(dirname,'LOG__',num2str(Tdown,'%04i'),'_Tmax_',num2str(Tup,'%04i'));
    filename1e = strcat(filename1,'.eps');
    filename1j = strcat(filename1,'.jpg');
    filename2e = strcat(filename2,'.eps');
    filename2j = strcat(filename2,'.jpg');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
    set(figure1,'WindowStyle','docked');
    axes1 = axes('Parent',figure1,'YTick',yl:dy:yu,...
        'YMinorTick','on',...
        'XTick',xl:dx:xu,...
        'XMinorTick','on',...
        'TickLength',[0.02 0.04],...
        'PlotBoxAspectRatio',[1 1 1],...
        'LineWidth',1.5,...
        'FontSize',16,...
        'FontName','Arial');    
    xlim(axes1,[xl xu]);
    ylim(axes1,[yl yu]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');    
    xlabel(xlab,'Interpreter','latex','FontSize',24,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);    
    ylabel(ylab,'Interpreter','latex','FontSize',24,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);    
    title(tit,'Interpreter','latex','FontSize',14,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 1.02 0]);    
    set(gcf,'renderer','painters')
    
    b = find(stack1(1,:) ~= 1e50);
    N = length(b);
    ydat = stack1(1,b);
    xdat = tbin(b);
    for i=1:N
        b1 = find(abs(xdat-xdat(i))<=100);
        ydat(i) = sum( stack1(1,b(b1)).*0.5.^(((xdat(b1)-xdat(i))./25).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(i))./25).^2) );
    end
    %errorbar(logtbin(b),stack1(1,b),stack1(2,b),stack1(3,b),'marker','none',...
    %    'linestyle','-','linewidth',4,'color','k');
    l(1) = plot(xdat,ydat,'marker','none',...
        'linestyle','-','linewidth',3,'color','k',...
        'DisplayName','${\rm sSFR}$');
    b = find(stack4(1,:) ~= 1e50);
    N = length(b);
    ydat = stack4(1,b);
    xdat = tbin(b);
    for i=1:N
        b1 = find(abs(xdat-xdat(i))<=100);
        ydat(i) = sum( stack4(1,b(b1)).*0.5.^(((xdat(b1)-xdat(i))./25).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(i))./25).^2) );
    end
    l(2) = plot(xdat,ydat,'marker','none',...
        'linestyle',':','linewidth',3,'color','r',...
        'DisplayName','${\rm stellar\:\:mass\:\:loss}$');
    b = find(stack3(1,:) ~= 1e50);
    N = length(b);
    ydat = stack3(1,b);
    xdat = tbin(b);
    for i=1:N
        b1 = find(abs(xdat-xdat(i))<=50);
        ydat(i) = sum( stack3(1,b(b1)).*0.5.^(((xdat(b1)-xdat(i))./25).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(i))./25).^2) );
    end
    l(3) = plot(xdat,ydat,'marker','none',...
        'linestyle','-.','linewidth',3,'color','b',...
        'DisplayName','${\rm gas\:\:outflow}$');
    b = find(stack2(1,:) ~= 1e50);
    N = length(b);
    ydat = stack2(1,b);
    xdat = tbin(b);
    for i=1:N
        b1 = find(abs(xdat-xdat(i))<=50);
        ydat(i) = sum( stack2(1,b(b1)).*0.5.^(((xdat(b1)-xdat(i))./25).^2) ) ./ sum( 0.5.^(((xdat(b1)-xdat(i))./25).^2) );
    end
    l(4) = plot(xdat,ydat,'marker','none',...
        'linestyle','--','linewidth',3,'color','g',...
        'DisplayName','${\rm gas\:\:inflow}$');
    
    legend1 = legend(gca,l(1:4));
    set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
        'Position',[0.7, 0.7, 0.16, 0.12],'FontSize',14,...
        'Interpreter','latex');
    legend boxoff
    print(gcf,'-depsc',filename1e);
    
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);    
    ylabel(ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);    
    title(tit,'Interpreter','latex','FontSize',14,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 1.02 0]);    
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1)+0.05 ti(2)+0.08 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    
    set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
        'Position',[0.6, 0.77, 0.16, 0.12],'FontSize',14,...
        'Interpreter','latex');
    legend boxoff
    print(gcf,'-djpeg',filename1j);
    close all
    fclose all
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
    set(figure1,'WindowStyle','docked');
    axes1 = axes('Parent',figure1,'YTick',log_yl:log_dy:log_yu,...
        'YMinorTick','on',...
        'XTick',log_xl:log_dx:log_xu,...
        'XMinorTick','on',...
        'TickLength',[0.02 0.04],...
        'PlotBoxAspectRatio',[1 1 1],...
        'LineWidth',1.5,...
        'FontSize',16,...
        'FontName','Arial');    
    xlim(axes1,[log_xl log_xu]);
    ylim(axes1,[log_yl log_yu]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');    
    xlabel(xlab_log,'Interpreter','latex','FontSize',24,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);    
    ylabel(ylab_log,'Interpreter','latex','FontSize',24,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);    
    title(tit,'Interpreter','latex','FontSize',14,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.45 1.02 0]);    
    set(gcf,'renderer','painters')
    
    b = find(log_stack1(1,:) ~= 1e50);
    N = length(b);
    ydat = log_stack1(1,b);
    xdat = logtbin(b);
    for i=1:N
        b1 = find(abs(10.^xdat-10^xdat(i))<=50);
        ydat(i) = sum( log_stack1(1,b(b1)).*0.5.^(((10.^xdat(b1)-10^xdat(i))./25).^2) ) ./ sum( 0.5.^(((10.^xdat(b1)-10^xdat(i))./25).^2) );
    end
    %errorbar(logtbin(b),stack1(1,b),stack1(2,b),stack1(3,b),'marker','none',...
    %    'linestyle','-','linewidth',4,'color','k');
    l(1) = plot(xdat,ydat,'marker','none',...
        'linestyle','-','linewidth',3,'color','k',...
        'DisplayName','${\rm sSFR}$');
    b = find(log_stack4(1,:) ~= 1e50);
    N = length(b);
    ydat = log_stack4(1,b);
    xdat = logtbin(b);
    for i=1:N
        b1 = find(abs(10.^xdat-10^xdat(i))<=50);
        ydat(i) = sum( log_stack4(1,b(b1)).*0.5.^(((10.^xdat(b1)-10^xdat(i))./25).^2) ) ./ sum( 0.5.^(((10.^xdat(b1)-10^xdat(i))./25).^2) );
    end
    l(2) = plot(xdat,ydat,'marker','none',...
        'linestyle',':','linewidth',3,'color','r',...
        'DisplayName','${\rm stellar\:\:mass\:\:loss}$');
    b = find(log_stack3(1,:) ~= 1e50);
    N = length(b);
    ydat = log_stack3(1,b);
    xdat = logtbin(b);
    for i=1:N
        b1 = find(abs(10.^xdat-10^xdat(i))<=50);
        ydat(i) = sum( log_stack3(1,b(b1)).*0.5.^(((10.^xdat(b1)-10^xdat(i))./25).^2) ) ./ sum( 0.5.^(((10.^xdat(b1)-10^xdat(i))./25).^2) );
    end
    l(3) = plot(xdat,ydat,'marker','none',...
        'linestyle','-.','linewidth',3,'color','b',...
        'DisplayName','${\rm gas\:\:outflow}$');
    b = find(log_stack2(1,:) ~= 1e50);
    N = length(b);
    ydat = log_stack2(1,b);
    xdat = logtbin(b);
    for i=1:N
        b1 = find(abs(10.^xdat-10^xdat(i))<=50);
        ydat(i) = sum( log_stack2(1,b(b1)).*0.5.^(((10.^xdat(b1)-10^xdat(i))./25).^2) ) ./ sum( 0.5.^(((10.^xdat(b1)-10^xdat(i))./25).^2) );
    end
    l(4) = plot(xdat,ydat,'marker','none',...
        'linestyle','--','linewidth',3,'color','g',...
        'DisplayName','${\rm gas\:\:inflow}$');
    
    legend1 = legend(gca,l(1:4));
    set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
        'Position',[0.7, 0.7, 0.16, 0.12],'FontSize',14,...
        'Interpreter','latex');
    legend boxoff    
    print(gcf,'-depsc',filename2e);
    
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);    
    ylabel(ylab_log2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);    
    title(tit,'Interpreter','latex','FontSize',14,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 1.02 0]);    
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1)+0.05 ti(2)+0.08 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    
    set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
        'Position',[0.6, 0.77, 0.16, 0.12],'FontSize',14,...
        'Interpreter','latex');
    legend boxoff
    print(gcf,'-djpeg',filename2j);
    close all
    fclose all
end

        