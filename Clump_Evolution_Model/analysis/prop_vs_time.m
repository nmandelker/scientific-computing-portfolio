function prop_vs_time(is, nis, ind, gal, prop, tmax_window, Mmax_window, zform_window)

% gal is the simulation index: V**
% gal==0 means stack all simulations together

%(prop = 1: ID)
%prop = 2: redshift
%prop = 3: Rc
%prop = 4: M_{gas}
%prop = 5: M_{star}
%prop = 6: M_{bar}
%prop = 7: f_{gas}
%prop = 8: f_{dm}
%prop = 9: {\Sigma}_{gas}
%prop = 10: {\Sigma}_{*}
%prop = 11: {\Sigma}_{bar}
%prop = 12: age
%prop = 13: z_{gas}
%prop = 14: z_{stars}
%prop = 15: SFR
%prop = 16: {\Sigma}_{SFR}
%prop = 17: sSFR
%prop = 18: t_{dep}
%prop = 19: d/Rd
%prop = 20: z/Hd
%prop = 21: d
%prop = 22: z
%prop = 23: mean_residual
%prop = 24: shape parameter
%prop = 25: dark matter contrast
%(prop = 26: Ex situ)
%prop = 28: tff
%prop = 29: tdyn_local
%prop = 30: tdyn_global
%prop = 45: time

dirname = './clump_evolution/';
mkdir(dirname);

if(gal~=0)
    dirname = strcat(dirname,'V',num2str(gal,'%02i'),'/');
    mkdir(dirname);
else
    dirname = strcat(dirname,'all/');
    mkdir(dirname);
end
if(prop==6)
    log_prop = 1;
    norm_stack = 1;
    yl = 6;
%     yl = 7.2;
    yu = 9.2;
    dy = 0.2;
    normyl = -1;
    normyu = 1;
    normdy = 0.2;
    dirname = strcat(dirname,'mass/');
    ylab = '${\rm log(}M\:{\rm [M_{\odot}])}$';
    ylab2 = '${\rm log(}\:\:M\:\:{\rm [M_{\odot}]\:)}$';
    norm_ylab = '${\rm log(}M\:{\rm [normalized])}$';
    norm_ylab2 = '${\rm log(}\:\:M\:\:{\rm [normalized]\:)}$';
elseif(prop==7)
    log_prop = 1;
    norm_stack = 1;
    yl = -4;
%     yl = -2;
    yu = 0;
    dy = 0.4;
%     dy = 0.2;
    normyl = -3;
    normyu = 1;
%     normyl = -1.4;
%     normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'fgas/');
    ylab = '${\rm log(}f_{\rm g})$';
    ylab2 = '${\rm log(\:\:}f_{\rm g}\:)$';
    norm_ylab = '${\rm log(}f_{\rm g}\:{\rm [normalized])}$';
    norm_ylab2 = '${\rm log(\:\:}f_{\rm g}\:\:{\rm [normalized]\:)}$';
elseif(prop==15)
    log_prop = 1;
    norm_stack = 1;
    yl = -4;
%     yl = -2;
    yu = 1;
    dy = 0.2;
    normyl = -1;
    normyu = 1;
%     normyl = -1.6;
%     normyu = 0.4;
    normdy = 0.2;
    dirname = strcat(dirname,'sfr/');
    ylab = '${\rm log(}SFR\:{\rm [M_{\odot}\:yr^{-1}])}$';
    ylab2 = '${\rm log(\:\:}SFR\:\:{\rm [M_{\odot}\:yr^{-1}\:])}$';
    norm_ylab = '${\rm log(}SFR\:{\rm [normalized])}$';
    norm_ylab2 = '${\rm log(\:\:}SFR\:\:{\rm [normalized])}$';
elseif(prop==21)
    log_prop = 1;
    norm_stack = 1;
    yl = -1;
    yu = 2;
%     yl = -0.4;
%     yu = 1.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
%     normyl = -1;
%     normyu = 0.4;
    normdy = 0.2;
    dirname = strcat(dirname,'dist/');
    ylab = '${\rm log(}d\:{\rm [kpc])}$';
    ylab2 = '${\rm log(\:\:}d\:\:{\rm [kpc\:])}$';
    norm_ylab = '${\rm log(}d\:{\rm [normalized])}$';
    norm_ylab2 = '${\rm log(\:\:}d\:\:{\rm [normalized])}$';
elseif(prop==22)
    log_prop = 1;
    norm_stack = 1;
    yl = -2;
    yu = 1;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'height/');
    ylab = '${\rm log(}h\:{\rm [kpc])}$';
    ylab2 = '${\rm log(\:\:}h\:\:{\rm [kpc\:])}$';
    norm_ylab = '${\rm log(}h\:{\rm [normalized])}$';
    norm_ylab2 = '${\rm log(\:\:}h\:\:{\rm [normalized])}$';
elseif(prop==49)
    log_prop = 1;
    norm_stack = 0;
    yl = 0;
    yu = 2;
    dy = 0.2;
    normyl = -1;
    normyu = 1;
    normdy = 0.2;
    dirname = strcat(dirname,'td_local_over_tff/');
    ylab = '${\rm log(}t_{\rm d,\,l}/t_{\rm ff})$';
    ylab2 = '${\rm log(\:\:}t_{\rm d,\,l}\:/t_{\rm ff}\:)$';
elseif(prop==50)
    log_prop = 1;
    norm_stack = 0;
    yl = 0;
    yu = 2;
    dy = 0.2;
    normyl = -1;
    normyu = 1;
    normdy = 0.2;
    dirname = strcat(dirname,'td_global_over_tff/');
    ylab = '${\rm log(}t_{\rm d,\,g}/t_{\rm ff})$';
    ylab2 = '${\rm log(\:\:}t_{\rm d,\,g}\:/t_{\rm ff}\:)$';
elseif(prop==65)
    log_prop = 1;
    norm_stack = 0;
    yl = -3;
    yu = 0;
    dy = 0.2;
    normyl = -1;
    normyu = 1;
    normdy = 0.2;
    dirname = strcat(dirname,'eps_ff/');
    ylab = '${\rm log(}\epsilon_{\rm ff})$';
    ylab2 = '${\rm log(\:\:}\epsilon_{\rm ff}\:)$';
elseif(prop==66)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_Rc_Vr_0/');
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';    
elseif(prop==67)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_Rc_Vr_Vesc/');
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';
elseif(prop==68)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_Rc_Vr_0_V_Vesc/');
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';
elseif(prop==69)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_avg_Vr_0/');
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';
elseif(prop==70)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_avg_Vr_Vesc/');
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';
elseif(prop==71)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_avg_Vr_0_V_Vesc/');
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';
elseif(prop==72)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_Frederic/');
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';
elseif(prop==73)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'eta_stars/');
    ylab = '${\rm log(}\eta_{\rm s})$';
    ylab2 = '${\rm log(\:\:}\eta_{\rm s}\:\:)$';
elseif(prop==74)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'eta_stars_net/');
    ylab = '${\rm log(}\eta_{\rm s})$';
    ylab2 = '${\rm log(\:\:}\eta_{\rm s}\:\:)$';
elseif(prop==80)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_avg_Vr_0/');
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==81)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Rc_Vr_Vesc_gas_cons/');
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==82)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Rc_Vr_0_V_Vesc_gas_cons/');
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==83)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_avg_Vr_Vesc_gas_cons/');
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==84)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_avg_Vr_0_V_Vesc_gas_cons/');
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==85)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Frederic_gas_cons/');
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==92)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Rc_Vr_Vesc_bar_cons/');
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==93)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Rc_Vr_0_V_Vesc_bar_cons/');
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==94)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_avg_Vr_Vesc_bar_cons/');
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==95)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_avg_Vr_0_V_Vesc_bar_cons/');
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==96)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Frederic_bar_cons/');
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
end
mkdir(dirname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col(:,1)  = [0 0 0];
col(:,2)  = [1 0 0];
col(:,3)  = [0 1 0];
col(:,4)  = [0 0 1];
col(:,5)  = [0 1 1];
col(:,6)  = [1 0 1];
col(:,7)  = [1 1 0];
col(:,8)  = [0.7 0.7 0.7];
col(:,9)  = [0.2 0.7 0.4];
col(:,10) = [0.7 0.4 0.2];
col(:,11)  = [0.1 0.2 0.8];
col(:,12)  = [0.8 0.1 0.5];
col(:,13) = [0.3 0.8 0.3];
Ncol = 13;

Tdown = tmax_window(1);
Tup   = tmax_window(2);

Mdown = Mmax_window(1);
Mup   = Mmax_window(2);

zdown = zform_window(1);
zup   = zform_window(2);

tit = strcat('$',num2str(Mdown,'%3.1f'),'<{\rm log}(M_{\rm c,\;50})<',num2str(Mup,'%3.1f'));
tit2 = strcat('$',num2str(Mdown,'%3.1f'),'<{\rm log}\:\:(M_{\rm c,\;50}\:\:)<',num2str(Mup,'%3.1f'));
dirname = strcat(dirname,num2str(Mdown,'%3.1f'),'_M50_',num2str(Mup,'%3.1f'),'/');

tit = strcat(tit,'\:\:\:\:',num2str(zdown,'%3.1f'),'<z_{\rm form}<',num2str(zup,'%3.1f'));
tit2 = strcat(tit2,'\:\:\:\:\:\:\:',num2str(zdown,'%3.1f'),'<z_{\rm form}\:<',num2str(zup,'%3.1f'));
dirname = strcat(dirname,num2str(zdown,'%3.1f'),'_zform_',num2str(zup,'%3.1f'),'/');

tit = strcat(tit,'\:\:\:\:',num2str(Tdown,'%4i'),'<t_{\rm c,\;max}<',num2str(Tup,'%4i'),'$');
tit2 = strcat(tit2,'\:\:\:\:\:\:\:',num2str(Tdown,'%4i'),'<t_{\rm c,\;max}\:\:<',num2str(Tup,'%4i'),'$');

dt = 20; %Myr
tmax = 4000;
tmin = 0;
tbin = (tmin+dt/2):dt:(tmax-dt/2);
Nt = length(tbin);

dlogt = 0.1; %Myr
logtmax = 4;
logtmin = -2;
logtbin = (logtmin+dlogt/2):dlogt:(logtmax-dlogt/2);
Nlogt = length(logtbin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(gal~=0)
    m = find(ind==gal);
    vec = (nis(m)+1):(nis(m+1)-1);
    vec2 = vec+1;
    minz = min(is(vec2,2))
    
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
        minz = min(is(vec2,2));
        
        b = find(is(vec,58)>=Tdown & is(vec,58)<Tup & is(vec,1)~=is(vec2,1));
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
    [0,0,0,0,length(start_index), length(end_index)]
    return
elseif( min(end_index-start_index)<=0 )
    [0,0,0,0,min(end_index-start_index)]
    return
else
    
    Nclump = length(start_index);
    clump_prop      = 1e50.*ones(Nclump,max(end_index-start_index));
    norm_clump_prop = 1e50.*ones(Nclump,max(end_index-start_index));
    clump_time      = 1e50.*ones(Nclump,max(end_index-start_index));
    clump_stack     = 1e50.*ones(3,Nlogt);
    last_snap       = 1e50.*ones(Nclump,1);
    [min(is(end_index(1:Nclump),2)), minz]
    for i=1:Nclump
        j = start_index(i);
        k = end_index(i);
        if(norm_stack == 1)
            b = find(abs(is(j:k,48)-50)<=10);
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
            avg_prop = mean(abs(is(j-1+b,prop)));
        else
            avg_prop = 1;
        end
        for n=j:k
            if(log_prop)
                clump_prop(i,n)      = log10(abs(is(n,prop)));
                norm_clump_prop(i,n) = log10(abs(is(n,prop)./avg_prop));
            else
                clump_prop(i,n) = is(n,prop);
                norm_clump_prop(i,n) = is(n,prop)./avg_prop;
            end
            clump_time(i,n) = log10(is(n,48));
        end
        if(is(k,2)==minz)
            if(log_prop)
                last_snap(i) = log10(abs(is(k,prop)));
            else
                last_snap(i) = is(k,prop);
            end
        end
    end
    for i=1:Nlogt
        b = find(clump_time>(logtbin(i)-dlogt/2) & clump_time<=(logtbin(i)+dlogt/2));
        temp_vec = sort(norm_clump_prop(b));
        Ntemp_vec = length(temp_vec);
        if(Ntemp_vec>=2)
            Nlow  = floor(Ntemp_vec/6) + 1;
            Nhigh = floor(5*Ntemp_vec/6);
            clump_stack(1,i) = median(temp_vec);
            clump_stack(2,i) = median(temp_vec) - temp_vec(Nlow);
            clump_stack(3,i) = temp_vec(Nhigh) - median(temp_vec);
        end
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dirname = strcat(dirname,'log_time/');
    mkdir(dirname);
    filename1 = strcat(dirname,'INDIVIDUAL__',num2str(Tdown,'%04i'),'_Tmax_',num2str(Tup,'%04i'));
    filename2 = strcat(dirname,'STACKED__',num2str(Tdown,'%04i'),'_Tmax_',num2str(Tup,'%04i'));
    filename1e = strcat(filename1,'.eps');
    filename1j = strcat(filename1,'.jpg');
    filename2e = strcat(filename2,'.eps');
    filename2j = strcat(filename2,'.jpg');
    
    xlab = '${\rm log(}t\:{\rm [Myr])}$';
    xlab2 = '${\rm log(}\:\:t\:\:{\rm [Myr])}$';
    xl = 0;
    xu = 3.6;
    dx = 0.4;
    
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
        'position',[0.45 1.02 0]);    
    set(gcf,'renderer','painters')
    
    for i=1:Nclump
        m = mod(i,10);
        m = m+1;
        if(mod(floor((i-1)/10),2)==0)
            lin = '-';
        else
            lin = '--';
        end
        b = find(clump_prop(i,:) ~= 1e50);
        if(log_prop)
            vec  = log10(abs(is(start_index(i):end_index(i),prop)));
        else
            vec  = is(start_index(i):end_index(i),prop);
        end
        time_vec = is(start_index(i):end_index(i),45);
        temp_vec = vec;
        Nvec = length(vec);
        for j=1:Nvec
            b = find(abs(time_vec-time_vec(j))<=15);
            vec(j) = mean(temp_vec(b));
        end
%         vec(1)        = ( temp_vec(1) + temp_vec(2) + temp_vec(3) ) / 3;
%         vec(2)        = ( temp_vec(1) + temp_vec(2) + temp_vec(3) + temp_vec(4) ) / 4;
%         vec(3:Nvec-2) = ( temp_vec(1:Nvec-4) + temp_vec(2:Nvec-3) + ...
%             temp_vec(3:Nvec-2) + temp_vec(4:Nvec-1) + temp_vec(5:Nvec) ) / 5;
%         vec(Nvec-1)   = ( temp_vec(Nvec-3) + temp_vec(Nvec-2) + temp_vec(Nvec-1) + temp_vec(Nvec) ) / 4;
%         vec(Nvec)     = ( temp_vec(Nvec-2) + temp_vec(Nvec-1) + temp_vec(Nvec) ) / 3;
        plot(log10(is(start_index(i):end_index(i),45)),vec,...
            'marker','none',...
            'linestyle',lin,'linewidth',2,'color',col(:,m));
        plot(log10(is(end_index(i),45)),last_snap(i),...
            'linestyle','none',...
            'marker','square','MarkerSize',10,'MarkerFaceColor',col(:,m),'MarkerEdgeColor',col(:,m));
    end
    print(gcf,'-depsc',filename1e);
    
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);    
    ylabel(ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);    
    title(tit2,'Interpreter','latex','FontSize',14,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.42 1.02 0]);    
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1)+0.05 ti(2)+0.08 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    
    print(gcf,'-djpeg',filename1j);
    close all
    fclose all
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
    set(figure1,'WindowStyle','docked');
    if(norm_stack == 0)
        normyl=yl;
        normdy=dy;
        normyu=yu;
        norm_ylab = ylab;
        norm_ylab2 = ylab2;
    end
    axes1 = axes('Parent',figure1,'YTick',normyl:normdy:normyu,...
        'YMinorTick','on',...
        'XTick',xl:dx:xu,...
        'XMinorTick','on',...
        'TickLength',[0.02 0.04],...
        'PlotBoxAspectRatio',[1 1 1],...
        'LineWidth',1.5,...
        'FontSize',16,...
        'FontName','Arial');    
    xlim(axes1,[xl xu]);
    ylim(axes1,[normyl normyu]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');    
    xlabel(xlab,'Interpreter','latex','FontSize',24,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);    
    ylabel(norm_ylab,'Interpreter','latex','FontSize',24,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);    
    title(tit,'Interpreter','latex','FontSize',14,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.45 1.02 0]);    
    set(gcf,'renderer','painters')
    
    b = find(clump_stack(1,:) ~= 1e50);
    errorbar(logtbin(b),clump_stack(1,b),clump_stack(2,b),clump_stack(3,b),'marker','none',...
        'linestyle','-','linewidth',4,'color','k');
    print(gcf,'-depsc',filename2e);
    
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);    
    ylabel(norm_ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);    
    title(tit2,'Interpreter','latex','FontSize',14,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.42 1.02 0]);    
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1)+0.05 ti(2)+0.08 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+3 pos(4)+ti(2)+ti(4)+3]);
    
    print(gcf,'-djpeg',filename2j);
    close all
    fclose all
end

        