clear

thin = 1;
gal = [7];
%gal = [1:17,19:35];

version = 10;
decomp = 5;
pc = 100;
average = 'volume';

Eind = 6:9;

xl = [0.1 0.1];
xu = [0.55 0.55];
dx = 0.05;
xtit = '$a$';
for k=1:length(Eind)
    if(Eind(k)==6 | Eind(k)==7)
        yl = [-1 -1];
        yu = [1 1];
        dy = 0.2;
        ytit = '$2{\vec {\sigma}}_{\rm c}\cdot{\vec {\sigma}}_{\rm s}\:\:/\:\sigma_{\rm tot}\:^2$';
    else
        yl = [0 0];
        yu = [4 4];
        dy = 0.5;
        ytit = '$E_{\rm sol}\:\:/\:E_{\rm comp}$';
    end
    for i=1:length(gal)
        if(thin==1)
            dirname0 = strcat('thin_timesteps/version_',num2str(version,'%02i'),'/');
        else
            dirname0 = strcat('version_',num2str(version,'%02i'),'/');
        end
        if(version < 9) then
            dirname = strcat(dirname0,'VELA_v2_',num2str(gal(i),'%02i'),'_',num2str(pc,'%03i'),'pc/');
        else
            dirname0 = strcat(dirname0,'decomp_',num2str(decomp,'%1i'),'/');
            if(decomp ==3)
                dirname0 = strcat(dirname0,num2str(pc,'%03i'),'pc_',average,'_average/');
            else
                dirname0 = strcat(dirname0,num2str(pc,'%03i'),'pc/');
            end
            dirname = strcat(dirname0,'VELA_v2_',num2str(gal(i),'%02i'),'/');
        end
        galname = strcat('V',num2str(gal(i),'%02i'));
        %input_name = strcat(dirname,'solenoidal_over_compressive.out');
        input_name = strcat(dirname,'compressive_over_solenoidal.out');
        if(Eind(k)==9)
            output_dir = strcat(dirname0,'sol_over_comp_figs_mass_weighted/');
        elseif(Eind(k)==8)
            output_dir = strcat(dirname0,'sol_over_comp_figs_volume_weighted/');
        elseif(Eind(k)==7)
            output_dir = strcat(dirname0,'dot_over_tot_figs_mass_weighted/');
        elseif(Eind(k)==6)
            output_dir = strcat(dirname0,'dot_over_tot_figs_volume_weighted/');
        end
        if(i==1)
            mkdir(output_dir)
        end
        output_name = strcat(output_dir,galname,'.jpg');
        x = load(input_name);
        nice_fig(xl,xu,yl,yu,dx,dy,xtit,ytit,galname,12,18,14,0,1,1)
        ydat = x(:,Eind(k));
        Ny = length(ydat);
        if(thin==1)
            temp1 = 0.25 .* ( ydat(1) + ydat(2) + ydat(3) + ydat(4) );
            temp2 = 0.25 .* ( ydat(Ny-3) + ydat(Ny-2) + ydat(Ny-1) + ydat(Ny) );
            temp3 = ( ydat(1) + ydat(2) + ydat(3) ) ./ 3;
            temp4 = ( ydat(Ny-2) + ydat(Ny-1) + ydat(Ny) ) ./ 3;
            ydat(3:Ny-2) = 0.2 .* ( ydat(1:Ny-4) + ydat(2:Ny-3) + ydat(3:Ny-2) + ydat(4:Ny-1) + ydat(5:Ny) );
            ydat(2) = temp1;
            ydat(Ny-1) = temp2;
            ydat(1) = temp3;
            ydat(Ny) = temp4;
        end
        if(thin==1)
            plot(x(:,1),ydat,'linestyle','-','color','r','linewidth',1);
        else
            plot(x(:,1),x(:,Eind(k)),'linestyle','-','color','k','marker','o',...
                'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',3);
        end
        saveas(gcf,output_name);
        close all
        fclose all
    end
end
    
