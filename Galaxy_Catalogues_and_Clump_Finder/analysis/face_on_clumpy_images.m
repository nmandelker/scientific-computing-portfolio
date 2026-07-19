function face_on_clumpy_images(gal_name, hh, ss, ch, cs)
% 'gal_name' should be 'VELA_v2_**' or 'VELA**'
% h, s are 0 or 1, if you want to plot zoomed-in Ha or zoomed-in stars respectively. 
% ch and cs are vectors with 3 entries: Lower and upper limits for color bar and tick-step.
% All tracers are plotted face on. 

dirname = strcat('./binary_grid_outputs/',gal_name,'/');
mkdir('./face_on_clumpy_images/');
mkdir(strcat('./face_on_clumpy_images/',gal_name,'/'));

filelist1 = dir(dirname);
ndir = length(filelist1)-2;
filelist = filelist1(3:ndir+2);
cmap = [0 0 0;0 0 0.25;0 0 0.5;0 0 0.75;0 0 1;0 0.07143 0.9929;0 0.1429 0.9857;0 0.2143 0.9786;...
        0 0.2857 0.9714;0 0.3571 0.9643;0 0.4286 0.9571;0 0.5 0.95;0 0.5714 0.9429;0 0.6429 0.9357;...
        0 0.7143 0.9286;0 0.7857 0.9214;0 0.8571 0.9143;0 0.9286 0.9071;0 1 0.9;0 1 0.825;0 1 0.75;...
        0 1 0.675;0 1 0.6;0 1 0.525;0 1 0.45;0 1 0.375;0 1 0.3;0 1 0.225;0 1 0.15;0 1 0.075;0 1 0;...
        0.15 1 0;0.3 1 0;0.45 1 0;0.6 1 0;0.75 1 0;0.9 1 0;0.92 1 0;0.94 1 0;0.96 1 0;0.98 1 0;1 1 0;...
        1 0.9352 0;1 0.8704 0;1 0.8056 0;1 0.7407 0;1 0.6759 0;1 0.6111 0;1 0.584 0;1 0.5568 0;1 0.5296 0;...
        1 0.5025 0;1 0.4753 0;1 0.3565 0;1 0.2377 0;1 0.1188 0;1 0 0;1 0.2 0.2;1 0.4 0.4;1 0.6 0.6;1 0.8 0.8;...
        1 0.8667 0.8667;1 0.9333 0.9333;1 1 1];
for i=1:ndir
    filename = strcat(dirname,filelist(i).name);
    filelist(i).name
    filelist(i).name(27)
    if(strcmp(filelist(i).name(27),'H'))
        if(hh==1)
            fid = fopen(filename, 'rb','ieee-be');
            ntemp = fread(fid,1, 'int')
            n = fread(fid,2, 'int')
            s = fread(fid,2, 'float64')
            Rd = fread(fid,1, 'float32')
            raw = fread(fid, n(1) * n(2), 'float32');
            rtemp = fread(fid,2, 'float32')
            smooth1 = fread(fid, n(1) * n(2), 'float32');
            rtemp = fread(fid,2, 'float32')
            smooth2 = fread(fid, n(1) * n(2), 'float32');
            rtemp = fread(fid,2, 'float32')
            diff = fread(fid, n(1) * n(2), 'float32');
            rtemp = fread(fid,2, 'float32')
            nclump = fread(fid, 1, 'int');
            idclump = fread(fid, nclump, 'int');
            new = fread(fid, nclump, 'int');
            es = fread(fid, nclump, 'int');
            xclump = fread(fid, nclump, 'float64');
            yclump = fread(fid, nclump, 'float64');
            rclump = fread(fid, nclump, 'float64');
            mclump = fread(fid, nclump, 'float64');
            shape_clump = fread(fid, nclump, 'float64');
            fclose(fid);
            clear raw smooth2 diff
        
            rclump(:) = rclump(:) + 2.*s(1)./n(1);
            smooth1=reshape(smooth1,[n(1),n(2)]);
            smooth1=permute(smooth1,[2 1]);
            b = find(smooth1<=0);
            smooth1(b) = 1e-10;
        
            x=linspace(-s(1),s(1),n(1));
            y=linspace(-s(2),s(2),n(2));
            if (mod(s(1),1)>0.93)
                xl1=ceil(s(1));
                yl1=ceil(s(2));
            else
                xl1=floor(s(1));
                yl1=floor(s(2));
            end
            if(xl1<8)
                d1=1;
            elseif(xl1==15 | xl1==18 | xl1==21)
                d1=3;
            elseif(xl1<16)
                d1=2;
            elseif(xl1>=16)
                d1=4;
            end
            cl = [ch(1) ch(2)];
            clab = '${\rm log(}\:\:\:\Sigma_{\rm cold}\:\:{\rm [M_{\odot}\:\:pc^{-2}\:\:]}\:\:\:\:{\rm )}$';
            xlab = '$x\:[kpc]$';
            ylab = '$y\:[kpc]$';
            filename2 = strcat('./face_on_clumpy_images/',gal_name,'/a0_',filename(length(filename)-6:length(filename)-4),'_Ha.jpg');
            filename3 = strcat('./face_on_clumpy_images/',gal_name,'/a0_',filename(length(filename)-6:length(filename)-4),'_Ha.eps');
        
            figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],'Colormap',cmap,'Visible','Off');
            set(figure1,'WindowStyle','docked');

        % Create axes
            axes1 = axes('Parent',figure1,'YTick',[-yl1:d1:yl1],...
            'XTick',[-xl1:d1:xl1],...
            'PlotBoxAspectRatio',[1 1 1],...
            'FontSize',16,...
            'FontName','Times',...
            'CLim',[cl(1) cl(2)]);
             box(axes1,'on');
             axis square
         
        % Create textbox
            aexp = str2num(strcat('0.',filename(length(filename)-6:length(filename)-4)));
            zed = 1./aexp - 1;
            annotation(figure1,'textbox',[0.1675 0.8605 0.11 0.0635],...
                'String',strcat('z=',num2str(zed,'%4.2f')),...
                'FontSize',16,...
                'FontName','Times',...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'BackgroundColor',[1 1 1]);
            ind = str2num(gal_name(length(gal_name)-1:length(gal_name)));
            annotation(figure1,'textbox',[0.8675 0.8605 0.11 0.0635],...
                'String',strcat('V',num2str(ind,'%02i')),...
                'FontSize',16,...
                'FontName','Times',...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'BackgroundColor',[1 1 1]);

            ylabel(ylab,'Interpreter','latex','FontSize',16,'FontName','Times','units','normalized',...
                'position',[-0.05,0.5,0.0]);
            xlabel(xlab,'Interpreter','latex','FontSize',16,'FontName','Times','units','normalized',...
                'position',[0.5,-0.065,0.0]);
            bar=colorbar('peer',axes1,'FontSize',16,'FontName','Times',...
            'CLim',[1 64],'YTick',[cl(1):ch(3):cl(2)]);
            set(get(bar,'Ylabel'),'String',clab,'Interpreter','latex','Fontsize',16,...
            'Rotation',270,'position',[10,(cl(1)+cl(2))/2,9.16],'FontName','Times');
             xlim([-xl1 xl1]);
             ylim([-xl1 xl1]);
%            xlim([-10 10]);
%            ylim([-10 10]);
            hold('all');
            surf(x,y,log10(smooth1),'Parent',axes1,'LineStyle','none');        
            set(gcf,'renderer','zbuffer')
        
             theta = linspace(0,2*pi,100);
%             x1 = Rd.*cos(theta);
%             y1 = Rd.*sin(theta);
%             z1 = 50.*x1./x1;
%             plot3(x1,y1,z1,'linewidth',2,'color','w','linestyle','--');                         

%             for j=1:nclump
%                 if( log10(mclump(j)) >= 7 )
%                     if(log10(mclump(j)) >= 9)
%                         col = 'k';
%                     elseif(log10(mclump(j)) >= 8)
%                         col = 'r';
%                     else
%                         col = 'w';
%                     end
%                     if(new(j)==0)
%                         x1 = xclump(j) + rclump(j).*cos(theta);
%                         y1 = yclump(j) + rclump(j).*sin(theta);
%                         z1 = 50.*x1./x1;
%                         plot3(x1,y1,z1,'linewidth',2,'color',col);
%                     elseif(new(j)==1)
%                         x1 = linspace(xclump(j)-rclump(j),xclump(j)+rclump(j),10);
%                         y1 = (yclump(j)+rclump(j)).*x1./x1;
%                         z1=50.*x1./x1;
%                         x2 = linspace(xclump(j)-rclump(j),xclump(j)+rclump(j),10);
%                         y2 = (yclump(j)-rclump(j)).*x2./x2;
%                         z2=50.*x2./x2;
%                         y3 = linspace(yclump(j)-rclump(j),yclump(j)+rclump(j),10);
%                         x3 = (xclump(j)-rclump(j)).*y3./y3;
%                         z3=50.*x3./x3;
%                         y4 = linspace(yclump(j)-rclump(j),yclump(j)+rclump(j),10);
%                         x4 = (xclump(j)+rclump(j)).*y4./y4;
%                         z4=50.*x4./x4;
%                         plot3(x1,y1,z1,'linewidth',2,'color',col);
%                         plot3(x2,y2,z2,'linewidth',2,'color',col);
%                         plot3(x3,y3,z3,'linewidth',2,'color',col);
%                         plot3(x4,y4,z4,'linewidth',2,'color',col);
%                     end
%                 end
%             end
%             clear xclump yclump rclump spec es idclump nclump smooth1
%             saveas(gcf,filename2);
%             close all
%             fclose all;
        end
    end
%%
    if(strcmp(filelist(i).name(27),'s'))
        if(ss==1)
            fid = fopen(filename, 'rb','ieee-be');
            ntemp = fread(fid,1, 'int')
            n = fread(fid,2, 'int')
            s = fread(fid,2, 'float64')
            Rd = fread(fid,1, 'float32')
            raw = fread(fid, n(1) * n(2), 'float32');
            rtemp = fread(fid,2, 'float32')
            smooth1 = fread(fid, n(1) * n(2), 'float32');
            rtemp = fread(fid,2, 'float32')
            smooth2 = fread(fid, n(1) * n(2), 'float32');
            rtemp = fread(fid,2, 'float32')
            diff = fread(fid, n(1) * n(2), 'float32');
            rtemp = fread(fid,2, 'float32')
            nclump = fread(fid, 1, 'int');
            idclump = fread(fid, nclump, 'int');
            new = fread(fid, nclump, 'int');
            es = fread(fid, nclump, 'int');
            xclump = fread(fid, nclump, 'float64');
            yclump = fread(fid, nclump, 'float64');
            rclump = fread(fid, nclump, 'float64');
            mclump = fread(fid, nclump, 'float64');
            shape_clump = fread(fid, nclump, 'float64');
            fclose(fid);
            clear raw smooth2 diff
        
            rclump(:) = rclump(:) + 2.*s(1)./n(1);
            smooth1=reshape(smooth1,[n(1),n(2)]);
            smooth1=permute(smooth1,[2 1]);
            b = find(smooth1<=0);
            smooth1(b) = 1e-10;
        
            x=linspace(-s(1),s(1),n(1));
            y=linspace(-s(2),s(2),n(2));
            if (mod(s(1),1)>0.93)
                xl1=ceil(s(1));
                yl1=ceil(s(2));
            else
                xl1=floor(s(1));
                yl1=floor(s(2));
            end
            if(xl1<8)
                d1=1;
            elseif(xl1==15 | xl1==18 | xl1==21)
                d1=3;
            elseif(xl1<16)
                d1=2;
            elseif(xl1>=16)
                d1=4;
            end
            cl = [cs(1) cs(2)];
            clab = '${\rm log(}\:\:\:\Sigma_{\rm *}\:\:{\rm [M_{\odot}\:\:pc^{-2}\:\:]}\:\:\:\:{\rm )}$';
            xlab = '$x\:[kpc]$';
            ylab = '$y\:[kpc]$';
            filename2 = strcat('./face_on_clumpy_images/',gal_name,'/a0_',filename(length(filename)-6:length(filename)-4),'_stars.jpg');
        
            figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],'Colormap',cmap,'Visible','Off');
            set(figure1,'WindowStyle','docked');

        % Create axes
            axes1 = axes('Parent',figure1,'YTick',[-yl1:d1:yl1],...
            'XTick',[-xl1:d1:xl1],...
            'PlotBoxAspectRatio',[1 1 1],...
            'FontSize',16,...
            'FontName','Times',...
            'CLim',[cl(1) cl(2)]);
             box(axes1,'on');
             axis square
         
       % Create textbox
            annotation(figure1,'textbox',[0.1675 0.8605 0.11 0.0635],...
                'String',strcat('a=0.',filename(length(filename)-6:length(filename)-4)),...
                'FontSize',16,...
                'FontName','Times',...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'BackgroundColor',[1 1 1]);

            ylabel(ylab,'Interpreter','latex','FontSize',16,'FontName','Times','units','normalized',...
                'position',[-0.05,0.5,0.0]);
            xlabel(xlab,'Interpreter','latex','FontSize',16,'FontName','Times','units','normalized',...
                'position',[0.5,-0.065,0.0]);
            bar=colorbar('peer',axes1,'FontSize',16,'FontName','Times',...
                'CLim',[1 64],'YTick',[cl(1):cs(3):cl(2)]);
            set(get(bar,'Ylabel'),'String',clab,'Interpreter','latex','Fontsize',16,...
                'Rotation',270,'position',[10,(cl(1)+cl(2))/2,9.16],'FontName','Times');
             xlim([-xl1 xl1]);
             ylim([-xl1 xl1]);
%            xlim([-10 10]);
%            ylim([-10 10]);
            hold('all');
            surf(x,y,log10(smooth1),'Parent',axes1,'LineStyle','none');        
            set(gcf,'renderer','zbuffer')
        
            theta = linspace(0,2*pi,100);
%             x1 = Rd.*cos(theta);
%             y1 = Rd.*sin(theta);
%             z1 = 50.*x1./x1;
%             plot3(x1,y1,z1,'linewidth',2,'color','w','linestyle','--');    
            
            for j=1:nclump
                if( log10(mclump(j)) >= 7 )
                    if(log10(mclump(j)) >= 9)
                        col = 'k';
                    elseif(log10(mclump(j)) >= 8)
                        col = 'r';
                    else
                        col = 'w';
                    end
                    if(new(j)==0)
                        x1 = xclump(j) + rclump(j).*cos(theta);
                        y1 = yclump(j) + rclump(j).*sin(theta);
                        z1 = 50.*x1./x1;
                        plot3(x1,y1,z1,'linewidth',2,'color',col);
                    elseif(new(j)==1)
                        x1 = linspace(xclump(j)-rclump(j),xclump(j)+rclump(j),10);
                        y1 = (yclump(j)+rclump(j)).*x1./x1;
                        z1=50.*x1./x1;
                        x2 = linspace(xclump(j)-rclump(j),xclump(j)+rclump(j),10);
                        y2 = (yclump(j)-rclump(j)).*x2./x2;
                        z2=50.*x2./x2;
                        y3 = linspace(yclump(j)-rclump(j),yclump(j)+rclump(j),10);
                        x3 = (xclump(j)-rclump(j)).*y3./y3;
                        z3=50.*x3./x3;
                        y4 = linspace(yclump(j)-rclump(j),yclump(j)+rclump(j),10);
                        x4 = (xclump(j)+rclump(j)).*y4./y4;
                        z4=50.*x4./x4;
                        plot3(x1,y1,z1,'linewidth',2,'color',col);
                        plot3(x2,y2,z2,'linewidth',2,'color',col);
                        plot3(x3,y3,z3,'linewidth',2,'color',col);
                        plot3(x4,y4,z4,'linewidth',2,'color',col);
                    end
                end
            end
            clear xclump yclump rclump spec es idclump nclump smooth1
            saveas(gcf,filename2);
            close all
            fclose all
        end
    end
end