function Tv_over_Tk = compressible_slab_growth_times(a_over_Rvir, lambda_over_a, Mh, Del)

Mc = sqrt(Del)*Mh;
Mtot = 1/(1/Mh+1/Mc);
[Mh, Del, Mtot]

Length = 100001;
dirname = strcat('./modes/Mh_',num2str(Mh,'%5.3f'),'_Del_',num2str(Del,'%3i'),'/')
filelist1 = dir(dirname);
filelist = filelist1(3:length(filelist1));
n = length(filelist)
nmode = n/6;

% uncomment to limit to first 20 reflected modes %
%nmode = min([nmode,10]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(2*nmode,Length);
ImP = zeros(2*nmode,Length);
Na = zeros(2*nmode,1);    
for j=1:n
    if(filelist(j).name(1)=='x' | filelist(j).name(5)=='i')
        n1=length(filelist(j).name);
        if(filelist(j).name(n1-7)=='_')
            n2 = str2num(filelist(j).name(n1-6));
        else
            n2 = str2num(filelist(j).name(n1-7:n1-6));
        end
        if(n2 <= nmode)            
            if(strcmp(filelist(j).name(n1-4),'p'))
                k = 2*n2 - 1;
            elseif(strcmp(filelist(j).name(n1-4),'s'))
                k = 2*n2;
            end
            a = load(strcat(dirname,filelist(j).name));
            Na(k) = min(Length,length(a));
            if(filelist(j).name(1)=='x')
                X(k,1:Na(k)) = abs(a(1:Na(k)));
            elseif(filelist(j).name(5)=='i')
                ImP(k,1:Na(k)) = abs(a(1:Na(k)));
            end            
        end
    end
end

growth_rates = zeros(1,2*nmode);
for i=1:2*nmode
    if(Na(i)>0)
        b = find(abs( lambda_over_a - 2*pi./X(i,1:Na(i)) ) == min(abs( lambda_over_a - 2*pi./X(i,1:Na(i)) )), 1, 'first');
        growth_rates(i) = ImP(i,b);
    end
end
max_growth_rate = max(growth_rates(1:2*nmode));
%Tkh_over_Tvir = (1/(2*pi*max_growth_rate)) * lambda_over_a * a_over_Rvir;
Tv_over_Tk = (2*pi*max_growth_rate) / (lambda_over_a) / (a_over_Rvir);