function launch_prop_vs_time(is, nis, ind, gal)

prop = [6 7];
Nprop = length(prop);

% mass(:,1) = [8.5 9.5];
% mass(:,2) = [8.0 8.5];
% mass(:,3) = [7.5 8.0];
% mass(:,4) = [7.0 7.5];

if(gal==7)
%     redshift(:,1) = [0.5 4.0];
%     redshift(:,2) = [0.5 1.0];
%     redshift(:,3) = [1.0 1.5];
%     redshift(:,4) = [1.5 2.0];
%     redshift(:,5) = [2.0 2.5];
%     redshift(:,6) = [2.5 3.0];
%     redshift(:,7) = [3.0 3.5];
%     redshift(:,8) = [3.5 4.0];

%     time(:,1) = [200, 5000];
%     time(:,2) = [100, 200];
%     time(:,3) = [50,  100];

%     mass(:,1) = [8.5 9.5];
%     mass(:,2) = [8.0 8.5];
%     mass(:,3) = [7.5 8.0];
%     mass(:,4) = [7.0 7.5];

    redshift(:,1) = [0.5 4.0];
    time(:,1) = [50, 5000];
    mass(:,1) = [8.3 9.5];
elseif(gal==8)
    redshift(:,1) = [1.0 2.0];
    redshift(:,2) = [1.0 1.5];
    redshift(:,3) = [1.5 2.0];

    time(:,1) = [200, 5000];
    time(:,2) = [100, 200];
    time(:,3) = [50,  100];

    mass(:,1) = [8.5 9.5];
    mass(:,2) = [8.0 8.5];
    mass(:,3) = [7.5 8.0];
    mass(:,4) = [7.0 7.5];
elseif(gal==19)
%     redshift(:,1) = [3.0 6.0];
%     redshift(:,2) = [3.0 3.5];
%     redshift(:,3) = [3.5 4.0];
%     redshift(:,4) = [4.0 4.5];
%     redshift(:,5) = [4.5 5.0];
%     redshift(:,6) = [5.0 5.5];
%     redshift(:,7) = [5.5 6.0];

%     time(:,1) = [100, 200];
%     time(:,2) = [50, 100];

%     mass(:,1) = [8.5 9.5];
%     mass(:,2) = [8.0 8.5];
%     mass(:,3) = [7.5 8.0];
%     mass(:,4) = [7.0 7.5];

    redshift(:,1) = [3.0 6.0];
    time(:,1) = [50, 5000];
    mass(:,1) = [8.0 9.5];
end
Nred = length(redshift(1,:));
Ntime = length(time(1,:));
Nmass = length(mass(1,:));

for i=1:Nprop
    for j=1:Nred
        for k=1:Nmass
            for l=1:Ntime
                prop_vs_time(is, nis, ind, gal, prop(i), time(:,l), mass(:,k), redshift(:,j))
            end
        end
    end
end