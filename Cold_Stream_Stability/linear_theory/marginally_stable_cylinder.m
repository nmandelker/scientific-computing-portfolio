function [Phi_marg, X_marg, X_new, Phi_new, F_min]=marginally_stable_cylinder(Del, Mh, m, Nind)

Mc = sqrt(Del)*Mh; % Mach number with respect to internal (cold) medium
Mtot = (sqrt(Del)/(1+sqrt(Del)))*Mh;
tit_mode = strcat('m=',num2str(m,'%1i'));

ind = 1:1:Nind;         % 'n' mode numbers
ind = 2.*ind + m - 0.5; % effective mode number for cylinder with 'm' symmetry mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(Mtot>1)
    0
    maxP = min([1/Mh, 1 - 1/Mc]);
    Phi = linspace(0.5*maxP,maxP,100000);
    N = length(Phi);
    Phi = Phi(1:N-1);
    N = N-1;
    
    A = 1-Mc^2.*((Phi-1).^2);
    B = 1-Mh^2.*(Phi.^2);
    C = ( Phi./(Phi-1) ).^2;
    Z = -(1/Del).*C.*sqrt(A./B);
    A2 = sqrt(Mc^2.*((Phi-1).^2)-1);
    
    K = zeros(1,N);
    Phi_MS = zeros(1,Nind);
    K_MS = zeros(1,Nind);
    for j=1:Nind
        K(1:N) = ( -atan(abs(Z(1:N))) + 0.5.*ind(j)*pi )./A2(1:N);
        for i=2:N-1
            if( (K(i)>max([K(i-1),K(i+1)]) | K(i)<min([K(i-1),K(i+1)])) & K(i)>=-0.0001 )
                Phi_MS(j) = Mh.*Phi(i);
                K_MS(j) = K(i);
            end
        end
    end
    clear A B C Z A2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Look for complex root near marginal stability
    1
    X_marg = [0 K_MS];
    X_new = zeros(1,Nind+1);
    X_new(1) = 0.001;
    for i=2:Nind+1
        X_new(i) = min(X_marg(i) + 0.001, 1.001*X_marg(i));
    end
    X_new = X_new';
    Phi_new = zeros(1,Nind+1);
    F_min = zeros(1,Nind+1);
    if(m==0)
        Phi_marg = [1, Phi_MS./Mh];
        Re_xlin = (X_new(1)./sqrt(2*Del)) .* sqrt( log( 1./(sqrt(abs(Mh^2-1))*X_new(1)) ) );
        Im_xlin = Re_xlin;
    elseif(m>=1)
        Phi_marg = [(Del+1i*sqrt(Del))./(1+Del), Phi_MS./Mh];
        Re_xlin = Del./(1+Del);
        Im_xlin = sqrt(Del)./(1+Del);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fundamental
    2
    Ngrid = 2000;
    [x,y] = meshgrid(linspace(-2*Re_xlin,2*Re_xlin,Ngrid),linspace(0.2*Im_xlin,2*Im_xlin,Ngrid));
    Phi_1 = Phi_marg(1) + x + 1i.*y;
    A = 1-Mc^2.*((Phi_1-1).^2);
    B = 1-Mh^2.*(Phi_1.^2);
    C = ( Phi_1./(Phi_1-1) ).^2;
    Z = -(1/Del).*C.*sqrt(A./B);
    A2 = X_new(1).*A;
    if(m==0)
        numer1 =  besseli(0,A2);
        numer2 = -besselk(1,A2);
        denom1 =  besselk(0,A2);
        denom2 =  besseli(1,A2);
    elseif(m>=1)
        numer1 =  besseli(m,A2);
        numer2 = -besselk(m+1,A2) + (m./(A2)).*besselk(m,A2);
        denom1 =  besselk(m,A2);
        denom2 =  besseli(m+1,A2) + (m./(A2)).*besseli(m,A2);
    end
    F = Z + (numer1.*numer2)./(denom1.*denom2);
    [c,j]=min(min(F));
    Phi_new(1) = Phi_1(j);
    F_min(1) = c;
    clear A B C Z F A2
    
    %reflected
    3
    [x,y] = meshgrid(linspace(0,0.01,Ngrid),linspace(0,0.02,Ngrid));
    for i=2:Nind+1
        i
        Phi_1 = Phi_marg(i) + x + 1i.*y;
        A = 1-Mc^2.*((Phi_1-1).^2);
        B = 1-Mh^2.*(Phi_1.^2);
        C = ( Phi_1./(Phi_1-1) ).^2;
        Z = -(1/Del).*C.*sqrt(A./B);
        A2 = X_new(i).*A;
        if(m==0)
            numer1 =  besseli(0,A2);
            numer2 = -besselk(1,A2);
            denom1 =  besselk(0,A2);
            denom2 =  besseli(1,A2);
        elseif(m>=1)
            numer1 =  besseli(m,A2);
            numer2 = -besselk(m+1,A2) + (m./(A2)).*besselk(m,A2);
            denom1 =  besselk(m,A2);
            denom2 =  besseli(m+1,A2) + (m./(A2)).*besseli(m,A2);
        end
        F = Z + (numer1.*numer2)./(denom1.*denom2);
        [c,j]=min(min(F));
        Phi_new(i) = Phi_1(j);
        F_min(i) = c;
    end
    clear A B C Z F A2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    4
    A = 1-Mc^2.*((Phi_marg-1).^2);
    B = 1-Mh^2.*(Phi_marg.^2);
    C = ( Phi_1./(Phi_marg-1) ).^2;
    Z = -(1/Del).*C.*sqrt(A./B);
    A2 = X_marg.*A;
    if(m==0)
        numer1 =  besseli(0,A2);
        numer2 = -besselk(1,A2);
        denom1 =  besselk(0,A2);
        denom2 =  besseli(1,A2);
    elseif(m>=1)
        numer1 =  besseli(m,A2);
        numer2 = -besselk(m+1,A2) + (m./(A2)).*besselk(m,A2);
        denom1 =  besselk(m,A2);
        denom2 =  besseli(m+1,A2) + (m./(A2)).*besseli(m,A2);
    end
    F = Z + (numer1.*numer2)./(denom1.*denom2);
    max(abs(F))
end