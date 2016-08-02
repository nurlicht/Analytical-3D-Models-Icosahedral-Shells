%
%
% Beginning of IcshdrSpctrm documentation
%
%
%
%% Tasks
% A. Coding icosahedral spectrum calculation over the theta-phi grid
%  * Icosahedral spectrum of a numerical solid icosahedron
%  * Higher-order icosahedral harmonics
%  * Icosahedral spectrum calculation with shell-related integrands
%  * Icosahedral spectrum (analytical) of a solid icosahedron
%  * Icosahedral spectrum (analytical) of an icosahedral shell

% B. Calculation of icosahedral spectrum
%  * Numerical (over theta-phi grid)
%  * Numerical (over modified coordinate w/ Jacobian factor)
%  * Analytical (over Theta)
%  * Analytical (over Phi)

%% Log
%  *_14.m:
%   f_l(r) 2D spectra
%
%  *_13.m:
%
%  *_12.m:
%  Seemingly successful implementation and monitoring of the progress.
%  Oversampling is critical, and the resolution is to be increased.
%  Possible nonuniform scaling of resolution in the spherical coordinate
%  corresponding to a uniform isotropic resolution in the Cartesian
%  coordinate (efficient tradeoff between accuracy/cost of computation)

%  *_11.m:
%  Implementation of the icosahedral grid with support (and a function!)
%  Meaningful spectra and relevance of different radii, but still issues

%  *_10.m:
%  Problem in Make_IcoVolFunc; inappropriate resampling of an icosahedron
%  on spherical shells AND/OR inappropriate "icosahedron" to begin with.
%  Examples of data/sample radii that do not make sense:
%  s:0.51232, d:0.6325 to 1.96, n(NaN):0
%  s:0.54302, d:0.6325 to 1.96, n(NaN):0
%  s:0.57462, d:0.6325 to 1.96, n(NaN):0
%  s:0.6071, d:0.6325 to 1.96, n(NaN):0
%  s:0.64049, d:0.6325 to 1.96, n(NaN):0
%  s:0.67476, d:0.6325 to 1.96, n(NaN):0
%  s:0.70993, d:0.6325 to 1.96, n(NaN):0
%  s:0.74599, d:0.6325 to 1.96, n(NaN):0
%  s:0.78294, d:0.6325 to 1.96, n(NaN):0
%  s:0.82079, d:0.6325 to 1.96, n(NaN):0
%  s:0.85953, d:0.6325 to 1.96, n(NaN):0
%  s:0.89916, d:0.6325 to 1.96, n(NaN):0
%  s:0.93969, d:0.6325 to 1.96, n(NaN):0
%  s:0.98111, d:0.6325 to 1.96, n(NaN):0
%  s:1.0234, d:0.6325 to 1.96, n(NaN):0
%  s:1.0666, d:0.6325 to 1.96, n(NaN):0
%  s:1.1107, d:0.6325 to 1.96, n(NaN):0
%  s:1.1557, d:0.6325 to 1.96, n(NaN):0
%  s:1.2016, d:0.6325 to 1.96, n(NaN):0
%  s:1.2484, d:0.6325 to 1.96, n(NaN):239
%  s:1.2961, d:0.6325 to 1.96, n(NaN):1251


%  *_9.m:
%  Converted to 1+2D model, but no room for 3D interpolation!
%
%  *_8.m:
%  Solid icosahedral function fed to the Icos spectrum; white spectra :(
%
%  *_7.m:
%  Resampling icosahedral volume on the spherical coordinate grid
%  (synthesizing coordinate-aligned function)
%
%  *_6.m:
%  Implementation of icosahedral shell/volume grids
%  Determination of 20 S2/20 rotation matrices
%
%  *_5.m:
%  Possible use of http://mathworld.wolfram.com/IcosahedralEquation.html
%  
%  *_4.m:
%  Test routine with arbitrary (symmetric) function implemented; odd!!
%  1/20 flat triangular region was used for spectral analysis
%
%  *_3.m:
%  Modularized
%
%  *_2.m:
%  Significant improvement by integrating over the entire S2 using the same
%  
%
%  *_1.m:
%  Acceptable Ico. spectrum of <IH_1,IH_2>
%  No specific improvement of the spectrum by use of quadrature int.
%  
%
%
% End of IcshdrSpctrm documentation
%
%
%% Original theta-phi formulation
%15-6: Implementation of the correct parameterization; meaningful yet
%questionable results; possibilities of numerical errors, anaytical
%singularities, ...


%% r'-Tau formulation
% Up to xxx14.m
%14:
%13: Specific values of range-resolution for both u and r implemented
%12_1: Documentation of direct (r',tau) parameterization
%12:   Integrated Legendre quadrature calculation and integration
%11_4: Integrated Legendre quadrature calculation and integration

%11_2: spherical integration on 1/20 of the sphere; different densities of
%theta,phi and r


%10: Best so far; complete retrieval of f_r and f_u
%9: A very good one!


%7: 2D Fourier transform of a triangle implemented



%Quadrature on a sphere
%http://link.springer.com/content/pdf/10.1007%2FBF00966954.pdf

%Range and sign of Tau (and Eta/w) to be controlled
%Correct v isualization possible with a change of sign (of the 'right' value)

%First implementation of integration of an icosahedral shell on a
%non-overlapping spherical domain; new domain of integration: INTERSECTION
%of the original integration (spehrical shell) domain and the object
%(icosahedral shell) domain

%Intersection domain parameterized with (r',Tau)
%r'=sqrt(r^2-cos(Alpha)^2), r: spherical coordinate variable
%Tau: Subset of [0,2*pi] with C3 symmetry

%r_Max=sin(Alpha)
%r_min=sqrt(cos(Alpha)^2+(r_Max/2)^2)

%r < cos(Alpha) --> Integral = 0
%0 < r' < r_min --> Tau: [pi/6,5*pi/6]
%r_min < r' < r_Max --> Tau: [pi/6+acos(r_min/r'),5*pi/6-acos(r_min/r')]
%r > 1 --> Integral = 0

%f_ico(r,theta,phi)=[N*q/(A*dt)]*[u(r-r_ico(theta,phi))-u(r-r_ico(theta,phi
%)-dt)]=[N*q/(A*dt)]*delta(r-r_ico(theta,phi))

%Possible improvement using dt(theta,phi)

%Also trying with the pulse function, not the delta function

%Integral=I(r)=[N*q/(A*dt)]*Int{delta(r-r_ico)*Jl*sin*Delta*dTau*dr'},
%where Int{delta(r-sqrt(r'^2+cos(Alpha)^2))}=r/sqrt(r^2-cos(Alpha)^2), so
%I(r)=[N*q/(A*dt)]*Int{Jl*sin*Delta*dTau} over [Tau_min(r) Tau_Max(r)]


%Integral = Int(f_ico(r_vector)*J_l(theta,phi)*sin(theta) d_theta d_phi) =
%Int(f_ico(r_vector)*J_l(theta,phi)*sin(theta) Delta d_r' d_Tau), where
%Delta is the determinant of the 2x2 Jacobian:

%Delta_11= d[Theta/r'] = [r'*cos(Theta)-0.6071*r*sin(Tau)]/(r^2*sin(Theta))

%Delta_12= d[Theta/Tau] = -0.6071*(r'/r)*cos(Tau)/sin(Theta)

%Delta_21= d[Phi/r'] = [sin(Theta)*sin(Phi)*cos(Alpha)^2/(r*r') -0.2836/r'
%-r*sin(Phi)*cos(Theta)*Delta_11] / [r*sin(Theta)*cos(Phi)]

%Delta_22= d[Phi/Tau] = [0.2836-r'*(0.8090*sin(Tau)+0.4671*cos(Tau))
%-r*sin(Phi)*cos(Theta)*Delta_11] / [r*sin(Theta)*cos(Phi)]

%Function value of icosahedron:
%Delta function over the (Eta,W) domain:
%f(r_vector) = delta(r - r_ico(theta,phi))
%IMPORTANT: This delta function should be appropriately converted to scaled
%delta function(s) by considering the Jacobian elements




%New functions
%f_ico=Sum{f 1/20 (R,r)} or Sum{f 1/60 (R,r)}

%f_T=Sum{ Pentamers + Hexamers}

%% Implemented functions


%%
function []=IndirectBessel()
    Initialize();
    IcshFaceMode=1; %1/0: Phi_min=0/-Omega/2
    [B,C,A,F,Omega,Theta_c,Alpha,Beta,Delta,Theta_F,Phi_min,Phi_Max]=...
        Load_IH_Tri(IcshFaceMode);
    
    
    %PlotIcshdrnFace(B,C,A);
    %return
    %hold on;
    
    %Integral grid (1)
    OS_Factor=2.5;  %2.5

    %Object (icosahedron) grid
    Radius1=0.8;
    Radius2=1;    %1/cos(Alpha)=1.2584
    R_Norm_Flag=1;
    Radius=Radius1;
    N_TriGrid=40;   %20
    Radius1_min=Radius1*cos(Alpha);
    Radius2_min=Radius2*cos(Alpha);
    Radius1_Max=Radius1*sqrt(cos(Alpha)^2+(sin(Alpha)/2)^2);
    Radius2_Max=Radius2*sqrt(cos(Alpha)^2+(sin(Alpha)/2)^2);
    if isequal(Radius1,Radius2)
        N_Radius=1;
        N_R=ceil(2*OS_Factor*N_Radius); %Integral grid (2)
        Radius_Res=(Radius2-Radius1_min)/(N_R-1);
    else
        N_Radius=2*ceil(N_TriGrid*(2/(sqrt(3)*tan(Alpha)))*...
            (1-Radius1_min/Radius2)/(1+Radius1_min/Radius2));
        N_R=ceil(2*OS_Factor*N_Radius); %Integral grid (2)
        Radius_Res=(Radius2-Radius1_min)/(N_R-1);
    end
    Radius_Res=(1+cos(Alpha))/(8*OS_Factor*N_TriGrid/sqrt(3)-1/tan(Alpha/2)^2);
    N_R=1+ceil((Radius2-Radius1_min)/Radius_Res);
    
    Plot_Grid_Flag=1;
    AxisOff_Flag=1;
    AppendSupport_Flag=1;
    IcoObjParam=struct('N_TriGrid',N_TriGrid,...
        'Radius_Res',Radius_Res,...
        'N_Radius',N_Radius,...
        'Radius1',Radius1,...
        'Radius2',Radius2,...
        'Radius1_min',Radius1_min,...
        'Radius2_min',Radius2_min,...
        'Radius',Radius,...
        'Plot_Grid_Flag',Plot_Grid_Flag,...
        'AxisOff_Flag',AxisOff_Flag,...
        'AppendSupport_Flag',AppendSupport_Flag,...
        'A',A,...
        'B',B,...
        'C',C,...
        'F',F);
    disp(['Radius1_min=' num2str(Radius1_min) ', Radius1_Max=' num2str(Radius1_Max) ', Radius1=' num2str(Radius1)]);
    disp(['Radius2_min=' num2str(Radius2_min) ', Radius2_Max=' num2str(Radius2_Max) ', Radius2=' num2str(Radius2)]);
    disp(['N_Radius=' num2str(N_Radius) ', Radius_Res=' num2str(Radius_Res)]);
    [IcoGridVol,IcoGrid,H0]=Form_IcosahedralGrids(IcoObjParam);
    disp(['Min IcoGridVol.R=' num2str(min(sqrt(IcoGridVol.x.^2+IcoGridVol.y.^2+IcoGridVol.z.^2)))]);
    
    %Integral grid (3)
    N_Theta=ceil(2*OS_Factor*N_TriGrid);
    N_Phi=N_Theta+1;

    %Object (icosahedral function f(r))
    [IcoFunc,FuncGrid]=Make_IcoVolFunc(IcoGridVol,N_R,N_Theta,N_Phi);
    
    
    
    disp(['Min FuncGrid.R=' num2str(min(FuncGrid.R_1D))]);
    Theta=FuncGrid.Theta_2D;
    Phi=FuncGrid.Phi_2D;
    R_1D=FuncGrid.R_1D;
    Rm=min(R_1D);
    RM=max(R_1D);

    %Indices of critical distances
    [~,Index]=min(abs(R_1D-Radius1_min));
    Radius1_min_Index=Index(1);
    [~,Index]=min(abs(R_1D-Radius2_min));
    Radius2_min_Index=Index(1);
    [~,Index]=min(abs(R_1D-Radius1_Max));
    Radius1_Max_Index=Index(1);
    [~,Index]=min(abs(R_1D-Radius2_Max));
    Radius2_Max_Index=Index(1);
    [~,Index]=min(abs(R_1D-Radius1));
    Radius1_Index=Index(1);
    [~,Index]=min(abs(R_1D-Radius2));
    Radius2_Index=Index(1);
    

    Integral_Mode=(min(size(Theta)) > 1);

    %Basis functions
    LMa=Load_LMa;
    l_Vector=unique(LMa(:,1));
    N_l=numel(l_Vector);
    %N_l=4;l_Vector=l_Vector(1:N_l);
    
    %Display and test of basis functions
    N_subplot=ceil(sqrt(N_R));
    N_subplot_l=ceil(sqrt(N_l));
    Plot_Flag=0;
    if Plot_Flag
        H=figure;
    else
        H=[];
    end
    [IH_Basis,ArrayAngleFlag]=Make_IHBasis(Theta,Phi,l_Vector,H);
    if Plot_Flag
        Test_IHBasis(H);
        %H2=figure;while 1;Test_IHBasis_Rand(H2);pause;end
        H2=figure;Test_IHBasis_Rand(H2);
    end
    
    Plot_RadialScan=0;
    if Plot_RadialScan
        x=sin(Theta).*cos(Phi);
        y=sin(Theta).*sin(Phi);
        z=cos(Theta);
        figure;
        for R_cntr=1:N_R
            scatter3(x,y,z,20,IcoFunc{R_cntr});daspect([1 1 1]);axis tight; axis off;title(['f(r_{' num2str(r_cntr) '}) @r=' num2str(R_Shell) ', d:' num2str(C) ' to ' num2str(D) ', n(NaN):' num2str(E)]);colorbar;view(-90,32);alpha(0.5);alpha('color');drawnow;
        end
    end
    
    save IcoFunc.mat IcoFunc;
    %return
            

    H3=figure;
    H4=figure;
    H5=figure;
    ISpectrumMtx=nan(N_l,N_R);
    for R_cntr=1:N_R
        %fr=squeeze(IcoFunc(R_cntr,:,:));
        fr=IcoFunc{R_cntr};
        ISpectrum_=Make_ISpectrum(fr);
        if R_Norm_Flag
            ISpectrum_=ISpectrum_/(Radius2-Radius1);
        end
        ISpectrumMtx(:,R_cntr)=ISpectrum_;
        F_Recon=SynthIcoFunc(ISpectrum_);

        figure(H3);
        subplot(N_subplot,N_subplot,R_cntr);
        bar(l_Vector(2:end),ISpectrum_(2:end),'b');
        xlim([min(l_Vector)-1 max(l_Vector)+1]);
        if isequal(R_cntr,N_R)
            title(['R_{' num2str(R_cntr) '}=' num2str(R_1D(R_cntr)) '(' num2str(Rm) '<R<' num2str(RM)]);
        else
            title(['R_{' num2str(R_cntr) '}=' num2str(R_1D(R_cntr))]);
        end
        drawnow;

        figure(H4);
        subplot(N_subplot,N_subplot,R_cntr);
        plot(fr(:),F_Recon(:),'b.');axis tight;
        if isequal(R_cntr,N_R)
            title(['R_{' num2str(R_cntr) '}=' num2str(R_1D(R_cntr)) '(' num2str(Rm) '<R<' num2str(RM)]);
        else
            title(['R_{' num2str(R_cntr) '}=' num2str(R_1D(R_cntr))]);
        end
        drawnow;
        
        figure(H5);
        for l_cntr=1:N_l
            subplot(N_subplot_l,N_subplot_l,l_cntr);
            Tempx=ISpectrumMtx(l_cntr,:);
            plot(R_1D,Tempx,'.--k');
            hold on;
            plot(   R_1D(Radius1_min_Index),Tempx(Radius1_min_Index),'b*',...
                    R_1D(Radius1_Max_Index),Tempx(Radius1_Max_Index),'b*',...
                    R_1D(Radius1_Index),Tempx(Radius1_Index),'b*',...
                    R_1D(Radius2_min_Index),Tempx(Radius2_min_Index),'ro',...
                    R_1D(Radius2_Max_Index),Tempx(Radius2_Max_Index),'ro',...
                    R_1D(Radius2_Index),Tempx(Radius2_Index),'ro');
            hold off;
            title(['l=' num2str(l_Vector(l_cntr))]);
            axis tight;
            xlabel('Radial position, R');
        end
        drawnow;
    end
    
    H6=figure;
    subplot(2,1,1);
    surf(R_1D,l_Vector,ISpectrumMtx);
    axis tight; xlabel('Radial Position');ylabel('Icosahedral order');colorbar;shading interp
    subplot(2,2,3);
    surf(R_1D,l_Vector(2:end),ISpectrumMtx(2:end,:));
    axis tight; xlabel('Radial Position');ylabel('Icosahedral order');colorbar;shading interp
    subplot(2,2,4);
    surf(R_1D,l_Vector(3:end),ISpectrumMtx(3:end,:));
    axis tight; xlabel('Radial Position');ylabel('Icosahedral order');colorbar;shading interp
    
    save All.mat;
    return;
    
        
            
            
    
    
    
    
    
    
    
    
    
    
    
    r_min=sqrt(cos(Alpha)^2+0*(sin(Alpha)/2)^2)+1e-16;
    r_boundary=sqrt(cos(Alpha)^2+1*(sin(Alpha)/2)^2);
    %r_Max=sqrt(cos(Alpha)^2+1*(sin(Alpha)/2)^2);
    r_Max=1-1e-20;
    
    
    d=zeros(10,1);
    d(1:4)=[sqrt(3)*tan(Omega/2)
        sin(Alpha)*(1+2*cos(Theta_c))
        2*(1-cos(Theta_c))
        2*sqrt(3)*sin(Theta_c)*sin(Omega/2)];
    d(5)=d(3)/d(4);
    d(6)=d(2)/d(4);
    
    
    N=201;
    Tau_1D=linspace(0,2*pi,N+1);
    Tau_1D=Tau_1D(1:N);
    Tau_1D=Tau_1D-mean(Tau_1D);
    r_1D=sqrt(cos(Alpha)^2+linspace(0,sin(Alpha),N).^2);
    r_1D(1)=cos(Alpha)+1e-8;
    r_1D(end)=1-1e-8;
    dTau_1D=([0 diff(Tau_1D)]+[diff(Tau_1D) 0])/2;
    [r,Tau]=meshgrid(r_1D,Tau_1D); %f(theta,r)
    [~,dTau]=meshgrid(r_1D,dTau_1D); %f(theta,r)
    
    
    Rho=real(sqrt(r+cos(Alpha)).*sqrt(r-cos(Alpha)));
    Correct_Small_rp_Flag=0;
    if Correct_Small_rp_Flag
        Threshold=1e-8;
        Index=(Rho < Threshold);
        Epsilon=r(Index)-cos(Alpha);
        Eta=Epsilon/(2*cos(Alpha));
        Rho(Index)=(sqrt(2*cos(Alpha))*sqrt(Epsilon)).*...
            (1+(1/2)*Eta-(1/8)*Eta.^2);
    end

    
    Phi=atan(d(1)*Rho.*cos(Tau)./(sin(Alpha)-Rho.*sin(Tau)));
    Theta=pi/2-atan((d(5)*tan(Tau)+d(6)./(Rho.*cos(Tau))).*sin(Phi));
    
    
    x=r.*sin(Theta).*cos(Phi+Omega/2);
    y=r.*sin(Theta).*sin(Phi+Omega/2);
    z=r.*cos(Theta);

    %Index=(Theta_ >= 0) & (Theta_ <= Theta_c);
    
    C_EMAN=[sin(Theta_c);0;cos(Theta_c)];    %Vertex
    A_EMAN=[sin(Theta_c)*cos(Omega);sin(Theta_c)*sin(Omega);cos(Theta_c)];   %Vertex
    Index=(x > 0) & (x < C_EMAN(1)) & (y > 0) & (y < A_EMAN(2)) & (z > C_EMAN(3)) ...
        & (z < B(3)) & (z < 1-(tan(Theta_c/2)/sin(Omega))*y) ...
        & (imag(Phi) == 0) & (Theta.*Phi ~= 0);
    
    %Phi=Phi(Index);
    %Theta=Theta(Index);
    %Rho=Rho(Index);
    %r=r(Index);
    %Tau=Tau(Index);
    dTau(~Index)=0;
    
    x=r.*sin(Theta).*cos(Phi);
    y=r.*sin(Theta).*sin(Phi);
    z=r.*cos(Theta);
    x(~Index)=nan;
    y(~Index)=nan;
    z(~Index)=nan;
    
    scatter3(x(:),y(:),z(:),40,ones(size(z(:))),'.');
    hold on
    PlotIcshdrnFace(B,C,A)
    hold off
    
    Delta_11=-(sin(Theta).^2.*sin(Phi)./(Rho.^2.*cos(Theta))).*(-d(6)+(sin(Alpha)/d(1))./(tan(Theta).*cos(Phi).^3));
    Delta_12=-((d(5)*Rho+d(6)*sin(Tau))+(Rho-sin(Alpha)*sin(Tau))./(d(1)*tan(Theta).*cos(Phi).^3)).*(sin(Theta).^2.*sin(Phi))./(Rho.*cos(Tau).^2);
    Delta_21=(sin(Alpha)/d(1)).*sin(Phi).^2./(Rho.^2.*cos(Tau).*cos(Phi).^4);
    Delta_22=sin(Phi).^2.*(Rho-sin(Alpha).*sin(Tau))./(d(1)*Rho.*cos(Tau).^2.*cos(Phi).^4);
    Delta=Delta_11.*Delta_22-Delta_12.*Delta_21;
    
    %min(Delta_11(:)),max(Delta_11(:)),return
    
    H6=figure;
    subplot(2,4,1)
    imagesc(Theta);
    daspect([1 1 1]);colorbar;xlabel('r');ylabel('\tau');axis tight;
    title('Theta');
    subplot(2,4,2)
    imagesc(Phi);
    daspect([1 1 1]);colorbar;xlabel('r');ylabel('\tau');axis tight;
    title('Phi');
    subplot(2,4,3)
    surf(log(Delta-min(Delta(:))+1));
    colorbar;xlabel('r');ylabel('\tau');axis tight;
    title('Delta');
    view(0,90);shading interp;
    subplot(2,4,5)
    surf(Delta_11);
    colorbar;xlabel('r');ylabel('\tau');axis tight;
    title('Delta_{11}');
    view(0,90);shading interp;
    subplot(2,4,6)
    surf(log(Delta_12-min(Delta_12(:))+1));
    colorbar;xlabel('r');ylabel('\tau');axis tight;
    title('Delta_{12}');
    view(0,90);shading interp;
    subplot(2,4,7)
    surf(log(Delta_21-min(Delta_21(:))+1));
    colorbar;xlabel('r');ylabel('\tau');axis tight;
    title('Delta_{21}');
    view(0,90);shading interp;
    subplot(2,4,8)
    surf(log(Delta_22-min(Delta_22(:))+1));
    colorbar;xlabel('r');ylabel('\tau');axis tight;
    title('Delta_{22}');
    view(0,90);shading interp;
    
    drawnow;
    
    %return
    H2=figure;
    H3=figure;
    H4=figure;
    H5=figure;
    N_subplot=ceil(sqrt(N_l));
    Recon=zeros(size(Theta));
    Energy=nan(N_l,1);
    for cntr=1:N_l
        Integrand0=IHrmncSngl(Theta,Phi,l_Vector(cntr));

        figure(H2);
        subplot(N_subplot,N_subplot,cntr)
        imagesc(log(Integrand0-min(Integrand0(:))+1));
        daspect([1 1 1]);colorbar;xlabel('r');ylabel('\tau');axis tight;
        title(['l=' num2str(l_Vector(cntr))]);
        drawnow;

        Integrand=Integrand0.*(20*r./Rho).*(sin(Theta)).*Delta;
        figure(H3);
        subplot(N_subplot,N_subplot,cntr)
        imagesc(log(Integrand0-min(Integrand0(:))+1));
        daspect([1 1 1]);colorbar;xlabel('r');ylabel('\tau');axis tight;
        title(['l=' num2str(l_Vector(cntr))]);
        drawnow;

        figure(H4);
        subplot(N_subplot,N_subplot,cntr)
        R_Profile=sum(Integrand.*dTau,1);
        plot(r(1,:),R_Profile);
        xlabel('r');ylabel(['f_{l=' num2str(l_Vector(cntr)) '}(r)']);axis tight;
        drawnow;
        
        Recon=Recon+Integrand0.*repmat(R_Profile,[N 1]);
        Energy(cntr)=R_Profile*R_Profile';
        
     
        %x=r.*sin(Theta).*cos(Phi);
        %y=r.*sin(Theta).*sin(Phi);
        %z=r.*cos(Theta);
        %x(~Index)=0;
        %y(~Index)=0;
        %z(~Index)=0;
        figure(H5);
        subplot(1,2,1)
        scatter3(x(:),y(:),z(:),20,Recon(:),'.');
        hold on
        plot3(F(1),F(2),F(3),'r*')
        plot3([0 F(1)],[0 F(2)],[0 F(3)],'r');
        plot3(A(1),A(2),A(3),'k*')
        plot3(B(1),B(2),B(3),'k*')
        plot3(C(1),C(2),C(3),'k*')
        plot3([B(1) A(1) C(1) B(1)],[B(2) A(2) C(2) B(2)],[B(3) A(3) C(3) B(3)],'k')
        axis tight;daspect([1 1 1]);xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
        title(['Normalized standard deviation, \sigma/\eta=' ...
            num2str(std(Recon(:))/mean(Recon(:)))]);
        view(148,32);
        hold off
        subplot(1,2,2)
        bar(l_Vector,Energy);
        ylabel('Icosahedral spectrum');
        xlabel('Icosahedral order');
        drawnow;
    end
    
return


    
    Theta_MaxPhi=atan(1./(1.6180.*cos(Omega/2-Phi)));
    CotGamma=Beta*cos(Phi-Omega/2);
    %CotGamma=Beta*cos(Phi);
    Gamma=atan(1./CotGamma);
    SinGammaOverR=sin(Gamma)./r;
    Theta_0=asin(sin(Gamma)./r)-Gamma;
    %Theta_0=(pi-asin(sin(Gamma)./r))-Gamma;
    Index_NoPhi=(imag(Theta_0) ~= 0) | (Theta_0 < 0);
    %Index_LowQ=(sum(Index_NoPhi,1) < N/10);
    %Index_NoPhi(Index_LowQ,:)=Index_NoPhi(Index_LowQ,:)-Index_NoPhi(Index_LowQ,:);
    Theta_0(Index_NoPhi)=0;
    %Theta_0(Theta_0 > Theta_MaxPhi)=0;
    %Phi_Max=-Omega/2+acos(sqrt(1./r.^2-1)*Beta);
    Phi_Max=acos(sqrt(1./r.^2-1)*Beta);
    
    
    dPhi=abs(Phi_1D(2)-Phi_1D(1));
    subplot(4,4,1)
    imagesc(Theta_MaxPhi)
    title('Theta_{MaxPhi}')
    daspect([1 1 1]);colorbar;xlabel('r');ylabel('\phi');axis tight;
    
    subplot(4,4,2)
    imagesc(CotGamma)
    title('CotGamma')
    daspect([1 1 1]);colorbar;xlabel('r');ylabel('\phi');axis tight;
    
    subplot(4,4,3)
    imagesc(Gamma)
    title('Gamma')
    daspect([1 1 1]);colorbar;xlabel('r');ylabel('\phi');axis tight;
    
    subplot(4,4,4)
    imagesc(SinGammaOverR)
    title(['SinGammaOverR; Max=' num2str(max(SinGammaOverR(:)))]);
    daspect([1 1 1]);colorbar;xlabel('r');ylabel('\phi');axis tight;

    subplot(4,4,5)
    imagesc(Theta_0)
    title('Theta_0')
    daspect([1 1 1]);colorbar;xlabel('r');ylabel('\phi');axis tight;

    subplot(4,4,6)
    x=r.*sin(Theta_0).*cos(Phi);
    y=r.*sin(Theta_0).*sin(Phi);
    z=r.*cos(Theta_0);
    scatter3(x(:),y(:),z(:),20,ones(size(z(:))),'.');
    hold on
    plot3(F(1),F(2),F(3),'r*')
    plot3([0 F(1)],[0 F(2)],[0 F(3)],'r');
    plot3(A(1),A(2),A(3),'k*')
    plot3(B(1),B(2),B(3),'k*')
    plot3(C(1),C(2),C(3),'k*')
    plot3([B(1) A(1) C(1) B(1)],[B(2) A(2) C(2) B(2)],[B(3) A(3) C(3) B(3)],'k')
    axis tight;daspect([1 1 1]);xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
    view(148,32);
    hold off
    daspect([1 1 1]);colorbar;xlabel('r');ylabel('\phi');axis tight;
    min(Theta_0(:)),max(Theta_0(:))
    min(Phi(:)),max(Phi(:))
    Theta_c
    Omega

    subplot(4,4,7)
    imagesc(Phi_Max)
    title('Phi_{Max}')
    daspect([1 1 1]);colorbar;xlabel('r');ylabel('\phi');axis tight;
return
    

    
    return;
    
    
    subplot(4,4,1)
    plot(Phi*180/pi,Theta_MaxPhi*180/pi);
    xlabel('Phi (\circ)');
    ylabel('Theta_{MaxPhi} (\circ)');
    axis tight;
    
    subplot(4,4,2)
    plot(Phi*180/pi,CotGamma);
    xlabel('Phi (\circ)');
    ylabel('CotGamma');
    %hold on
    %plot([min(Phi) max(Phi)]*180/pi,sqrt(1/r_min^2-1)*[1 1],'r');
    %plot([min(Phi) max(Phi)]*180/pi,sqrt(1/r_Max^2-1)*[1 1],'k');
    %hold off
    axis tight;
    ylim([0 1])
    
    subplot(4,4,3)
    plot(Phi*180/pi,1./sqrt(1+CotGamma.^2));
    xlabel('Phi (\circ)');
    ylabel('SinGamma');
    hold on
    hold off
    axis tight;
    

    subplot(4,4,4)
    plot(Phi*180/pi,Gamma*180/pi);
    xlabel('Phi (\circ)');
    ylabel('Gamma');
    axis tight;
    
    subplot(4,4,5)
    plot(Phi*180/pi,sin(Gamma)/r_min);
    xlabel('Phi (\circ)');
    ylabel('sin(Gamma)/r_{min}');
    axis tight;
    
    subplot(4,4,6)
    plot(Phi*180/pi,sin(Gamma)/r_Max);
    xlabel('Phi (\circ)');
    ylabel('sin(Gamma)/r_{Max}');
    axis tight;
    
    
    return
    
    subplot(4,4,6)
    plot(Phi*180/pi,Sin_Theta0PlusGamma_rmin);
    xlabel('Phi (\circ)');
    ylabel('Sin_Theta0PlusGamma_rmin');
    axis tight;
    
    subplot(4,4,7)
    plot(Phi*180/pi,Sin_Theta0PlusGamma_rMax);
    xlabel('Phi (\circ)');
    ylabel('Sin_Theta0PlusGamma_rMax');
    axis tight;
    
    subplot(4,4,8)
    plot(Phi*180/pi,Theta_0_rmin*180/pi);
    xlabel('Phi (\circ)');
    ylabel('Theta_0rmin');
    axis tight;
    
    subplot(4,4,9)
    plot(Phi*180/pi,Theta_0_rMax*180/pi);
    xlabel('Phi (\circ)');
    ylabel('Theta_0rMax');
    axis tight;
    
    
    return;
    
    
    

    %Theta and Phi angles of a 1/60 unit (ABF=1/3 of ABC)
    Theta_1D=linspace(0,Theta_c,N);
    Phi_1D=linspace(0,pi/5,N);
    [Theta,Phi]=meshgrid(Theta_1D,Phi_1D);
    Index=(1 >= tan(Theta)*1.6180.*cos(Omega-Phi));
    Theta=Theta(Index);
    Phi=Phi(Index);

    %Equation of a 1/60 icosahedron: r=f(theta,phi)
    r=1./(cos(Theta)+Beta*cos(Phi-Omega/2).*sin(Theta));

    R0=single(Rotate_Tri_to_EMAN());
    R=IcosSymRot();
    

    %____________________
    %theta-phi formulation
    N=100;
    Theta_F=real(acos(F(3)/norm(F)));
    Theta_1D=linspace(0,Theta_c,N);
    Phi_1D=linspace(0,pi/5,N);
    [Theta,Phi]=meshgrid(Theta_1D,Phi_1D);
    Index=(1 >= tan(Theta)*1.6180.*cos(Omega-Phi));
    Theta=Theta(Index);
    Phi=Phi(Index);
    Gamma=(pi/2)-atan(Beta*cos(Phi-Omega/2));
    %r=1./(cos(Theta)+Beta*cos(Phi-Omega/2).*sin(Theta));
    [r,~]=meshgrid(linspace(cos(Alpha),1,N),Phi_1D);
    r=r(Index);
    Theta_DeltaRoot=real((asin(sin(Gamma)./r)-Gamma));
    plot(Phi(:),Theta_DeltaRoot(:),'.');
    return;
    
    [x,y,z]=Polar2Cartesian(r,Theta,Phi);
    hold on
    plot3(F(1),F(2),F(3),'r*')
    plot3([0 F(1)],[0 F(2)],[0 F(3)],'r');
    plot3(A(1),A(2),A(3),'k*')
    plot3(B(1),B(2),B(3),'k*')
    plot3(C(1),C(2),C(3),'k*')
    plot3([B(1) A(1) C(1) B(1)],[B(2) A(2) C(2) B(2)],[B(3) A(3) C(3) B(3)],'k')
    plot3(x,y,z,'b.');
    axis tight;
    daspect([1 1 1]);
    xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
    %view(108,20)
    view(148,32);
    hold off
    
    return;
    
    
    


    
    %____________________
    %r'-Tau formulation
    
    
    r_sampling=0.004;
    rp_min=sin(Alpha)/2;
    rp_Max=sin(Alpha);
    r_min=cos(Alpha);
    r_Max=sqrt(cos(Alpha)^2+rp_Max^2)-4*r_sampling;
    r_Vector=single(r_min:r_sampling:r_Max);
    N_r=numel(r_Vector);
    N_Tau0=20;
    N_Tau=3*N_Tau0;
    N_TauHalf=N_Tau/2;
    %Tau=linspace(pi/6,5*pi/6,N_Tau)';
    Tau=linspace(0,2*pi,N_Tau+1)';
    Tau=Tau(1:N_Tau);
    
    LMa=Load_LMa;
    %l_Vector=unique(LMa(:,1));
    N_l=numel(l_Vector);
    I_r_l=nan(N_r,N_l,'single');
    Energy=nan(N_l,1,'single');

    %Ico_Recon=zeros(N_r,N_Tau,'single');
    
    Plot_Ico_1_Flag=1;
    for l_cntr=1:N_l
        if Plot_Ico_1_Flag
            X=nan(60*N_Tau,N_r,'single');
            Y=nan(60*N_Tau,N_r,'single');
            Z=nan(60*N_Tau,N_r,'single');
        end
        l=l_Vector(l_cntr);
        for r_cntr=1:N_r
            r=r_Vector(r_cntr);
            rp=sqrt(r.^2-cos(Alpha)^2);
            if rp <= rp_min
                Tau=linspace(pi/2,7*pi/6,N_Tau)';
            else
                Tau_HR=acos(rp_min/rp);
                Tau=[linspace(pi/2,5*pi/6-Tau_HR,N_TauHalf),...
                    linspace(5*pi/6+Tau_HR,7*pi/6,N_TauHalf)]';
            end
            %Tau=[Tau;Tau+2*pi/3;Tau+4*pi/3];
            %Tau=unique(Tau);
            

            I_r_l(r_cntr,l_cntr)=sum(Integrand_t_Tau(Tau));
            
            subplot(2,2,1)
            %imagesc(abs(I_r_l));
            surf(double(log(abs(I_r_l'))));view(3);shading interp;
            axis tight;
            ylabel('{\bfIndex} of icosahedral order, l');xlabel('r index');
            view(80,50);

            if Plot_Ico_1_Flag

                subplot(2,2,2)
                hold on
                if rp < rp_min
                    plot3(Temp0.x,Temp0.y,Temp0.z,'b.');
                else
                    plot3(Temp0.x,Temp0.y,Temp0.z,'c.');
                end
                if r_cntr == 1
                    plot3(F(1),F(2),F(3),'r*')
                    plot3([0 F(1)],[0 F(2)],[0 F(3)],'r');
                    plot3(A(1),A(2),A(3),'k*')
                    plot3(B(1),B(2),B(3),'k*')
                    plot3(C(1),C(2),C(3),'k*')
                    plot3([B(1) A(1) C(1) B(1)],[B(2) A(2) C(2) B(2)],[B(3) A(3) C(3) B(3)],'k')
                end
                hold off
                axis tight;
                daspect([1 1 1]);
                xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
                title(['l=' num2str(l)]);
                %view(108,20)
                view(148,32);

                subplot(2,2,4)
                hold on
                if rp < rp_min
                    plot3(Temp.x,Temp.y,Temp.z,'b.');
                else
                    plot3(Temp.x,Temp.y,Temp.z,'c.');
                end
                if r_cntr == 1
                    plot3(F(1),F(2),F(3),'r*')
                    plot3([0 F(1)],[0 F(2)],[0 F(3)],'r');
                    plot3(A(1),A(2),A(3),'k*')
                    plot3(B(1),B(2),B(3),'k*')
                    plot3(C(1),C(2),C(3),'k*')
                    plot3([B(1) A(1) C(1) B(1)],[B(2) A(2) C(2) B(2)],[B(3) A(3) C(3) B(3)],'k')
                end
                hold off
                axis tight;
                daspect([1 1 1]);
                xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
                title(['l=' num2str(l)]);
                %view(108,20)
                view(148,32);

                drawnow;
            end
        end
        
        if Plot_Ico_1_Flag
            if l_cntr < N_l
                subplot(2,2,2)
                cla;
                subplot(2,2,4)
                cla;
                drawnow;
            end
        end
        
        Temp_b=I_r_l(:,l_cntr);
        Energy(l_cntr)=sum(Temp_b(Temp_b == Temp_b).^2);
        subplot(2,2,3)
        bar(l_Vector,Energy);
        axis tight;
        xlabel('Icosahedral order, l');ylabel('Energy');

        
        
        drawnow;
    end
    return;
            
        
    
    
    I_r=nan(N_r,1,'single');
    Integrands=nan(N_r,N_Tau,'single');
    l=0;
    for cntr=1:N_r
        r=r_Vector(cntr);
        rp=sqrt(r.^2-cos(Alpha)^2);
        Temp=Integrand_t_Tau(Tau);
        Integrands(cntr,:)=Temp';
        I_r(cntr)=sum(Temp);

        subplot(1,2,1)
        imagesc(Integrands);
        xlabel('Tau index');ylabel('r index')
        axis tight;
        subplot(1,2,2)
        plot(r_Vector,I_r,'b--.');
        axis tight;
    end
    
    function Integrand=Integrand_t_Tau(Tau)
        Theta=real(acos((0.6315+0.6071*rp.*sin(Tau))./r));
        Phi=real(asin(...
            (0.2836+rp.*(0.8090*cos(Tau)-0.4671*sin(Tau)))./(r.*sin(Theta))...
            ));

        RsinThetacosPhi=r.*sin(Theta).*cos(Phi);
        %Delta_11:d(Theta)/d(r')
        Delta_11=(rp.*cos(Theta)-0.6071*r.*sin(Tau))./(r.^2.*sin(Theta));
        %Delta_12: d(Theta)/d(Tau)
        Delta_12=-0.6071*(cos(Tau)./sin(Theta)).*(rp./r);
        %Delta_21: d(Phi)/d(r')
        Delta_21=(sin(Theta).*sin(Phi)*cos(Alpha)^2./(r.*rp)-0.2836./rp-...
            r.*sin(Phi).*cos(Theta).*Delta_11)./RsinThetacosPhi;
        %Delta_22: d(Phi)/d(Tau)
        Delta_22=(0.2836-rp.*(0.8090*sin(Tau)+0.4671*cos(Tau))...
            -r.*Delta_11.*cos(Theta).*sin(Phi))./RsinThetacosPhi;
        Delta=Delta_11.*Delta_22-Delta_12.*Delta_21;

        Integrand=(r/rp)*IHrmncSngl(Theta,Phi,l).*sin(Theta).*Delta;
        
        if Plot_Ico_1_Flag
            Temp0=struct(...
                'x',r*sin(Theta).*cos(Phi),...
                'y',r*sin(Theta).*sin(Phi),...
                'z',r*cos(Theta));
            Temp=struct('x',[],'y',[],'z',[]);
            for cntr=1:60
                Temp_a=Rotate_Coordinate(Temp0,single(R{cntr}));
                Temp.x=[Temp.x;Temp_a.x];        
                Temp.y=[Temp.y;Temp_a.y];        
                Temp.z=[Temp.z;Temp_a.z];        
            end
            X(:,r_cntr)=Temp.x;
            Y(:,r_cntr)=Temp.y;
            Z(:,r_cntr)=Temp.z;
        end
    end
    return
    
    
    
    
    
    
    %H=Int{f_ico(r,theta,phi) J_l(theta,phi) sin(theta) d_theta d_phi}=
    %Int{f_ico(r,theta(r',tau),phi(r',tau)) Jl(theta,phi) sin(theta)
    %Jacobian{(theta,phi)/(r',tau)} d_r' d_Tau}
    
    %r_min=sin(Alpha)/2=0.3035

    %1) r < cos(Alpha) OR r > 1 --> H(r)=0

    %2) cos(Alpha) < r < sqrt[r_min^2+cos(Alpha)^2]
    %Tau: U_n [pi/6+n*2*pi/3,5*pi/6+n*2*pi/3]
    %n from {0,1,2} for 3 asymmetric units and n=0 for 1 asymmetric unit

    %3) sqrt[r_min^2+cos(Alpha)^2] < r < 1
    %Tau: U_n [pi/6+X+n*2*pi/3,5*pi/6-X+n*2*pi/3], X=acos[r_min./rp]
    %n from {0,1,2} for 3 asymmetric units and n=0 for 1 asymmetric unit
    
    
    
    
    %rp=sqrt(r.^2-cos(Alpha).^2);
    %Theta=acos((0.6315+0.6071*rp.*sin(Tau))./r);
    %Phi=asin(...
    %    (0.2836+rp.*(0.8090*cos(Tau)-0.4671*sin(Tau)))./(r.*sin(Theta))...
    %    );
    
    %RsinThetacosPhi=r.*sin(Theta).*cos(Phi);
    
    %Delta_11: d(Theta)/d(r')
    %Delta_11=(rp.*cos(Theta)-0.6071*r.*sin(Tau))./(r.^2.*sin(Theta));

    %Delta_12: d(Theta)/d(Tau)
    %Delta_12=-0.6071*(cos(Tau)./sin(Theta)).*(rp./r);
    
    %Delta_21: d(Phi)/d(r')
    %Delta_21=(sin(Theta).*sin(Phi)*cos(Alpha)^2./(r.*rp)-0.2836./rp-...
    %    r.*sin(Phi).*cos(Theta).*Delta_11)./RsinThetacosPhi;

    %Delta_22: d(Phi)/d(Tau)
    %Delta_22=(0.2836-rp.*(0.8090*sin(Tau)+0.4671*cos(Tau))...
    %    -r.*Delta_11.*cos(Theta).*sin(Phi))./RsinThetacosPhi;
    
    %Delta=Delta_11.*Delta_22-Delta_12.*Delta_21;

    %clear Delta_11;
    %clear Delta_12;
    %clear Delta_21;
    %clear Delta_22;
    %clear RsinThetacosPhi;
    
    %I_l(r)=(r./rp)*Int{J_l(Theta,Phi).*sin(Theta).*Delta d_Tau}, where
    %Tau: [Tau_min(r'),Tau_max(r')]
    %Tau: U_n [pi/6+X+n*2*pi/3,5*pi/6-X+n*2*pi/3], X=acos[r_min./rp]
    %n from {0,1,2} for 3 asymmetric units and n=0 for 1 asymmetric unit
    
    %Test: Numerical orthonormality of icosahedral harmonics
    
    
    %_______________________
    return;
    

    WB=waitbar(0,'Calculating coordinate matrices');drawnow;
    a=single(1);
    N=single(20);
    OverSampling=single(1.5);
    N_R=fix(100*N);
    Spectral_Sampling=1/(4*a*OverSampling*2*pi);
    Spectral_Range=Spectral_Sampling*(N_R-1);
    u_R_1D=single(0:Spectral_Sampling:Spectral_Range);
    
    legendre_rule(num2str(N),'0',num2str(Theta_c),['./WTheta_N' num2str(N)]);
    legendre_rule(num2str(N),'0',num2str(Omega),['./WPhi_N' num2str(N)]);
    pause(1);
    W_Theta_1D=importdata(['./WTheta_N' num2str(N) '_w.txt']);
    u_Theta_1D=importdata(['./WTheta_N' num2str(N) '_x.txt']);
    W_Phi_1D=importdata(['./WPhi_N' num2str(N) '_w.txt']);
    u_Phi_1D=importdata(['./WPhi_N' num2str(N) '_x.txt']);
    %u_Theta_1D=single(linspace(0,Theta_c,N));
    %u_Phi_1D=single(linspace(0,Omega,N));
    
    
    [u_R,u_Theta,u_Phi]=ndgrid(u_R_1D,u_Theta_1D,u_Phi_1D);
    [~,W_Theta,W_Phi]=ndgrid(u_R_1D,W_Theta_1D,W_Phi_1D);
    Index_null=(tan(u_Theta) > tan(Theta_c)*sin(Omega)./(1-cos(Omega)+sin(Omega)*cos(u_Phi)));
    [u_x,u_y,u_z]=Polar2Cartesian(u_R,u_Theta,u_Phi);
    r_sampling=(1/Spectral_Range)*1.1;
    r_range=1.8*a;
    r=single(0:r_sampling:r_range);
    N_r=numel(r);
    Window_Cutoff=0.80;
    Window=ones(N_R,1,'single');
    N1=single(ceil(Window_Cutoff*N_R));
    Index=N1:N_R;
    Window(Index)=cos((single(pi/2)/(N_R-N1))*(Index-N1)).^1;
    
    u=struct('x',u_x,'y',u_y,'z',u_z);
    clear u_x;
    clear u_y;
    clear u_z;

    WB=waitbar(0,WB,'Calculating rotation matrices');drawnow;
    R0=single(Rotate_Tri_to_EMAN());
    R=IcosSymRot(R0);
    WB=waitbar(0,WB,'Calculating the diffraction volume');drawnow;
    DV=u.x*single(0);
    for cntr=1:60;
        DV=DV+Fourier_Tri(Rotate_Coordinate(u,single(R{cntr})));
        WB=waitbar(cntr/60,WB,'Calculating the diffraction volume');drawnow;
    end
    DV(Index_null)=0;
    close(WB);drawnow;
    
    [J_l,l,N_l]=Icosahedral_Harmonics(u_Theta,u_Phi);
    clear u_Phi;
    
    f_l_u=cell(N_l,1);
    Integrand_1=sin(u_Theta).*DV.*W_Theta.*W_Phi;
    
    %clear DV;
    clear u_Theta;
    
    N3=ceil(N_l/3);
    Energy_u=zeros(N_l,1,'single');
    Recon_Error=zeros(N_l,1,'single');
    DV_Recon=DV*single(0);
    H1=figure;
    for cntr=1:N_l
        %Temp=sum(sum(J_l{cntr}.*Integrand_1,2),3);
        Temp=sum(sum(J_l{cntr}.*Integrand_1,2),3);
        Temp(Temp ~= Temp)=0;
        f_l_u{cntr}=Temp;
        Energy_u(cntr)=sum(abs(Temp).^2);
        DV_Recon=DV_Recon+(-(1i))^l(cntr)*J_l{cntr}.*repmat(Temp,[1 N N]);
        Index=(DV==DV);
        Recon_Error(cntr)=norm(DV(Index)-DV_Recon(Index))/norm(DV(Index));
        subplot(3,N3,cntr);
        plot(u_R_1D,real(Temp),'b');
        hold on
        plot(u_R_1D,imag(Temp),'r');
        title(['f_{l=' num2str(l(cntr)) '}(u)']);
        legend('Real[f]','Imag[f]');
        drawnow;
    end
    
    subplot(3,N3,N_l+1);
    semilogy(Energy_u/max(Energy_u),'*--')
    title('Icosahedral spectrum of an icosahedral shell');
    subplot(3,N3,N_l+2);
    semilogy(Recon_Error,'*--')
    title('Reconctruction Error vs. icosahedral order');
    
    
    
    
    f_l_r=cell(N_l,1);
    Energy_r=zeros(N_l,1,'single');
    H2=figure;
    for cntr=1:N_l
        Integrand_1=u_R_1D'.^2.*f_l_u{cntr}.*Window;
        Temp=zeros(N_r,1);
        for cntr2=1:numel(r)
            Temp(cntr2)=sum(Integrand_1(:).*...
                Spherical_Bessel(l(cntr),2*pi*r(cntr2)*u_R_1D'));
        end
        Energy_r(cntr)=sum(abs(Temp).^2);
        f_l_r{cntr}=Temp;
        subplot(3,N3,cntr);
        plot(r,real(Temp),'b');
        hold on
        plot(r,imag(Temp),'r');
        title(['f_{l=' num2str(l(cntr)) '}(r)']);
        legend('Real[f]','Imag[f]');
        drawnow;
    end
    subplot(3,N3,N_l+1);
    semilogy(Energy_r/max(Energy_r),'*--')
    title('Icosahedral spectrum of an icosahedral shell');
    return;
    
    
    
    
    IsoValue=[5 10 20 40]*median(DV(:));
    M=numel(IsoValue);
    figure;
    for cntr=1:M
        subplot(1,M,cntr);
        Plot_Iso_Surface(u.x,u.y,u.z,DV,IsoValue(cntr));
        view(-46,-54)
    end
    drawnow;
    
    
    
    
    N=100;
    x_1D=linspace(0,C(1),N);
    y_1D=linspace(0,A(2),N);
    z_1D=linspace(C(3),B(3),N);
    [x,y,z]=meshgrid(x_1D,y_1D,z_1D);
    Index=(([x(:)-F(1),y(:)-F(2),z(:)-F(3)]*F).^2 < 1e-5) & (z(:) <= 1-(tan(Theta_c/2)/sin(Omega))*y(:));
    x=x(Index);y=y(Index);z=z(Index);

    
    scatter3(x(:),y(:),z(:),40,ones(size(z(:))),'.');
    hold on
    plot3(F(1),F(2),F(3),'r*')
    plot3([0 F(1)],[0 F(2)],[0 F(3)],'r');
    plot3(A(1),A(2),A(3),'k*')
    plot3(B(1),B(2),B(3),'k*')
    plot3(C(1),C(2),C(3),'k*')
    plot3([B(1) A(1) C(1) B(1)],[B(2) A(2) C(2) B(2)],[B(3) A(3) C(3) B(3)],'k')
    axis tight;daspect([1 1 1]);xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
    view(148,32);
    hold off
    daspect([1 1 1]);colorbar;xlabel('x');ylabel('y');zlabel('z');axis tight;
    
    
    return;
    
    function I=Integral(f)
        switch Integral_Mode
            case 1
                I=mean(mean(f,1),2);
            case 2
                I=sum(sum(f.*W_Theta.*W_Phi,1),2);
        end
    end
    function Test_IHBasis_Rand(H)
        %W([5 8 11 13])=0;
        W=2*(rand(N_l,1)-0.5);
        F_Test=SynthIcoFunc(W);
        ISpectrum=Make_ISpectrum(F_Test);
        F_Recon=SynthIcoFunc(ISpectrum);
        
        figure(H);
        subplot(2,2,1);
        Plot_Shell(Theta,Phi,F_Test,nan);
        colorbar;        
        subplot(2,2,2);
        Plot_Shell(Theta,Phi,F_Recon,nan);
        colorbar;        
        subplot(2,4,5);
        bar(l_Vector,real(ISpectrum));xlim([min(l_Vector(1:N_l))-1 max(l_Vector(1:N_l))+1])
        subplot(2,4,6);
        bar(l_Vector,W);xlim([min(l_Vector(1:N_l))-1 max(l_Vector(1:N_l))+1])
        subplot(2,4,7);
        plot(W,real(ISpectrum),'*');
        xlabel('Input Weight');ylabel('Output Weight');
        subplot(2,4,8);
        plot(F_Test(:),F_Recon(:),'*');
        daspect([1 1 1]);xlabel('Input Function');ylabel('Reconstructed Function');xlim(1.5*[-1 1]);ylim(1.5*[-1 1]);
        drawnow;
    end
    function F=SynthIcoFunc(W)
        if size(W,2) > 1
            W=W';
        end
        if ArrayAngleFlag
            F=IH_Basis*W;
        else
            F=zeros(N_Phi,N_Theta);
            for Cntr=1:N_l
                F=F+W(Cntr)*IH_Basis(:,:,Cntr)*(-1)^(l_Vector(Cntr));
            end
        end
        %F=real(F)+imag(F);
        %F=F/F(1);
    end
    function Test_IHBasis(H)
        if ~numel(H)
            H=figure;
        end
        for input_cntr=1:N_l
            if ArrayAngleFlag
                ISpectrum=Make_ISpectrum(IH_Basis(:,input_cntr));
            else
                ISpectrum=Make_ISpectrum(squeeze(IH_Basis(:,:,input_cntr)));
            end
            figure(H);
            if N_subplot > 1
                subplot(N_subplot,2*N_subplot,input_cntr+N_subplot^2);
            else
                subplot(2,N_l,input_cntr+N_l);
            end
            bar(l_Vector,abs(ISpectrum));xlim([min(l_Vector(1:N_l))-1 max(l_Vector(1:N_l))+1])
        end
        drawnow;
    end
    function Plot_IcoBoundaries(x,y)
    end

    function ISpectrum=Make_ISpectrum(f)
        Integrand_=f.*sin(Theta);
        ISpectrum=zeros(N_l,1);
        for cntr=1:N_l
            if ArrayAngleFlag
                ISpectrum(cntr)=Integral(Integrand_.*IH_Basis(:,cntr));
            else
                ISpectrum(cntr)=Integral(Integrand_.*squeeze(IH_Basis(:,:,cntr)));
            end
            if mod(l_Vector(cntr),2)
                %ISpectrum(cntr)=-ISpectrum(cntr);
            end
        end
    end
end
function ISpectrum=Make_ISpectrumGlbl(f,Theta,IH_Basis,N_l,ArrayAngleFlag)
    Integrand_=f.*sin(Theta);
    ISpectrum=zeros(N_l,1);
    for cntr=1:N_l
        if ArrayAngleFlag
            ISpectrum(cntr)=Integral(Integrand_.*IH_Basis(:,cntr));
        else
            %ISpectrum(cntr)=Integral(Integrand_.*squeeze(IH_Basis(:,:,cntr)));
            ISpectrum(cntr)=mean(mean(Integrand_.*IH_Basis,1),2);
        end
    end
end

function IcoGridVol=Make_IcoVolGrid(B,C,A,Radius1,Radius2,N_TriGrid,N_R)
    if isequal(Radius1,Radius2)
        IcoGridVol=Make_IcoShellGrid(B,C,A,Radius1,N_TriGrid);
    else
        Radius=linspace(Radius1,Radius2,N_R);
        Shell=Make_IcoShellGrid(B,C,A,1,N_TriGrid);
        N0=numel(Shell.x);
        Zero=zeros(N_R*N0,1);
        X=Zero;
        Y=Zero;
        Z=Zero;
        for cntr=1:N_R
            X((1:N0)+(cntr-1)*N0)=Radius(cntr)*Shell.x;
            Y((1:N0)+(cntr-1)*N0)=Radius(cntr)*Shell.y;
            Z((1:N0)+(cntr-1)*N0)=Radius(cntr)*Shell.z;
        end
        IcoGridVol=struct('x',X,'y',Y,'z',Z);
    end
    %numel(IcoGridVol.x)
    IcoGridVol=RemoveRedundancy(IcoGridVol);
    %numel(IcoGridVol.x)
end
function [IH_Basis,ArrayAngleFlag]=Make_IHBasis(Theta,Phi,l_Vector,H)
    N_l=numel(l_Vector);
    IH_Basis=zeros(size(Theta,1),size(Theta,2),N_l);
    if isequal(min(size(Theta)),1)
        ArrayAngleFlag=1;
    else
        ArrayAngleFlag=0;
    end
    
    if numel(H)
        Plot_Flag=1;
        N_subplot=ceil(sqrt(N_l));
    else
        Plot_Flag=0;
    end
    for l_cntr=1:N_l
        l=l_Vector(l_cntr);
        Temp=IHrmncSngl(Theta,Phi,l);
        if ArrayAngleFlag
            IH_Basis(:,l_cntr)=Temp;
        else
            IH_Basis(:,:,l_cntr)=Temp;
        end
        if Plot_Flag
            if N_subplot > 1
                subplot(N_subplot,2*N_subplot,l_cntr);
            else
                subplot(2,N_l,l_cntr);
            end
            Plot_Shell(Theta,Phi,Temp,l);
            drawnow;
        end
    end
end

function [B,C,A,F,Omega,Theta_c,Alpha,Beta,Delta,Theta_F,Phi_min,Phi_Max]=...
    Load_IH_Tri(IcshFaceMode)
    %Constant Angles
    Omega=2*pi/5;
    Theta_c=acos(cos(Omega)/(1-cos(Omega)));
    Alpha=acos(1/(sqrt(3)*tan(Omega/2)));%37.38 deg or 0.6524 rad
    %tan(Alpha), tan(Theta_c/2)/cos(Omega/2)
        %r=1./(cos(Theta)+Beta*cos(Phi-Omega/2).*sin(Theta));

    Beta=tan(Theta_c/2)/cos(Omega/2);   %0.7639
    Delta=atan(1/Beta); %52.62 deg (91.84 rad); sin(Delta)=0.7947=cos(Alpha)
    %Alpha+Delta=pi/2
    if IcshFaceMode %1/0: Phi_min=0/-Omega/2
        Phi_min=0;
    else
        Phi_min=-Omega/2;
    end
    Phi_Max=Phi_min+Omega;
    
    %Constant Points
    B=[0;0;1];  %Vertex
    C=[sin(Theta_c)*cos(Phi_min);sin(Theta_c)*sin(Phi_min);cos(Theta_c)];    %Vertex
    A=[sin(Theta_c)*cos(Phi_Max);sin(Theta_c)*sin(Phi_Max);cos(Theta_c)];   %Vertex
    F=(A+B+C)/3;    %Center of gravity of face (corner of an 1/60 unit)
    Theta_F=(180/pi)*acos(F(3)/norm(F)); %37.3774 Alpha!
    %PlotIcshdrnFace(B,C,A)
end

function Plot_Shell(varargin)
    switch nargin
        case 5
            x=varargin{1};
            y=varargin{2};
            z=varargin{3};
            V=varargin{4};
            l=varargin{5};
        case 4
            Theta=varargin{1};
            Phi=varargin{2};
            V=varargin{3};
            l=varargin{4};
            x=sin(Theta).*cos(Phi);
            y=sin(Theta).*sin(Phi);
            z=cos(Theta);
    end
    if ~isreal(V)
        V=real((-1i)*V);
    end
    if min(size(x)) > 1
        surf(x,y,z,double(squeeze(V)));
    else
        scatter3(x,y,z,20,double(squeeze(V)),'.');
    end
    shading interp;
    %xlabel('x');
    %ylabel('y');
    %zlabel('z');
    if ~isnan(l)
        title(['l=' num2str(l)]);
    else
        title(['{\bfAll} icosahedral orders']);
    end
    axis tight;daspect([1 1 1]);xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
    view(148,32);
    %axis off
end

function PlotIcshdrnFace(B,C,A)
    IsHold=ishold;
    hold on
    F=(A+B+C)/3;
    plot3(F(1),F(2),F(3),'r*')
    plot3([0 F(1)],[0 F(2)],[0 F(3)],'r');
    plot3(A(1),A(2),A(3),'k*')
    plot3(B(1),B(2),B(3),'k*')
    plot3(C(1),C(2),C(3),'k*')
    plot3([B(1) A(1) C(1) B(1)],[B(2) A(2) C(2) B(2)],[B(3) A(3) C(3) B(3)],'k')
    axis tight;daspect([1 1 1]);xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
    view(148,32);
    %view(100,0);
    daspect([1 1 1]);colorbar;xlabel('x');ylabel('y');zlabel('z');
    if ~IsHold
        hold off
    end
end
function R0=Rotate_Tri_to_EMAN()
    b=sqrt(3)/2;
    Bp=[0;b;0];
    Cp=[-0.5;0;0];
    Ap=[0.5;0;0];
    R0=Rotate_to_EMAN(Bp,Cp,Ap);
end
function R=IcosSymRot20(varargin)
    if nargin
        R60=IcosSymRot(varargin{1});
    else
        R60=IcosSymRot();
    end
    Indices=[1:5,8:10,14:15,19:20,23,25,30:35];
    R=R60(Indices);
end
function R=IcosSymRot(varargin)
    %Ref: page 35 of http://docs.lib.purdue.edu/ecetr/203
    
    Beta=atan(2);
    Gamma=(2*pi)/5;

    S=Rot_z(Gamma);
    U=Rot_y(Beta);
    P=Rot_y(pi);
    T=U*S*inv(U);

    R=cell(60,1);

    for cntr=1:5
        R{cntr}=S^(cntr-1);
    end
    R{6}=S*T;
    R{7}=T*R{6};
    R{8}=T*R{7};
    R{9}=inv(T)*R{6};
    R{10}=inv(T)*R{9};
    for cntr=11:20
        R{cntr}=S*R{cntr-5};
    end
    for cntr=21:25
        R{cntr}=inv(S)*R{cntr-15};
    end
    for cntr=26:30
        R{cntr}=inv(S)*R{cntr-5};
    end
    for cntr=31:60
        R{cntr}=P*R{cntr-30};
    end
    
    if nargin
        ROT=varargin{1}*inv(R{1});
        for cntr=1:60
            R{cntr}=ROT*R{cntr};
        end
    end
end

function R=Rot_z(x)
    C=cos(x);
    S=sin(x);
    R=[ C    -S   0
        S    C    0
        0    0    1];
end
function R=Rot_y(x)
    C=cos(x);
    S=sin(x);
    R=[ C    0    S
        0    1    0
        -S   0    C];
end
function F=Fourier_Tri(u)
    %u is frequency and not radian frequency
    %u is normalized: u<-- u*a
    Sqrt3=single(sqrt(3));
    

    %with a
    %F=((1i*a/(2*pi))./u.y).*...
    %    (exp(-(1i*Sqrt3*pi*a)*u.y).*sinc(a*(Sqrt_3*u.y-u.x))-sinc(a*u.x));

    %normalized w/o a
    F=((1i/(2*pi))./u.y).*...
        (exp(-(1i*Sqrt3*pi)*u.y).*sinc(Sqrt3*u.y-u.x)-sinc(u.x));

    %only real part
    %Temp=Sqrt3*u.y;
    %F=sinc(Temp).*sinc(Temp-u.x);
end
function v=Rotate_Coordinate(u,R)
    [M,N,P]=size(u.x);
    Temp=R'*[u.x(:)';u.y(:)';u.z(:)'];
    if min(M,(min(N,P))) == 1
        v=struct('x',Temp(1,:)','y',Temp(2,:)','z',Temp(3,:)');
    else
        v=struct('x',reshape(Temp(1,:),[M N P]),...
            'y',reshape(Temp(2,:),[M N P]),...
            'z',reshape(Temp(3,:),[M N P]));
    end
end
function R=Rotate_to_EMAN(Bp,Cp,Ap)
    Omega=2*pi/5;
    Theta_c=acos(cos(Omega)/(1-cos(Omega)));
    Alpha=acos(1/(sqrt(3)*tan(Omega/2)));

    B=[0;0;1];
    C=[sin(Theta_c);0;cos(Theta_c)];
    A=[sin(Theta_c)*cos(Omega);sin(Theta_c)*sin(Omega);cos(Theta_c)];
    F=(A+B+C)/3;

    Fp=(Ap+Bp+Cp)/3;

    BCA_p_centroid=[Bp-Fp,Cp-Fp,Ap-Fp];
    BCA_centroid=[B-F,C-F,A-F]*(norm(Bp-Ap)/norm(B-A));

    M=[Bp-Fp,Cp-Fp,Ap-Fp]*[B-F,C-F,A-F]';
    %R=[Bp-Fp,Cp-Fp,Ap-Fp]/[B-F,C-F,A-F];
    [U,~,V]=svd(BCA_centroid*BCA_p_centroid');
    R=U*V';
    [V,~]=eig(R);
    for cntr=1:3
        if isreal(V(:,cntr))
            n=V(:,cntr);
        end
    end
    Theta=(180/pi)*acos((trace(R)-1)/2);
end
function Plot_Iso_Surface(x,y,z,DV,IsoValue)
    %p = patch(shrinkfaces(isosurface(x,y,z,DV,IsoValue),0.8));
    p = patch(isosurface(x,y,z,DV,IsoValue));
    isonormals(x,y,z,DV,p);
    set(p,'FaceColor','magenta','EdgeColor','none');
    daspect([1,1,1])
    view(3); axis tight
    
    Color={[1 0 0],[0 0 1]};
    for cntr=1:2
        h{cntr} = light;
        set(h{cntr},'Color',Color{cntr});
        lightangle(0,90+cntr*180);
    end
    
    %for az = -50:10:50
    %    lightangle(h,az,30);
    %end
    %light('Position',[1 3 2]);
    %light('Position',[-3 -1 3]);
%    material shiny
    
%    Q=3;
%    for cntr=1:Q
        %camlight(cntr*360/Q,40);%
        %camlight(cntr*360/Q,-90+cntr*180/Q);
    %end
    %camlight right;
    %lighting phong
    xlim([min(x(:)) max(x(:))])
    ylim([min(y(:)) max(y(:))])
    zlim([min(z(:)) max(z(:))])
end




%Current "Main"
function []=Icosahedron_Object_1()
    Initialize();
    
    
    
        
    %b=sqrt(3)/2;
    %Bp=[0;b;0];
    %Cp=[-0.5;0;0];
    %Ap=[0.5;0;0];
    %R=Rotate_to_EMAN(Bp,Cp,Ap);
    
    
    %Icosahedron=Load_Icosahedron();
    %Ico_Grid=Make_Icosahedral_Grid(Icosahedron);
    %Plot_Icosahedron(Icosahedron,Ico_Grid);
    %figure
    %subplot(121)
    %Plot_Vertex_Edges(V_parallel);
    %subplot(122)
    %Plot_Vertex_Edges(V_EMAN)
end

















%_______________________________
% Icosahedral harmonic stuff
%_______________________________


function [] = IH_Rot_3()
    Initialize();
    M=1024;
    N_1D_Cart=100;
    Mesh_GridBar_Flag=0;
    [x,y]=meshgrid(linspace(-1,1,M));
    r_norm=sqrt(x.^2+y.^2);
    Hamming1=hamming(M);
    Hamming2=Hamming1'*Hamming1;

    q=Scaled_q();
    q_Matrix=[q.X(:),q.Y(:),zeros(numel(q.X),1)]';
    %Load_Ex_Data=0;
    %WB=waitbar(0,'Loading {\bfmeasured} snapshots');drawnow;
    %[q,q_Matrix,Filter,Snapshots,N_Snapshots]=Load_Experiment(Load_Ex_Data);
    %close(WB);drawnow;

    %WB=waitbar(0,'Calculating the mask and outer indices');drawnow;
    %[Filter,OuterInd,Filter2]=Make_Mask(M);
    %close(WB);drawnow;
    
    %Plot_Log(Snapshots{1});return;

    %WB=waitbar(0,'Calculating I_q indices/weights and Blur');drawnow;
    %[Iq_Index,kn_factor,kn_factor_1D,Blur]=Find_Iq_Pixel_Indices();
    %kn_factor_Index=(kn_factor ~= 0);
    %kn_factor=kn_factor(kn_factor_Index);
    %close(WB);drawnow;
    
    %Icosahedral harmonics
    %N=80;
    %[Theta,Phi]=SphGrid(N);
    %global Zero_DV;
    %N_1D_Cart=30;
    Zero_DV=zeros(N_1D_Cart,N_1D_Cart,N_1D_Cart);
    WB=waitbar(0,'Calculating the 3D grid in real- and Fourier spaces');drawnow;
    [u_R,u_Theta,u_Phi,u_x,u_y,u_z,OS_1D,u_x_2D,u_y_2D,u_1D,Tri]=u_Grid(N_1D_Cart,q,Mesh_GridBar_Flag);

    WB=waitbar(0,WB,'Calculating basis functions');drawnow;
    [Basis_3D,N_total,n,l]=Make_3D_Basis(u_R,u_Theta,u_Phi);
    %[Iq_Index,kn_factor,kn_factor_1D,Blur]=Find_Iq_Pixel_Indices();
    Temp=max([u_x(:);u_y(:);u_z(:)])*(1/2);
    Range=[min(u_x(:)) max(u_x(:))];
    clear Temp;
    
 
    %WB=waitbar(0,WB,'Simulating rotations of icosahedral harmonics');drawnow;
    close(WB);drawnow;
    %Original snapshot
    %load Weights;
    %Weights([1 numel(Weights/2)])=8*max(abs(Weights));

    IH1_Max=1;
    IH2_Max=1;
    RB1_Max=1;
    RB2_Max=1;
    DVPhi_Max=(pi/5);
    DVTheta_Max=(pi/180)*63.43;

    DVPhi_Default=DVPhi_Max/2;
    DVTheta_Default=DVTheta_Max/2;

    DVTheta_Level=DVTheta_Default;
    DVPhi_Level=DVPhi_Default;

    %Coefficients of first- and higher-orders (1 for 0-order)
    IH_Default=[0.28+(1i)*0; 0+(1i)*0];
    RB_Default=[0+(1i)*0; 0+(1i)*0];
    Load_Defaults();

    H_Sns=figure;

    %Two Euler angles of rotation
    Ui{1}=uicontrol(H_Sns,'style','slider','Min',0,'Max',DVTheta_Max,'Value',...
        DVTheta_Level,'Callback',@DVThetaHandle,'Units','norm',...
        'Position',[0.04 0.83 0.02 0.15]);
    Ui{2}=uicontrol(H_Sns,'style','slider','Min',0,'Max',DVPhi_Max,'Value',...
        DVPhi_Level,'Callback',@DVPhiHandle,'Units','norm',...
        'Position',[0.09 0.83 0.02 0.15]);

    %Complex weights of icosahedral harmonics
    Ui{3}=uicontrol(H_Sns,'style','slider','Min',-IH1_Max,'Max',IH1_Max,'Value',...
        IHR_Level(1),'Callback',@IH1RHandle,'Units','norm',...
        'Position',[0.04 0.65 0.02 0.15]);
    Ui{4}=uicontrol(H_Sns,'style','slider','Min',-IH1_Max,'Max',IH1_Max,'Value',...
        IHI_Level(1),'Callback',@IH1IHandle,'Units','norm',...
        'Position',[0.09 0.65 0.02 0.15]);
    Ui{5}=uicontrol(H_Sns,'style','slider','Min',-IH2_Max,'Max',IH2_Max,'Value',...
        IHR_Level(2),'Callback',@IH2RHandle,'Units','norm',...
        'Position',[0.14 0.65 0.02 0.15]);
    Ui{6}=uicontrol(H_Sns,'style','slider','Min',-IH2_Max,'Max',IH2_Max,'Value',...
        IHI_Level(2),'Callback',@IH2IHandle,'Units','norm',...
        'Position',[0.19 0.65 0.02 0.15]);

    %Complex weights of radial basis functions
    Ui{7}=uicontrol(H_Sns,'style','slider','Min',-RB1_Max,'Max',RB1_Max,'Value',...
        RBR_Level(1),'Callback',@RB1RHandle,'Units','norm',...
        'Position',[0.04 0.47 0.02 0.15]);
    Ui{8}=uicontrol(H_Sns,'style','slider','Min',-RB1_Max,'Max',RB1_Max,'Value',...
        RBI_Level(1),'Callback',@RB1IHandle,'Units','norm',...
        'Position',[0.09 0.47 0.02 0.15]);
    Ui{9}=uicontrol(H_Sns,'style','slider','Min',-RB2_Max,'Max',RB2_Max,'Value',...
        RBR_Level(2),'Callback',@RB2RHandle,'Units','norm',...
        'Position',[0.14 0.47 0.02 0.15]);
    Ui{10}=uicontrol(H_Sns,'style','slider','Min',-RB2_Max,'Max',RB2_Max,'Value',...
        RBI_Level(2),'Callback',@RB2IHandle,'Units','norm',...
        'Position',[0.19 0.47 0.02 0.15]);

    %Defining "End" button
    Ui{11}=uicontrol(H_Sns,'style','pushbutton','string','End',...
        'Callback',@FinishFilt,'Units','norm','Position',[0.04 0.24 0.04 0.02]);
    
    %Defining the "Load Default" button
    Ui{12}=uicontrol(H_Sns,'style','pushbutton','string','Defaults',...
        'Callback',@LoadDefaultsHandle,'Units','norm','Position',[0.09 0.24 0.04 0.02]);
    
    %Labels
    Ui{13}=uicontrol(H_Sns,'Style','text','Units','norm',...
        'Position',[0.005 0.90 0.032 0.012],'String','Rotation');
    Ui{14}=uicontrol(H_Sns,'Style','text','Units','norm',...
        'Position',[0.005 0.72 0.032 0.012],'String','Icosahedral');
    Ui{15}=uicontrol(H_Sns,'Style','text','Units','norm',...
        'Position',[0.005 0.54 0.032 0.012],'String','Radial');
   


    %Weights=[1 10+(1i)*2];
    DV_Interp=Make_DV(Weights);
    Draw_Snapshot([0,DVTheta_Level,DVPhi_Level])
    Update_Orientation_Graph(DVTheta_Level,DVPhi_Level);


    

    %Infinite loop (interrupts taking over)
    FiltBlur_Flag=1;
    while FiltBlur_Flag
        drawnow;
        pause(1e-4);
    end
    %End of Main program
    
    
    
    
    
    
    function []=Load_Defaults()
        DVPhi_Level=DVPhi_Default;
        DVTheta_Level=DVTheta_Default;

        IHR_Level=real(IH_Default);
        IHI_Level=imag(IH_Default);
        RBR_Level=real(RB_Default);
        RBI_Level=imag(RB_Default);

        Update_Weights();
    end

    function []=Update_GUI()
        set(Ui{1},'Value',DVTheta_Level);
        set(Ui{2},'Value',DVPhi_Level);
        set(Ui{3},'Value',IHR_Level(1));
        set(Ui{4},'Value',IHI_Level(1));
        set(Ui{5},'Value',IHR_Level(2));
        set(Ui{6},'Value',IHI_Level(2));
        set(Ui{7},'Value',RBR_Level(1));
        set(Ui{8},'Value',RBI_Level(1));
        set(Ui{9},'Value',RBR_Level(2));
        set(Ui{10},'Value',RBI_Level(2));
    end

    function []=Update_Weights()
        IH_Level=IHR_Level+(1i)*IHI_Level;
        RB_Level=RBR_Level+(1i)*RBI_Level;
        Weights=[1;RB_Level]*transpose([1;IH_Level]);
        Weights=Weights(:);
    end

    function Update_Orientation_Graph(DVTheta,DVPhi)
        Gamma=36;
        Theta_c=63.43;
        subplot(9,8,1);
        hold off
        plot([0 Gamma],[0 0],'b');
        hold on
        plot([0 0],[Theta_c 0],'b');
        plot([Gamma Gamma],[0 (180/pi)*atan(0.618./cos((pi/180)*Gamma))],'b');
        Phi=linspace(0,Gamma,200);
        Theta=(180/pi)*atan(0.618./cos((pi/180)*(2*Gamma-Phi)));
        plot(Phi,Theta,'b');
        plot(DVPhi*(180/pi),DVTheta*(180/pi),'r*');
        hold off;
        daspect([1 1 1]);
        xlim([0 Gamma]);
        ylim([0 Theta_c]);
        xlabel('Euler angle \phi');
        ylabel('Euler angle \theta');
        box off;
        %axis off;
    end

    function []=Draw_Snapshot(ELR)
        subplot(1,2,2);
        PSI=ELR(1);
        THETA=ELR(2);
        PHI=ELR(3);
        Plot_Log(abs(Section_DV(Euler2RotMat(PSI,THETA,PHI))));
        title(['(\theta,\phi)=(' num2str(fix(THETA*(180/pi))) ...
            '\circ,' num2str(fix(PHI*(180/pi))) '\circ)' ...
            ', Icosahedral spectrum: I_0=1' ...
            ', I_1=' num2str(IHR_Level(1)) '+' num2str(IHI_Level(1)) 'i' ...
            ', I_2=' num2str(IHR_Level(2)) '+' num2str(IHI_Level(2)) 'i' ...
            ', Radial spectrum: R_0=1' ...
            ', R_1=' num2str(RBR_Level(1)) '+' num2str(RBI_Level(1)) 'i' ...
            ', R_2=' num2str(RBR_Level(2)) '+' num2str(RBI_Level(2)) 'i' ...
            ]);
        drawnow;
    end
    
    %DV rotations
    function DVThetaHandle(hObj,Event)
        DVTheta_Level=get(hObj,'Val');
        Draw_Snapshot([0,DVTheta_Level,DVPhi_Level]);
        Update_Orientation_Graph(DVTheta_Level,DVPhi_Level);
        drawnow;
    end
    function DVPhiHandle(hObj,Event)
        DVPhi_Level=get(hObj,'Val');
        Draw_Snapshot([0,DVTheta_Level,DVPhi_Level])
        Update_Orientation_Graph(DVTheta_Level,DVPhi_Level);
        drawnow;
    end

    %Icosahedral harmonics
    function IH1RHandle(hObj,Event)
        IHR_Level(1)=get(hObj,'Val');
        Update_Weights();
        DV_Interp=Make_DV(Weights);
        Draw_Snapshot([0,DVTheta_Level,DVPhi_Level])
        drawnow;
    end
    function IH1IHandle(hObj,Event)
        IHI_Level(1)=get(hObj,'Val');
        Update_Weights();
        DV_Interp=Make_DV(Weights);
        Draw_Snapshot([0,DVTheta_Level,DVPhi_Level])
        drawnow;
    end
    function IH2RHandle(hObj,Event)
        IHR_Level(2)=get(hObj,'Val');
        Update_Weights();
        DV_Interp=Make_DV(Weights);
        Draw_Snapshot([0,DVTheta_Level,DVPhi_Level])
        drawnow;
    end
    function IH2IHandle(hObj,Event)
        IHI_Level(2)=get(hObj,'Val');
        Update_Weights();
        DV_Interp=Make_DV(Weights);
        Draw_Snapshot([0,DVTheta_Level,DVPhi_Level])
        drawnow;
    end

    %Radial basis
    function RB1RHandle(hObj,Event)
        RBR_Level(1)=get(hObj,'Val');
        Update_Weights();
        DV_Interp=Make_DV(Weights);
        Draw_Snapshot([0,DVTheta_Level,DVPhi_Level])
        drawnow;
    end
    function RB1IHandle(hObj,Event)
        RBI_Level(1)=get(hObj,'Val');
        Update_Weights();
        DV_Interp=Make_DV(Weights);
        Draw_Snapshot([0,DVTheta_Level,DVPhi_Level])
        drawnow;
    end
    function RB2RHandle(hObj,Event)
        RBR_Level(2)=get(hObj,'Val');
        Update_Weights();
        DV_Interp=Make_DV(Weights);
        Draw_Snapshot([0,DVTheta_Level,DVPhi_Level])
        drawnow;
    end
    function RB2IHandle(hObj,Event)
        RBI_Level(2)=get(hObj,'Val');
        Update_Weights();
        DV_Interp=Make_DV(Weights);
        Draw_Snapshot([0,DVTheta_Level,DVPhi_Level])
        drawnow;
    end

    %Miscellaneous
    function LoadDefaultsHandle(hObj,Event)
        Load_Defaults();
        Update_Weights();
        Update_GUI();
        DV_Interp=Make_DV(Weights);
        Draw_Snapshot([0,DVTheta_Level,DVPhi_Level])
    end
    function FinishFilt(hObj,Event)
        FiltBlur_Flag=0;
        for cntr=1:numel(Ui)
            delete(Ui{cntr});
        end
    end

    
    
    








    %Distortions
    function BlurHandle(hObj,Event)
        Blur_Level=get(hObj,'Val');
        Plot_All_Menu(AddBlur(I0,Blur_Level),2)        
        Plot_All_Menu(AddOffset(AddAddNoise(AddMulNoise(AddBlur(I0,...
            Blur_Level),MulNoise_Level),AddNoise_Level),Offset_Level),6);
        drawnow;
    end
    function AddNoiseHandle(hObj,Event)
        AddNoise_Level=get(hObj,'Val')-1;
        Plot_All_Menu(AddAddNoise(I0,AddNoise_Level),3);
        Plot_All_Menu(AddOffset(AddAddNoise(AddMulNoise(AddBlur(I0,...
            Blur_Level),MulNoise_Level),AddNoise_Level),Offset_Level),6);
        drawnow;
    end
    function MulNoiseHandle(hObj,Event)
        MulNoise_Level=get(hObj,'Val')-1;
        Plot_All_Menu(AddMulNoise(I0,MulNoise_Level),4);
        Plot_All_Menu(AddOffset(AddAddNoise(AddMulNoise(AddBlur(I0,...
            Blur_Level),MulNoise_Level),AddNoise_Level),Offset_Level),6);
        drawnow;
    end
    function OffsetHandle(hObj,Event)
        Offset_Level=get(hObj,'Val')-1;
        Plot_All_Menu(AddOffset(I0,Offset_Level),5);
        Plot_All_Menu(AddOffset(AddAddNoise(AddMulNoise(AddBlur(I0,...
            Blur_Level),MulNoise_Level),AddNoise_Level),Offset_Level),6);
        drawnow;
    end

    %Rotations
    function DVPsiHandle(hObj,Event)
        DVPsi_Level=get(hObj,'Val');
        I0=abs(Section_DV(Rot2R(Rots(1:4,cntr)))).^2;
        
        Plot_All_Menu(I0,1)        
        Plot_All_Menu(AddBlur(I0,Blur_Level),2)        
        Plot_All_Menu(AddAddNoise(I0,AddNoise_Level),3);
        Plot_All_Menu(AddMulNoise(I0,MulNoise_Level),4);
        Plot_All_Menu(AddOffset(I0,Offset_Level),5);
        Plot_All_Menu(AddOffset(AddAddNoise(AddMulNoise(AddBlur(I0,...
            Blur_Level),MulNoise_Level),AddNoise_Level),Offset_Level),6);
        drawnow;
    end



    function I=AddBlur(I0,BL)
        I=abs(FFT2(FFT2(I0,1).*exp(-r_norm.^2/(BL/100)^2),-1)+Offset_Level);
    end
    function I=AddAddNoise(I0,NL)
        I=abs(I0+(NL/20)*randn(M,M));
    end
    function I=AddMulNoise(I0,NL)
        I=abs(I0.*(1+(NL/2)*randn(M,M)));
    end
    function I=AddOffset(I0,OL)
        I=abs(I0+OL*max(abs(I0(:)))/1e2);
    end
    function Plot_All_Menu(I,Column_Index)
        switch Column_Index
            case 1
                Title='Original';
            case 2
                Title='Blurred';
            case 3
                Title='Noisy (Additive)';
            case 4
                Title='Noisy (Multiplicative)';
            case 5
                Title='Offset';
            case 6
                Title='All distortions';
        end
        for cntr=1:SubPlot(1)
            Cntr=(cntr-1)*SubPlot(2)+Column_Index;
            subplot(SubPlot(1),SubPlot(2),Cntr);
            switch cntr
                case 1
                    Plot_Log(I,Title);
                case 2
                    %Image_2_Iq(I.*kn_factor,Iq_Index);
                    J=Image_2_Iq(I,Iq_Index);
                    semilogy(J);
                    if isequal(Column_Index,1)
                        Original{2}=J;
                    else
                        hold off;
                        hold on;
                        semilogy(Original{2},'r');
                        hold off;
                    end
                    legend('I_q','Location','Northeast');
                    axis tight;
                    box on;
                case 3
                    H=xcov(J);
                    H=H/max(H)+0.2;
                    %H=H(1:M/2);
                    semilogy(H);
                    if isequal(Column_Index,1)
                        Original{3}=H;
                    else
                        hold off;
                        hold on;
                        semilogy(Original{3},'r');
                        hold off;
                    end
                    legend('Cov[I_q]','Location','Northeast')
                    axis tight;
                    box on;
                case 4
                    K=Image_2_Iq(abs(FFT2(I,1)),Iq_Index);
                    K=K/max(K);
                    loglog(K);
                    if isequal(Column_Index,1)
                        Original{4}=K;
                    else
                        hold off;
                        hold on;
                        loglog(Original{4},'r');
                        hold off;
                    end
                    legend('J_d (J=|FT2[I]|)','Location','Southwest')
                    axis tight;
                    box on;
                case 5
                    K=Image_2_Iq(abs(FFT2(I-mean(I(:)),1)),Iq_Index);
                    K=K/max(K);
                    loglog(K);
                    if isequal(Column_Index,1)
                        Original{5}=K;
                    else
                        hold off;
                        hold on;
                        loglog(Original{5},'r');
                        hold off;
                    end
                    legend('J_d (J=|FT2[I-I_{mean}]|','Location','Southwest')
                    axis tight;
                    box on;
                case 6
                    Plot_Log(FastCorr2(I-mean(I(:))));
            end
            if isequal(cntr,1)
                title(Title)
            elseif ~isequal(cntr,6)
                ylim([min(Original{cntr}) max(Original{cntr})]);
            end
        end
    end
    function [Rot,Theta,Phi]=RotationMatrix(N)
        Theta=acos(linspace(1,-1,N+1));
        Theta=Theta(1:N);
        Phi=(2*pi)*linspace(0,1,N+1);
        Phi=Phi(1:N);
        [Theta,Phi]=meshgrid(Theta,Phi);
        [Theta,Phi]=IrrRotZone(Theta,Phi);


        Theta=53*(pi/180);Phi=0;
        
        
        N_IRZ=numel(Theta);
        Rot=cell(N_IRZ,1);
        for cntr=1:N_IRZ
            Rot{cntr}=Euler2RotMat(pi/2,Theta(cntr),Phi(cntr));
        end
    end
    function [Theta,Phi]=IrrRotZone(Theta,Phi)
        m=5;
        Omega=(2*pi/m);
        Phi_Max=0.5*Omega;
        Theta_Max=atan(1./(1.6180*cos(Omega-Phi)));
        Index=(0 <= Phi) & (Phi < Phi_Max) & (Theta > 0) & (Theta < Theta_Max);
        Theta=[Theta(Index);0];
        Phi=[Phi(Index);0];
    end
    function Rot=Euler2RotMat(Psi,Theta,Phi)
        Rot=Rot_z(Phi)*Rot_y(Theta)*Rot_z(Psi);
    end
    function []=SaveStructure()
        DV=Make_DV(Weights);
        FitStructure=struct();
        FitStructure.DV=DV;
        FitStructure.DV=DV;
        FitStructure.Basis_3D=Basis_3D;
        FitStructure.u_x=u_x;
        FitStructure.u_y=u_y;
        FitStructure.u_z=u_z;
        FitStructure.q_Matrix=q_Matrix;
        FitStructure.q=q;
        FitStructure.OuterInd=OuterInd;
        FitStructure.Filter=Filter;
        FitStructure.Iq_Index=Iq_Index;
        FitStructure.kn_factor=kn_factor;
        FitStructure.Snapshots=Snapshots;
        FitStructure.N_total=N_total;
        FitStructure.N_1D_Cart=N_1D_Cart;
        save FitStructure.mat FitStructure;pause;
    end
    function H=PlotIRot()
        H=figure;
        N_1D_plot=ceil(sqrt(N_Rot));
        for plot_cntr=1:N_Rot
            subplot(N_1D_plot,N_1D_plot,plot_cntr);
            Image=I_Rot{plot_cntr};
            Image=abs(FFT2(FFT2(Image,1).*Blur,-1));
            Plot_Log(Image);
            title(['\theta=' num2str(round((180/pi)*Theta_(plot_cntr))) '\circ, ' ...
                '\phi=' num2str(round((180/pi)*Phi_(plot_cntr))) '\circ']);
            %title(['Rotation index: ' num2str(plot_cntr)]);
        end
        drawnow;
    end
    function Param=OptimizeParameters()
        Opt_Loop_Flag=1;
        U_opt=uicontrol(H_comp,'style','pushbutton',...
            'string','End','Callback',@Finish_Opt);
        %<Levenberg>
        %[Options]=Opt_param();
        %UB=ones(size(Parameters));LB=-UB;
        %UB=[];LB=[];
        %[Parameters,ResNorm1,Residual1,ExitFlag1]= ...
        %    lsqnonlin(@(Param)OR_Func(...
        %    Param,N_total,N_1D_Cart,Basis_3D,Snapshots,H_comp,u_x,u_y,u_z,q_Matrix),...
        %    Parameters,LB,UB,Options);
        % [Weights,Rots]=Split_Opt_Param(Parameters,N_total);
         %disp(['ResNorm1 = ' num2str(ResNorm1)])
         %disp(['Weigths (final) = ' num2str(Weights')])

        %Plot_3D(Theta,Phi,J_l{cntr},l(cntr),cntr);
        %</Levenberg>
        %Parameters=Initialize_Parameters(N_Snapshots,N_total);
        %<CMA>
        Covariance_Plot_Flag=0;
        if Covariance_Plot_Flag
            clear global Stat;
            global Stat;
        end
        clear global Dim;
        global Dim;
        Dim=N_total+4*N_Snapshots;
        if (numel(l) == 1)
            Dim=N_total;
        end
        clear global Lambda;
        global Lambda;
        Lambda=40;
        N_Generations=300;
        Generations=1:N_Generations;
        Fitness=zeros(Lambda,1)';
        Best_Fitness=zeros(N_Generations,1)*nan;
        Mid_Fitness=zeros(N_Generations,1)*nan;
        Worst_Fitness=zeros(N_Generations,1)*nan;
        AllParameters=zeros(Dim,N_Generations);
        genes=Initialization('myoptions.m');
        %return;
        %genes=NextGeneration(rand(Lambda,1)');
        Base=1;
        g_cntr=1;
        while (g_cntr <= N_Generations) & Opt_Loop_Flag
            Param=Gene2Param(genes);  %input parameters
            for p_cntr=1:Lambda
                %2*genes(:,p_cntr)-1);  %input parameters
                %Fitness(p_cntr)=;      %output fitness
                %Plot_Flag=isequal(p_cntr,Lambda);
                Plot_Flag=0;
                Fitness(p_cntr)=1-OR_Func(Param(:,p_cntr),Plot_Flag); %output fitness
            end
            Fitness=Fitness/Base;
            [~,Ind]=sort(Fitness,'descend');
            Best_Fitness(g_cntr)=Fitness(Ind(1));
            Worst_Fitness(g_cntr)=Fitness(Ind(end));
            Mid_Fitness(g_cntr)=median(Fitness);
            [~,Median_Index]=sort(abs(Mid_Fitness(g_cntr)-Fitness),'ascend');
            %AllParameters(:,g_cntr)=Param(:,Median_Index(1));
            AllParameters(:,g_cntr)=Param(:,Ind(1));

            Plot_Flag=1;
            %OR_Func(Param(:,Median_Index(1)),Plot_Flag); %Just to reproduce/plot mid-fitness
            OR_Func(Param(:,Ind(1)),Plot_Flag); %Just to reproduce/plot mid-fitness
            subplot(432);
            plot(Generations,Best_Fitness,'b.-',...
                Generations,Mid_Fitness,'k.-',Generations,Worst_Fitness,'r.-');
            xlabel('Generation Index');
            ylabel('Fitness');
            legend('Best','Median','Worst','Location','Southeast');

            if Covariance_Plot_Flag
                Temp_=Stat(end);
                subplot(4,6,9);
                imagesc(Temp_.sig2);title('Covariance');daspect([1 1 1]);axis off;
                subplot(4,6,10);
                bar(Temp_.sig1(end:-1:1)/max(abs(Temp_.sig1)));
                title('Eigenvalues');axis tight;ylim([0 1]);axis off;
            end

            drawnow

            genes=NextGeneration(Fitness);
            g_cntr=g_cntr+1;
            
            %global stat options;
            %save stat.mat stat;
            %saveEvolutionStatus(stat, options);
            %opts = myrestart
            %load(opt.loadStatusFileName);
        end
        AllParameters=AllParameters(:,1:(g_cntr-1));
        save(['Parameters_' mfilename '.dat'],'AllParameters','-ascii');
        disp(['Parameters saved as: "Parameters_' mfilename '.dat"']);

        %Mean_Fitness=mean(Fitness);
        %</CMA>
        function Param=Gene2Param(Gene)
            Param=2*Gene-1;  %input parameters
        end
        function []=Finish_Opt(hObject,Event)
            Opt_Loop_Flag=0;
            disp(['Optimization interrupted; ' num2str(g_cntr) ...
                ' generations; Best fitness = ' ...
                num2str(max(Best_Fitness))]);
            delete(U_opt);
        end
    end
    function Error=OR_Func(Param,Plot_Flag)
        [Weights,Rots]=Split_Opt_Param(Param,N_total);
        DV=abs(Make_DV(Weights));
        %disp(['Weigths = ' num2str(Weights')])
        %disp(['Rots = ' num2str(Rots')])
        Error=1;
        Blur_Flag=0;
        for cntr=1:N_Snapshots
            I=Snapshots{cntr};
            %I_recon=Section_DV(DV,Rot2R(Rots(:,cntr)));
            %I_recon=abs(Section_DV(DV,Rot2R([1;0;0;0]))).^2;
            %I_recon=abs(FFT2(FFT2(I_recon,1).*Blur,-1));
            %I_recon=I_recon.*(1+1e-2*randn(M,M));
            %I_recon=Offset_Filter(I_recon,OuterInd,Filter);
            %I_recon=Offset_Filter(I_recon,OuterInd,Filter);
            %I_recon=medfilt2(I_recon,[3 3]);
            I_recon_raw=abs(Section_DV(Rot2R(Rots(1:4,cntr)))).^2;
            %I_recon_raw=abs(Section_DV(DV,Rot2R([1;0;0;0]))).^2;
            I_recon=Offset_Filter(I_recon_raw,OuterInd,Filter);
            %I_recon=I_recon-mean(I_recon(:));
            %I_recon=I_recon*(norm(I-mean(I(:)))/norm(I_recon(:)));
            %I_recon=I_recon+std(std(I(OuterInd)))-std(std(I_recon(OuterInd)));
            %I_recon=I_recon-min(I_recon(:));
            if Blur_Flag
                I=abs(FFT2(FFT2(I,1).*Blur,-1));
                I_recon=abs(FFT2(FFT2(I_recon,1).*Blur,-1));
            end
            if ~any(I_recon)
                disp('All-zero I_recon');
            elseif isnan(I_recon)
                disp('I_recon returning NaN');
            elseif isinf(I_recon)
                disp('I_recon returning Inf');
            end
            if isnan(I)
                disp('I returning NaN');
            elseif isinf(I)
                disp('I returning Inf');
            end
            if Plot_Flag
                Plot_Recon(I,I_recon,I_recon_raw,Weights);
                Plot_Proj([])
                %Plot_3D([])
            end
            Fitness_Method=4;
            switch Fitness_Method
                case 1  % Normalized norm of 2D covariance
                    Temp1=kn_factor.*(I-mean(I(:)));
                    Temp2=kn_factor.*(I_recon-mean(I_recon(:)));
                    Error=1-Norm(FastCorr2(Temp1,Temp2))/sqrt(...
                        Norm(FastCorr2(Temp1))*Norm(FastCorr2(Temp2)));
                case 2  % Normalized norm of 2D FT's
                    Temp1=Hamming2.*fft2(Hamming2.*(I-mean(I(:))));
                    Temp2=Hamming2.*fft2(Hamming2.*(I_recon-mean(I_recon(:))));
                    Error=acos(abs(Temp1(:)'*Temp2(:))/sqrt(...
                        abs(Temp1(:)'*Temp1(:))*abs(Temp2(:)'*Temp2(:))))/(pi/2);
                case 3  % Normalized norm of 2D FT's & 1D xcov (sum of acos)
                    mm=1;
                    lmm=1;
                   
                    Temp1=Hamming2.*fft2(Hamming2.*(I-mean(I(:))));
                    Temp2=Hamming2.*fft2(Hamming2.*(I_recon-mean(I_recon(:))));
                    Temp3=Image_2_Iq((I-mean(I(:))).*kn_factor,Iq_Index);
                    Temp4=Image_2_Iq((I_recon-mean(I_recon(:))).*kn_factor,Iq_Index);
                    Temp3=Temp3-mean(Temp3);
                    Temp4=Temp4-mean(Temp4);
                    Error=(mm/(pi/2))*real(acos(abs(Temp1(:)'*Temp2(:))/sqrt(...
                        abs(Temp1(:)'*Temp1(:))*abs(Temp2(:)'*Temp2(:)))))+...
                        real(acos(abs(Temp3'*Temp4)/sqrt((Temp3'*Temp3)*...
                        (Temp4'*Temp4))))*(lmm/(pi/2));
                case 4  % Normalized norm of 2D FT's & 1D xcov (acos of product)
                    mm=1;
                    lmm=1;
                    Temp1=Hamming2.*fft2(Hamming2.*(I-mean(I(:))));
                    Temp2=Hamming2.*fft2(Hamming2.*(I_recon-mean(I_recon(:))));
                    Temp3=Image_2_Iq((I-mean(I(:))).*kn_factor,Iq_Index);
                    Temp4=Image_2_Iq((I_recon-mean(I_recon(:))).*kn_factor,Iq_Index);
                    Temp3=Temp3-mean(Temp3);
                    Temp4=Temp4-mean(Temp4);
                    %Error=Error*(2/pi)*real(acos(...
                    Error=Error*...
                        (abs(Temp1(:)'*Temp2(:))/sqrt(...
                        abs(Temp1(:)'*Temp1(:))*abs(Temp2(:)'*Temp2(:))))^mm*...
                        (abs(Temp3'*Temp4)/sqrt((Temp3'*Temp3)*(Temp4'*Temp4)))^lmm;
                case 5  % Difference of 1D xcov's
                    %mm=1;
                    %lmm=1;
                    Temp1=Hamming2.*fft2(Hamming2.*(I-mean(I(:))));
                    Temp2=Hamming2.*fft2(Hamming2.*(I_recon-mean(I_recon(:))));
                    Temp1=Temp1/norm(Temp1(:));
                    Temp2=Temp2/norm(Temp2(:));
                    %Temp3=Image_2_Iq((I-mean(I(:))).*kn_factor,Iq_Index);
                    %Temp4=Image_2_Iq((I_recon-mean(I_recon(:))).*kn_factor,Iq_Index);
                    Temp3=Image_2_Iq(I,Iq_Index);
                    Temp4=Image_2_Iq(I_recon,Iq_Index);
                    Temp3=Temp3-sum(Temp3(end-(0:20)));
                    Temp4=Temp4-sum(Temp4(end-(0:20)));
                    Temp3=real(xcov(Temp3.*kn_factor_1D));
                    Temp4=real(xcov(Temp4.*kn_factor_1D));
                    Error=(norm(Temp3-Temp4)/norm(Temp3))*...
                        (norm(Temp1(:)-Temp2(:))/norm(Temp1(:)));
                    %Error=(2/pi)*real(acos(...
                    %    (abs(Temp1(:)'*Temp2(:))/sqrt(...
                    %    abs(Temp1(:)'*Temp1(:))*abs(Temp2(:)'*Temp2(:))))^mm*...
                    %    (abs(Temp3'*Temp4)/sqrt((Temp3'*Temp3)*(Temp4'*Temp4)))^lmm));
                case 6  % Normalized norm of 2D FT's (acos)
                    Temp1=Hamming2.*fft2(Hamming2.*(I-mean(I(:))));
                    Temp2=Hamming2.*fft2(Hamming2.*(I_recon-mean(I_recon(:))));
                    Temp1=Temp1-mean(Temp1(:));
                    Temp2=Temp2-mean(Temp2(:));
                    Error=(2/pi)*real(acos(abs(Temp1(:)'*Temp2(:))/sqrt(...
                        abs(Temp1(:)'*Temp1(:))*abs(Temp2(:)'*Temp2(:)))));
            end
            
            %Temp1=kn_factor(:).*(I(:)-mean(I(:)));
            %Temp2=kn_factor(:).*(I_recon(:)-mean(I_recon(:)));
            %Error=Error+norm(kn_factor(:).*(I_recon(:)/norm(I_recon(:))-I(:)/norm(I(:))));
            %Temp1=Image_2_Iq((I-mean(I(:))).*kn_factor,Iq_Index);
            %Temp2=Image_2_Iq((I_recon-mean(I_recon(:))).*kn_factor,Iq_Index);
            %Error=Error+(1-Temp1'*Temp2/(norm(Temp1)*norm(Temp2)));
            %Error=Error+(1-Temp1(:)'*Temp2(:)/(norm(Temp1(:))*norm(Temp2(:))));
            %Error=Error+1-norm(xcorr(Temp1,Temp2))/sqrt(...
            %    norm(xcorr(Temp1))*norm(xcorr(Temp2)));
            
            %Temp1=kn_factor.*I;
            %Temp2=kn_factor.*I_recon;
            %Error=Error+(1/pi)*acos(Norm(FastCorr2(Temp1,Temp2))/(Norm(FastCorr2(Temp1))*Norm(FastCorr2(Temp2))));
            %Error=norm(abs(Temp1/norm(Temp1)-Temp2/norm(Temp2)));
            if isnan(Error)
                disp('Objective function returning NaN');
            elseif isinf(Error)
                disp('Objective function returning Inf');
            end
        end
        Error=Error^(1/N_Snapshots);
        Error=(2/pi)*real(acos(Error));
        %Error=Error/N_Snapshots;
        %disp(['norm(Error) = ' num2str(norm(Error))])
        function Plot_Proj(Title)
            subplot(235);
            hold off
            %slice(u_x,u_y,u_z,DV,u_1D(end),u_1D(end),u_1D(1),'cubic');
            imagesc(squeeze(max(DV)));
            %view(3);
            daspect([1 1 1]);
            %caxis([min(DV(:)) max(DV(:))]);
            axis off;
            
            return
            %imagesc(squeeze(max(DV)));daspect([1 1 1]);axis off;return
            %scatter3(u_x(:),u_y(:),u_z(:),20,DV(:),'filled');return;
            p=patch(isosurface(u_x,u_y,u_z,DV));
            isonormals(u_x,u_y,u_z,DV,p);
            set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
            daspect([1 1 1]); view(3); camlight; lighting phong; axis image
            title(Title);
            xlabel('x(nm)');
            ylabel('y(nm)');
            zlabel('z(nm)');
            xlim(Range)
            ylim(Range)
            zlim(Range)
            %colorbar
        end
        function Plot_3D(Title)
            subplot(235);
            hold off
            DV_=DV;
            DV_(DV/max(DV(:)) < 1e-2)=nan;
            scatter3(u_x(:),u_y(:),u_z(:),20,DV_(:),'filled');return;
            return;
            p=patch(isosurface(u_x,u_y,u_z,DV_));
            isonormals(u_x,u_y,u_z,DV,p);
            set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
            daspect([1 1 1]); view(3); camlight; lighting phong; axis image
            title(Title);
            xlabel('x(nm)');
            ylabel('y(nm)');
            zlabel('z(nm)');
            xlim(Range)
            ylim(Range)
            zlim(Range)
            %colorbar
        end
    end
    function DV=Make_DV(Weights)
        %figure;plot(Weights,'--*')
        %Zero_DV=zeros(N_1D_Cart,N_1D_Cart,N_1D_Cart);
        %DV=zeros(N_1D_Cart,N_1D_Cart,N_1D_Cart);
        DV=Zero_DV;
        if ~any(Weights(:))
            disp('All-zero Weights (beginning of Make_DV)');
        elseif isnan(Weights(:))
            disp('Weights returning NaN');
        elseif isinf(Weights(:))
            disp('Weights returning Inf');
        end
        for cntr_q=1:numel(Weights)
            DV=DV+Basis_3D{cntr_q}*Weights(cntr_q);
        end
        if ~Mesh_GridBar_Flag
            DV=griddedInterpolant({u_1D,u_1D,u_1D},DV,'spline');
        end
    end
    function SingleHit=Section_DV(R)
        qR=R*q_Matrix;
        %if ~any(DV(:))
        %    disp('All-zero DV');
        %elseif isnan(DV(:))
        %    disp('DV returning NaN');
        %elseif isinf(DV(:))
        %    disp('DV returning Inf');
        %end
        %SingleHitComplex=interp3(u_x,u_y,u_z,DV,...
        %    qR(1,:),qR(2,:),qR(3,:),'linear',0);
        %DV_Interp=griddedInterpolant({u_1D,u_1D,u_1D},DV,'linear');
        if ~Mesh_GridBar_Flag
            SingleHitComplex=DV_Interp(qR(1,:),qR(2,:),qR(3,:));
            SingleHitComplex(SingleHitComplex ~= SingleHitComplex)=0;
        else
            SingleHitComplex=interp3(u_x,u_y,u_z,DV_Interp,...
                qR(1,:),qR(2,:),qR(3,:),'linear',0);
        end
        SingleHit=reshape(SingleHitComplex,[M M]);
        
        
        %From common arc
        %Intersection of images and the common arc
        %u = 1 + (CA_Nu-Temp_0_Nu(1))/(Temp_0_Nu(mx)-Temp_0_Nu(1))*(ncols-1);
        %v = 1 + (CA_Miu-Temp_0_Miu(1))/(Temp_0_Miu(my)-Temp_0_Miu(1))*(nrows-1);
        %uout = (u<.5)|(u>=ncols+.5);
        %anyuout = any(uout(:));
        %if anyuout, u(uout) = 1; end
        %vout = (v<.5)|(v>=nrows+.5);
        %anyvout = any(vout(:));
        %if anyvout, v(vout) = 1; end
        %Index=round(v)+(round(u)-1)*nrows;
        
        
        
        
    end
    function Plot_Recon(I,I_recon,I_recon_raw,Weights)
        figure(H_comp)

        subplot(231);
        hold off;
        [Min,Max]=Plot_Log(I);
        %colorbar;
        title(['Measured snapshot; Min = ' num2str(Min) ', Max = ' num2str(Max)]);

        subplot(233);
        hold off;
        [Min,Max]=Plot_Log(I_recon);
        %colorbar;
        title(['Reconstructed snapshot; Min = ' num2str(Min) ', Max = ' num2str(Max)]);

        subplot(4,3,7);
        hold off;
        Temp1=Image_2_Iq(I.*kn_factor,Iq_Index);
        Temp2=Image_2_Iq(I_recon.*kn_factor,Iq_Index);
        Temp1=Temp1/norm(Temp1(:));
        Temp2=Temp2/norm(Temp2(:));
        semilogy(Temp1,'b');
        hold on;
        semilogy(Temp2,'r');
        xlabel('Transverse scattering vector, q');
        legend('Measured','Reconstructed','Location','Southeast');
        hold off
        box on;
        title('w/ weighting')
        axis tight
        Temp3=[Temp1(:);Temp2(:)];
        ylim([ max(min(Temp3),max(Temp3)/1e4) max(Temp3)])

        subplot(4,3,10);
        hold off;
        Temp1=Image_2_Iq(I,Iq_Index);
        Temp2=Image_2_Iq(I_recon,Iq_Index);
        Temp1=Temp1/norm(Temp1(:));
        Temp2=Temp2/norm(Temp2(:));
        semilogy(Temp1,'b');
        hold on;
        semilogy(Temp2,'r');
        xlabel('Transverse scattering vector, q');
        legend('Measured','Reconstructed','Location','Southeast');
        hold off
        box on;
        title('w/o weighting')
        axis tight
        Temp3=[Temp1(:);Temp2(:)];
        ylim([ max(min(Temp3),max(Temp3)/1e4) max(Temp3)])
        

        subplot(4,3,9);
        hold off
        Temp1=Image_2_Iq(I-mean(I(:)),Iq_Index);
        Temp2=Image_2_Iq(I_recon-mean(I_recon(:)),Iq_Index);
        Temp=abs(fftshift(fft(fftshift(Temp1))));
        plot(Temp/norm(Temp(:)),'b');
        Temp=abs(fftshift(fft(fftshift(Temp2))));
        hold on;
        plot(Temp/norm(Temp(:)),'r');
        xlabel('Displacement, \deltad');
        legend('Measured','Reconstructed');
        hold off
        box on;
        title('Fourier transform');
        %xlim(250*[1 -1]+[1 512]);
        axis tight
        xlim(256+(100/2)*[-1 1]);

        subplot(4,3,12);
        hold off
        Temp=real(xcorr(Temp1));
        plot(Temp/norm(Temp(:)),'b');
        hold on;
        Temp=real(xcorr(Temp2));
        plot(Temp/norm(Temp(:)),'r');
        xlabel('Shift of transverse scattering vector, \deltaq');
        legend('Measured','Reconstructed');
        hold off
        box on;
        title('1D Covariance function');
        axis tight

        subplot(4,3,5);
        cla;
        hold off;
        [Min,Max]=Plot_Log(I_recon_raw.*Filter2);
        %colorbar;
        axis off;
        %text(0.05,0,['A(q)=\Sigma\lambda_n\psi_n(q) with \lambda_n: ' num2str(Weights')]);

        %drawnow;
    end
    function [Iq_Index,kn_factor,kn_factor_1D,Blur]=Find_Iq_Pixel_Indices()
        MM=1024;
        M_Half=MM/2;
        x_1D=linspace(-1,1,MM);
        Pixel=x_1D(2)-x_1D(1);
        %[X,Y]=meshgrid(x_1D);
        %r_norm=sqrt(X.^2+Y.^2);
        Iq_Index=cell(M_Half,1);
        for cntr_=1:M_Half
           Index_=find(abs(r_norm/Pixel - cntr_) <= 0.5);
           Iq_Index{cntr_}=Index_;
        end

        %x=k_n
        x1=0.60;
        x2=0.70;
        P=4;
        a=sqrt(1+(P/x1)^2);
        b=asin(1/a)/(x1-x2);
        c=-b*x2;
        %a*sin(b*x1+c)
        %a*sin(b*x2+c)
        
        kn_factor=(r_norm/x1).^P;
        Index__=(r_norm >= x1) & (r_norm >= x1);
        kn_factor(Index__)=a*sin(b*r_norm(Index__)+c);
        kn_factor(r_norm >= x2)=0;
        
        %figure;imagesc(kn_factor);colorbar;pause    
        Blur=exp(-r_norm.^2/0.6^2);
        %figure;imagesc(Blur);colorbar;pause
        kn_factor_1D=Image_2_Iq(kn_factor,Iq_Index);
    end
end

function [Min,Max]=Plot_Log(varargin)
    I=varargin{1};
    %I=abs(I);
    DR=50;
    Min=min(I(:));
    Max=max(I(:));
    I=log(1/DR+(abs(I)-Min)/(Max-Min));
    imagesc(I);
    daspect([1 1 1]);
    axis off;
    %colormap('hot');
    caxis(1.2*[min(I(:)) max(I(:))]);
    if nargin > 1
        title(varargin{2});
    end
    xlim(mean(xlim)+0.6*(xlim-mean(xlim)))
    ylim(mean(ylim)+0.6*(ylim-mean(ylim)))
end
    
   
function [Weights,Rots]=Split_Opt_Param(Parameters,N_total)
    Weights=Parameters(1:N_total);
    if numel(Parameters) > N_total
        Rots=Wrap(Parameters((N_total+1):end));
    else
        Rots=[1;0;0;0];
    end
    if ~any(Parameters)
        disp('All-zero Parameters (beginning of Split_Opt_Param)');
    elseif isnan(Parameters)
        disp('Parameters returning NaN');
    elseif isinf(Parameters)
        disp('Parameters returning Inf');
    end
    if ~any(Weights)
        disp('All-zero Weights (beginning of Split_Opt_Param)');
    elseif isnan(Weights)
        disp('Weights returning NaN');
    elseif isinf(Weights)
        disp('Weights returning Inf');
    end
end

function Rots=Wrap(Rots_)
    N_Snapshots=numel(Rots_)/4;
    Rots=zeros(4,N_Snapshots);
    for cntr=1:4
        Rots(cntr,1:N_Snapshots)=...
        Rots_((cntr-1)*N_Snapshots+(1:N_Snapshots));
    end
end






function Iq=Image_2_Iq(I,Iq_Index)
    MM=1024;
    M_Half=MM/2;
    Iq=zeros(M_Half,1);
    for cntr_=1:M_Half
       Iq(cntr_)=mean(mean(I(Iq_Index{cntr_})));
    end
end

  
function Rots=Unwrap(Rots_)
    Rots=[Rots_(1,:),Rots_(2,:),Rots_(3,:),Rots_(4,:)]';
end
    
function R=Rot2R(Rot)
    Threshold=1e-4;
    Norm=norm(Rot);
    if Norm > Threshold
        Rot=Rot/Norm;
    else
        Rot=[1;0;0;0];
    end
    a=Rot(1);
    b=Rot(2);
    c=Rot(3);
    d=Rot(4);
    R=[a^2+b^2-c^2-d^2, 2*(b*c-a*d), 2*(b*d+a*c); ...
        2*(b*c+a*d), a^2-b^2+c^2-d^2, 2*(c*d-a*b); ...
        2*(b*d-a*c), 2*(c*d+a*b), a^2-b^2-c^2+d^2];
end

function Parameters=Initialize_Parameters(N_Snapshots,N_total)
    Weights=Initialize_Weights(N_total);
    Rots=Initialize_Rots(N_Snapshots);
    Parameters=[Weights;Rots];

    function Weights=Initialize_Weights(N)
        %Weights=2*rand(N,1)-1;

        
        
        Weights=ones(N,1);
    end

    function Rots=Initialize_Rots(N)
        Rots_=2*rand(4,N)-1;
        Rots=Unwrap(Rots_);
        
        
        
        %Rots=[1;0;0;0];
        %Rots=[];
        
        
    end
end

function [Basis_3D,N_total,n,l]=Make_3D_Basis(u_R,u_Theta,u_Phi)
    [J_l,l,N_l]=Icosahedral_Harmonics(u_Theta,u_Phi);
    
    
    
    
    %New Icosahedron
    %Icosahedron=Load_Icosahedron();
    %Ico_Grid=Make_Icosahedral_Grid(Icosahedron);
    
    
    
    
    
    %New
    Size=0.8;
    x=Size*u_R.*sin(u_Theta).*cos(u_Phi);
    y=Size*u_R.*sin(u_Theta).*sin(u_Phi);
    z=Size*u_R.*cos(u_Theta);
    
    %x=x(:);
    %y=y(:);
    %z=z(:);
    
    
    Omega=(2*pi)/5;
    Theta_c=acos(cos(Omega)/(1-cos(Omega)));

    Thickness=0.1;
    Threshold_=Thickness/2;
    
    b=[0;0;1];
    c=[sin(Theta_c);0;cos(Theta_c)];
    a=[sin(Theta_c)*cos(Omega);sin(Theta_c)*sin(Omega);cos(Theta_c)];
    da2=sum([x(:)-a(1),y(:)-a(2),z(:)-a(3)].^2,2);
    db2=sum([x(:)-b(1),y(:)-b(2),z(:)-b(3)].^2,2);
    dc2=sum([x(:)-c(1),y(:)-c(2),z(:)-c(3)].^2,2);

    Index_=(y > sin(Theta_c)*sin(Omega));
    Index_=Index_ | (x < y/tan(Omega)) | (x > (y/tan(Omega)+sin(Theta_c)));
    Index_=Index_ | (z < cos(Theta_c));
    Index_ = Index_ | (abs(z - (1-(1-cos(Theta_c))*...
        (y*(1-cos(Omega))/(sin(Theta_c)*sin(Omega))+x/sin(Theta_c)))) > Threshold_);
    
    %Index_ = Index_ | (da2 < dc2) | (da2 < db2);

    V=zeros(size(u_R));
    V(Index_)=1;
    V((da2 < dc2) | (da2 < db2))=1;
    
    
    
    
    % /New
    
    
    
    
    
    
    N_n=1;
    n=1;
    N_total=N_n*N_l;
    figure;
    Basis_3D=cell(N_total,1);
    N1D=size(u_R,3);
    
    Recon=zeros(size(u_Theta));
    Index_Surf=(u_R(:) == max(u_R(:)));
    Energy=zeros(N_total,1);
    for cntr_=1:N_total
         %if (l(cntr_) ~= 15) & (l(cntr_) ~= 21) & (l(cntr_) ~= 0)
         if (l(cntr_) ~= 0)
            Temp=J_l{cntr_};
            Radial=sum(sum(J_l{cntr_}.*V.*sin(u_Theta),3),2);
            Recon=Recon+repmat(Radial,[1 N1D N1D]).*Temp;
            Basis_3D{cntr_}=Radial;
            Energy(cntr_)=Radial'*Radial;

            subplot(2,N_total+1,cntr_);
            plot(Radial);
            %title(num2str(cntr_))
            subplot(2,N_total+1,cntr_+N_total+1);
            scatter3(x(Index_Surf),y(Index_Surf),z(Index_Surf),20,Temp(Index_Surf));
            title(num2str(l(cntr_)))
            daspect([1 1 1])
            axis off
         end
    end
    save Recon.mat Recon;
    
    
    Recon=Recon(:);
    Index_Surf=(Recon > 0.10*max(Recon));
    subplot(2,N_total+1,N_total+1);
    scatter3(x(Index_Surf),y(Index_Surf),z(Index_Surf),20,Recon((Index_Surf)));
    title('Reconstructed Object')
    daspect([1 1 1])
    axis off
    
    subplot(2,N_total+1,2*(N_total+1));
    bar(l,Energy);
    xlim([min(l) max(l)])
    title('Icosahedral Spectrum')
    
    pause;
    

    
    return;
    
    
    
    
    [R_n,n,N_n]=Radial_Harmonics(u_R);
    N_total=N_n*N_l;
    Basis_3D=cell(N_total,1);
    Cntr=1;
    for cntr1=1:N_l
        for cntr2=1:N_n
            Basis_3D{Cntr}=J_l{cntr1}.*R_n{cntr2};
            Cntr=Cntr+1;
        end
    end
    
    
    
    

    function [R_n,n,N_n]=Radial_Harmonics(u_R)
        N_n=3;
        n=(1:N_n)-1;
        n=n+1;
        %n=[0 2];
        %n=1;N_n=numel(n);
        
        

        %q_Max*R=l_Max, q_Max=2*pi/d (d: resolution d), d/R=2*pi/l_max ~1/5
        Q=11.5*pi*u_R/max(u_R(:));
        R_n=cell(N_n,1);
        for cntr=1:N_n
            R_n{cntr}=Spherical_Bessel(n(cntr),Q)./Q;
            %R_n{cntr}=Spherical_Bessel(n(cntr),Q);
            %Temp=R_n{cntr};H=figure;plot(u_R(:),Temp(:));pause;close(H);drawnow;
        end
    end
end
function [J_l,l,N_l]=Icosahedral_Harmonics(Theta,Phi)
    %LMa=Load_LMa;
    %l_Vector=unique(LMa(:,1));
    %N_l=numel(l_Vector);






    %N_l=3;
    %l_Vector(1:N_l)
    %l=6p+10q
    %l'=6p+10q+15



    l=zeros(N_l,1);
    J_l=cell(N_l,1);
    WB=waitbar(0,'Calculating Icosahedral harmonics');drawnow;
    for cntr=1:N_l
        l(cntr)=l_Vector(cntr);
        J_l{cntr}=IHrmncSngl(Theta,Phi,l(cntr));
        WB=waitbar(cntr/N_l,WB,'Calculating Icosahedral harmonics');drawnow;
    end
    close(WB);drawnow;
    %H=figure;Temp=J_l{1};plot(Temp(:));pause;close(H);drawnow;
end
function j_n=Spherical_Bessel(Nu,r)
    Threshold=1e-20;
    j_n=sqrt(pi/2)*besselj(Nu+0.5,r)./sqrt(r);
    Index=(abs(r) < Threshold);
    if Nu & numel(Index)
        j_n(Index)=0;
    else
        j_n(Index)=1;
    end
end


function [q,q_Matrix,Filter,Snapshots,N_Snapshots]=Load_Experiment(Load_Ex_Data)
    
    %Scaled scattering vector
    q=Scaled_q();
    q_Matrix=[q.X(:),q.Y(:),zeros(numel(q.X),1)]';
    
    %[~,~,~,Blur]=Find_Iq_Pixel_Indices();
    
    %Parameters for basic pre-processing
    M=1024;
    Mid=round(M/2);
    Zero_M2=zeros(M,M);
    %camera offsets
    L1=5;   % camera offset (gap half size)
    L2=6;   % camera offset (gap half size)
    % Indices w/wo gap
    I1_Index=[(L1+1):Mid (Mid+1):(M-L2)];
    I2_Index=[1:(Mid-L1) (Mid+1+L2):M];
    [Filter,OuterInd,Filter2]=Make_Mask(M);

    %Mask
    %load Filter.mat;

    Path='./';
    Dataset_Index=3;
    if ~Load_Ex_Data
        Snapshots=[];
        N_Snapshots=0;
        return;
    end
    switch Dataset_Index
        case 1
            I1=importdata([Path 'S1_100_071_Raw.mat']);
            I2=Zero_M2; % Creation of the gap
            I2(I2_Index,1:M)=I1(I1_Index,1:M);
            Snapshots={Offset_Filter(I2,OuterInd,Filter)};
            N_Snapshots=numel(Snapshots);
        case 2
            load([Path 'Single_Hits']);
            Indices=[3 6 7 9];
            N_Snapshots=numel(Indices);
            Snapshots=cell(N,1);
            for cntr=1:N
                I1=Single_Hits{cntr};
                I2=Zero_M2; % Creation of the gap
                I2(I2_Index,1:M)=I1(I1_Index,1:M);
                Snapshots{cntr}=Offset_Filter(I2,OuterInd,Filter);
            end
        case 3
            Snapshots=importdata([Path 'I_Er_St_1.mat']);
            N_Snapshots=numel(Snapshots);
            L1=[6 8 2 0 0 4 -2];
            L2=16;
            for cntr=1:N_Snapshots
                I1=Snapshots{cntr};
                I2=Zero_M2; % Creation of the gap
                if L1(cntr) > 0
                    I2(end-M+(1:M),(L1(cntr)+1):M)=...
                    I1(end-round(L2/2)-M+(1:M),1:(M-L1(cntr)));
                else
                    I2(end-M+(1:M),1:(M+L1(cntr)))=...
                    I1(end-round(L2/2)-M+(1:M),(-L1(cntr)+1):M);
                end
                Snapshots{cntr}=Offset_Filter(I2,OuterInd,Filter);
            end
    end
    %Nr=numel(Snapshots);[x,y]=meshgrid(linspace(-1,1,M));r=sqrt(x.^2+y.^2);N_=ceil(sqrt(Nr));for cntr_=1:Nr;Temp=Snapshots{cntr_};Temp(r < 0.10)=0;Min=min(Temp(:));Max=max(Temp(:));Temp=log(1/50+(Temp-Min)/(Max-Min));subplot(N_,N_,cntr_);imagesc(Temp);caxis(1.2*caxis);daspect([1 1 1]);axis off;end;pause
end

function I2=Offset_Filter(I2,OuterInd,Filter)
    I2=abs((I2-mean(mean(I2(OuterInd))))).*Filter;
end
    
function [W,OuterInd,W2]=Make_Mask(M)
    %Window
    [x,y]=meshgrid(linspace(-1,1,M));
    r=sqrt(x.^2+y.^2);
    Theta=atan2(y,x);

    W=ones(size(r));
    
    a=0.065;
    b=0.135;
    m=0.5/(b-a);n=0.5-m*a;
    Index_23=(r>a & r<b);
    %W(Index_23)=cos(pi*(m*r(Index_23)+n)).^10;
    W(r <= a)=0;
    
    W2=W;
    

    a=0.90;
    b=0.80;
    m=0.5/(b-a);n=0.5-m*a;
    Index_23=(r>b & r<a);
    %W(Index_23)=cos(pi*(m*r(Index_23)+n)).^10;
    W(r >= a)=0;
    
    
    %M=size(I1,1);
    Mid=round(M/2);
    %Angle parameters
    dtheta=5*(pi/180);
    dtheta2=22*(pi/180);
    dtheta3=6*(pi/180);
    dtheta4=6*(pi/180);
    Pi_2=pi/2;
            
    Margin=(pi/180)*8;
    %Minus the CCW angle (because of meshgrid rather than ndgrid)
    Theta_mid=[-Pi_2; 0; Pi_2; pi; -pi];
    Theta_HalfRange=[dtheta2; dtheta4; dtheta; dtheta3; dtheta3];
    Theta_min=Theta_mid-Theta_HalfRange;
    Theta_max=Theta_mid+Theta_HalfRange;
    
    Index_pi=isequal(Theta_mid,pi);
    Theta_max(Index_pi)=pi;
    Theta_min(Index_pi)=pi-dtheta3;
    Index_Mpi=isequal(Theta_mid,-pi);
    Theta_max(Index_Mpi)=pi+dtheta3;
    Theta_min(Index_Mpi)=-pi;
    
    Theta_HalfRange(Index_pi)=Theta_HalfRange(Index_pi)*2;
    
    %Index_Arcs=(abs(Theta + Pi_2) < dtheta2) | (abs(Theta -0) < dtheta) ...
    %    | (abs(Theta - Pi_2) < dtheta) | (abs(Theta -pi) < dthetaHalf) ...
    %    | (abs(Theta + pi) < dthetaHalf);

    Index_Arcs=~isequal(Theta,Theta);
    for cntr=1:numel(Theta_mid)
        Index_Arcs=Index_Arcs | (abs(Theta-Theta_mid(cntr)) < Theta_HalfRange(cntr));
        %W=W.*Smooth_Arc(Theta,Theta_min(cntr),Theta_max(cntr),Margin);
    end
    
    
    
    

    %Outer region
    %OuterInd=(r > 1.01);    
    %Low-Q region
    LowQInd=(r < (50/(M/2)));
    HighQInd=(r > 0.8);
    % Jet
    Jet_Span=1:M;           
    Jet_Index=Mid+(-2:10);
    JetInd=isnan(x);
    JetInd(Jet_Index,Jet_Span)=1;
    %W(Index_Arcs | LowQInd | HighQInd)=0;
    W(Index_Arcs | LowQInd)=0;

    %OuterInd=(r > 1.3) & (W == 1);
    %586:1024, 333:336 -- 231:447,354:355 -- 691:1024,99 -- 555:1024,899
    W(333:336,586:M)=0;
    W(354:355,231:447)=0;
    W(99,691:M)=0;
    W(899,555:M)=0;
    W(179,287:371)=0;

    
    %OuterInd=(r > 0.84) & (W == 1);
    OuterInd=(r ~= r);
    OuterInd(30:70,30:70)=1;
    %W(HighQInd)=0;
    
    %W=isreal(r);
    %W=(r == r);
    %W(497:521,:)=0;
    %W(:,508:518)=0;
    %W(321:327,M/2+(1:M/2))=0;
    %W(344:345,1:M/2)=0;
    %W(r > 0.855)=0;
    %W(r < 0.065)=0;
    %W((abs(Theta + Pi_2) < dtheta))=0;
    %subplot(133);imagesc(W);daspect([1 1 1]);colorbar;pause;
end

function W=Smooth_Arc(Theta,Theta_min,Theta_max,Margin)
    W=ones(size(Theta));
    Factor=pi/(2*Margin);

    Theta_min_=Theta_min-Margin;
    if (Theta_min_ < -pi)
        Theta_min_=Theta_min_+(2*pi);
    end
    Index=(Theta_min_ <= Theta) & (Theta_min >= Theta);
    W(Index)=sin(Factor*(Theta(Index)-Theta_min)).^2;

    Theta_max_=Theta_max+Margin;
    if (Theta_max_ > pi)
        Theta_max_=Theta_max_-(2*pi);
    end
    Index=(Theta_max_ >= Theta) & (Theta_max <= Theta);
    W(Index)=sin(Factor*(Theta(Index)-Theta_max)).^2;
end


function [x,y,z]=Polar2Cart(Theta,Phi,J_l,Mode)
    switch Mode
        case 1
            r=1+0*abs(J_l);
        case 2
            r=1+0.1*abs(J_l);
        case 3
            r=1+0.1*J_l;
    end
    x=r.*cos(Phi).*sin(Theta);
    y=r.*sin(Phi).*sin(Theta);
    z=r.*cos(Theta);
end

function LMa=Load_LMa
    LMa=[0   0   1;
        6   0   0.531085;
        6   5   0.847318;
        10  0   0.265539;
        10  5   -0.846143;
        10  10  0.462094;
        12  0   0.454749;
        12  5   0.469992;
        12  10  0.756513;
        15  0   0;
        15  5   -0.730479;
        15  10  0.553390;
        15  15  0.400200;
        16  0   0.334300;
        16  5   -0.493693;
        16  10  -0.634406;
        16  15  0.491975;
        18  0   0.399497;
        18  5   0.450611;
        18  10  0.360958;
        18  15  0.712083;
        20  0   0.077539;
        20  5   -0.460748;
        20  10  0.747888;
        20  15  -0.231074;
        20  20  0.411056;
        21  0   0;
        21  5   -0.686874;
        21  10  -0.165198;
        21  15  0.589866;
        21  20  0.391117;
        22  0   0.374046;
        22  5   -0.305627;
        22  10  -0.478054;
        22  15  -0.524717;
        22  20  0.512658;
        24  0   0.349347;
        24  5   0.459541;
        24  10  0.312605;
        24  15  0.323072;
        24  20  0.681676;
        25  0   0;
        25  5   -0.287213;
        25  10  0.690434;
        25  15  -0.517311;
        25  20  0.082394;
        25  25  0.407935;
        26  0   0.134096;
        26  5  -0.516658;
        26  10  0.334011;
        26  15  0.585483;
        26  20  -0.295937;
        26  25  0.416113;
        27  0   0;
        27  5   -0.626928;
        27  10  -0.332517;
        27  15  0.039428;
        27  20  0.590696;
        27  25  0.381993;
        28  0   0.396335;
        28  5   -0.161426;
        28  10  -0.387561;
        28  15  -0.426225;
        28  20  -0.454516;
        28  25  0.527640;
        30  0   0.302281;
        30  5   0.454268;
        30  10  0.354171;
        30  15  0.198039;
        30  20  0.331042;
        30  25  0.653810;
        30  30  0.033905;
        ];
end

function J_l=IHrmncSngl(Theta,Phi,l)
    %Y_(l,m)=sqrt(2)*P_(m,l)(cos(Theta),'norm')*cos(m*Phi)
    %P_(m={0,1,...l},x) <-- legendre(l,x)
    %J_(l)(Theta,Phi)=Sum_(m){a_(l,m)S_(l,m)(Theta,Phi)}
    
    LMa=Load_LMa;
    Indices=find(LMa(:,1) == l);
    m_Vector=single(LMa(Indices,2));
    a_Vector=single(LMa(Indices,3));
    Legendre=legendre(l,cos(Theta),'norm');
    J_l=zeros(size(Theta),'single');
    Sqrt_2=sqrt(2);
    for cntr=1:numel(Indices)
        m=m_Vector(cntr);
        a=a_Vector(cntr);
        if l
            if ndims(Legendre) == 4 %3D arrays of theta and phi
                Temp=squeeze(Legendre(m+1,:,:,:));
            elseif ndims(Legendre) == 3 %vectors of theta and phi
                Temp=squeeze(Legendre(m+1,:,:));
            elseif ndims(Legendre) == 2 %vectors of theta and phi
                Temp=squeeze(Legendre(m+1,:))';
            else
                size(Legendre)
            end
        else
            if min(size(Phi)) > 1
                Temp=Legendre;
            else
                Temp=Legendre';
            end
        end

        %Y_lm=sqrt(2)*Temp.*cos(m*Phi);
        %J_l=J_l+Y_lm*a;

        %J_l=J_l+Temp.*cos(m*Phi)*a;
        if mod(l,2)
            %J_l=J_l+Temp.*sin(m*Phi)*(2*a*(1i));
            J_l=J_l+Temp.*sin(m*Phi)*(2*a);
        else
            J_l=J_l+Temp.*cos(m*Phi)*(2*a);
        end
    end
    J_l=real(J_l)*Sqrt_2;
end


function Y=SphHrmnc(Theta,Phi,LMa)
    N=size(LMa,1);
    L=LMa(:,1);
    M=LMa(:,2);
    Y=cell(N,1);
    for cntr=1:N
        Y{cntr}=SphHrmncSngl(Theta,Phi,L(cntr),M(cntr));
    end
end

function Y_lm=SphHrmncSngl(Theta,Phi,l,m)
    %Y_(l,m)=sqrt(2)*P_(m,l)(cos(Theta),'norm')*cos(m*Phi)
    %P_(m={0,1,...l},x) <-- legendre(l,x)
    Legendre=legendre(l,cos(Theta),'norm');
    if isequal(size(Theta,1),1)
        Legendre=Legendre(m+1,:);
    elseif isequal(size(Theta,2),1)
        Legendre=Legendre(m+1,:)';
    else
        Legendre=squeeze(Legendre(m+1,:,:));
    end
    Y_lm=sqrt(2)*Legendre.*cos(m*Phi);
end

function [Theta,Phi]=SphGrid(N)
    Phi_1D=(2*pi)*linspace(0,1,N+1);
    Phi_1D=Phi_1D(1:N);
    Phi_1D=Phi_1D-mean(Phi_1D);
    Theta_1D=(pi)*linspace(0,1,N);
    [Theta,Phi]=meshgrid(Theta_1D,Phi_1D);
end

function []=Initialize()
    close all force;
    clear all;
    clc;
    drawnow;
end



function q=Scaled_q(varargin)
    M=1024;
    switch nargin
        case 0
            E_eV=520;   %Energy in eV
            dPixel_m=75e-6;       %Pixel picth in m
            zD_m=740e-3;  %Distance to camera in m

            Lambda_nm=1240/E_eV;   % in m
            Q0_nm_inv=1/Lambda_nm;
            dq=(dPixel_m/zD_m)*Q0_nm_inv;
            Dq=M*dq;
    end
    Grid_1D=Dq*linspace(-0.5,0.5,M);

    q.x=[min(Grid_1D(:)) max(Grid_1D(:))];
    q.y=[max(Grid_1D(:)) min(Grid_1D(:))];
    q.Step=dq;
    q.Span=Dq;

    q.Corr.x=q.x*(2*M-1)/(M);
    q.Corr.y=q.y*(2*M-1)/(M);
    q.Corr.Step=dq*(2*M-1)/(M);
    q.Corr.Span=Dq*(2*M-1)/(M);

    q.Inverse.Span=M/(2*Dq);
    q.Inverse.Step=1/Dq;
    q.Inverse.x=[-q.Inverse.Span/2 q.Inverse.Span/2];
    q.Inverse.y=[q.Inverse.Span/2 -q.Inverse.Span/2];
    
    [q.X,q.Y]=meshgrid(Grid_1D);
end

function [u_R,u_Theta,u_Phi,u_x,u_y,u_z,OS_1D,u_x_2D,u_y_2D,u_1D,Tri]=u_Grid(N,q,Mesh_GridBar_Flag)
    M=size(q.X,1);
    Span_Half=q.Span/2;
    
    
    
    Span_Half=2;
    
    
    
    
    
    OS_1D=M/N;
    u_1D=linspace(-Span_Half,Span_Half,N)';
    [u_x_2D,u_y_2D]=meshgrid(u_1D,u_1D);
    
    Polar_CartesianBar_Flag=1;
    if Polar_CartesianBar_Flag
        u_R_1D=linspace(Span_Half*0.3,Span_Half,N);
        %u_Theta_1D=acos(linspace(1,-1,N));
        u_Theta_1D=linspace(0,pi,N);
        u_Phi_1D=linspace(0,2*pi,N+1);u_Phi_1D=u_Phi_1D(1:N);
        [u_R,u_Theta,u_Phi]=ndgrid(u_R_1D,u_Theta_1D,u_Phi_1D);
        [u_x,u_y,u_z]=Polar2Cartesian(u_R,u_Theta,u_Phi);
    else
        [u_R,u_Theta,u_Phi]=Cartesian2Polar(u_x,u_y,u_z);
        if Mesh_GridBar_Flag
            [u_x,u_y,u_z]=meshgrid(u_1D,u_1D,u_1D);
        else
            [u_x,u_y,u_z]=ndgrid(u_1D,u_1D,u_1D);
        end
        [u_R,u_Theta,u_Phi]=Cartesian2Polar(u_x,u_y,u_z);
    end
    Tri=DelaunayTri(u_1D,u_1D,u_1D);
    
end
function [u_R,u_Theta,u_Phi]=Cartesian2Polar(u_x,u_y,u_z)
    u_R=sqrt(u_x.^2+u_y.^2+u_z.^2);
    u_Theta=real(acos(u_z./u_R));
    u_Phi=atan2(u_y,u_x);
end
function [u_x,u_y,u_z]=Polar2Cartesian(u_R,u_Theta,u_Phi)
    u_x=u_R.*sin(u_Theta).*cos(u_Phi);
    u_y=u_R.*sin(u_Theta).*sin(u_Phi);
    u_z=u_R.*cos(u_Theta);
end

function Corr2=FastCorr2(varargin)
    %Corr2=xcorr2(I);return;
    I=varargin{1};
    N=size(I,1);
    J=zeros(2*N-1,2*N-1);
    J(1:N,1:N)=I;
    switch nargin
        case 1
            Corr2=real(fftshift(ifft2(abs(fft2(fftshift(J))).^2)));
        case 2
            I2=varargin{2};
            J2=zeros(2*N-1,2*N-1);
            J2(1:N,1:N)=I2;
            Corr2=real(fftshift(ifft2(abs(fft2(fftshift(J)).*fft2(fftshift(J2))))));
    end
end

function J=FFT2(I,n)
    switch n
        case 1
            J=fftshift(fft2(ifftshift(I)));
        case -1
            J=ifftshift(ifft2(fftshift(I)));
    end
end

function N=Norm(I)
    N=norm(I(:));
end


function Euler=IIRZ(N)
    Gamma=pi/5;
    Theta_c=(pi/180)*63.43;
    A=[sin(Theta_c).*cos(Gamma);sin(Theta_c).*sin(Gamma);cos(Theta_c)];
    B=[0;0;1];
    C=[sin(Theta_c);0;cos(Theta_c)];
    CA_Mid=(C+A)/2;


    Theta_1D=acos(linspace(1,cos(Theta_c),N));
    Phi_1D=(Gamma/2)*linspace(0,1,N);
    %Index=randperm(N);
    [Phi,Theta]=meshgrid(Phi_1D,Theta_1D);
    Euler=zeros(3,numel(Theta));
    Euler(2,:)=Theta(:)';
    Euler(3,:)=Phi(:)';
    Index=(Theta(:) > atan(0.618./cos(Gamma-Phi(:))));
    %Theta=Theta(Index);
    %Phi=Phi(Index);
    Euler(1:3,Index)=nan;
    Plot_Flag=1;
    if Plot_Flag
        H=figure;
        [x,y,z]=sphere(N);
        surf(x,y,z);
        shading interp;
        colormap('hot');
        caxis(-1+[0 1e-8]);
        hold on;
        X=sin(Euler(2,:)).*cos(Euler(3,:));
        Y=sin(Euler(2,:)).*sin(Euler(3,:));
        Z=cos(Euler(2,:));
        plot3(X,Y,Z,'*');
        hold off;
        daspect([1 1 1]);
        axis off;
        view(180-72,16);
        pause;close(H);drawnow;
    end
end






%_______________________
%Solid icosahedron stuff
%_______________________

function [IcoGridVol,IcoGrid,H0]=Form_IcosahedralGrids(IcoObjParam)
    N_TriGrid=IcoObjParam.N_TriGrid;
    Radius1=IcoObjParam.Radius1;
    Radius2=IcoObjParam.Radius2;
    Radius1_min=IcoObjParam.Radius1_min;
    Radius=IcoObjParam.Radius;
    Radius_Res=IcoObjParam.Radius_Res;
    N_Radius=IcoObjParam.N_Radius;
    Plot_Grid_Flag=IcoObjParam.Plot_Grid_Flag;
    AxisOff_Flag=IcoObjParam.AxisOff_Flag;
    AppendSupport_Flag=IcoObjParam.AppendSupport_Flag;
    A=IcoObjParam.A;
    B=IcoObjParam.B;
    C=IcoObjParam.C;
    F=IcoObjParam.F;


    Ico=Make_TriangleGrid(B,C,A,Radius,N_TriGrid);
    IcoGrid=Make_IcoShellGrid(B,C,A,Radius,N_TriGrid);
    IcoGridVol=Make_IcoVolGrid(B,C,A,Radius1,Radius2,N_TriGrid,N_Radius);
    IcoGridVolSupportIn=struct('x',[],'y',[],'z',[]);
    IcoGridVolSupportOut=struct('x',[],'y',[],'z',[]);
    if AppendSupport_Flag
        N_R_S=2;
        IcoGridVolSupportIn=Make_IcoVolGrid(B,C,A,...
            max(0,Radius1-Radius_Res*N_R_S),Radius1-Radius_Res,N_TriGrid,N_R_S);
        IcoGridVolSupportOut=Make_IcoVolGrid(B,C,A,...
            Radius2+Radius_Res,Radius2+Radius_Res*N_R_S,N_TriGrid,N_R_S);
    end
    IcoGridVol=Append_Grids(IcoGridVol,IcoGridVolSupportIn,IcoGridVolSupportOut);
    if Plot_Grid_Flag
        H0=figure;
        subplot(1,4,1);hold on;plot3(A(1),A(2),A(3),'*');plot3(B(1),B(2),B(3),'*');plot3(C(1),C(2),C(3),'*');plot3(F(1),F(2),F(3),'r*');plot3([0 F(1)],[0 F(2)],[0 F(3)],'r');plot3(B(1)-A(1),B(2)-A(2),B(3)-A(3),'k');plot3([B(1) A(1) C(1)],[B(2) A(2) C(2)],[B(3) A(3) C(3)],'k');plot3([B(1) A(1) C(1) B(1)],[B(2) A(2) C(2) B(2)],[B(3) A(3) C(3) B(3)],'k');plot3(0,0,0,'+r');view(148,32);daspect([1 1 1]);xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);hold off;title('S^2/20 flat boundary');camva(4.6);if AxisOff_Flag; axis off; else xlabel('x');ylabel('y');zlabel('z');end;alpha(0.5);alpha('color');
        subplot(1,4,2);plot3(Ico.x,Ico.y,Ico.z,'.k');daspect([1 1 1]);view(148,32);xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);camva(4.6);title('S^2/20 flat grid');if AxisOff_Flag; axis off; else xlabel('x');ylabel('y');zlabel('z');end;alpha(0.5);alpha('color');
        subplot(1,4,3);plot3(IcoGrid.x,IcoGrid.y,IcoGrid.z,'.k');daspect([1 1 1]);view(148,32);xlim(Radius2*[-1 1]);ylim(Radius2*[-1 1]);zlim(Radius2*[-1 1]);camva(4.6);title('S^2 flat grid');if AxisOff_Flag; axis off; else xlabel('x');ylabel('y');zlabel('z');end;alpha(0.5);alpha('color');
        %subplot(1,4,3);patch(IcoGrid.x,IcoGrid.y,IcoGrid.z,1);daspect([1 1 1]);view(148,32);xlim(Radius2*[-1 1]);ylim(Radius2*[-1 1]);zlim(Radius2*[-1 1]);camva(4.6);title('S^2 flat grid');if AxisOff_Flag; axis off; else xlabel('x');ylabel('y');zlabel('z');end;
        subplot(1,4,4);plot3(IcoGridVol.x,IcoGridVol.y,IcoGridVol.z,'.k');daspect([1 1 1]);view(148,32);xlim(Radius2*[-1 1]);ylim(Radius2*[-1 1]);zlim(Radius2*[-1 1]);camva(4.6);title('S^2xR (thick shell) 3D grid');if AxisOff_Flag; axis off; else xlabel('x');ylabel('y');zlabel('z');end;alpha(0.5);alpha('color');
    else
        H0=[];
    end
end
function Grid=Append_Grids(MainGrid,GridIn,GridOut)
    Factor=1;   %1/nan
    Grid=struct(...
        'x',[GridIn.x;MainGrid.x;GridOut.x],...
        'y',[GridIn.y;MainGrid.y;GridOut.y],...
        'z',[GridIn.z;MainGrid.z;GridOut.z],...
        'F',[Factor*zeros(size(GridIn.z));ones(size(MainGrid.z));...
            Factor*zeros(size(GridOut.z))]);
        %size(MainGrid.x)
        %size(Grid.x)
end
function []=Superimpose_Triangle_Glob()
    hold on;
    plot3(A(1),A(2),A(3),'*')
    plot3(B(1),B(2),B(3),'*')
    plot3(C(1),C(2),C(3),'*')
    plot3(F(1),F(2),F(3),'r*')
    plot3([0 F(1)],[0 F(2)],[0 F(3)],'r');
    plot3(B(1)-A(1),B(2)-A(2),B(3)-A(3),'k')
    plot3([B(1) A(1) C(1)],[B(2) A(2) C(2)],[B(3) A(3) C(3)],'k')
    plot3([B(1) A(1) C(1) B(1)],[B(2) A(2) C(2) B(2)],[B(3) A(3) C(3) B(3)],'k')
    plot3(Bp(1),Bp(2),Bp(3),'*k')
    plot3(Ap(1),Ap(2),Ap(3),'*k')
    plot3(Cp(1),Cp(2),Cp(3),'*k')
    plot3(0,0,0,'+r')
    Adjust_Axes();
end

function [IcoFunc,FuncGrid]=Make_IcoVolFunc(IcoGridVol,N_R,N_Theta,N_Phi)
    Margin=0;
    OverSampling=1;
    N_S2=N_Theta*N_Phi;
    N_Total=N_R*N_S2;

    R_IcoGridVol=sqrt(IcoGridVol.x(:).^2+IcoGridVol.y(:).^2+IcoGridVol.z(:).^2);
    R_min_data=min(R_IcoGridVol);   %R1*cos(Alpha)
    R_Max_data=max(R_IcoGridVol);   %R2
    R_min=R_min_data*(1-Margin);
    R_Max=R_Max_data*OverSampling;
    clear R_IcoGridVol;
    
    %FuncGrid coordinate
    R_1D=linspace(R_min,R_Max,N_R);
    Theta_1D=linspace(0,pi,N_Theta);
    Phi_1D=linspace(0,2*pi,N_Phi);
    [Theta_2D,Phi_2D]=meshgrid(Theta_1D,Phi_1D);
    Xn=sin(Theta_2D).*cos(Phi_2D);
    Yn=sin(Theta_2D).*sin(Phi_2D);
    Zn=cos(Theta_2D);
    %R_XYZn=min(sqrt(Xn.^2+Yn.^2+Zn.^2));
    FuncGrid=struct('R_1D',R_1D,'Theta_2D',Theta_2D,'Phi_2D',Phi_2D);

    %Model making w/ Object data
    Model=TriScatteredInterp(IcoGridVol.x,IcoGridVol.y,IcoGridVol.z,...
        IcoGridVol.F,'nearest');
        %ones(size(IcoGridVol.x)),'natural');

    %Sampling model in the FuncGrid coordinate
    IcoFunc=cell(N_R,1);
    X=zeros(N_Total,1);
    Y=zeros(N_Total,1);
    Z=zeros(N_Total,1);
    Plot_Flag=1;
    if Plot_Flag
        H=figure;
    else
        WB=waitbar(0,'Calculating f{\bf_r}(\theta,\phi) functions');
    end
    for r_cntr=1:N_R
        R_Shell=R_1D(r_cntr);
        x=R_Shell*Xn;
        y=R_Shell*Yn;
        z=R_Shell*Zn;
        X((1:N_S2)+(r_cntr-1)*N_S2)=x;
        Y((1:N_S2)+(r_cntr-1)*N_S2)=y;
        Z((1:N_S2)+(r_cntr-1)*N_S2)=z;

        Temp=Model(x,y,z);

        %A=min(sqrt(x(:).^2+y(:).^2+z(:).^2));
        C=R_min_data;
        D=R_Max_data;
        E=sum(isnan(Temp(:)));
        %disp(['sm:' num2str(A) ', dm:' num2str(C) ', sM:' num2str(B) ', dM:' num2str(D)])
        %disp(['s:' num2str(A) ', R_XYZn:' num2str(R_XYZn) ', d:' num2str(C) ' to ' num2str(D) ', n(NaN):' num2str(E)])
        disp(['f(r_' num2str(r_cntr) ') @r=' num2str(R_Shell) ', d:' num2str(C) ' to ' num2str(D) ', n(NaN):' num2str(E)])
        
        Index=~isnan(Temp);
        if Plot_Flag 
            %imagesc(Temp);daspect([1 1 1]);axis tight; axis off;title(num2str(r_cntr));colorbar;pause;
            scatter3(x(Index),y(Index),z(Index),20,Temp(Index));daspect([1 1 1]);axis tight; axis off;title(['f(r_{' num2str(r_cntr) '}) @r=' num2str(R_Shell) ', d:' num2str(C) ' to ' num2str(D) ', n(NaN):' num2str(E)]);colorbar;view(-90,32);alpha(0.5);alpha('color');drawnow;
        else
            WB=waitbar(r_cntr/N_R,WB);drawnow;
        end
        Temp(~Index)=0;
        IcoFunc{r_cntr}=Temp;
    end
    if ~Plot_Flag
        close(WB);drawnow;
    end
end
function Ico_Grid=Make_IcoShellGrid(B,C,A,Radius,N_TriGrid)
    X=Make_TriangleGrid(B,C,A,Radius,N_TriGrid);
    N_Face=numel(X.x);
    Zero=zeros(20*N_Face,1);
    Ico_Grid=struct('x',Zero,'y',Zero,'z',Zero);
    R=IcosSymRot20();   %Only for EMAN configuration
    for Cntr=1:20
        Temp=R{Cntr};
        B_=Temp*B;
        C_=Temp*C;
        A_=Temp*A;
        X=Make_TriangleGrid(B_,C_,A_,Radius,N_TriGrid);
        Ico_Grid.x((1:N_Face)+(Cntr-1)*N_Face)=X.x;
        Ico_Grid.y((1:N_Face)+(Cntr-1)*N_Face)=X.y;
        Ico_Grid.z((1:N_Face)+(Cntr-1)*N_Face)=X.z;
    end
    Ico_Grid=RemoveRedundancy(Ico_Grid);
end    
function Grid=RemoveRedundancy(Grid)
    if size(Grid.x,1) > 1
        Temp=[Grid.x,Grid.y,Grid.z];
        Temp=unique(Temp,'rows');
        Grid.x=Temp(:,1);
        Grid.y=Temp(:,2);
        Grid.z=Temp(:,3);
    else
        Temp=[Grid.x;Grid.y;Grid.z]';
        Temp=unique(Temp,'rows');
        Grid.x=Temp(1,:);
        Grid.y=Temp(2,:);
        Grid.z=Temp(3,:);
    end
end
function X=Make_TriangleGrid(B,C,A,Radius,N)
    [Lambda,Mu]=meshgrid(Radius*linspace(0,1,N));
    Index=(Lambda+ Mu <= 1);
    Lambda=Lambda(Index);
    Mu=Mu(Index);
    X=cell(3,1);
    for cntr_=1:3
        X{cntr_}=B(cntr_)+Lambda*(A(cntr_)-B(cntr_))+Mu*(C(cntr_)-B(cntr_));
    end
    X=cell2struct(X,{'x','y','z'},1);
end

function [V_parallel,V_EMAN,E_parallel,E_EMAN]=Full_Icosahedron()
    P=(1+sqrt(5))/2;
    Q=sqrt(1+P^2);

    Vx=[0  0  0  0  1  1 -1 -1  P  P -P -P]/Q;
    Vy=[1  1 -1 -1  P -P  P -P  0  0  0  0]/Q;
    Vz=[P -P  P -P  0  0  0  0  1  -1 1 -1]/Q;
    
    
    V_parallel=struct('x',Vx,'y',Vy,'z',Vz);

    %Rotation matrix calculation
    B=[0;0;1];C=[0.8944;0;0.4472];A=[0.2764;0.8507;0.4472];BCA=[B,C,A];
    Bp=[0;1;P]/Q;Cp=[0;-1;P]/Q;Ap=[P;0;1]/Q;BCA_p=[Bp,Cp,Ap];           
    %R=(BCA*inv(BCA_p));
    f=0;
    c=0.8944*Q/(2*P);
    a=(0.2764*Q-c)/P;
    d=(0.8507*Q-f)/P;
    i=1.4472*Q/(2*P);
    g=(0.4472*Q-i)/P;
    b=-c*P;
    e=-f*P;
    h=Q-i*P;
    R=[ a,b,c;
        d,e,f;
        g,h,i];
    [U,~,V] = svd(R);
    R=U*V';
    
    for cntr=1:12
        Temp=R*[Vx(cntr);Vy(cntr);Vz(cntr)];
        Vx(cntr)=Temp(1);
        Vy(cntr)=Temp(2);
        Vz(cntr)=Temp(3);
    end
    V_EMAN=struct('x',Vx,'y',Vy,'z',Vz);
    E_parallel=Find_Edges(V_parallel);
    E_EMAN=Find_Edges(V_EMAN);
end
function []=Plot_Vertex(V)
    plot3(V.x,V.y,V.z,'*');
    daspect([1 1 1]);
    box on;
    xlabel('x');ylabel('y');zlabel('z');
    xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
end
function []=Plot_Edges(V,E)
    for cntr1=1:12
        for cntr2=1:5
            plot3([V.x(E(cntr1,cntr2)) V.x(cntr1)],...
                [V.y(E(cntr1,cntr2)) V.y(cntr1)],...
                [V.z(E(cntr1,cntr2)) V.z(cntr1)]);
        end
    end
end
function Plot_Vertex_Edges(V)
    IsHold=ishold();
    Plot_Vertex(V);
    hold on;
    Plot_Edges(V,Find_Edges(V));
    if ~IsHold
        hold off;
    end
end
function Icosahedron=Load_Icosahedron()
    [V_parallel,V_EMAN,E_parallel,E_EMAN]=Full_Icosahedron();
    Omega=2*pi/5;
    Theta_c=acos(cos(Omega)/(1-cos(Omega)));
    Alpha=acos(1/(sqrt(3)*tan(Omega/2)));

    A=[sin(Theta_c)*cos(Omega);sin(Theta_c)*sin(Omega);cos(Theta_c)];
    B=[0;0;1];
    C=[sin(Theta_c);0;cos(Theta_c)];
    F=(A+B+C)/3;

    CA=A-C;
    BF=F-B;
    
    Bp=(A+C)/2;
    Ap=(B+C)/2;
    Cp=(A+B)/2;

    r_inscribed=norm(F);
    rp_Max=norm(BF);
    rp_min=rp_Max/2;

    w_Max=0;
    w_min=-0.5;
    Eta_Max=1;
    Eta_min=0;

    r_Max=norm([r_inscribed rp_Max]);
    r_min=norm([r_inscribed rp_min]);
    Tau_min=pi/2;
    Tau_Max=Tau_min+2*pi/3;
    
    Icosahedron=Var2Struct(Omega,Theta_c,Alpha,B,C,A,F,CA,BF,Bp,Ap,Cp,...
        r_inscribed,rp_Max,rp_min,w_Max,w_min,Eta_Max,Eta_min,...
        r_Max,r_min,Tau_Max,Tau_min,V_parallel,V_EMAN,E_parallel,E_EMAN);

end
function Edges=Find_Edges(V)
    Edges=zeros(12,5);
    for cntr1=1:12
        Distance=(V.x-V.x(cntr1)).^2+(V.y-V.y(cntr1)).^2+(V.z-V.z(cntr1)).^2;
        [~,I]=sort(Distance,'ascend');
        Edges(cntr1,:)=I(1+(1:5));
    end
end

function Icosahedron=Var2Struct(Omega,Theta_c,Alpha,B,C,A,F,CA,BF,...
    Bp,Ap,Cp,r_inscribed,rp_Max,rp_min,w_Max,w_min,Eta_Max,Eta_min,...
        r_Max,r_min,Tau_Max,Tau_min,V_parallel,V_EMAN,E_parallel,E_EMAN)


    Icosahedron=struct(...
        'Omega',Omega,...
        'Theta_c',Theta_c,...
        'Alpha',Alpha,...
        'B',B,...
        'C',C,...
        'A',A,...
        'F',F,...
        'CA',CA,...
        'BF',BF,...
        'Bp',Bp,...
        'Cp',Cp,...
        'Ap',Ap,...
        'r_inscribed',r_inscribed,...
        'rp_Max',rp_Max,...
        'rp_min',rp_min,...
        'w_Max',w_Max,...
        'w_min',w_min,...
        'Eta_Max',Eta_Max,...
        'Eta_min',Eta_min,...
        'r_Max',r_Max,...
        'r_min',r_min,...
        'Tau_Max',Tau_Max,...
        'Tau_min',Tau_min,...
        'V_parallel',V_parallel,...
        'V_EMAN',V_EMAN,...
        'E_parallel',E_parallel,...
        'E_EMAN',E_EMAN...
        );
end
function [Omega,Theta_c,Alpha,B,C,A,F,CA,BF,Bp,Ap,Cp...
    r_inscribed,rp_Max,rp_min,w_Max,w_min,Eta_Max,Eta_min,...
    r_Max,r_min,Tau_Max,Tau_min,V_parallel,V_EMAN,E_parallel,E_EMAN]=...
    Struct2Var(Icosahedron)

        Omega=Icosahedron.Omega;
        Theta_c=Icosahedron.Theta_c;
        Alpha=Icosahedron.Alpha;
        B=Icosahedron.B;
        C=Icosahedron.C;
        A=Icosahedron.A;
        F=Icosahedron.F;
        CA=Icosahedron.CA;
        BF=Icosahedron.BF;
        Bp=Icosahedron.Bp;
        Cp=Icosahedron.Cp;
        Ap=Icosahedron.Ap;
        r_inscribed=Icosahedron.r_inscribed;
        rp_Max=Icosahedron.rp_Max;
        rp_min=Icosahedron.rp_min;
        w_Max=Icosahedron.w_Max;
        w_min=Icosahedron.w_min;
        Eta_Max=Icosahedron.Eta_Max;
        Eta_min=Icosahedron.Eta_min;
        r_Max=Icosahedron.r_Max;
        r_min=Icosahedron.r_min;
        Tau_Max=Icosahedron.Tau_Max;
        Tau_min=Icosahedron.Tau_min;
        V_parallel=Icosahedron.V_parallel;
        V_EMAN=Icosahedron.V_EMAN;
        E_parallel=Icosahedron.E_parallel;
        E_EMAN=Icosahedron.E_EMAN;
end

function Ico_Grid=Make_Icosahedral_Grid(Icosahedron)
    [Omega,Theta_c,Alpha,B,C,A,F,CA,BF,Bp,Ap,Cp,r_inscribed,rp_Max,rp_min,...
        w_Max,w_min,Eta_Max,Eta_min,r_Max,r_min,Tau_Max,Tau_min,...
        V_parallel,V_EMAN,E_parallel,E_EMAN]=Struct2Var(Icosahedron);
    N=100;
    w=linspace(w_min,w_Max,N);
    Eta=linspace(Eta_min,Eta_Max,N+1);
    Eta=Eta(1:N);
    [w,Eta]=meshgrid(w,Eta);
    Index=(Eta > 2*abs(w)) & (Eta <= (2/3)*(1-w));
    w=w(Index);
    Eta=Eta(Index);
    clear Index;
    Grid_w_Eta=cell(3,1);
    for cntr=1:3
        Grid_w_Eta{cntr}=B(cntr)+w(:)*(CA(cntr))+Eta(:)*(Bp(cntr)-B(cntr));
    end
    
    Tau=linspace(Tau_min,Tau_Max,N);
    rp=linspace(0,rp_Max,N);
    [Tau,rp]=meshgrid(Tau,rp);
    Index=(rp/rp_min <= 1) | (...
        (rp_min./rp > (cos(Tau-5*pi/6))) & 1 ...
       );
    Tau=Tau(Index);
    rp=rp(Index);
    clear Index;
    Grid_rp_Tau=cell(3,1);
    for cntr=1:3
        Grid_rp_Tau{cntr}=F(cntr)+rp(:).*cos(Tau(:))*(CA(cntr)/norm(CA))+...
            rp(:).*sin(Tau(:))*(-BF(cntr)/norm(BF));
    end
    Ico_Grid={Grid_w_Eta,Grid_rp_Tau,Tau,rp,w,Eta,numel(Eta),numel(Tau)};
end
function Ico_Grid=Plot_Icosahedron(Icosahedron,Ico_Grid)
    [Omega,Theta_c,Alpha,B,C,A,F,CA,BF,Bp,Ap,Cp,r_inscribed,rp_Max,rp_min,...
        w_Max,w_min,Eta_Max,Eta_min,r_Max,r_min,Tau_Max,Tau_min,...
        V_parallel,V_EMAN,E_parallel,E_EMAN]=Struct2Var(Icosahedron);
    
    Grid_w_Eta=Ico_Grid{1};
    Grid_rp_Tau=Ico_Grid{2};
    

    figure;
    subplot(221);
    plot3(Grid_w_Eta{1},Grid_w_Eta{2},Grid_w_Eta{3},'k.')
    Superimpose_Triangle()
    subplot(222);
    plot3(Grid_rp_Tau{1},Grid_rp_Tau{2},Grid_rp_Tau{3},'k.')
    Superimpose_Triangle()
    subplot(223);
    plot3(Grid_w_Eta{1},Grid_w_Eta{2},Grid_w_Eta{3},'k.')
    Superimpose_Icosahedron(V_EMAN)
    subplot(224);
    plot3(Grid_rp_Tau{1},Grid_rp_Tau{2},Grid_rp_Tau{3},'k.')
    Superimpose_Icosahedron(V_EMAN)

    function []=Superimpose_Triangle()
        hold on;
        plot3(A(1),A(2),A(3),'*')
        plot3(B(1),B(2),B(3),'*')
        plot3(C(1),C(2),C(3),'*')
        plot3(F(1),F(2),F(3),'r*')
        plot3([0 F(1)],[0 F(2)],[0 F(3)],'r');
        plot3(B(1)-A(1),B(2)-A(2),B(3)-A(3),'k')
        plot3([B(1) A(1) C(1)],[B(2) A(2) C(2)],[B(3) A(3) C(3)],'k')
        plot3([B(1) A(1) C(1) B(1)],[B(2) A(2) C(2) B(2)],[B(3) A(3) C(3) B(3)],'k')
        plot3(Bp(1),Bp(2),Bp(3),'*k')
        plot3(Ap(1),Ap(2),Ap(3),'*k')
        plot3(Cp(1),Cp(2),Cp(3),'*k')
        plot3(0,0,0,'+r')
        Adjust_Axes();
    end
    function []=Superimpose_Icosahedron(V)
        hold on;
        plot3(F(1),F(2),F(3),'r*')
        plot3([0 F(1)],[0 F(2)],[0 F(3)],'r');
        plot3(A(1),A(2),A(3),'*')
        plot3(B(1),B(2),B(3),'*')
        plot3(C(1),C(2),C(3),'*')
        plot3(Bp(1),Bp(2),Bp(3),'*k')
        plot3(Ap(1),Ap(2),Ap(3),'*k')
        plot3(Cp(1),Cp(2),Cp(3),'*k')
        plot3(0,0,0,'+r')
        Plot_Edges(V,Find_Edges(V));
        %Plot_Vertex_Edges(V)
        Adjust_Axes();
    end
    function Adjust_Axes()
        view(155,50)
        daspect([1 1 1])
        xlim([-1 1])
        ylim([-1 1])
        zlim([-1 1])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        set(gca,'Xdir','normal')
        set(gca,'Zdir','normal')
        box on;
    end
end



%% Legendre Quadrature
function legendre_rule ( order, a, b, filename )

%*****************************************************************************80
%
%% LEGENDRE_RULE generates a Gauss-Legendre rule.
%
%  Discussion:
%
%    This program computes a standard Gauss-Legendre quadrature rule
%    and writes it to a file.
%
%    The user specifies:
%    * the ORDER (number of points) in the rule;
%    * A, the left endpoint;
%    * B, the right endpoint;
%    * FILENAME, the root name of the output files.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 February 2010
%
%  Author:
%
%    John Burkardt
%
  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'LEGENDRE_RULE\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Compute a Gauss-Legendre rule for approximating\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '    Integral ( A <= x <= B ) f(x) dx\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  of order ORDER.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  The user specifies ORDER, A, B, and FILENAME.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  ORDER is the number of points;\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  A is the left endpoint;\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  B is the right endpoint;\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  FILENAME is used to generate 3 files:\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '    filename_w.txt - the weight file\n' );
  fprintf ( 1, '    filename_x.txt - the abscissa file.\n' );
  fprintf ( 1, '    filename_r.txt - the region file.\n' );
%
%  Initialize the parameters.
%
  alpha = 0.0;
  beta = 0.0;
%
%  Get ORDER.
%
  if ( 1 <= nargin )
    order = str2num ( order );
  else
    order = input ( '  Enter the rule order ORDER.' );
  end
%
%  Get A.
%
  if ( 2 <= nargin )
    a = str2num ( a );
  else
    a = input ( '  Enter A.' );
  end
%
%  Get B.
%
  if ( 3 <= nargin )
    b = str2num ( b );
  else
    b = input ( '  Enter B.' );
  end
%
%  Get FILENAME.
%
  if ( 4 <= nargin )

  else
    fprintf ( 1,  '\n' );
    fprintf ( 1,  '  FILENAME specifies the ''root name'' of the quadrature files).\n' );
    filename = input ( '  Enter the value of FILENAME as a quoted string:' );
  end
%
%  Input summary.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, '  ORDER = %d\n', order );
  fprintf ( 1, '  A = %f\n', a );
  fprintf ( 1, '  B = %f\n', b );
  fprintf ( 1, '  FILENAME = "%s".\n', filename );
%
%  Construct the rule.
%
  kind = 1;
  [ x, w ] = cgqf ( order, kind, alpha, beta, a, b );
%
%  Write the rule.
%
  r = [ a, b ]';
  rule_write ( order, filename, x, w, r );
%
%  Terminate.
%
  fprintf ( 1,  '\n' );
  fprintf ( 1,  'LEGENDRE_RULE:\n' );
  fprintf ( 1,  '  Normal end of execution.\n' );
  fprintf ( 1,  '\n' );
  timestamp ( );

  return
end
function [ t, wts ] = cdgqf ( nt, kind, alpha, beta )

%*****************************************************************************80
%
%% CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
%
%  Discussion:
%
%    This routine computes all the knots and weights of a Gauss quadrature
%    formula with a classical weight function with default values for A and B,
%    and only simple knots.
%
%    There are no moments checks and no printing is done.
%
%    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 January 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Parameters:
%
%    Input, integer NT, the number of knots.
%
%    Input, integer KIND, the rule.
%    1, Legendre,             (a,b)       1.0
%    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
%    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
%    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
%    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
%    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
%    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
%    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
%
%    Input, real ALPHA, the value of Alpha, if needed.
%
%    Input, real BETA, the value of Beta, if needed.
%
%    Output, real T(NT), the knots.
%
%    Output, real WTS(NT), the weights.
%
  parchk ( kind, 2 * nt, alpha, beta );
%
%  Get the Jacobi matrix and zero-th moment.
%
  [ aj, bj, zemu ] = class_matrix ( kind, nt, alpha, beta );
%
%  Compute the knots and weights.
%
  [ t, wts ] = sgqf ( nt, aj, bj, zemu );

  return
end
function [ t, wts ] = cgqf ( nt, kind, alpha, beta, a, b )

%*****************************************************************************80
%
%% CGQF computes knots and weights of a Gauss quadrature formula.
%
%  Discussion:
%
%    The user may specify the interval (A,B).
%
%    Only simple knots are produced.
%
%    The user may request that the routine print the knots and weights,
%    and perform a moment check.
%
%    Use routine EIQFS to evaluate this quadrature formula.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 February 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Parameters:
%
%    Input, integer NT, the number of knots.
%
%    Input, integer KIND, the rule.
%    1, Legendre,             (a,b)       1.0
%    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
%    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
%    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
%    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
%    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
%    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
%    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
%    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
%
%    Input, real ALPHA, the value of Alpha, if needed.
%
%    Input, real BETA, the value of Beta, if needed.
%
%    Input, real A, B, the interval endpoints.
%
%    Output, real T(NT), the knots.
%
%    Output, real WTS(NT), the weights.
%

%
%  Compute the Gauss quadrature formula for default values of A and B.
%
  [ t, wts ] = cdgqf ( nt, kind, alpha, beta );
%
%  All knots have multiplicity = 1.
%
  mlt = zeros(nt,1);
  mlt(1:nt) = 1;
%
%  NDX(I) = I.
%
  ndx = ( 1 : nt );
%
%  Scale the quadrature rule.
%
  [ t, wts ] = scqf ( nt, t, mlt, wts, nt, ndx, kind, alpha, beta, a, b );

  return
end
function [ aj, bj, zemu ] = class_matrix ( kind, m, alpha, beta )

%*****************************************************************************80
%
%% CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
%
%  Discussion:
%
%    This routine computes the diagonal AJ and subdiagonal BJ
%    elements of the order M tridiagonal symmetric Jacobi matrix
%    associated with the polynomials orthogonal with respect to
%    the weight function specified by KIND.
%
%    For weight functions 1-7, M elements are defined in BJ even
%    though only M-1 are needed.  For weight function 8, BJ(M) is
%    set to zero.
%
%    The zero-th moment of the weight function is returned in ZEMU.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 January 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Parameters:
%
%    Input, integer KIND, the rule.
%    1, Legendre,             (a,b)       1.0
%    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
%    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
%    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
%    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
%    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
%    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
%    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
%
%    Input, integer M, the order of the Jacobi matrix.
%
%    Input, real ALPHA, the value of Alpha, if needed.
%
%    Input, real BETA, the value of Beta, if needed.
%
%    Output, real AJ(M), BJ(M), the diagonal and subdiagonal
%    of the Jacobi matrix.
%
%    Output, real ZEMU, the zero-th moment.
%
  temp = eps;

  parchk ( kind, 2 * m - 1, alpha, beta );

  temp2 = 0.5;

  if ( 500.0 * temp < abs ( ( gamma ( temp2 ) )^2 - pi ) )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'CLASS - Fatal error!\n' );
    fprintf ( 1, '  Gamma function does not match machine parameters.\n' );
    error ( 'CLASS - Fatal error!' );
  end

  bj = zeros(m,1);
  aj = zeros(m,1);

  if ( kind == 1 )

    ab = 0.0;

    zemu = 2.0 / ( ab + 1.0 );

    aj(1:m) = 0.0;

    for i = 1 : m
      abi = i + ab * mod ( i, 2 );
      abj = 2 * i + ab;
      bj(i) = abi * abi / ( abj * abj - 1.0 );
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 2 )

    zemu = pi;

    aj(1:m) = 0.0;

    bj(1) =  sqrt ( 0.5 );
    bj(2:m) = 0.5;

  elseif ( kind == 3 )

    ab = alpha * 2.0;
    zemu = 2.0^( ab + 1.0 ) * gamma ( alpha + 1.0 )^2 ...
      / gamma ( ab + 2.0 );

    aj(1:m) = 0.0;
    bj(1) = 1.0 / ( 2.0 * alpha + 3.0 );
    for i = 2 : m
      bj(i) = i * ( i + ab ) / ( 4.0 * ( i + alpha )^2 - 1.0 );
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 4 )

    ab = alpha + beta;
    abi = 2.0 + ab;
    zemu = 2.0^( ab + 1.0 ) * gamma ( alpha + 1.0 ) ...
      * gamma ( beta + 1.0 ) / gamma ( abi );
    aj(1) = ( beta - alpha ) / abi;
    bj(1) = 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) ...
      / ( ( abi + 1.0 ) * abi * abi );
    a2b2 = beta * beta - alpha * alpha;

    for i = 2 : m
      abi = 2.0 * i + ab;
      aj(i) = a2b2 / ( ( abi - 2.0 ) * abi );
      abi = abi^2;
      bj(i) = 4.0 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) ...
        / ( ( abi - 1.0 ) * abi );
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 5 )

    zemu = gamma ( alpha + 1.0 );

    for i = 1 : m
      aj(i) = 2.0 * i - 1.0 + alpha;
      bj(i) = i * ( i + alpha );
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 6 )

    zemu = gamma ( ( alpha + 1.0 ) / 2.0 );

    aj(1:m) = 0.0;

    for i = 1 : m
      bj(i) = ( i + alpha * mod ( i, 2 ) ) / 2.0;
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 7 )

    ab = alpha;
    zemu = 2.0 / ( ab + 1.0 );

    aj(1:m) = 0.0;

    for i = 1 : m
      abi = i + ab * mod(i,2);
      abj = 2 * i + ab;
      bj(i) = abi * abi / ( abj * abj - 1.0 );
    end
    bj(1:m) =  sqrt ( bj(1:m) );

  elseif ( kind == 8 )

    ab = alpha + beta;
    zemu = gamma ( alpha + 1.0 ) * gamma ( - ( ab + 1.0 ) ) ...
      / gamma ( - beta );
    apone = alpha + 1.0;
    aba = ab * apone;
    aj(1) = - apone / ( ab + 2.0 );
    bj(1) = - aj(1) * ( beta + 1.0 ) / ( ab + 2.0 ) / ( ab + 3.0 );
    for i = 2 : m
      abti = ab + 2.0 * i;
      aj(i) = aba + 2.0 * ( ab + i ) * ( i - 1 );
      aj(i) = - aj(i) / abti / ( abti - 2.0 );
    end

    for i = 2 : m - 1
      abti = ab + 2.0 * i;
      bj(i) = i * ( alpha + i ) / ( abti - 1.0 ) * ( beta + i ) ...
        / ( abti^2 ) * ( ab + i ) / ( abti + 1.0 );
    end

    bj(m) = 0.0;
    bj(1:m) =  sqrt ( bj(1:m) );

  end

  return
end
function [ d, z ] = imtqlx ( n, d, e, z )

%*****************************************************************************80
%
%% IMTQLX diagonalizes a symmetric tridiagonal matrix.
%
%  Discussion:
%
%    This routine is a slightly modified version of the EISPACK routine to
%    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
%
%    The authors thank the authors of EISPACK for permission to use this
%    routine.
%
%    It has been modified to produce the product Q' * Z, where Z is an input
%    vector and Q is the orthogonal matrix diagonalizing the input matrix.
%    The changes consist (essentialy) of applying the orthogonal transformations
%    directly to Z as they are generated.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 January 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%    Roger Martin, James Wilkinson,
%    The Implicit QL Algorithm,
%    Numerische Mathematik,
%    Volume 12, Number 5, December 1968, pages 377-383.
%
%  Parameters:
%
%    Input, integer N, the order of the matrix.
%
%    Input, real D(N), the diagonal entries of the matrix.
%
%    Input, real E(N), the subdiagonal entries of the
%    matrix, in entries E(1) through E(N-1). 
%
%    Input, real Z(N), a vector to be operated on.
%
%    Output, real D(N), the diagonal entries of the diagonalized matrix.
%
%    Output, real Z(N), the value of Q' * Z, where Q is the matrix that 
%    diagonalizes the input symmetric tridiagonal matrix.
%
  itn = 30;

  prec = eps;

  if ( n == 1 )
    return
  end

  e(n) = 0.0;

  for l = 1 : n

    j = 0;

    while ( 1 )

      for m = l : n

        if ( m == n )
          break
        end

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) )
          break
        end

      end

      p = d(l);

      if ( m == l )
        break
      end

      if ( j == itn )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'IMTQLX - Fatal error!\n' );
        fprintf ( 1, '  Iteration limit exceeded.\n' );
        error ( 'IMTQLX - Fatal error!' );
      end

      j = j + 1;
      g = ( d(l+1) - p ) / ( 2.0 * e(l) );
      r =  sqrt ( g * g + 1.0 );
      g = d(m) - p + e(l) / ( g + r8_sign ( g ) * abs ( r ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ii = 1 : mml

        i = m - ii;
        f = s * e(i);
        b = c * e(i);

        if ( abs ( f ) >= abs ( g ) )
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e(i+1) = f * r;
          s = 1.0 / r;
          c = c * s;
        else
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e(i+1) = g * r;
          c = 1.0 / r;
          s = s * c;
        end

        g = d(i+1) - p;
        r = ( d(i) - g ) * s + 2.0 * c * b;
        p = s * r;
        d(i+1) = g + p;
        g = c * r - b;
        f = z(i+1);
        z(i+1) = s * z(i) + c * f;
        z(i) = c * z(i) - s * f;

      end

      d(l) = d(l) - p;
      e(l) = g;
      e(m) = 0.0;

    end

  end

  for ii = 2 : n

     i = ii - 1;
     k = i;
     p = d(i);

     for j = ii : n
       if ( d(j) < p )
         k = j;
         p = d(j);
       end
     end

     if ( k ~= i )
       d(k) = d(i);
       d(i) = p;
       p = z(i);
       z(i) = z(k);
       z(k) = p;
     end

  end

  return
end
function parchk ( kind, m, alpha, beta )

%*****************************************************************************80
%
%% PARCHK checks parameters ALPHA and BETA for classical weight functions.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 January 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Parameters:
%
%    Input, integer KIND, the rule.
%    1, Legendre,             (a,b)       1.0
%    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
%    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
%    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
%    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
%    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
%    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
%    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
%
%    Input, integer M, the order of the highest moment to
%    be calculated.  This value is only needed when KIND = 8.
%
%    Input, real ALPHA, BETA, the parameters, if required
%    by the value of KIND.
%
  if ( kind <= 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'PARCHK - Fatal error!\n' );
    fprintf ( 1, '  KIND <= 0.\n' );
    error ( 'PARCHK - Fatal error!' );
  end
%
%  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
%
  if ( 3 <= kind && alpha <= -1.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'PARCHK - Fatal error!\n' );
    fprintf ( 1, '  3 <= KIND and ALPHA <= -1.\n' );
    error ( 'PARCHK - Fatal error!' );
  end
%
%  Check BETA for Jacobi.
%
  if ( kind == 4 && beta <= -1.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'PARCHK - Fatal error!\n' );
    fprintf ( 1, '  KIND == 4 and BETA <= -1.0.\n' );
    error ( 'PARCHK - Fatal error!' );
  end
%
%  Check ALPHA and BETA for rational.
%
  if ( kind == 8 )
    tmp = alpha + beta + m + 1.0;
    if ( 0.0 <= tmp || tmp <= beta )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'PARCHK - Fatal error!\n' );
      fprintf ( 1, '  KIND == 8 but condition on ALPHA and BETA fails.\n' );
      error ( 'PARCHK - Fatal error!' );
    end
  end

  return
end
function value = r8_sign ( x )

%*****************************************************************************80
%
%% R8_SIGN returns the sign of an R8.
%
%  Discussion:
%
%    The value is +1 if the number is positive or zero, and it is -1 otherwise.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    21 March 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the number whose sign is desired.
%
%    Output, real VALUE, the sign of X.
%
  if ( 0 <= x )
    value = +1.0;
  else
    value = -1.0;
  end

  return
end
function r8mat_write ( output_filename, m, n, table )

%*****************************************************************************80
%
%% R8MAT_WRITE writes an R8MAT file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 August 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string OUTPUT_FILENAME, the output filename.
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of points.
%
%    Input, real TABLE(M,N), the points.
%

%
%  Open the file.
%
  output_unit = fopen ( output_filename, 'wt' );

  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8MAT_WRITE - Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'R8MAT_WRITE - Error!' );
  end
%
%  Write the data.
%
%  For smaller data files, and less precision, try:
%
%     fprintf ( output_unit, '  %14.6f', table(i,j) );
%
  for j = 1 : n
    for i = 1 : m
      fprintf ( output_unit, '  %24.16f', table(i,j) );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  Close the file.
%
  fclose ( output_unit );

  return
end
function rule_write ( order, filename, x, w, r )

%*****************************************************************************80
%
%% RULE_WRITE writes a quadrature rule to a file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    18 February 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer ORDER, the order of the rule.
%
%    Input, string FILENAME, specifies the output files.
%    write files 'filename_w.txt', 'filename_x.txt', 'filename_r.txt' defining 
%    weights, abscissas, and region.
%
%    Input, real X(ORDER), the abscissas.
%
%    Input, real W(ORDER), the weights.
%
%    Input, real R(2), the region.
%
  filename_x = strcat ( filename, '_x.txt' );
  filename_w = strcat ( filename, '_w.txt' );
  filename_r = strcat ( filename, '_r.txt' );

  fprintf ( 1, '\n' );
  fprintf ( 1,'  Creating quadrature files.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  "Root" file name is   "%s".\n', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Weight file will be   "%s".\n', filename_w );
  fprintf ( 1, '  Abscissa file will be "%s".\n', filename_x );
  fprintf ( 1, '  Region file will be   "%s".\n', filename_r );

  r8mat_write ( filename_w, 1, order, w' );
  r8mat_write ( filename_x, 1, order, x' );
  r8mat_write ( filename_r, 1, 2,     r' );

  return
end
function [ t, wts ] = scqf ( nt, t, mlt, wts, nwts, ndx, kind, alpha, ...
  beta, a, b )

%*****************************************************************************80
%
%% SCQF scales a quadrature formula to a nonstandard interval.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 February 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Parameters:
%
%    Input, integer NT, the number of knots.
%
%    Input, real T(NT), the original knots.
%
%    Input, integer MLT(NT), the multiplicity of the knots.
%
%    Input, real WTS(NWTS), the weights.
%
%    Input, integer NWTS, the number of weights.
%
%    Input, integer NDX(NT), used to index the array WTS.
%    For more details see the comments in CAWIQ.
%
%    Input, integer KIND, the rule.
%    1, Legendre,             (a,b)       1.0
%    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
%    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
%    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
%    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
%    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
%    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
%    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
%    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
%
%    Input, real ALPHA, the value of Alpha, if needed.
%
%    Input, real BETA, the value of Beta, if needed.
%
%    Input, real A, B, the interval endpoints.
%
%    Output, real T(NT), the scaled knots.
%
%    Output, real WTS(NWTS), the scaled weights.
%
  temp = eps;

  parchk ( kind, 1, alpha, beta )

  if ( kind == 1 )

    al = 0.0;
    be = 0.0;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  elseif ( kind == 2 )

    al = -0.5;
    be = -0.5;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  elseif ( kind == 3 )

    al = alpha;
    be = alpha;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  elseif ( kind == 4 )

    al = alpha;
    be = beta;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  elseif ( kind == 5 )

    if ( b <= 0.0 )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  B <= 0.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = a;
    slp = 1.0 / b;
    al = alpha;
    be = 0.0;

  elseif ( kind == 6 )

    if ( b <= 0.0 )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  B <= 0.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = a;
    slp = 1.0 / sqrt ( b );
    al = alpha;
    be = 0.0;

  elseif ( kind == 7 )

    al = alpha;
    be = 0.0;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  elseif ( kind == 8 )

    if ( a + b <= 0.0 )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  A + B <= 0.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = a;
    slp = a + b;
    al = alpha;
    be = beta;

  elseif ( kind == 9 )

    al = 0.5;
    be = 0.5;

    if ( abs ( b - a ) <= temp )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'SCQF - Fatal error!\n' );
      fprintf ( 1, '  |B - A| too small.\n' );
      fprintf ( 1, '  A = %f\n', a );
      fprintf ( 1, '  B = %f\n', b );
      error ( 'SCQF - Fatal error!' );
    end

    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;

  end

  p = slp^( al + be + 1.0 );

  for k = 1 : nt

    t(k) = shft + slp * t(k);
    l = abs ( ndx(k) );

    if ( l ~= 0 )
      tmp = p;
      for i = l : l + mlt(k) - 1
        wts(i) = wts(i) * tmp;
        tmp = tmp * slp;
      end
    end

  end

  return
end
function [ t, wts ] = sgqf ( nt, aj, bj, zemu )

%*****************************************************************************80
%
%% SGQF computes knots and weights of a Gauss Quadrature formula.
%
%  Discussion:
%
%    This routine computes all the knots and weights of a Gauss quadrature
%    formula with simple knots from the Jacobi matrix and the zero-th
%    moment of the weight function, using the Golub-Welsch technique.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    12 February 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Parameters:
%
%    Input, integer NT, the number of knots.
%
%    Input, real AJ(NT), the diagonal of the Jacobi matrix.
%
%    Input, real BJ(NT), the subdiagonal of the Jacobi
%    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
%
%    Input, real ZEMU, the zero-th moment of the weight function.
%
%    Output, real T(NT), the knots.
%
%    Output, real WTS(NT), the weights.
%

%
%  Exit if the zero-th moment is not positive.
%
  if ( zemu <= 0.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'SGQF - Fatal error!\n' );
    fprintf ( 1, '  ZEMU <= 0.\n' );
    error ( 'SGQF - Fatal error!' );
  end
%
%  Set up vectors for IMTQLX.
%
  wts = zeros ( nt, 1 );

  wts(1) = sqrt ( zemu );
  wts(2:nt) = 0.0;
%
%  Diagonalize the Jacobi matrix.
%
  [ t, wts ] = imtqlx ( nt, aj, bj, wts );

  wts(1:nt) = wts(1:nt).^2;

  return
end
function timestamp ( )

%*****************************************************************************80
%
%% TIMESTAMP prints the current YMDHMS date as a timestamp.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 February 2003
%
%  Author:
%
%    John Burkardt
%
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );
  fprintf ( 1, '%s\n', s );

  return
end
