%% Implemented functions


%%
function []=icosahedralSpectrumOfIcosahedralShell()
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

    H=[];
    [IH_Basis,ArrayAngleFlag]=Make_IHBasis(Theta,Phi,l_Vector,H);
    
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
    ISpectrumMtx=nan(N_l,N_R);
    for R_cntr=1:N_R
        %fr=squeeze(IcoFunc(R_cntr,:,:));
        fr=IcoFunc{R_cntr};
        ISpectrum_=Make_ISpectrum(fr);
        if R_Norm_Flag
            ISpectrum_=ISpectrum_/(Radius2-Radius1);
        end
        ISpectrumMtx(:,R_cntr)=ISpectrum_;
 
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
    
    H5=figure;
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

    
    function I=Integral(f)
        switch Integral_Mode
            case 1
                I=mean(mean(f,1),2);
            case 2
                I=sum(sum(f.*W_Theta.*W_Phi,1),2);
        end
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

function Edges=Find_Edges(V)
    Edges=zeros(12,5);
    for cntr1=1:12
        Distance=(V.x-V.x(cntr1)).^2+(V.y-V.y(cntr1)).^2+(V.z-V.z(cntr1)).^2;
        [~,I]=sort(Distance,'ascend');
        Edges(cntr1,:)=I(1+(1:5));
    end
end
