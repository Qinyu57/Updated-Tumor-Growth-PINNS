% The main solver. Please see the document attached for usage and examples.
% Xu'an Dou dxa@pku.edu.cn 2020.9.6
function [X,Y,rho,p,cell_rho]=solver_2D(varargin)
%% Parsing the Input
   pa = inputParser;
   pa.CaseSensitive=false;   validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
   
   defaultTime = 1;validScalarNonNeg = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
   addParameter(pa,'T',defaultTime,validScalarNonNeg);
   
   defaultM=3;
   addParameter(pa,'m',defaultM,validScalarPosNum);

   defaultDx = 0.1;
   addParameter(pa,'dx',defaultDx,validScalarPosNum);
   
   defaultDt=0.005;
   addParameter(pa,'dt',defaultDt,validScalarPosNum);

   defaultInitialFun=@(X,Y)(0.9*((X.^2+Y.^2)<0.3));
   defaultGfun=@(p,X,Y)(0*p+1);
   addParameter(pa,'InitialFun',defaultInitialFun,@(x) isa(x,'function_handle'));
   addParameter(pa,'GrowthFun',defaultGfun,@(x) isa(x,'function_handle'));
   
   defaultOutVelocity_x=@(X,Y)(X*0);
   defaultOutVelocity_y=@(X,Y)(Y*0);
   addParameter(pa,'Evx',defaultOutVelocity_x,@(x) isa(x,'function_handle'));
   addParameter(pa,'Evy',defaultOutVelocity_y,@(x) isa(x,'function_handle'));
   
   defaultIfVideo=false;
   defaultVideoName='defaultName.avi';
   addParameter(pa,'IfVideo',defaultIfVideo,@(x) isa(x,'logical'));
   addParameter(pa,'VideoName',defaultVideoName,@(x) isa(x,'char'));
   
   default_a=-3;
   default_b=3;
   addParameter(pa,'a',default_a,@(x)isscalar(x));
   addParameter(pa,'b',default_b,@(x)isscalar(x));

   default_shottime=[];% should be a increasing positive sequence
   addParameter(pa,'shottime',default_shottime);% 


   parse(pa, varargin{:});

%%
m=pa.Results.m;
T=pa.Results.T;
dx=pa.Results.dx; dy=dx; dt=pa.Results.dt; 
funG=pa.Results.GrowthFun; VideoName=pa.Results.VideoName; IfVideo=pa.Results.IfVideo;

a=pa.Results.a;b=pa.Results.b; 
Ntime=floor(T/dt);
Nx=length(a:dx:b)-1;Ny=length(a:dy:b)-1;

[X,Y]=meshgrid(a:dx:b,a:dy:b);
u_e=pa.Results.Evx(X,Y);u_e=u_e';u_e=u_e(2:end-1,2:end-1);
v_e=pa.Results.Evy(X,Y);v_e=v_e';v_e=v_e(2:end-1,2:end-1);
%% shot time
shot_time=pa.Results.shottime;
num_shot=length(shot_time); cell_rho=cell(1,num_shot);
cnt_shot=0;

%% initial value
rho=pa.Results.InitialFun(X,Y);%0.9*((sqrt(X.^2+Y.^2)-0.5-sin(4*atan(Y./X))/2)<0);
%surf(X,Y,rho)
rho=rho';

p=rho.^(m-1)*m/(m-1);
nx=Nx-1;ny=Ny-1; %u=zeros(nx,ny);
u=-(p(3:end,2:end-1)-p(1:end-2,2:end-1))/(2*dx); % Why we need u and v here?
v=-(p(2:end-1,3:end)-p(2:end-1,1:end-2))/(2*dy);

if IfVideo
    cnt_frame=0;
    viedo = VideoWriter(VideoName);
    viedo.FrameRate=15;
    open(viedo)
end
%% Declare some global variables used in solving the diffusion matrix(maybe useless)
    U_m=zeros(nx,ny);
    V_m=U_m;uv_pp=U_m;uv_pm=U_m;uv_mm=U_m;uv_mp=U_m;vu_pp=U_m;vu_pm=U_m;vu_mm=U_m;vu_mp=U_m;
    V_up=U_m;V_low=U_m;   U_up=U_m;U_low=U_m; b_sol=zeros(nx*ny*2,1);

%% one step time evolution
for i=1:Ntime
    if cnt_shot<num_shot&&abs((i-1)*dt-shot_time(cnt_shot+1))<dt
        cell_rho{cnt_shot+1}=rho;
        cnt_shot=cnt_shot+1;
    end
    if IfVideo&&mod(i,10)==1%Some code for recording a movie
        cnt_frame=cnt_frame+1;
        plot_frame;
        title(sprintf('t=%g',i*dt))
        M(cnt_frame)=getframe(gcf);
    end
    % May add some modifaction here, recording more data, nutrients, update the
    % external velocity field......
    %% Solve the equation from n dt to (n+1)dt   
    % Most details are written in some nested functions.
    % Updating Growth Rate
    G=funG(p',X,Y);G=G';
    % Solve the diffusion 
        Build_Matrix;%Build A and b for Ax=b;
        [sol,flag]=gmres(@fA,b_sol, 20);% using matrix-free gmres, with 10 restarts
    % handling the failure of GMRES. maybe move to a nested function? maybe
    if flag~=0
        if flag==1
            warning('gmres fails to converge within maximum iterations, try more......')
            if IfVideo
                cnt_frame=cnt_frame+1;
                plot_frame;             
                title(sprintf('t=%g',i*dt))
                M(cnt_frame)=getframe(gcf);
            end
            [sol,flag]=gmres(@fA,b_sol,400);
            if flag==0
                disp('success')
            else
                if IfVideo
                    writeVideo(viedo,M);
                    close(viedo);
                end
                error('GMRES fails')
            end
        else
            if IfVideo
                writeVideo(viedo,M);
                close(viedo);
            end
        error('GMRES fails')% maybe use a flag in output instead of simply 'error'?
        end
    end
    %From big vectors to matrix
    u_star=sol(1:nx*ny);v_star=sol(nx*ny+1:end);u_star=reshape(u_star,nx,ny);v_star=reshape(v_star,nx,ny);
    % Convection and Correction
    Convec_Cor;
end
% Transforming the coordinate to be consistent with meshgrid
    rho=rho';
    p=p';
if IfVideo
    writeVideo(viedo,M);
    close(viedo);
end
    % Build A and b for Ax=b. 
    % Actually build auxiliary matrixs for evaluating Ax with given x in
    % function handle fA
    function Build_Matrix
    %rho_{i+1/2,j} etc i,j from 1 to n Value on stagger grid
    r_ip=((rho(3:end,2:end-1)+rho(2:end-1,2:end-1))/2).^(m-2);
    r_im=((rho(1:end-2,2:end-1)+rho(2:end-1,2:end-1))/2).^(m-2);
    
    r_jp=((rho(2:end-1,3:end)+rho(2:end-1,2:end-1))/2).^(m-2);
    r_jm=((rho(2:end-1,1:end-2)+rho(2:end-1,2:end-1))/2).^(m-2);

    %Auu
    % main diag  nx*ny
    U_m=1+m*(dt/dx/dx)*rho(2:end-1,2:end-1).*(r_ip+r_im);
    % uppper diag (nx-1)*ny
    U_up=-m*(dt/dx/dx)*r_ip.*rho(3:end,2:end-1);
    U_up=U_up(1:end-1,:);    
    % lower diag (nx-1)*ny
    U_low=-m*(dt/dx/dx)*r_im.*rho(1:end-2,2:end-1);
    U_low=U_low(2:end,:);    
    %Avv
    % main diag  nx*ny
    V_m=1+m*(dt/dy/dy)*rho(2:end-1,2:end-1).*(r_jp+r_jm);
    % uppper diag (nx)*(ny-1)
    V_up=-m*(dt/dy/dy)*r_jp.*rho(2:end-1,3:end); 
    V_up=V_up(:,1:end-1);    
    % lower diag (nx)*(ny-1)
    V_low=-m*(dt/dy/dy)*r_jm.*rho(2:end-1,1:end-2); 
    V_low=V_low(:,2:end);
    
    % Auv & Avv
    lam=m*dt/dx/dy/4;
    %Auv
    uv_pp=-lam*(rho(3:end,2:end-1).^(m-2)).*rho(3:end,3:end);%i+1,j+1
    uv_pp=uv_pp(1:end-1,1:end-1);
    
    uv_pm=lam*(rho(3:end,2:end-1).^(m-2)).*rho(3:end,1:end-2);%i+1,j-1
    uv_pm=uv_pm(1:end-1,2:end);
   
    uv_mp=lam*(rho(1:end-2,2:end-1).^(m-2)).*rho(1:end-2,3:end);%i-1,j+1
    uv_mm=-lam*(rho(1:end-2,2:end-1).^(m-2)).*rho(1:end-2,1:end-2);%i-1,j-1    
    uv_mp=uv_mp(2:end,1:end-1); uv_mm=uv_mm(2:end,2:end);

    %Avu
    vu_pp=-lam*(rho(2:end-1,3:end).^(m-2)).*rho(3:end,3:end);%i+1,j+1
    vu_pm=lam*(rho(2:end-1,3:end).^(m-2)).*rho(3:end,1:end-2);%i+1,j-1
    vu_pp=vu_pp(1:end-1,1:end-1); vu_pm=vu_pm(1:end-1,2:end);
 
    vu_mp=lam*(rho(2:end-1,1:end-2).^(m-2)).*rho(1:end-2,3:end);%i-1,j+1
    vu_mm=-lam*(rho(2:end-1,1:end-2).^(m-2)).*rho(1:end-2,1:end-2);%i-1,j-1    
    vu_mp=vu_mp(2:end,1:end-1); vu_mm=vu_mm(2:end,2:end);
    
    %% Build_Vec b
        rhoG=rho.^(m-1).*G;
        bu=u-m*(dt/dx/2)*(rhoG(3:end,2:end-1)-rhoG(1:end-2,2:end-1));
        bv=v-m*(dt/dy/2)*(rhoG(2:end-1,3:end)-rhoG(2:end-1,1:end-2));
    % For External Velocity
        bu=bu-(U_m-1).*u_e;
        bu(1:end-1,:)=bu(1:end-1,:)-U_up.*u_e(2:end,:);
        bu(2:end,:)=bu(2:end,:)-U_low.*u_e(1:end-1,:);   
    
        bu(1:end-1,1:end-1)=bu(1:end-1,1:end-1)-uv_pp.*v_e(2:end,2:end);
        bu(1:end-1,2:end)= bu(1:end-1,2:end)-uv_pm.*v_e(2:end,1:end-1);
        bu(2:end,1:end-1)=bu(2:end,1:end-1)-uv_mp.*v_e(1:end-1,2:end);
        bu(2:end,2:end)= bu(2:end,2:end)-uv_mm.*v_e(1:end-1,1:end-1);



        bv=bv-(V_m-1).*v_e;
        bv(:,1:end-1)=bv(:,1:end-1)-V_up.*v_e(:,2:end);
        bv(:,2:end)=bv(:,2:end)-V_low.*v_e(:,1:end-1);
    
    bv(1:end-1,1:end-1)=bv(1:end-1,1:end-1)-vu_pp.*u_e(2:end,2:end);
    bv(1:end-1,2:end)= bv(1:end-1,2:end)-vu_pm.*u_e(2:end,1:end-1);
    bv(2:end,1:end-1)=bv(2:end,1:end-1)-vu_mp.*u_e(1:end-1,2:end);
    bv(2:end,2:end)=bv(2:end,2:end)-vu_mm.*u_e(1:end-1,1:end-1);

    
        bu=reshape(bu,nx*ny,1); bv=reshape(bv,nx*ny,1);
        b_sol=[bu;bv];  
    end
    %evaluate Ax with given x for gmres
    function ax=fA(x_tmp)
    %turn big vectors to matrixs
    u_tmp=reshape(x_tmp(1:nx*ny),nx,ny);
    v_tmp=reshape(x_tmp(nx*ny+1:end),nx,ny);
    ax_u=u_tmp*0;ax_v=ax_u;
    
    ax_u=U_m.*u_tmp;
    ax_u(1:end-1,:)=ax_u(1:end-1,:)+U_up.*u_tmp(2:end,:);
    ax_u(2:end,:)=ax_u(2:end,:)+U_low.*u_tmp(1:end-1,:);   
    
    ax_u(1:end-1,1:end-1)=ax_u(1:end-1,1:end-1)+uv_pp.*v_tmp(2:end,2:end);
    ax_u(1:end-1,2:end)= ax_u(1:end-1,2:end)+uv_pm.*v_tmp(2:end,1:end-1);
    ax_u(2:end,1:end-1)=ax_u(2:end,1:end-1)+uv_mp.*v_tmp(1:end-1,2:end);
    ax_u(2:end,2:end)= ax_u(2:end,2:end)+uv_mm.*v_tmp(1:end-1,1:end-1);



    ax_v=V_m.*v_tmp;
    ax_v(:,1:end-1)=ax_v(:,1:end-1)+V_up.*v_tmp(:,2:end);
    ax_v(:,2:end)=ax_v(:,2:end)+V_low.*v_tmp(:,1:end-1);
    
    ax_v(1:end-1,1:end-1)=ax_v(1:end-1,1:end-1)+vu_pp.*u_tmp(2:end,2:end);
    ax_v(1:end-1,2:end)= ax_v(1:end-1,2:end)+vu_pm.*u_tmp(2:end,1:end-1);
    ax_v(2:end,1:end-1)=ax_v(2:end,1:end-1)+vu_mp.*u_tmp(1:end-1,2:end);
    ax_v(2:end,2:end)=ax_v(2:end,2:end)+vu_mm.*u_tmp(1:end-1,1:end-1);
    
    %turn back
    ax_u=reshape(ax_u,nx*ny,1);
    ax_v=reshape(ax_v,nx*ny,1);
    
    ax=[ax_u;ax_v];

    end
    function Convec_Cor% Convection and Correction
   %% Get the value of u in staggered grid
     u_sg=(u_star(2:end,:)+u_star(1:end-1,:)+u_e(2:end,:)+u_e(1:end-1,:))/2;
    v_sg=(v_star(:,2:end)+v_star(:,1:end-1)+v_e(:,2:end)+v_e(:,1:end-1))/2;
    %% Convection
    diff_Fx=getdFx(rho,u_sg);%call functions to compute the differences in flux
    diff_Fy=getdFy(rho,v_sg);
    F_total=dt/dx*diff_Fx+dt/dy*diff_Fy;    
    rho(2:end-1,2:end-1)=(rho(2:end-1,2:end-1)-F_total)./(1-dt*G(2:end-1,2:end-1));% implicit in rho, explicit in G
    %% Correction
    p=rho.^(m-1)*m/(m-1);
    u=-(p(3:end,2:end-1)-p(1:end-2,2:end-1))/(2*dx);
    v=-(p(2:end-1,3:end)-p(2:end-1,1:end-2))/(2*dy);

    end
    function plot_frame
%         subplot(2,1,1)
%         surf(X,Y,p')
%         subplot(2,1,2)
        s=surf(X,Y,rho');
        s.EdgeColor = 'none';
        view(2)

    end
 
%Compute the difference of flux in x direction
function rst=getdFx(rho,u_star) % difference in x, fixed y
rho=rho(:,2:end-1);% 0,1,...Nx
diff_l=(rho(2:end-1,:)-rho(1:end-2,:));%rho_{i}-rho_{i-1} i=1,...Nx-1
diff_r=(rho(3:end,:)-rho(2:end-1,:));%rho_{i+1}-rho_{i} i=1,...Nx-1
   % We omit the stepsize since it is eliminated in the denominator and numerator.
   parial_rho=min(diff_l,diff_r).*(diff_l>0).*(diff_r>0)+max(diff_l,diff_r).*(diff_l<0).*(diff_r<0);
   %rho^l_{i+1/2} i=1,...Nx-2
   rho_l=rho(2:end-2,:)+1/2*parial_rho(1:end-1,:);
   %rho^r_{i+1/2} i=1,...Nx-2
   rho_r=rho(3:end-1,:)-1/2*parial_rho(2:end,:);
   %F_{i+1/2} i=1,...Nx-2
   F=1/2*(rho_l.*u_star+rho_r.*u_star-abs(u_star).*(rho_r-rho_l));
   %Add zero rows
   F(end+1,:)=0;F(end+1,:)=0; F(2:end,:)=F(1:end-1,:); F(1,:)=0;
   rst=F(2:end,:)-F(1:end-1,:); 
end
%Compute the difference of flux in y direction
function rst=getdFy(rho,v_star) % difference in y, fixed x
rho=rho(2:end-1,:);% 0,1,...Ny
diff_l=(rho(:,2:end-1)-rho(:,1:end-2));%rho_{i}-rho_{i-1} i=1,...Ny-1
diff_r=(rho(:,3:end)-rho(:,2:end-1));%rho_{i+1}-rho_{i} i=1,...Ny-1
   % We omit the stepsize since it is eliminated in the denominator and numerator.
   parial_rho=min(diff_l,diff_r).*(diff_l>0).*(diff_r>0)+max(diff_l,diff_r).*(diff_l<0).*(diff_r<0);
   %rho^l_{i+1/2} i=1,...Ny-2
   rho_l=rho(:,2:end-2)+1/2*parial_rho(:,1:end-1);
   %rho^r_{i+1/2} i=1,...Ny-2
   rho_r=rho(:,3:end-1)-1/2*parial_rho(:,2:end);
   %F_{i+1/2} i=1,...Ny-2
   F=1/2*(rho_l.*v_star+rho_r.*v_star-abs(v_star).*(rho_r-rho_l));
   %Add zero columns
   F(:,end+1)=0;F(:,end+1)=0; F(:,2:end)=F(:,1:end-1); F(:,1)=0;
   rst=F(:,2:end)-F(:,1:end-1); 
end
end
