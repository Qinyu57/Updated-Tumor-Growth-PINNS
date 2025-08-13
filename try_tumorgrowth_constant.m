% Test program 4
%  Spatial Dependence in growth rate

initialFun=@(X,Y)(0.34*((sqrt(X.^2+Y.^2)-0.5)<0));
gf=@(p,X,Y)(p*0+2.2+(-0)*sin(sqrt(X.^2 + Y.^2)));

%%
[X,Y,rho,p]=solver_2D('T',1,'GrowthFun',gf,'dt',0.005/2,'InitialFun',initialFun,'m',4,'IfVideo',true,'videoname','demo4.avi');
save('demo4.mat')
%[X,Y,rho,p]=solver_2D('m',40,'T',3,'IfVideo',true,'videoname','bvm=40p.avi')
%%
load('demo4.mat')
figure(2)
subplot(3,2,5)
surf(X,Y,rho)
view(90,0)
title('t=1')
subplot(3,2,6)
surf(X,Y,rho)
view(2)
title('t=1')
subplot(3,2,1)
surf(X,Y,1*((X.^2+Y.^2)<0.5))
view(90,0)
title('t=0')
subplot(3,2,2)
surf(X,Y,1*((X.^2+Y.^2)<0.5))
view(2)
title('t=0')
%%
tic
[X,Y,rho,p]=solver_2D('T',0.875,'GrowthFun',gf,'dt',0.005/2);
toc

subplot(3,2,3)
surf(X,Y,rho)
view(90,0)
title('t=0.875')
subplot(3,2,4)
surf(X,Y,rho)
view(2)
title('t=0.875')

set(gcf,'unit','centimeters','position',[10 5 14 42]);

print(['demo4_1','.eps'],'-depsc');
disp(rho)



%[X,Y,rho,p]=solver_pro('T',6,'Evx',Velocity_x,'IfVideo',true,'videoname','demo3.avi');
