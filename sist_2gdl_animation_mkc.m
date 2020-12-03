function sist_2gdl_animation_mkc(t_final, freq)

%freq = 1;

%t_final = 10;
   
global  a b LT0 ns hs H

%--- dados geometricos
H = .2;
a = .2;
b = .2;
LT0 = .8;
hs = .05;
ns = 6;

%----- dados do sistema
m1 = 10;
m2 = 5;

k1 = 20000;
k2 = 15000;

c1 = 10.0;
c2 = 10.0;

%---- tempo
t1=1
t2=2
F0 = 100;
Force_sin = [ F0
               0];
          
              
%  freq = 2.5;
omega = 2. * pi * freq;

C =[ c1+c2     -c2
      -c2     c2]
      
%C =[  0    0
%      0   c2+c3 ];

K =[ k1+k2   -k2
      -k2    k2]

M = [ m1  0. 
      0.  m2]
%-------       
[V,D] = eig(K,M) ;

for i=1:2
   omega_n(i) =  sqrt(D(i,i));
   f_n(i)     = omega_n(i)  / (2*pi);
   V(:,i)     = V(:,i) / V(1,i);

end
%----
disp('autovetores e autovalores')

V
omega_n
f_n

%------------
%- estados
%------------
       
A = zeros(4,4);

A(1:2,3:4) = [1 0; 
              0 1]; 

A(3:4,1:2) = - K / M;
A(3:4,3:4) = - C / M;


B = zeros(4,1);

B(3:4,1) = inv(M) *  Force_sin ;


% -----------------
% time vector
fs = 512;
t_step = 1/fs;


tspan = 0.:t_step: (t_final-t_step);
n_pontos = length(tspan)

%----------

y0 = [ 0.0; 0.0; 0.0; 0.0]
%----------

opts = odeset('MaxStep',0.01,'RelTol',0.1,'Stats','on');

[t,y] = ode45(@df_2gdl_AB,tspan,y0, opts , A,B,omega);


dxdt = y(:,3);%+3*omega*cos(omega*t);
%  dfghs

%------------
figure(2) ;
clf;
set(gcf,'Units','pixels')


%  set(gcf,'Position',[200 200 500 600])
hold on
%------------
% x
% set(gca,'Position',[0.075 0.6 .85 .35])
hold on
plot(t,y(:,1),'r')
plot(t,y(:,2),'b')
xlabel('t (s)','FontSize',12);
ylabel('x(t) (m/s)','FontSize',12,'interpreter','tex');
grid on

%------------------------------------
%--- plotagem
%------------------------------------
%ff=figure(1, 'position', [ 100 200 1200 600]);
figure(1);
clf

% create empty toolbar
%tf = uitoolbar (ff);
%% create a 19x19x3 black square
%img=zeros(19,19,3);
%% add pushtool button to toolbar
%b = uipushtool (tf, "cdata", img);
%------------

subplot(2,1,1)
hold on
xlim( [0. 2*(LT0+2.*a)] )
ylim( [-(a/2.) (a/2.) ] )  
axis ("equal");

plot( [0 0], [-b/2 b/2] , '--r');
plot( [0 2.75*LT0], [-b/2 -b/2] , '--r');
plot( [LT0+a/2 LT0+a/2], [0 b ] , '-.r');
plot( [2*LT0+a+a/2 2*LT0+a+a/2], [0 b ] , '-.r');

%------------------------------------
%--- first shot  --------------------
%------------------------------------
shot = update_k_m(y(1,1:2));

hs1 = plot(shot.spr1x,shot.spr1y, 'ko-','MarkerSize', 2);
hb1 = plot(shot.bl1x,shot.bl1y, 'x-r');  

hs2 = plot(shot.spr2x,shot.spr2y, 'ko-','MarkerSize', 2);
hb2 = plot(shot.bl2x,shot.bl2y, 'x-b');    

if(c1>0) 
   hd1 = plot(shot.dp1x,shot.dp1y, 's-g');   
end
if(c2>0)
   hd2 = plot(shot.dp2x,shot.dp2y, 's-g');   
end


ht = text(2*LT0+3*a,1.5*b,[ 't=' num2str(t(1)) ] , 'fontsize', 12);
hfr = text(0,1.5*b,['frequencia: ' num2str(freq) ' hz'], 'fontsize', 12 );
text(0,1.25*b,['fn_1: ' num2str(f_n(1)) ' Hz'] , 'fontsize', 12);
text(0,1.00*b,['fn_2: ' num2str(f_n(2)) ' Hz'] , 'fontsize', 12);

text(LT0+a/2,1.5*b,[ 'm1=' num2str(m1) ' kg'], 'fontsize', 12 );
text(2*LT0+a+a/2,1.5*b,[ 'm2=' num2str(m2) ' kg'], 'fontsize', 12 );

text(0.5*LT0    ,1.5*b,[ 'k1=' num2str(k1) ' N/m'] , 'fontsize', 12);
text(1.5*LT0+a/2,1.5*b,[ 'k2=' num2str(k2) ' N/m'] , 'fontsize', 12);
text(0.5*LT0    ,1.25*b,[ 'c1=' num2str(c1) ' Ns/m'] , 'fontsize', 12);
text(1.5*LT0+a/2,1.25*b,[ 'c2=' num2str(c2) ' Ns/m'] , 'fontsize', 12);

axis off

%----
subplot(2,1,2)

hx1 = plot(t,y(:,1),'r')
hold on
hx2 = plot(t,y(:,2),'b')

ylim( [-H H] )
hold on
xlim( [0. t_final] )
ylabel('x(t)')

   
for i=2:4:length(t)
      
   set(ht, 'string',[ 't=' num2str(t(i))] );
   
   fr = omega * t(i) / 2/pi;
   set(hfr, 'string',[ 'frequencia: ' num2str(fr) ' hz' ] );
   
   
   shot = update_k_m(y(i,1:2));  
   
   set(hs1, 'XData', shot.spr1x) % set new y-data in plot handle
   set(hs1, 'YData', shot.spr1y) % set new y-data in plot handle   
   
   set(hb1, 'XData', shot.bl1x) % set new y-data in plot handle
   set(hb1, 'YData', shot.bl1y) % set new y-data in plot handle
   
      
   if(c1>0) 
      set(hd1, 'XData', shot.dp1x) % set new y-data in plot handle
      set(hd1, 'YData', shot.dp1y) % set new y-data in plot handle
   end
   
   set(hs2, 'XData', shot.spr2x) % set new y-data in plot handle
   set(hs2, 'YData', shot.spr2y) % set new y-data in plot handle   
   
   set(hb2, 'XData', shot.bl2x) % set new y-data in plot handle
   set(hb2, 'YData', shot.bl2y) % set new y-data in plot handle
      
   if(c2>0)
      set(hd2, 'XData', shot.dp2x) % set new y-data in plot handle
      set(hd2, 'YData', shot.dp2y) % set new y-data in plot handle
   end

   
   set(hx1, 'XData', t(1:i)) % set new y-data in plot handle
   set(hx1, 'YData', y(1:i,1)) % set new y-data in plot handle
   
   set(hx2, 'XData', t(1:i)) % set new y-data in plot handle
   set(hx2, 'YData', y(1:i,2)) % set new y-data in plot handle
   
     
   %pause(t_step);
   drawnow;
   
   
end


%print('-dpng', 'animation_km.png', '-S640,480');
   
endfunction


%function [springx, springy, blockx, blocky, sp2x, sp]=update_k_m(x)
function shot=update_k_m(x)
      
global  a b LT0 ns hs H

   fsh = x(1,1)/H;
   ys = b/4;
   yd = -b/4;
      
   LT1 = LT0 + x(1,1);
      
   L = LT1 - 2* hs;  
   
   %--- mola 1
   x0 = 0;
   [springx springy] = spring(x0,ys,L,hs,ns,fsh);
   %--- dashpot 1
   x0 = 0;
   [dpx dpy] = dashpot(x0,yd,L,hs);      
   %--- bloco 1
   x0 = LT1;
   [blockx blocky] = block(x0,a,b);
   
   shot=[];
   shot.spr1x = springx;
   shot.spr1y = springy;
   shot.bl1x = blockx;
   shot.bl1y = blocky;
   shot.dp1x = dpx;
   shot.dp1y = dpy;
  
   LT2 = LT0 + x(1,2);
    
   L = LT2 - 2* hs;
   
   ds = L / (2*ns-2);
    
   %--- mola 2
   x0 = blockx(4) ;
   [springx springy] = spring(x0,ys,L,hs,ns,fsh);
   %--- dashpot 2
   x0 =  blockx(4);
   [dpx dpy] = dashpot(x0,yd,L,hs);
      
   %--- bloco 2
   x0 = LT1+LT2+a;
   [blockx blocky] = block(x0,a,b);
   
   
   shot.spr2x = springx;
   shot.spr2y = springy;
   shot.bl2x  = blockx;
   shot.bl2y  = blocky;
   shot.dp2x = dpx;
   shot.dp2y = dpy;
   
   
endfunction

function [springx springy] = spring(x0,ys,L,hs,ns,fsh)
   
   dhs = -hs/20 * fsh;
    
   ds = L / (2*ns-1);
   
   
   springx = zeros(2*ns+2,1);   
   springy = ones(2*ns+2,1)*ys;   
   springx(1) = x0;
   springx(2) = x0+hs    ; springy(2) = ys;
   for i=3:2:2*ns
      springx(i) = springx(i-1,1) + ds;
      springy(i) = ys + (hs/2+dhs) ;      
      
      springx(i+1) = springx(i,1) + ds;
      springy(i+1) = ys -(hs/2+dhs) ;
   end
   
   springx(2*ns+1,1) = x0+hs+L; 
   springx(2*ns+2,1) = x0+hs+L+hs  ; 
   %springx(2*ns+1,1) = x0+L+ds+hs  ;  
   
endfunction
%------------------------------------------------------------------------------|

function [x y] = dashpot(x0,ys,L,hs)
   
   dp = L-2*hs;
   
   x = zeros(8,1);   
   y = ones(8,1)*ys;   
   x(1) = x0;
   x(2) = x0+2*hs    ;  
   x(3) = x0+2*hs    ; y(3) = ys-hs/2;
   x(4) = x0+2*hs+dp  ; y(4) = ys-hs/2;
   x(5) = x0+2*hs+dp  ; y(5) = ys+hs/2;   
   x(6) = x0+2*hs    ; y(6) = ys+hs/2; 
   x(7) = x0+2*hs    ; 
   x(8) = x0+L+2*hs; 
   
   
endfunction
%------------------------------------------------------------------------------|
function [blockx blocky] = block(x0,a,b)
   
   blockx = zeros(6,1);  
   blocky = zeros(6,1);  
   blockx(1) = x0;
   blockx(2) = x0;   blocky(2) = -b/2;
   blockx(3) = x0+a; blocky(3) = -b/2;
   blockx(4) = x0+a; blocky(4) =  b/2;
   blockx(5) = x0  ; blocky(5) =  b/2;
   blockx(6) = x0  ;  
   
endfunction
%------------------------------------------------------------------------------|
% function that solves a system with a general movement of the base passing 
% by a bump
%------------------------------------------------------------------------------|
function dydx = df_2gdl_AB(t,y,A,B,omega)

force = cos(omega * t *  t);

dydx = zeros(4,1);      
      
dydx =   A * y +   B *  force;

endfunction
