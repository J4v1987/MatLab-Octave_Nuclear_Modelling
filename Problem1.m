%Student Laboratories
%Problem 1: Application of the Participant-Spectator-Model

%Set project directory.
if (!strcmp("/home/javier/Documents/Octave/Nuclear Modelling/Problem1",pwd()))
  cd "/home/javier/Documents/Octave/Nuclear Modelling/Problem1"
endif

%Define input particle analysis constants
R0=1.2;                  % SI: [fm].
Ap=12;                 % Projectile particle mass number in absolute value, e.g. Ap(C)=12[u]→12.
At=197;                    % Target particle mass number in absolute value, e.g. At(Au)=197[u]→197.
Ep=85;                    % Impact parameter engagement percentage, 100: just disengaged, 0: center to center
plotIt=true;            % Set to true if you want to see the analysis 3D plot (may take longer to process).

%Define algorithm operation constants
N=10000;                  % Number of Monte Carlo random points to be tried in space.
delineants=15;         % Number of points to delineate the 2D shadows of each particle.
del_thkns=8;             % Thickness of delinenant points.
k=0;

%Define input variables
%It is assumed both particles are above XY plane along Z(+) in euclidean 3D space (Eⁿ)
Rp=R0*power(Ap,1/3);                     % Projectile particle radius
Rt=R0*power(At,1/3);                       % Target particle radius
b=(Rp+Rt)*Ep*power(100,-1) ;        % Impact parameter scaled to ensure engagement of particles
z_off_t=Rt;                                           % Target particle vertical 'z' offset, such that center coordinates are x=0, y=0, z=radius(Rt)
z_off_p=Rt+b;                                     % Projectile particle vertical 'z' offset, such that center coordinates are x=0, y=0, z=target diameter+projectile radius
P_top=b+Rt+Rp;                                 % Projectile top point (along z axis in Eⁿ)
P_bot=P_top-(2*Rp);                         % Projectile bottom point (along z axis in Eⁿ)
T_top=2*Rt;                                        % Target top point (along z axis in Eⁿ)
T_bot=0;                                              % Target bottom point (along z axis in Eⁿ)
%Test plots, depicting the external 3D boundaries of projectile and test particles.
  %Test target particle shadow on YZ & XZ in Eⁿ.
  y_t = linspace(-Rt,Rt,delineants);
  %Create top half-circle vector Eⁿ : YZ
  z_pos_t = (sqrt(Rt.^2 - y_t.^2)+z_off_t);
  %Create bottom half-circle vector Eⁿ : YZ
  z_neg_t = -(sqrt(Rt.^2 - y_t.^2)-z_off_t);
  %Create top half-circle vector Eⁿ : XY
  x_pos_t = sqrt(Rt.^2 - y_t.^2);
  %Create bottom half-circle vector Eⁿ : XY
  x_neg_t = -sqrt(Rt.^2 - y_t.^2);

 %Delineate 2D shadows of both particles onto YZ, XZ, XY
 if (plotIt==true)
    for i=1:length(y_t);
      %Plot upper half of target particle, YZ.
      plot3(0,y_t(i),z_pos_t(i),'r.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
      %Plot lower half of target particle, YZ.
      plot3(0,y_t(i),z_neg_t(i),'r.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
      %Plot upper half of target particle, XY.
      plot3(x_pos_t(i),y_t(i),z_off_t,'r.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
      %Plot lower half of target particle, XY.
      plot3(x_neg_t(i),y_t(i),z_off_t,'r.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
     endfor
    %Test projectile particle shadow on YZ in Eⁿ.
    y_p = linspace(-Rp,Rp,delineants);
    %Create top half-circle vector, Eⁿ : YZ
    z_pos_p = (sqrt(Rp.^2 - y_p.^2)+z_off_p);
    %Create bottom half-circle vector, Eⁿ : YZ
    z_neg_p = -(sqrt(Rp.^2 - y_p.^2)-z_off_p);
    %Create top half-circle vector, Eⁿ : XY.
    x_pos_p = sqrt(Rp.^2 - y_p.^2);
    %Create bottom half-circle vector, Eⁿ : XY.
    x_neg_p = -sqrt(Rp.^2 - y_p.^2);
    for i=1:length(y_t);
      %Plot upper half of projectile particle, YZ.
      plot3(0,y_p(i),z_pos_p(i),'b.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
      %Plot lower half of projectile particle, YZ.
      plot3(0,y_p(i),z_neg_p(i),'b.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
      %Plot upper half of projectile particle, XY.
      plot3(x_pos_p(i),y_p(i),z_off_p,'b.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
      %Plot lower half of projectile particle, XY.
      plot3(x_neg_p(i),y_p(i),z_off_p,'b.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
     endfor

   %Test target particle shadow on XZ in Eⁿ.
    x_t=y_t;
    for i=1:length(y_t);
      %Plot upper half of target particle.
      plot3(x_t(i),0,z_pos_t(i),'r.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
      %Plot lower half of target particle.
      plot3(x_t(i),0,z_neg_t(i),'r.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
     endfor


   %Test projectile particle shadow circle on XZ in Eⁿ.
    x_p=y_p;
    for i=1:length(y_p);
      %Plot upper half of projectile particle.
      plot3(x_p(i),0,z_pos_p(i),'b.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
      %Plot lower half of projectile particle.
      plot3(x_p(i),0,z_neg_p(i),'b.', 'MarkerSize', del_thkns, 'LineWidth', 2);
      hold on;
     endfor

  %Delineate top and bottom trajectories of each particle on XZ
  if (Rt>=Rp)
  for i=1:length(x_t);
    %Plot tangent top trajectory line of projectile particle.
    plot3(x_t(i),0,P_top,'g.', 'MarkerSize', del_thkns, 'LineWidth', 2);
    hold on;
    %Plot tangent bottom trajectory line of projectile particle.
    plot3(x_t(i),0,P_bot,'g.', 'MarkerSize', del_thkns, 'LineWidth', 2);
    hold on;
    %Plot tangent top trajectory line of target particle.
    plot3(x_t(i),0,T_top,'m.', 'MarkerSize', del_thkns, 'LineWidth', 2);
    hold on;
    %Plot tangent bottom trajectory line of target particle.
    plot3(x_t(i),0,T_bot,'m.', 'MarkerSize', del_thkns, 'LineWidth', 2);
    hold on;
  endfor
elseif (Rt<Rp)
  for i=1:length(x_p);
    %Plot tangent top trajectory line of projectile particle.
    plot3(x_p(i),0,P_top,'g.', 'MarkerSize', del_thkns, 'LineWidth', 2);
    hold on;
    %Plot tangent bottom trajectory line of projectile particle.
    plot3(x_p(i),0,P_bot,'g.', 'MarkerSize', del_thkns, 'LineWidth', 2);
    hold on;
    %Plot tangent top trajectory line of target particle.
    plot3(x_p(i),0,T_top,'m.', 'MarkerSize', del_thkns, 'LineWidth', 2);
    hold on;
    %Plot tangent top trajectory line of target particle.
    plot3(x_p(i),0,T_bot,'m.', 'MarkerSize', del_thkns, 'LineWidth', 2);
    hold on;
  endfor
endif


endif
% Create random points in a prismatic volume tightly enclosing both particles: target and projectile.
for i=1:N
  % Create a single random point on every i iteration, initially normalized from 0-1 and later to fit in the scope of the particles hosting prism.
  if ((P_bot<T_top)&&(Rp<=Rt))
    xp=(2*Rt)*(rand()-0.5);
    yp=(2*Rt)*(rand()-0.5);
    if (P_top<=T_top)
      zp=(T_top-T_bot)*rand();
    elseif (P_top>T_top)
      zp=(P_top-T_bot)*rand();
    endif
  elseif (Rp>Rt)
    xp=(2*Rp)*(rand()-0.5);
    yp=(2*Rp)*(rand()-0.5);
    if (P_bot>T_bot)
      zp=(P_top-T_bot)*rand();
    elseif (P_bot<=T_bot)
      zp=((Rp*(2*rand()-1)))+(Rt+b);
    endif
  endif

  if ((zp>P_bot )&&(zp<T_top)&&((yp.^2+(zp-z_off_p).^2)<Rp.^2)&&((yp.^2+(zp-z_off_t).^2)<Rt.^2)&&((xp.^2+(zp-z_off_t).^2)<Rt.^2)&&(T_bot<P_bot))
    k++;
    if (plotIt==true)
      plot3(xp, yp, zp, 'k.', 'MarkerSize', del_thkns/2, 'LineWidth', 2);
      hold on;
    endif
  elseif ((Rp>Rt)&&(T_bot<P_bot))
    if ((zp>P_bot )&&(zp<T_top)&&((yp.^2+(zp-z_off_p).^2)<Rp.^2)&&((yp.^2+(zp-z_off_t).^2)<Rt.^2)&&((xp.^2+(zp-z_off_t).^2)<Rt.^2))
      k++;
      if (plotIt==true)
        plot3(xp, yp, zp, 'k.', 'MarkerSize', del_thkns/2, 'LineWidth', 2);
        hold on;
      endif
    endif
  elseif ((Rp>Rt)&&(T_bot>=P_bot)) %The projectile is larger than target, and runs through the entirety of the target.
    %if the point (xp,yp,zp) is within the projectile shadows on YZ, XZ, and XY, then plot and count it.
    if (((yp.^2+(zp-z_off_t).^2)<Rt.^2)&&((xp.^2+(zp-z_off_t).^2)<Rt.^2)&&((xp.^2+yp.^2)<Rt.^2))
      k++;
      if (plotIt==true)
        plot3(xp, yp, zp, 'k.', 'MarkerSize', del_thkns/2, 'LineWidth', 2);
        hold on;
      endif
    endif
  endif;
  i++;
  plot3(xp,yp,zp,'c.', 'MarkerSize', del_thkns/2, 'LineWidth', 2) %if needed be, fill the prism with cyan random points to confirm MC volume zone distribution
endfor;

if (plotIt==true)
  axis equal;
endif
