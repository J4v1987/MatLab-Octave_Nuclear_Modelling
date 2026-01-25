%Student Laboratories
%Problem 1: Application of the Participant-Spectator-Model

%Ensure Octave starts fresh
clear;
close all;

%Set project directory.
if (!strcmp("/home/javier/Documents/Octave/Nuclear Modelling/Problem1",pwd()))
  cd "/home/javier/Documents/Octave/Nuclear Modelling/Problem1"
endif

%[USER INPUTS]---------------------------------------------------------------------------
%Define input particle analysis constants
R0=1.2;  % in Femtometres [fm], in consistency with "Introductory Nuclear Physics", by Krane, Kenneth S. ISBN: 0-471-80553-X, page 49.                                                                                                                                                              % SI: [fm].
Ap=input("Projectile atomic mass [u] or [g⋅mol⁻¹]: ");                                                                               % Projectile particle mass number in absolute value, e.g. Ap(C)=12[u]→12.
At=input("Target atomic mass [u] or [g⋅mol⁻¹]: ");                                                                                     % Target particle mass number in absolute value, e.g. At(Au)=197[u]→197.
Ep=input("Projectile-target impact parameter scaling factor (100: just disengaged, 0: center to center): ");        % Impact parameter engagement percentage, 100: just disengaged, 0: center to center


%Define algorithm operation constants
N=input("Number Monte Carlo trials (e.g. 500): ");    % Number of Monte Carlo random points to be tried in space.
delineants=15;                                                                   % Number of points to delineate the 2D shadows of each particle.
del_thkns=8;                                                                       % Thickness of delinenant points.
j=0;                                                                                       % Counts the number of points found within the prism tightly containing both projectile and target particles.
k=0;                                                                                      % Counts the number of points found in the intersection through-all of projectile and target.

%Ask user if they wish to view shadow particles plot
plotIt=yes_or_no("Do you want to see the analysis 3D plot? (may take longer to process)");                                                               % Set to true if you want to see the analysis 3D plot (may take longer to process).

%Debug variables
if (plotIt==1)
  plotPrismPoints=1;                                                                    % If set to true, the algorithm will skip the if conditioned blocks tied to this boolean (e.g. points inside MC prism).
else
  plotPrismPoints=0;
endif

%[END OF USER INPUTS]---------------------------------------------------------------------------

disp(char(10)) ;
disp(["PROCESSING..."]) ;
disp(char(10)) ;

%Define input variables
%It is assumed both particles are above XY plane along Z(+) in euclidean 3D space (Eⁿ). Projectile particle center is always at or above target particle center.
Rp=R0*power(Ap,1/3);                     % Projectile particle radius, in consistency with "Introductory Nuclear Physics", by Krane, Kenneth S. ISBN: 0-471-80553-X, page 48.
Rt=R0*power(At,1/3);                       % Target particle radius, in consistency with "Introductory Nuclear Physics", by Krane, Kenneth S. ISBN: 0-471-80553-X, page 48.
b=(Rp+Rt)*Ep*power(100,-1) ;        % Impact parameter scaled to ensure engagement of particles
z_off_t=Rt;                                           % Target particle vertical 'z' offset, such that center coordinates are x=0, y=0, z=radius(Rt)
z_off_p=Rt+b;                                     % Projectile particle vertical 'z' offset, such that center coordinates are x=0, y=0, z=target diameter+projectile radius
P_top=b+Rt+Rp;                                 % Projectile top point (along z axis in Eⁿ)
P_bot=P_top-(2*Rp);                         % Projectile bottom point (along z axis in Eⁿ)
T_top=2*Rt;                                        % Target top point (along z axis in Eⁿ)
T_bot=0;                                              % Target bottom point (along z axis in Eⁿ)

%Define 2D-arrays describing the shadows of both particles onto YZ, XZ, XY in Eⁿ
%Target particle shadow, input arrays.
y_t = linspace(-Rt,Rt,delineants);
x_t=y_t;
tp_shadow_list=[];
for i=1:length(y_t); % Define the target 2D shadow arrays onto YZ, XZ, XY
  tp_shadow_list(end+1, :) = [0, y_t(i), sqrt(Rt.^2 - y_t(i).^2)+z_off_t]; %Create top half-circle vector Eⁿ : YZ
  tp_shadow_list(end+1, :) = [0, y_t(i), -(sqrt(Rt.^2 - y_t(i).^2)-z_off_t)]; %Create bottom half-circle vector Eⁿ : YZ
  tp_shadow_list(end+1, :) = [sqrt(Rt.^2 - y_t(i).^2),y_t(i),z_off_t]; %Create top half-circle vector Eⁿ : XY
  tp_shadow_list(end+1, :) = [-sqrt(Rt.^2 - y_t(i).^2),y_t(i),z_off_t]; %Create bottom half-circle vector Eⁿ : XY
  tp_shadow_list(end+1, :) = [x_t(i),0,sqrt(Rt.^2 - y_t(i).^2)+z_off_t]; %Create top half-circle vector Eⁿ : XZ
  tp_shadow_list(end+1, :) = [x_t(i),0,-(sqrt(Rt.^2 - y_t(i).^2)-z_off_t)]; %Create bottom half-circle vector Eⁿ : XZ
  tp_shadow_list(end+1, :) = [x_t(i),0,T_top] ;% Tangent top trajectory line of target particle.
  tp_shadow_list(end+1, :) = [x_t(i),0,T_bot] ;% Tangent bottom trajectory line of target particle.
endfor

%Projectile particle shadow, input arrays.
y_p = linspace(-Rp,Rp,delineants);
x_p=y_p;
pp_shadow_list=[];
for i=1:length(y_t); % Define the projectile 2D shadow arrays onto YZ, XZ, XY
  pp_shadow_list(end+1, :) = [0,y_p(i),sqrt(Rp.^2 - y_p(i).^2)+z_off_p]; %Create top half-circle vector, Eⁿ : YZ
  pp_shadow_list(end+1, :) = [0,y_p(i),-(sqrt(Rp.^2 - y_p(i).^2)-z_off_p)]; %Create bottom half-circle vector, Eⁿ : YZ
  pp_shadow_list(end+1, :) = [sqrt(Rp.^2 - y_p(i).^2),y_p(i),z_off_p]; %Create top half-circle vector, Eⁿ : XY.
  pp_shadow_list(end+1, :) = [-sqrt(Rp.^2 - y_p(i).^2),y_p(i),z_off_p]; %Create bottom half-circle vector, Eⁿ : XY.
  pp_shadow_list(end+1, :) = [x_p(i),0,sqrt(Rp.^2 - y_p(i).^2)+z_off_p]; %Create top half-circle vector, Eⁿ : XZ.
  pp_shadow_list(end+1, :) = [x_p(i),0,-(sqrt(Rp.^2 - y_p(i).^2)-z_off_p)]; %Create bottom half-circle vector, Eⁿ : XZ.
  pp_shadow_list(end+1, :) = [x_p(i),0,P_top] ;% Tangent top trajectory line of projectile particle.
  pp_shadow_list(end+1, :) = [x_p(i),0,P_bot] ;% Tangent bottom trajectory line of projectile particle.
endfor

% Create random points in a prismatic volume tightly enclosing both particles: target and projectile.
prismPointsList=[];
for i=1:N
  % Create a single random point on every i iteration, initially normalized from 0-1 and later to fit in the scope of the particles hosting prism.
  if ((P_bot<T_top)&&(Rp<=Rt))
    xp=(2*Rt)*(rand()-0.5); %let a random unitary distribution 0-1 be centered at Eⁿ to represent a target particle unitary diameter, and then scaled to match target diameter (2*Rt)
    yp=(2*Rt)*(rand()-0.5);
    if (P_top<=T_top)
      zp=(T_top-T_bot)*rand();
    elseif (P_top>T_top)
      zp=(P_top-T_bot)*rand();
    endif
  elseif ((P_bot<T_top)&&(Rp>Rt))
    xp=(2*Rp)*(rand()-0.5);
    yp=(2*Rp)*(rand()-0.5);
    if (P_bot>T_bot)
      zp=(P_top-T_bot)*rand();
    elseif (P_bot<=T_bot)
      zp=((rand()-0.5)*(P_top-P_bot))+z_off_p;
    endif
  endif
  prismPointsList(end+1, :) = [xp,yp,zp];
endfor

%Select only the projectile nucleons and the projectile participant nucleons
projectileNucleons=[];
targetNucleons=[];
projectileParticipants=[];
targetParticipants=[];
for i=1:length(prismPointsList)
  tmpX=prismPointsList(i,1);
  tmpY=prismPointsList(i,2);
  tmpZ=prismPointsList(i,3);

  if( ((((tmpX).^2)+((tmpY).^2)+((tmpZ-z_off_p).^2))<Rp.^2) )
    projectileNucleons(end+1,:)=[tmpX,tmpY,tmpZ];
  endif
  if( ((((tmpX).^2)+((tmpY).^2)+((tmpZ-z_off_t).^2))<Rt.^2) )
    targetNucleons(end+1,:)=[tmpX,tmpY,tmpZ];
  endif
  if( ((((tmpY).^2)+(tmpZ-z_off_t).^2)<Rt.^2)&&((((tmpX).^2)+((tmpY).^2)+((tmpZ-z_off_p).^2))<Rp.^2) )
    projectileParticipants(end+1, :)=[tmpX,tmpY,tmpZ];
  endif
  if( ((((tmpY).^2)+(tmpZ-z_off_p).^2)<Rp.^2)&&((((tmpX).^2)+((tmpY).^2)+((tmpZ-z_off_t).^2))<Rt.^2) )
    targetParticipants(end+1, :)=[tmpX,tmpY,tmpZ];
  endif
endfor

%Fit plot to equal proportions in all axes and add labels.
if (plotIt==1)
  plot_fig=figure("visible", "off");
  plot_ax=axes("parent", plot_fig);

  % Projectile contour plot
  projectilePlot=plot3(plot_ax,pp_shadow_list(:,1),pp_shadow_list(:,2),pp_shadow_list(:,3),'b.', 'MarkerSize', del_thkns, 'LineWidth', 2);
  hold(plot_ax, 'on');

  % Target contour plot
  targetPlot=plot3(plot_ax,tp_shadow_list(:,1),tp_shadow_list(:,2),tp_shadow_list(:,3),'r.', 'MarkerSize', del_thkns, 'LineWidth', 2);
  hold(plot_ax, 'on');

  % Prism Points Plot
  %prismPointsPlot=plot3(plot_ax,prismPointsList(:,1),prismPointsList(:,2),prismPointsList(:,3),'k.', 'MarkerSize', del_thkns, 'LineWidth', 0.5, 'DisplayName', 'Acknowledgement Volume');

  % Projectile points plot
  %projectileNucleonsPlot=plot3(plot_ax,projectileNucleons(:,1),projectileNucleons(:,2),projectileNucleons(:,3),'k.', 'MarkerSize', del_thkns, 'LineWidth', 0.5);
  %hold(plot_ax, 'on');

  % Target points plot
  %targetNucleonsPlot=plot3(plot_ax,targetNucleons(:,1),targetNucleons(:,2),targetNucleons(:,3),'k.', 'MarkerSize', del_thkns, 'LineWidth', 0.5);
  %hold(plot_ax, 'on');

  % Projectile participant points plot
  projectileParticipantsPlot=plot3(plot_ax,projectileParticipants(:,1),projectileParticipants(:,2),projectileParticipants(:,3),'c.', 'MarkerSize', del_thkns, 'LineWidth', 0.5);
  %hold(plot_ax, 'on');

  % Target participant points plot
  targetParticipantsPlot=plot3(plot_ax,targetParticipants(:,1),targetParticipants(:,2),targetParticipants(:,3),'m.', 'MarkerSize', del_thkns, 'LineWidth', 0.5);

 % Axis and Title setup
  xlabel(plot_ax, 'X-axis');
  ylabel(plot_ax, 'Y-axis');
  zlabel(plot_ax, 'Z-axis');
  title(plot_ax, 'Acknowledgement volume');
  axis(plot_ax, 'equal');
  %view(plot_ax, 3);

  %view(plot_ax, 0, 90); %XY view (trajectory view, top)
  view(plot_ax, 155, 20); %pseudo-dimetric
  %view(plot_ax, 0, 0); %XZ view (trajectory view, side)
  %view(plot_ax, 90, 0); %YZ view (impact parameter view)
  %view(plot_ax, 0, 90); %XY view (trajectory view, top)

  grid(plot_ax, 'on');

  drawnow();

  % Save the interactive figure (the only one that will support toggling)
  savefig(plot_fig, "02-prism points.ofig");
  % Save the static image for the report (it will show everything ON by default)
  set(plot_fig, "paperpositionmode", "auto");
  print(plot_fig, "02-prism points.png", "-dpng");
  close(plot_fig);
endif

%Results report
disp("[00. ATOMIC MASSES]");
disp(sprintf("Ap: Projectile particle atomic mass: %d [u] ", Ap));
disp(sprintf("At: Target particle atomic mass: %d [u] ", At));
disp(" ");
disp("[01. NUCLEI RADII & MAXIMUM IMPACT PARAMETER]");
disp("∵ R=R₀⋅A^1/3 (Kenneth S Krane. ISBN: 0-471-80553-X, page 49.)");
disp(sprintf("∧ R₀= %.2f [fm] (Kenneth S. Krane. ISBN: 0-471-80553-X, page 49.) ", R0));
disp(sprintf("∧ A= %.2f [u] (user input, script variable Ap)", Ap));
disp(sprintf("⟹ R(Ap)=Rp= %.2f [fm]", (R0*(Ap^(1/3)))));
disp("");
disp("∵ R=R₀⋅A^1/3 (Kenneth S. Krane. ISBN: 0-471-80553-X, page 49.)");
disp(sprintf("∧ R₀= %.2f [fm] (Kenneth S. Krane. ISBN: 0-471-80553-X, page 49.) ", R0));
disp(sprintf("∧ A= %.2f [u] (user input, script variable Ap)", At));
disp(sprintf("⟹ R(At)=Rt= %.2f [fm]", (R0*(At^(1/3)))));
disp("");
disp("∵ ""...the impact parameter (the perpendicular distance from the center of the target nucleon to the line of flight of the incident nucleon). See: Kenneth S. Krane. ISBN: 0-471-80553-X, page 86. ")
disp("→b=Rp+Rt (maximum impact parameter, b, equals the projectile radius, Rp, plus the target radius, Rt.");
disp(sprintf("⟹ b=%.2f [fm] max.", Rp+Rt));
disp("");
disp(sprintf("Note: based in user prompt, b has been scaled to %.2f %%, where 0%% denotes center-to-center impact, whereas 100%% denotes maximum impact parameter b (projectile particle just grazed target particle).", Ep));
disp("");
disp(sprintf("The effective impact parameter in the analysis is: b = %.2f [fm]", b));
disp("");
disp("[02. PARTICIPANT-SPECTATOR FRAMEWORK]");
disp("While considering the ""Nuclear Fireball Model"" (see: https://doi.org/10.1103/PhysRevLett.37.1202), where the projectile particle collides through-all in straight trajectory line with the target particle, it then follows that ""...The swept-out nucleons from the projectile and target are called the fireball..."", where the fireball would constitute the participant nucleons, while the spectators would be the remainder nucleons.");
disp("");
disp("Let there be a right-handed Cartesian coordinate system in standard Euclidean space as per ISO 80000-2:2009(E), Eⁿ:");
disp([sprintf("The target particle of radius Rt= %.2f,", Rt), sprintf(" is envisioned still, with center at (0,0,%.2f)", z_off_t)]);
disp([sprintf("The projectile particle of radius Rp= %.2f,", Rp), sprintf(" is envisioned displacing along Eⁿ x-axis. For Monte-Carlo analysis purposes it is envisioned to be in an instant of time with its center at (0,0,%.2f).", z_off_p)]);
disp(sprintf("This array allowed to randomly scatter N=%.2f points within a prismatic volume tightly enveloping both the projectile and the target particles.", N));
disp("[03. VOLUMETRIC ANALYSIS]");
disp("ANALYTIC VOLUMES:");
prismaticVolume=(2*(2*max(Rp,Rt)))*(Rt+b+Rp); %Refers to the analytical volume containing all Monte Carlo N points. Obtained by choosing the largest particle radius, double it to obtain a prism base side. Double again to obtain prism base area. Multiply by the prism height, equivalent to the product of Rt+b+Rp.
disp(sprintf("Monte Carlo test volume: %.2f [fm³]", prismaticVolume));
projectileVolume=(4/3)*pi()*(Rp.^3); % see: https://mathworld.wolfram.com/Sphere.html
disp(sprintf("Projectile volume: %.2f [fm³]", projectileVolume));
targetVolume=(4/3)*pi()*(Rt.^3); % see: https://mathworld.wolfram.com/Sphere.html
disp(sprintf("Target volume: %.2f [fm³]", projectileVolume));
disp("");
disp("MONTE CARLO POINTS:");
disp(sprintf("Monte Carlo test points (encompassing the volumes of projectile and target), N: %.2f [points]", N));
projectileSequesteredPoints=length(projectileNucleons); %from all Monte Carlo points N found in prismaticVolume, save the number of points found within the volume of the projectile.
disp(sprintf("Projectile monte carlo points (assumed to be proportional to all projectile nucleons): %.2f [points]", projectileSequesteredPoints));
targetSequesteredPoints=length(projectileNucleons); %from all Monte Carlo points N found in prismaticVolume, save the number of points found within the volume of the target.
disp(sprintf("Target monte carlo points (assumed to be proportional to all target nucleons): %.2f [points]", targetSequesteredPoints));
disp("");
disp("PARTICIPANT POINTS:");
projectileParticipantPoints=length(projectileParticipants); %from all Monte Carlo points N found in prismaticVolume, save the number of participant projectile points.
disp(sprintf("Projectile participant monte carlo points: %.2f [points]", projectileParticipantPoints));
targetParticipantPoints=length(targetParticipants); %from all Monte Carlo points N found in prismaticVolume, save the number of participant target points.
disp(sprintf("Target participant monte carlo points: %.2f [points]", targetParticipantPoints));
disp("");
disp(sprintf("PARTICIPANT FRACTION (Impact parameter, b=%.2f [fm] implicit): ", b));
disp("");
disp("Let:");
disp("Participant fraction: Nₚ");
disp(sprintf("Target participant Monte Carlo points: TPP=%.2f [points]",targetParticipantPoints));
disp(sprintf("Projectile participant Monte Carlo points: PPP=%.2f [points]",projectileParticipantPoints));
disp(sprintf("Projectile total Monte Carlo points: PMP=%.2f [points]",projectileSequesteredPoints));
disp(sprintf("Target total Monte Carlo points: TMP=%.2f [points]",targetSequesteredPoints));
disp("");
disp("It follows that:");
disp(sprintf("∵ ""...The participant fraction is the ratio of the number of participants to the number of all nucleons..."" (C. Hartnack, ""Student Laboratories, Master SNEAM, Application of the Participant-Spectator-Model, Problem 1"") "));
Np=(targetParticipantPoints+projectileParticipantPoints)/(projectileSequesteredPoints+targetSequesteredPoints)
disp(sprintf("→ Nₚ = (TPP+PPP)/(PMP+TMP)", Np));
disp("");
disp(sprintf("⟹ Nₚ = %.2f", Np));

%Ensure Octave workspace is clean for other parallel projects
%clear;
%close all;



