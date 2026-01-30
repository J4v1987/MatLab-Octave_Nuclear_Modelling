%Student Laboratories
%Problem 2:

%Ensure Octave starts fresh
clear;
close all;

%Set project directory.
if (!strcmp("/home/javier/Documents/Octave/Nuclear Modelling/Problem1",pwd()))
  cd "/home/javier/Documents/Octave/Nuclear Modelling/Problem1"
endif

%[USER INPUTS]---------------------------------------------------------------------------
cylinderRadius=input("Enter cylinder radius in [mm]: ");
cylinderHeight=input("Enter cylinder height in [mm]: ");
meanFreePath=input("Enter the medium mean free path [mm]: ");
N=input("Number of Monte Carlo trials per iteration (e.g. 10000): ");    % Number of Monte Carlo random points to be tried in space.
stat_N=input("Number of Monte Carlo iterations (e.g. 100): ");
plotIt=yes_or_no("Do you want to see the analysis 3D plot? (may take longer to process)");
%[END OF USER INPUTS]---------------------------------------------------------------------------

disp(char(10)) ;
disp(["PROCESSING..."]) ;
disp(char(10)) ;

%Algorithm constants:
delineantPoints=cylinderHeight*round(2/exp(1));%3+round(30/cylinderHeight);%round(cylinderHeight*(4/5));%15;
delineantThickness=8;

%Trigonometric constants:
octantAngle=1/sqrt(2);

%Construction of algorithm variables:

%Form a matrix of p⁰(x,y,z) points delineating the contour of the cylinder bases, cylinderBaseProfile.
cylinderBaseProfileX = linspace(-cylinderRadius,cylinderRadius,15+round(2/(cylinderRadius+3)));%*(round(cylinderHeight/(15+exp(3)))));
cylinderBaseProfile = [];
for i=1:length(cylinderBaseProfileX); % Define the target 2D shadow arrays onto YZ, XZ, XY
  cylinderBaseProfile(end+1, :) = [cylinderBaseProfileX(i), sqrt(cylinderRadius^2-cylinderBaseProfileX(i)^2), 0]; %Create top half-circle vector in Eⁿ : XY
  cylinderBaseProfile(end+1, :) = [cylinderBaseProfileX(i), -sqrt(cylinderRadius^2-cylinderBaseProfileX(i)^2), 0]; %Create bottom half-circle vector in Eⁿ : XY
  cylinderBaseProfile(end+1, :) = [cylinderBaseProfileX(i), sqrt(cylinderRadius^2-cylinderBaseProfileX(i)^2), cylinderHeight]; %Create top half-circle vector in Eⁿ : XY+cylinderHeight
  cylinderBaseProfile(end+1, :) = [cylinderBaseProfileX(i), -sqrt(cylinderRadius^2-cylinderBaseProfileX(i)^2), cylinderHeight]; %Create bottom half-circle vector in Eⁿ : XY+cylinderHeight
endfor

cylinderHeightFaceProfile = linspace(0,cylinderHeight,delineantPoints);
cylinderHeightFace = [];
for i=1:length(cylinderHeightFaceProfile);
  cylinderHeightFace(end+1,:)=[cylinderRadius,0,cylinderHeightFaceProfile(i)];
  cylinderHeightFace(end+1,:)=[-cylinderRadius,0,cylinderHeightFaceProfile(i)];
  cylinderHeightFace(end+1,:)=[0,cylinderRadius,cylinderHeightFaceProfile(i)];
  cylinderHeightFace(end+1,:)=[0,-cylinderRadius,cylinderHeightFaceProfile(i)];
  cylinderHeightFace(end+1,:)=[cylinderRadius*octantAngle,cylinderRadius*octantAngle,cylinderHeightFaceProfile(i)];
  cylinderHeightFace(end+1,:)=[-cylinderRadius*octantAngle,-cylinderRadius*octantAngle,cylinderHeightFaceProfile(i)];
  cylinderHeightFace(end+1,:)=[cylinderRadius*octantAngle,-cylinderRadius*octantAngle,cylinderHeightFaceProfile(i)];
  cylinderHeightFace(end+1,:)=[-cylinderRadius*octantAngle,cylinderRadius*octantAngle,cylinderHeightFaceProfile(i)];
endfor

listOfDetectedPoints=[];
for i=1:100
  %Form a matrix of random p¹(x,y,z) points filling the cylinder volume, volumePoints.
  volumePoints=[];
  r1=0;
  phi1=0;
  px1=0;
  py1=0;
  pz1=0;
  for i=1:N;
    r1=cylinderRadius*sqrt(rand());
    phi1=2*pi()*rand();
    px1=r1*cos(phi1);
    py1=r1*sin(phi1);
    pz1=cylinderHeight*rand();
    volumePoints(end+1,:)=[px1,py1,pz1];
  endfor

  %Form a matrix of random p²(x,y,z) points weighting the volumetric distribution of p¹(x,y,z) points to meet the essential β⁻-decay irradiation distribution pattern in 3D space.
  radiationPoints=[];
  r2 = 0;
  phi2 = 0;
  theta2 = 0;
  px2 = 0;
  py2 = 0;
  pz2 = 0;
  for i = 1:length(volumePoints)
    r2=-meanFreePath*reallog(rand());
    phi2 = 2*pi()*rand();
    theta2 = pi()*rand();
    px2 = r2*cos(phi2)+volumePoints(i,1);
    py2 = r2*sin(phi2)+volumePoints(i,2);
    pz2 = r2*cos(theta2)+volumePoints(i,3);
    radiationPoints(end+1,:)=[px2,py2,pz2];
  endfor

  %Form a matrix of random p³(x,y,z) points, located above z=cylinderHeight, and within the extrusion scope along z-axis of x²+y²=cylinderRadius².
  detectedPoints=[];
  px3=0;
  py3=0;
  pz3=0;
  for i=1:length(radiationPoints)
    px3=radiationPoints(i,1);
    py3=radiationPoints(i,2);
    pz3=radiationPoints(i,3);
    if((((px3)^2+(py3)^2)<((cylinderRadius)^2)) && (pz3>cylinderHeight))
      detectedPoints(end+1,:)=[px3,py3,pz3];
    endif
  endfor
  listOfDetectedPoints(end+1)=length(detectedPoints);
endfor

%[STATISTICAL ANALYSIS]----------------------------------------------------------------------------------------

%Mean of detected points, μ(Nₑ)
%See: context & formula 2.6.1, in: https://stats.libretexts.org/Bookshelves/Introductory_Statistics/Introductory_Statistics_1e_(OpenStax)/02%3A_Descriptive_Statistics/2.06%3A_Measures_of_the_Center_of_the_Data
sumDetectedPoints=0;
for i=1:length(listOfDetectedPoints)
  sumDetectedPoints+=listOfDetectedPoints(i);
endfor
meanDetectedPoints=sumDetectedPoints/length(listOfDetectedPoints);

%Standard deviation of detected points, σ(Nₑ)
%See: formula 2.8.3, in: https://stats.libretexts.org/Bookshelves/Introductory_Statistics/Introductory_Statistics_1e_(OpenStax)/02%3A_Descriptive_Statistics/2.08%3A_Measures_of_the_Spread_of_the_Data
sumOfSquaredDifferences=0;
for i=1:length(listOfDetectedPoints)
  sumOfSquaredDifferences+=(listOfDetectedPoints(i)-meanDetectedPoints)^2;
endfor
standardDeviationOfDetectedPoints=(sumOfSquaredDifferences/length(listOfDetectedPoints))^(1/2);

%Error propagation
%From the quadrature rule: "...Relative errors add in quadrature." (Marek Gierliński. ISBN: 9781119106890. "Understanding Statistical Error : A Primer for Biologists", Chapter 7.3 Multiple variables.)
% Let:
% μᵢ: mean of detected e⁻ points (or, script variable: meanDetectedPoints).
% μⱼ: mean of total emitted e⁻ points (or, user input constant, N).
% σᵢ: standard deviation of detected e⁻ points (or, script variable: standardDeviationOfDetectedPoints).
% σⱼ: standard deviation of total emitted e⁻ points.
% σₑ: standard deviation of the ratio of detected e⁻ points and total e⁻ points (or, script variable: propagatedStandardDeviation).
% ηₑ: ratio of detected e⁻ points and total e⁻ points (or, script varaible: detectionEfficiency).
%
% ∵σₑ/ηₑ=sqrt(((σᵢ/μᵢ)^2)+((σⱼ/μⱼ)^2))
% ∧σⱼ=0 (user input constant, N)
% →σₑ/ηₑ=sqrt(((σᵢ/μᵢ)^2)+((0/μⱼ)^2))
% →σₑ=σᵢ⋅(ηₑ/μᵢ)
% ∧ηₑ=μᵢ/μⱼ
% ⇒σₑ=σᵢ⋅(1/μⱼ)
detectionEfficiency=meanDetectedPoints/N;
propagatedStandardDeviation=standardDeviationOfDetectedPoints*(1/N);

%[PLOTTER SECTION, LAST MONTE CARLO ITERATION ONLY]-------------------------------------
%Plot Target monte carlo points
if (plotIt==1)
  plot_fig=figure("visible", "off");
  plot_ax=axes("parent", plot_fig);

  % Cylinder contour plots of side faces
  projectilePlot=plot3(cylinderBaseProfile(:,1),cylinderBaseProfile(:,2),cylinderBaseProfile(:,3),'r.', 'MarkerSize', delineantThickness, 'LineWidth', 2);
  hold(plot_ax, 'on');

  % Cylinder contour plot of height face
  targetPlot=plot3(plot_ax,cylinderHeightFace(:,1),cylinderHeightFace(:,2),cylinderHeightFace(:,3),'r.', 'MarkerSize', delineantThickness, 'LineWidth', 2);
  hold(plot_ax, 'on');

  % Cylinder volume points
  %targetPlot=plot3(plot_ax,volumePoints(:,1),volumePoints(:,2),volumePoints(:,3),'m.', 'MarkerSize', delineantThickness/4, 'LineWidth', 2);
  %hold(plot_ax, 'on');

  % Radiation points
  targetPlot=plot3(plot_ax,radiationPoints(:,1),radiationPoints(:,2),radiationPoints(:,3), '.', 'LineStyle', 'none', 'Color', [0, 0.5, 0], 'MarkerSize', delineantThickness/2, 'LineWidth', 2);
  hold(plot_ax, 'on');

  % Detected points
  targetPlot=plot3(plot_ax,detectedPoints(:,1),detectedPoints(:,2),detectedPoints(:,3),'b.', 'MarkerSize', delineantThickness/2, 'LineWidth', 2);
  %hold(plot_ax, 'on');

 % Axis and Title setup
  xlabel(plot_ax, 'X-axis');
  ylabel(plot_ax, 'Y-axis');
  zlabel(plot_ax, 'Z-axis');
  title(plot_ax, 'Cylinder Plot');
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
  savefig(plot_fig, "Problem2.ofig");
  % Save the static image for the report (it will show everything ON by default)
  set(plot_fig, "paperpositionmode", "auto");
  print(plot_fig, "Problem2.png", "-dpng");
  close(plot_fig);
endif

%[TERMINAL REPORT SECTION]---------------------------------------------------------------------------------

disp("[00. KEY INPUT PARAMETERS & VARIABLES]");
disp(sprintf("Number of scattered electrons from β⁻-decay (user input), N: %.2f [e⁻]", N));
disp(sprintf("Cylinder radius, ρ: %d [mm] ", cylinderRadius));
disp(sprintf("Cylinder height (entirely filled), H: %d [mm] ", cylinderRadius));
analyticalVolume= (pi()*(cylinderRadius^2))*cylinderHeight;% see: https://mathworld.wolfram.com/Cylinder.html
disp(sprintf("Analytical cylinder volume, Vₐ: %.5f [mm] (∝π⋅ρ²⋅H, see: https://mathworld.wolfram.com/Cylinder.html) ", analyticalVolume));
disp(" ");
disp("[01. DETECTION EFFICIENCY]");
disp(sprintf("Mean of the number of detected electrons from β⁻-decay, μ(Nₑ): %.5f±%.5f[e⁻] ", meanDetectedPoints,standardDeviationOfDetectedPoints));
disp("∵ ""Detection efficiency, ηₑ, is understood as the mean of the number of detected electrons from β⁻-decay, μ(Nₑ), compared to the total number of emitted electrons, N."" ");
disp("→ ηₑ = μ(Nₑ)/N");
disp(sprintf("→ ηₑ = (%.5f[e⁻])/(%.5f[e⁻])", meanDetectedPoints, N));
disp(sprintf("⟹ ηₑ = %.5f±%.5f", detectionEfficiency, propagatedStandardDeviation));

