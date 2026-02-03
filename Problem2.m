%Student Laboratories
%Problem 2:

%Ensure Octave starts fresh
clear;
close all;

%Runtime variables
fs = 8000;              % Sampling frequency (Hz)
t = 0:1/fs:0.15;        % Duration of each note (0.15 seconds)

%Set project directory.
if (!strcmp("/home/javier/Documents/Octave/Nuclear Modelling/Problem1",pwd()))
  cd "/home/javier/Documents/Octave/Nuclear Modelling/Problem1"
endif

%[USER INPUTS]---------------------------------------------------------------------------
cylinderRadius=input("Enter cylinder radius in [mm], ρ: ");
cylinderHeight=input("Enter cylinder height in [mm], H: ");
meanFreePath=input("Enter the medium mean free path [mm], λ: ");
N=input("Number of Monte Carlo trials per iteration (e.g. 10000), N: ");    % Number of Monte Carlo random points to be tried in space.
stat_N=input("Number of Monte Carlo iterations (e.g. 100), Iₘ꜀: ");
plotIt=yes_or_no("Do you want to generate a 3D plot? (may take longer to process)");
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

%Form a matrix of random p¹(x,y,z) points filling the cylinder volume, volumePoints.
listOfDetectedPoints=[];
for i=1:stat_N
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
  %
  % DERIVATION OF ISOTROPIC EMISSION PARAMETERS (ISO 80000-2:2009(E))
  % Goal: Find the probability distribution for the polar angle θ (from pole-to-pole)
  % to ensure uniform area coverage on a sphere (Isotropy).
  %
  % 1. Differential Area Element (dA):
  %    dA = (R ⋅ dθ) ⋅ [(R ⋅ sinθ) ⋅ dφ]
  %    dA = R² ⋅ sinθ ⋅ dθ ⋅ dφ
  %    The probability density function (PDF) is thus proportional to sinθ.
  %
  % 2. Normalization of PDF (Weighting factor W):
  %    Let P(θ) = ∫ W * sinθ dθ over the interval [0, π].
  %    ∵ 1 = W * ∫[0,π] sinθ dθ
  %    → 1 = W * [-cosθ][0,π]
  %    → 1 = W * [(-(-1)) - (-1)]
  %    → 1 = W * (2)
  %    ∴ W = 1/2
  %
  % 3. Cumulative Distribution Function (CDF):
  %    Let u be a uniform random number u ∈ [0, 1].
  %    ∵ u = ∫[0,θ] (1/2) * sin(t) dt
  %    → u = (1/2) * [-cos(t)][0,θ]
  %    → u = (1/2) * [1 - cosθ]
  %    → 2u = 1 - cosθ
  %    ∴ cosθ = 1 - 2u
  %    (Note: Since u is uniform [0,1], 1-2u is statistically identical to 2u-1)
  %
  % 4. Pythagorean Identity for sinθ:
  %    ∵sin²θ + cos²θ = 1 (Pythagorean theorem, see: https://mathworld.wolfram.com/PythagoreanTheorem.html)
  %    ∴ sinθ = sqrt(1 - cos²θ)
  radiationPoints=[];
  r2 = 0;
  phi2 = 0;
  sin_theta2 = 0;
  cos_theta2=0;
  px2 = 0;
  py2 = 0;
  pz2 = 0;
  for i = 1:length(volumePoints)
    r2=-meanFreePath*reallog(rand());
    phi2=2*pi()*rand();
    cos_theta2=1-(2*rand());
    sin_theta2=sqrt(1-(cos_theta2^2));
    px2=(r2*cos(phi2)*sin_theta2)+volumePoints(i,1);
    py2=(r2*sin(phi2)*sin_theta2)+volumePoints(i,2);
    pz2=(r2*cos_theta2)+volumePoints(i,3);
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

%ERROR PROPAGATION
%
%From the cuadrature principle: "...Relative errors add in quadrature." (Marek Gierliński. ISBN: 9781119106890. "Understanding Statistical Error : A Primer for Biologists", Chapter 7.3 Multiple variables.)
% Let:
% μ(Nₑ): mean of detected e⁻ points (or, script variable: meanDetectedPoints).
% μⱼ: mean of total emitted e⁻ points (or, user input constant, N).
% σₑ: standard deviation of detected e⁻ points (or, script variable: standardDeviationOfDetectedPoints).
% σⱼ: standard deviation of total emitted e⁻ points.
% σₚ: standard deviation of the ratio of detected e⁻ points and total e⁻ points (or, script variable: propagatedStandardDeviation).
% ηₑ: ratio of detected e⁻ points and total e⁻ points (or, script varaible: detectionEfficiency).
%
%∵ σₚ/ηₑ=sqrt(((σₑ/μ(Nₑ))^2)+((σⱼ/μⱼ)^2)))  (the Cuadrature Principle)
%∧ σⱼ=0 	(because of user input constant, N)
% →σₚ/ηₑ=sqrt(((σₑ/μ(Nₑ))^2)+((0/μⱼ)^2))
% →σₚ=σₑ⋅(ηₑ/μ(Nₑ))
% ∧ηₑ=μ(Nₑ)/μⱼ
% ⇒σₚ=σₑ⋅(1/μⱼ)
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
  savefig(plot_fig, "Problem2-3Dplot.ofig");
  % Save the static image for the report (it will show everything ON by default)
  set(plot_fig, "paperpositionmode", "auto");
  print(plot_fig, "Problem2-3Dplot.png", "-dpng");
  close(plot_fig);
endif

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
  targetPlot=plot3(plot_ax,volumePoints(:,1),volumePoints(:,2),volumePoints(:,3),'m.', 'MarkerSize', delineantThickness/4, 'LineWidth', 2);
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
  savefig(plot_fig, "Problem2-MCpoints.ofig");
  % Save the static image for the report (it will show everything ON by default)
  set(plot_fig, "paperpositionmode", "auto");
  print(plot_fig, "Problem2-MCpoints.png", "-dpng");
  close(plot_fig);
endif

%[TERMINAL REPORT SECTION]---------------------------------------------------------------------------------
disp("[00. KEY INPUT PARAMETERS & VARIABLES]");
disp(sprintf("Number of scattered electrons from β⁻-decay (user input), N: %.2f [e⁻]", N));
disp(sprintf("Cylinder radius, ρ: %d [mm] ", cylinderRadius));
disp(sprintf("Cylinder height (entirely filled), H: %d [mm] ", cylinderHeight));
analyticalVolume= (pi()*(cylinderRadius^2))*cylinderHeight;% see: https://mathworld.wolfram.com/Cylinder.html
disp(sprintf("Analytical cylinder volume, Vₐ: %.5f [mm³] (∝π⋅ρ²⋅H, see: https://mathworld.wolfram.com/Cylinder.html) ", analyticalVolume));
disp(" ");
disp("[01. DETECTION EFFICIENCY]");
disp(sprintf("Number of Monte Carlo iterations, Iₘ꜀: %.2f [iterations]", stat_N));
disp(sprintf("Mean of the number of detected electrons from β⁻-decay, μ(Nₑ): %.5f±%.5f[e⁻] ", meanDetectedPoints,standardDeviationOfDetectedPoints));
disp(" ");
disp("∵ ""Detection efficiency"", ηₑ, is understood as the mean of the number of detected electrons from β⁻-decay, μ(Nₑ), compared to the total number of emitted electrons, N. ");
disp("→ ηₑ = μ(Nₑ)/N");
disp(sprintf("→ ηₑ = (%.5f[e⁻])/(%.5f[e⁻])", meanDetectedPoints, N));
disp(sprintf("⟹ ηₑ = %.5f±%.5f", detectionEfficiency, propagatedStandardDeviation));
disp(" ");
disp("∵ ""Detection yield"", Y, is understood as the arithmetic product of detection efficiency, ηₑ, and analytical cylinder volume, Vₐ. ");
disp("→ Y = ηₑ⋅Vₐ");
detectionYield=detectionEfficiency*analyticalVolume;
disp(sprintf("→ Y = (%.5f[e⁻])⋅(%.5f[mm³])", detectionEfficiency, analyticalVolume));
disp(sprintf("⟹ Y = %.5f [mm³]", detectionYield));

%[EXTERNAL RAW DATA]
listOfResults=[N,stat_N,meanFreePath,cylinderRadius,cylinderHeight,analyticalVolume,meanDetectedPoints,standardDeviationOfDetectedPoints,detectionEfficiency,propagatedStandardDeviation,detectionYield];
dlmwrite('Problem2.csv',listOfResults,'-append', 'delimiter', ',', 'newline', 'pc');

% [END OF SCRIPT NOTIFICATION MELODY]
notes = [523.25, 659.25, 783.99, 1046.50];
melody = [];

for f = notes
    % Create a sine wave for the frequency 'f' with a quick fade-out
    wave = sin(2 * pi * f * t) .* exp(-4 * t);
    melody = [melody, wave];
end

sound(melody, fs);
