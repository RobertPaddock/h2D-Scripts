%Meshing script. Attempts to mesh using 4 layers, with a split in the CH
%region. If this fails or the split occurs too close to the Ice/CH
%interface, then it reverts to a layer mesh. The split is chosen so that
%the percentage mass difference of the CH end zone equals that of the ice
%layer. Changes now made for higher pressures

%Figures==1 sets manual mode - specify the capsule dimensions below, and
%figures will be plotted for the meshing.
close all
clear all

QuartzType = 1; %Quartz type, 1= alpha, 2 = amorphous

% %Select Values for Radii to vary over - original meshing, with 1mm diameter
% %target
% AxialLength1 = 0.02;
% AxialLength2 = 0.05-AxialLength1;
% DensityAxial = 1;
% 
% Axial1Thickness1 = 0.001;
% Axial1Thickness2 = 0.0015;
% Axial2Thickness1 = Axial1Thickness2*1.1;
% Axial2Thickness2 = 0.005;

% 500um diameter target
AxialLength1 = 0.01;
AxialLength2 = 0.025-AxialLength1;
AxialLength2 = 0.025-AxialLength1;
DensityAxial = 1;

Axial1Thickness1 = 0.001;
Axial1Thickness2 = 0.0015;
Axial2Thickness1 = Axial1Thickness2*1.1;
Axial2Thickness2 = 0.001;
% 
% %500um diameter target
% AxialLength1 = 0.0175;
% % AxialLength2 = 0.025-AxialLength1;
% AxialLength2 = 0.035-AxialLength1;
% DensityAxial = 1;
% 
% Axial1Thickness1 = 0.001;
% Axial1Thickness2 = 0.005;
% Axial2Thickness1 = Axial1Thickness2*1.1;
% Axial2Thickness2 = 0.001;
% 
% %500um diameter target
% AxialLength1 = 0.0175;
% % AxialLength2 = 0.025-AxialLength1;
% AxialLength2 = 0.035-AxialLength1;
% DensityAxial = 1;
% 
% Axial1Thickness1 = 0.003;
% Axial1Thickness2 = 0.005;
% Axial2Thickness1 = Axial1Thickness2*1.1;
% Axial2Thickness2 = 0.003;

% 
% %500um diameter target - higher res
% AxialLength1 = 0.01;
% AxialLength2 = 0.025-AxialLength1;
% DensityAxial = 1;
% 
% Axial1Thickness1 = 0.0003;
% Axial1Thickness2 = 0.0005;
% Axial2Thickness1 = Axial1Thickness2*1.1;
% Axial2Thickness2 = 0.0003;


% %440um diameter target - lower res
% AxialLength1 = 0.011;
% AxialLength2 = 0.022-AxialLength1;
% DensityAxial = 1;
% 
% Axial1Thickness1 = 0.001;
% Axial1Thickness2 = 0.0015;
% Axial2Thickness1 = Axial1Thickness2;
% Axial2Thickness2 = 0.001;

% %500um diameter target - extreme feathering
% AxialLength1 = 0.023;
% AxialLength2 = 0.025-AxialLength1;
% DensityAxial = 1;
% 
% Axial1Thickness1 = 0.001;
% Axial1Thickness2 = 0.0002;
% Axial2Thickness1 = Axial1Thickness2*1.1;
% Axial2Thickness2 = 0.00001;


% %420um diameter target - higher res
% AxialLength1 = 0.01;
% AxialLength2 = 0.021-AxialLength1;
% DensityAxial = 1;
% 
% Axial1Thickness1 = 0.0003;
% Axial1Thickness2 = 0.0005;
% Axial2Thickness1 = Axial1Thickness2*1.1;
% Axial2Thickness2 = 0.0003;


% %420um diameter target - higher res
% AxialLength1 = 0.01;
% AxialLength2 = 0.021-AxialLength1;
% DensityAxial = 1;
% 
% Axial1Thickness1 = 0.001;
% Axial1Thickness2 = 0.0015;
% Axial2Thickness1 = Axial1Thickness2*1.1;
% Axial2Thickness2 = 0.001;

% %New test meshing - just run these values directly without the script.
% AxialLength1 = 0.021;
% AxialLength2 = 0.009;
% RatioAxial1 = 1;
% RatioAxial2 =1;
% LayersAxial1 = 14;
% LayersAxial2 = 3;

[RatioAxial1, LayersAxial1] = ThicknessSolver(0, AxialLength1, Axial1Thickness1, Axial1Thickness2);
[RatioAxial2, LayersAxial2] = ThicknessSolver(0, AxialLength2, Axial2Thickness1, Axial2Thickness2);

LayersAxial1 = round(LayersAxial1);
LayersAxial2 = round(LayersAxial2);

[Axial, Density] = BoundarySolver(AxialLength1, AxialLength2, RatioAxial1, RatioAxial2, LayersAxial1,  LayersAxial2, DensityAxial)
[ZoneMass, MassDiff] = MassDifference(Axial, Density);

AxialThicknessesOverall = diff(Axial);



%Calculate number of layers, mass difference
NumLayers = LayersAxial1+LayersAxial2
MaxDiff = max(abs(MassDiff(1:end)));

%Plot if in manual mode

figure
plot(MassDiff)
ylabel('Percentage Mass Difference')
yyaxis right;
plot(ZoneMass, ':');
ylabel('Zone Mass')
title('Zoning Plot')
xlabel('Zone')



% Testing - not part of code
% [Axial, Density] = BoundarySolver(0.021, 0.009, 1, 1, 14,  3, 1)
% [ZoneMass, MassDiff] = MassDifference(Axial, Density);
% 
% AxialThicknessesOverall = diff(Axial);
% 
% 
% 
% %Calculate number of layers, mass difference
% NumLayers = LayersAxial1+LayersAxial2
% MaxDiff = max(abs(MassDiff(1:end)));
% 
% %Plot if in manual mode
% 
% figure
% plot(MassDiff)
% ylabel('Percentage Mass Difference')
% yyaxis right;
% plot(ZoneMass, ':');
% ylabel('Zone Mass')
% title('Zoning Plot')
% xlabel('Zone')














%Function to recreate the Hyades mesh function - given the mesh boundaries
%j and radii to stretch between r and the ratio, it will create the mesh
%coordinates.
function Radii = RatioIncrement(j1, j2, r1, r2, Ratio)
Index = (1+j1-j1):(j2-j1);
Increment = Ratio.^(Index);
Radius = cumsum(Increment);
Scaling = (r2-r1)/Radius(end);
Radii = [r1+(Radius*Scaling)];
end

%Given the radii and density, will calculate the zone masses and return
%the percentage mass difference
function [ZoneMass, MassDiff] = MassDifference(Radii, Density)
Vol = pi() * (Radii.^2);
ZoneMass = diff(Vol).*Density;
MassDiff = 100*(diff(ZoneMass)./ZoneMass(2:end));
end

%Optimisation function. Given the Ratios, Radii and Layers, will construct
%the mesh and calculate the mass difference using above functions. This
%solves for 3 layer systems
function [Radii, Density] = BoundarySolver(AxialLength1, AxialLength2, RatioAxial1, RatioAxial2, LayersAxial1,  LayersAxial2, DensityAxial)

BoundaryAxial1 = LayersAxial1+1;
BoundaryAxial2 = BoundaryAxial1+LayersAxial2;

Axial1zonepositions = RatioIncrement(1, BoundaryAxial1, 0, AxialLength1, RatioAxial1);
Axial2zonepositions = RatioIncrement(BoundaryAxial1, BoundaryAxial2, AxialLength1, AxialLength1+AxialLength2, RatioAxial2);

Density = [ DensityAxial*ones(1, LayersAxial1), DensityAxial*ones(1, (LayersAxial2))];

Radii = [0, Axial1zonepositions, Axial2zonepositions];


end

%Optimisation function. Given the Ratios, Radii and Layers, will construct
%the mesh and calculate the mass difference using above functions. This
%solves for 4 layers sytems where the split is in the CH.
function [Radii, Density] = BoundarySolver4Layer(RatioVapour, RatioIce, RatioCH1,RatioCH2, LayersVapour,  LayersIce, LayersCH1,LayersCH2, RadiusVapour, RadiusIce,RadiusSplit, RadiusCH, DensityVapour, DensityIce, DensityCH)

BoundaryVapour = LayersVapour+1;
BoundaryIce = BoundaryVapour + LayersIce;
BoundaryCH1 = BoundaryIce + LayersCH1;
BoundaryCH2 = BoundaryCH1 + LayersCH2;

VapourRadii = RatioIncrement(1, BoundaryVapour, 0, RadiusVapour, RatioVapour);
IceRadii = RatioIncrement(BoundaryVapour, BoundaryIce, RadiusVapour, RadiusIce, RatioIce);
CH1Radii = RatioIncrement(BoundaryIce, BoundaryCH1, RadiusIce, RadiusSplit, RatioCH1);
CH2Radii = RatioIncrement(BoundaryCH1, BoundaryCH2, RadiusSplit, RadiusCH, RatioCH2);

Density = [DensityVapour*ones(1, LayersVapour), DensityIce*ones(1, (LayersIce)), DensityCH*ones(1, LayersCH1), DensityCH*ones(1, LayersCH2)];

Radii = [0, VapourRadii, IceRadii, CH1Radii, CH2Radii];

end


%Optimisation function. Given the Ratios, Radii and Layers, will construct
%the mesh and calculate the mass difference using above functions. This
%solves for 4 layers sytems where the split is in the ice.
function [Radii, Density] = BoundarySolver4Layer2(RatioVapour, RatioIce1, RatioIce2,RatioCH, LayersVapour,  LayersIce1, LayersIce2,LayersCH, RadiusVapour, RadiusIce1,RadiusIce, RadiusCH, DensityVapour, DensityIce, DensityCH)

BoundaryVapour = LayersVapour+1;
BoundaryIce1 = BoundaryVapour + LayersIce1;
BoundaryIce = BoundaryIce1 + LayersIce2;
BoundaryCH = BoundaryIce + LayersCH;

VapourRadii = RatioIncrement(1, BoundaryVapour, 0, RadiusVapour, RatioVapour);
Ice1Radii = RatioIncrement(BoundaryVapour, BoundaryIce1, RadiusVapour, RadiusIce1, RatioIce1);
IceRadii = RatioIncrement(BoundaryIce1, BoundaryIce, RadiusIce1, RadiusIce, RatioIce2);
CHRadii = RatioIncrement(BoundaryIce, BoundaryCH, RadiusIce, RadiusCH, RatioCH);

Density = [DensityVapour*ones(1, LayersVapour), DensityIce*ones(1, (LayersIce1)), DensityIce*ones(1, LayersIce2), DensityCH*ones(1, LayersCH)];

Radii = [0, VapourRadii, Ice1Radii, IceRadii, CHRadii];

end


%Function to solve for optimal ratio and number of layers. Given the radii
%and desired thickness at either end of the region, this function will
%solve the two simultaneous equations to determine the optimal number of
%layers and ratio to use.
function [Ratio, Layers] = ThicknessSolver(LowerRadius, UpperRadius, LowerThickness, UpperThickness)
%Turn off numerical solver warning message
warning('off','symbolic:solve:FallbackToNumerical');

%Solve analytic formula to find ratio and layers for ice region .
syms a n
eqn1 = log(1 - (1-a)*(UpperRadius-LowerRadius)/LowerThickness)/log(a) == n;
eqn2 = a == (UpperThickness/LowerThickness)^(1/(n-1));
sol = solve([eqn1, eqn2], [a, n], 'Real', true);
aSol = sol.a;
nSol = sol.n;
Ratio=double(aSol);
Layers=double(nSol);

end

