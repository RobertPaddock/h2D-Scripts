%Meshing script. Attempts to mesh using 4 layers, with a split in the CH
%region. If this fails or the split occurs too close to the Ice/CH
%interface, then it reverts to a layer mesh. The split is chosen so that
%the percentage mass difference of the CH end zone equals that of the ice
%layer. Changes now made for higher pressures

%Figures==1 sets manual mode - specify the capsule dimensions below, and
%figures will be plotted for the meshing.
%close all
%clear all

QuartzType = 1; %Quartz type, 1= alpha, 2 = amorphous

%Select Values for Radii to vary over
CHLength = 0.0040; %Originally 0.002
GoldLength = 0.0003;
QuartzLength = 0.004;
FoamLength = 0.0040; %Originally 0.004, reduced to 0.02

%Insert Densities of the materials (if not specified in another script).
if QuartzType == 1
    DensityQuartz = 2.6; %Alpha quartz / silica density
else
    DensityQuartz = 2.2; %Amorphous quartz density
end
DensityFoam=2.530000e-01;
DensityCH = 1.040000e+00;
DensityGold=19.3;

% % Original mesh for 10% difference - gives around 3%
% % Front thickness is set, and back should be low too if want to resolve
% % details. Because foam is low density, if the quartz thickness 2 is too
% % high, then foam regions will need to be large and struggle to find
% % solution - so need to reduce quartz quite a lot too.
% % CHThickness1 = 0.00001; %Front res is 0.1um
% % CHThickness2 = 0.00005;
% % GoldThickness1 = CHThickness2*DensityCH/DensityGold;
% % GoldThickness2 = GoldThickness1*4;
% % QuartzThickness1 = GoldThickness2*DensityGold/DensityQuartz;
% % QuartzThickness2 = QuartzThickness1/5;
% % FoamThickness1 = QuartzThickness2*DensityQuartz/DensityFoam; %For thick foam
% % FoamThickness1 = 0.00001*DensityCH/DensityFoam; %For thin foam
% % QuartzThickness2 = FoamThickness1*DensityFoam/DensityQuartz; %For thin foam
% % FoamThickness2 = 0.00001;

% %6% mass difference
% CHThickness1 = 0.00001; %Front res is 0.1um
% CHThickness2 = 0.0002;
% GoldThickness1 = CHThickness2*DensityCH/DensityGold;
% GoldThickness2 = GoldThickness1*2;
% QuartzThickness1 = GoldThickness2*DensityGold/DensityQuartz;
% QuartzThickness2 = QuartzThickness1/7;
% FoamThickness1 = QuartzThickness2*DensityQuartz/DensityFoam; %For thick foam
% %FoamThickness1 = 0.00001*DensityCH/DensityFoam; %For thin foam
% %QuartzThickness2 = FoamThickness1*DensityFoam/DensityQuartz; %For thin foam
% FoamThickness2 = 0.00001;
% 
% %23th November - High precision meshing
% CHThickness1 = 0.00001; %Front res is 0.1um
% CHThickness2 = 0.0004;
% GoldThickness1 = CHThickness2*DensityCH/DensityGold;
% GoldThickness2 = GoldThickness1*2;
% QuartzThickness1 = GoldThickness2*DensityGold/DensityQuartz;
% QuartzThickness2 = QuartzThickness1/9;
% FoamThickness1 = QuartzThickness2*DensityQuartz/DensityFoam; %For thick foam
% %FoamThickness1 = 0.00001*DensityCH/DensityFoam; %For thin foam
% %QuartzThickness2 = FoamThickness1*DensityFoam/DensityQuartz; %For thin foam
% FoamThickness2 = 0.0001;


% %19th November - High precision meshing
% CHThickness1 = 0.00001; %Front res is 0.1um
% CHThickness2 = 0.0004;
% GoldThickness1 = CHThickness2*DensityCH/DensityGold;
% GoldThickness2 = GoldThickness1*3;
% QuartzThickness1 = GoldThickness2*DensityGold/DensityQuartz;
% QuartzThickness2 = QuartzThickness1/9;
% FoamThickness1 = QuartzThickness2*DensityQuartz/DensityFoam; %For thick foam
% %FoamThickness1 = 0.00001*DensityCH/DensityFoam; %For thin foam
% %QuartzThickness2 = FoamThickness1*DensityFoam/DensityQuartz; %For thin foam
% FoamThickness2 = 0.0001;

% %19th November - Low precision meshing
% CHThickness1 = 0.00001; %Front res is 0.1um
% CHThickness2 = 0.0006;
% GoldThickness1 = CHThickness2*DensityCH/DensityGold;
% GoldThickness2 = GoldThickness1*3;
% QuartzThickness1 = GoldThickness2*DensityGold/DensityQuartz;
% QuartzThickness2 = QuartzThickness1/9;
% FoamThickness1 = QuartzThickness2*DensityQuartz/DensityFoam; %For thick foam
% %FoamThickness1 = 0.00001*DensityCH/DensityFoam; %For thin foam
% %QuartzThickness2 = FoamThickness1*DensityFoam/DensityQuartz; %For thin foam
% FoamThickness2 = 0.0001;




[RatioCH, LayersCH] = ThicknessSolver(0, CHLength, CHThickness1, CHThickness2);
[RatioGold, LayersGold] = ThicknessSolver(0, GoldLength, GoldThickness1, GoldThickness2);
[RatioQuartz, LayersQuartz] = ThicknessSolver(0, QuartzLength, QuartzThickness1, QuartzThickness2);
[RatioFoam, LayersFoam] = ThicknessSolver(0, FoamLength, FoamThickness1, FoamThickness2);

LayersCH = round(LayersCH);
LayersGold = round(LayersGold);
LayersQuartz = round(LayersQuartz);
LayersFoam = round(LayersFoam);

% LayersGold = round(GoldLength/GoldThickness1);
% RatioGold=1;
% LayersQuartz = round(QuartzLength/QuartzThickness1); %comment out for thin foam
% RatioQuartz=1;%comment out for thin foam

[Radii, Density] = BoundarySolver(CHLength, GoldLength, QuartzLength, FoamLength, RatioCH, RatioGold, RatioQuartz, RatioFoam, LayersCH,  LayersGold, LayersQuartz, LayersFoam, DensityCH, DensityGold, DensityQuartz, DensityFoam)
[ZoneMass, MassDiff] = MassDifference(Radii, Density);




%Calculate number of layers, mass difference
NumLayers = LayersCH+LayersGold+LayersQuartz+LayersFoam;
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
ZoneMass = diff(Radii).*Density;
MassDiff = 100*(diff(ZoneMass)./ZoneMass(2:end));
end

%Optimisation function. Given the Ratios, Radii and Layers, will construct
%the mesh and calculate the mass difference using above functions. This
%solves for 3 layer systems
function [Radii, Density] = BoundarySolver(CHLength, GoldLength, AlLength, FoamLength, RatioCH, RatioGold, RatioAl, RatioFoam, LayersCH,  LayersGold, LayersAl, LayersFoam, DensityCH, DensityGold, DensityAl, DensityFoam)

BoundaryCH = LayersCH+1;
BoundaryGold = BoundaryCH+LayersGold;
BoundaryAl = BoundaryGold+LayersAl;
BoundaryFoam = BoundaryAl + LayersFoam;

CHzonepositions = RatioIncrement(1, BoundaryCH, 0, CHLength, RatioCH);
Goldzonepositions = RatioIncrement(BoundaryCH, BoundaryGold, CHLength, CHLength+GoldLength, RatioGold);
Alzonepositions = RatioIncrement(BoundaryGold, BoundaryAl, CHLength+GoldLength, CHLength+GoldLength+AlLength, RatioAl);
Foamzonepositions = RatioIncrement(BoundaryAl, BoundaryFoam, CHLength+GoldLength+AlLength, CHLength+GoldLength+AlLength+FoamLength, RatioFoam);

Density = [ DensityCH*ones(1, LayersCH), DensityGold*ones(1, (LayersGold)), DensityAl*ones(1, (LayersAl)),  DensityFoam*ones(1, (LayersFoam))];

Radii = [0, CHzonepositions, Goldzonepositions, Alzonepositions, Foamzonepositions];


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

