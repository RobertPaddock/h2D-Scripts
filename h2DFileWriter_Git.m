%Script to write Hyades input decks. Can write for three or four layers
%(specified by Regions) and three or 4 pulses.

RadialZones=30;
RadialLength = 0.03;

Pulse1Power = 100; %100 for full, 25 for quarter
LaserEndTime = 8E-9;
RiseTime = 2E-10;
ChangeTimeAt = [0 8]*10^-9;
TimeSpacing = [1E-11]; 
GaussPulse = 0;

%Output timing
%ChangeTimeAt = [0 16 22]*10^-9;
%TimeSpacing = [1E-10 1E-11];
Time=[0];
for i=1:(length(ChangeTimeAt)-1)
    Time = [Time, (ChangeTimeAt(i)+TimeSpacing(i)):TimeSpacing(i):ChangeTimeAt(i+1)];
end


%Densities of the three materials
IonTempAll = 1.551000e-06;

%Calculates relevant boundaries from layers
BoundaryCH = LayersCH+1;
BoundaryGold = BoundaryCH+LayersGold;
BoundaryQuartz = BoundaryGold+LayersQuartz;
BoundaryFoam = BoundaryQuartz + LayersFoam;
BoundaryAxial1=LayersAxial1+1;
BoundaryAxial2 = BoundaryAxial1+LayersAxial2;


%Converts Pulse Power in TW into erg, and reduces to account for cross beam
%energy transfer
Pulse1erg = Pulse1Power * 10^19 * pi() * (0.02^2);
if exist('Pulse2Time','var')==1
Pulse2erg = Pulse2Power * 0.8*10^19;
end

%Print the file
filename = sprintf('File');
filepath = sprintf(['MultiFile/', filename, '.inf']);
fileID = fopen(filepath,'w');

fprintf(fileID,'%s\r\nc\r\n', filename);

%K zones are spaced in Z, and l zones are spaced in R
% Mesh k1 l1 k2 l2 r1 z1 r2 z2
% Ratio k1 l1 k2 l2 rat
% if k2=k1, thenpoints beween l1 and l2 (i.e. R points) are mapped
% according to rat - so we want l1=l2
% Spacing k1 l1 k2 l2 k3 l3

% CH ablator layer
fprintf(fileID,'mesh 1 1 %.0f %.0f 0. 0. %.6e %.6e \r\n', BoundaryCH, BoundaryAxial1, AxialLength1, CHLength);
fprintf(fileID,'ratio 1 1 %.0f 1 %.6e \r\n', BoundaryCH, RatioCH);
fprintf(fileID,'spacing 1 1 %.0f %.0f 0 1 \r\n', BoundaryCH, BoundaryAxial1);
fprintf(fileID,'ratio 1 1 1 %.0f %.6e \r\n', BoundaryAxial1, RatioAxial1);
fprintf(fileID,'spacing 1 1 %.0f %.0f 1 0 \r\n', BoundaryCH, BoundaryAxial1);

fprintf(fileID,'mesh 1 %.0f %.0f %.0f %.6e 0. %.6e %.6e \r\n', BoundaryAxial1, BoundaryCH, BoundaryAxial2, AxialLength1, AxialLength1+AxialLength2, CHLength);
fprintf(fileID,'ratio 1 %.0f %.0f %.0f %.6e \r\n', BoundaryAxial1, BoundaryCH, BoundaryAxial1, RatioCH);
fprintf(fileID,'spacing 1 %.0f %.0f %.0f 0 1 \r\n', BoundaryAxial1, BoundaryCH, BoundaryAxial2);
fprintf(fileID,'ratio 1 %.0f 1 %.0f %.6e \r\n', BoundaryAxial1, BoundaryAxial2, RatioAxial2);
fprintf(fileID,'spacing 1 %.0f %.0f %.0f 1 0 \r\n', BoundaryAxial1, BoundaryCH, BoundaryAxial2);

fprintf(fileID, 'region 2 2 %.0f %.0f 1 %.3e \r\n', BoundaryCH, BoundaryAxial2, DensityCH);
fprintf(fileID, 'material 1 1. 2.014  0.5 \r\n');
fprintf(fileID, 'material 1 6. 12.012 0.5 \r\n');
fprintf(fileID, 'ioniz 1 4 \r\n');
fprintf(fileID, 'eos /work3/clf/rad_hydro/hyades/EOS-Opacity/SESAME//eos_32.dat 1 \r\nc\r\n');


% Gold layer
fprintf(fileID,'mesh %.0f 1 %.0f %.0f 0. %.6e %.6e %.6e \r\n', BoundaryCH, BoundaryGold, BoundaryAxial1, CHLength, AxialLength1, CHLength+GoldLength);
fprintf(fileID,'ratio %.0f 1 %.0f 1 %.6e \r\n', BoundaryCH, BoundaryGold, RatioGold);
fprintf(fileID,'spacing %.0f 1 %.0f %.0f 0 1 \r\n', BoundaryCH, BoundaryGold, BoundaryAxial1);
fprintf(fileID,'ratio %.0f 1 %.0f %.0f %.6e \r\n', BoundaryCH, BoundaryCH, BoundaryAxial1, RatioAxial1);
fprintf(fileID,'spacing %.0f 1 %.0f %.0f 1 0 \r\n', BoundaryCH, BoundaryGold, BoundaryAxial1);

fprintf(fileID,'mesh %.0f %.0f %.0f %.0f %.6e %.6e %.6e %.6e \r\n', BoundaryCH, BoundaryAxial1, BoundaryGold, BoundaryAxial2, AxialLength1, CHLength, AxialLength1+AxialLength2, CHLength+GoldLength);
fprintf(fileID,'ratio %.0f %.0f %.0f %.0f %.6e \r\n', BoundaryCH, BoundaryAxial1, BoundaryGold, BoundaryAxial1, RatioGold);
fprintf(fileID,'spacing %.0f %.0f %.0f %.0f 0 1 \r\n', BoundaryCH, BoundaryAxial1, BoundaryGold, BoundaryAxial2);
fprintf(fileID,'ratio %.0f %.0f %.0f %.0f %.6e \r\n', BoundaryCH, BoundaryAxial1, BoundaryCH, BoundaryAxial2, RatioAxial2);
fprintf(fileID,'spacing %.0f %.0f %.0f %.0f 1 0 \r\n', BoundaryCH, BoundaryAxial1, BoundaryGold, BoundaryAxial2);

fprintf(fileID, 'region %.0f 2 %.0f %.0f 2 %.3e \r\n', (BoundaryCH+1), BoundaryGold, BoundaryAxial2, DensityGold);
fprintf(fileID, 'material 2 au \r\n');
fprintf(fileID, 'ioniz 2 4 \r\n');
fprintf(fileID, 'eos /work3/clf/rad_hydro/hyades/EOS-Opacity/SESAME//eos_51.dat 2 \r\nc\r\n');

%Quartz layer
fprintf(fileID,'mesh %.0f 1 %.0f %.0f 0. %.6e %.6e %.6e \r\n', BoundaryGold, BoundaryQuartz, BoundaryAxial1, CHLength+GoldLength, AxialLength1, CHLength+GoldLength+QuartzLength);
fprintf(fileID,'ratio %.0f 1 %.0f 1 %.6e \r\n', BoundaryGold, BoundaryQuartz, RatioQuartz);
fprintf(fileID,'spacing %.0f 1 %.0f %.0f 0 1 \r\n', BoundaryGold, BoundaryQuartz, BoundaryAxial1);
fprintf(fileID,'ratio %.0f 1 %.0f %.0f %.6e \r\n', BoundaryGold, BoundaryGold, BoundaryAxial1, RatioAxial1);
fprintf(fileID,'spacing %.0f 1 %.0f %.0f 1 0 \r\n', BoundaryGold, BoundaryQuartz, BoundaryAxial1);

fprintf(fileID,'mesh %.0f %.0f %.0f %.0f %.6e %.6e %.6e %.6e \r\n', BoundaryGold, BoundaryAxial1, BoundaryQuartz, BoundaryAxial2, AxialLength1, CHLength+GoldLength, AxialLength1+AxialLength2, CHLength+GoldLength+QuartzLength);
fprintf(fileID,'ratio %.0f %.0f %.0f %.0f %.6e \r\n', BoundaryGold, BoundaryAxial1, BoundaryQuartz, BoundaryAxial1, RatioQuartz);
fprintf(fileID,'spacing %.0f %.0f %.0f %.0f 0 1 \r\n', BoundaryGold, BoundaryAxial1, BoundaryQuartz, BoundaryAxial2);
fprintf(fileID,'ratio %.0f %.0f %.0f %.0f %.6e \r\n', BoundaryGold, BoundaryAxial1, BoundaryGold, BoundaryAxial2, RatioAxial2);
fprintf(fileID,'spacing %.0f %.0f %.0f %.0f 1 0 \r\n', BoundaryGold, BoundaryAxial1, BoundaryQuartz, BoundaryAxial2);

fprintf(fileID, 'region %.0f 2 %.0f %.0f 3 %.3e \r\n', (BoundaryGold+1), BoundaryQuartz, BoundaryAxial2, DensityQuartz);
fprintf(fileID, 'material 3 sio2 \r\n');
fprintf(fileID, 'ioniz 3 4 \r\n');
if QuartzType == 1
    fprintf(fileID, 'eos /work3/clf/rad_hydro/hyades/EOS-Opacity/SESAME//eos_24.dat 3 \r\nc\r\n'); %Silica EOS
else
fprintf(fileID, 'eos /work3/clf/rad_hydro/hyades/EOS-Opacity/SESAME//eos_22.dat 3 \r\nc\r\n'); %Quartz EOS
end

%Foam Layer
fprintf(fileID,'mesh %.0f 1 %.0f %.0f 0. %.6e %.6e %.6e \r\n', BoundaryQuartz, BoundaryFoam, BoundaryAxial1, CHLength+GoldLength+QuartzLength, AxialLength1, CHLength+GoldLength+QuartzLength+FoamLength);
fprintf(fileID,'ratio %.0f 1 %.0f 1 %.6e \r\n', BoundaryQuartz, BoundaryFoam, RatioFoam);
fprintf(fileID,'spacing %.0f 1 %.0f %.0f 0 1 \r\n', BoundaryQuartz, BoundaryFoam, BoundaryAxial1);
fprintf(fileID,'ratio %.0f 1 %.0f %.0f %.6e \r\n', BoundaryQuartz, BoundaryQuartz, BoundaryAxial1, RatioAxial1);
fprintf(fileID,'spacing %.0f 1 %.0f %.0f 1 0 \r\n', BoundaryQuartz, BoundaryFoam, BoundaryAxial1);

fprintf(fileID,'mesh %.0f %.0f %.0f %.0f %.6e %.6e %.6e %.6e \r\n', BoundaryQuartz, BoundaryAxial1, BoundaryFoam, BoundaryAxial2, AxialLength1, CHLength+GoldLength+QuartzLength, AxialLength1+AxialLength2, CHLength+GoldLength+QuartzLength+FoamLength);
fprintf(fileID,'ratio %.0f %.0f %.0f %.0f %.6e \r\n', BoundaryQuartz, BoundaryAxial1, BoundaryFoam, BoundaryAxial1, RatioFoam);
fprintf(fileID,'spacing %.0f %.0f %.0f %.0f 0 1 \r\n', BoundaryQuartz, BoundaryAxial1, BoundaryFoam, BoundaryAxial2);
fprintf(fileID,'ratio %.0f %.0f %.0f %.0f %.6e \r\n', BoundaryQuartz, BoundaryAxial1, BoundaryQuartz, BoundaryAxial2, RatioAxial2);
fprintf(fileID,'spacing %.0f %.0f %.0f %.0f 1 0 \r\n', BoundaryQuartz, BoundaryAxial1, BoundaryFoam, BoundaryAxial2);

fprintf(fileID, 'region %.0f 2 %.0f %.0f 4 %.3e \r\n', (BoundaryQuartz+1), BoundaryFoam, BoundaryAxial2, DensityFoam);
fprintf(fileID, 'material 4 1. 2.014  0.5 \r\n');
fprintf(fileID, 'material 4 6. 12.012 0.5 \r\n');
fprintf(fileID, 'ioniz 4 4 \r\n');
fprintf(fileID, 'eos /work3/clf/rad_hydro/hyades/EOS-Opacity/SESAME//eos_32.dat 4 \r\nc\r\n');

% fprintf(fileID,'ratio 1 1 1 %.0f %.6e \r\n', BoundaryAxial1, RatioAxial1);
% fprintf(fileID,'ratio 1 %.0f 1 %.0f %.6e \r\n', BoundaryAxial1, BoundaryAxial2, RatioAxial2);
% fprintf(fileID,'spacing 1 1 %.0f %.0f 1 0 \r\n', BoundaryFoam, BoundaryAxial1);
% fprintf(fileID,'spacing 1 %.0f %.0f %.0f 1 0 \r\n', BoundaryAxial1, BoundaryFoam, BoundaryAxial2);

fprintf(fileID, 'eosxtrp  1  1  2  1  2 \r\n');
fprintf(fileID, 'eosxtrp  2  1  2  1  2 \r\n');
fprintf(fileID, 'eosxtrp  3  1  2  1  2 \r\n');
fprintf(fileID, 'eosxtrp  4  1  2  1  2 \r\nc\r\n');


fprintf(fileID, 'source laser 0.527 \r\n');
fprintf(fileID, 'raybndl 100 0.2783 175.882 \r\n');
    
fprintf(fileID, 'tv 0. 0. \r\n');
fprintf(fileID, 'tv %.3e %.3e \r\n', RiseTime, Pulse1erg);
fprintf(fileID, 'tv %.3e %.3e \r\n', (LaserEndTime), Pulse1erg);
fprintf(fileID, 'tv %.3e 0. \r\nc\r\n', (LaserEndTime+RiseTime));
 
fprintf(fileID, 'group  0 20 0.03 1.0 \r\n');
fprintf(fileID, 'group 20 50 1.00 5.0 \r\n');
fprintf(fileID, 'group 50 70 5.00 300.0 \r\nc\r\n');

fprintf(fileID, 'pparray rho te ti tr pres ur uz deplas xmass xlsint vol bpeprd bpeprdr bpedep dene eion eelc \r\nc\r\n');
fprintf(fileID, 'parm xlibam 1.0 \r\n');
fprintf(fileID, 'parm flxlem 0.050 \r\n');
fprintf(fileID, 'parm qstimx 2.e-4 \r\n');
fprintf(fileID, 'parm lrdtrn 1 \r\n');
fprintf(fileID, 'parm lrstrt 1 \r\n');
fprintf(fileID, 'parm irdtrn 2 \r\n');
fprintf(fileID, 'parm jlrsgmx 1000 \r\n');
fprintf(fileID, 'parm nstop 1e8 \r\n');
fprintf(fileID, 'parm dt 1e-15  \r\n');
fprintf(fileID, 'parm dtmin 1e-25 \r\n');
fprintf(fileID, 'parm rpltdt 1.0e-13 \r\n');
fprintf(fileID, 'parm postdt 1.0e-13 \r\n');

for i=2:1:(length(ChangeTimeAt)-1)
    fprintf(fileID, 'change %.4e postdt %.4e \n', ChangeTimeAt(i), TimeSpacing(i));
end
fprintf(fileID,'parm postdt %.4e \n', TimeSpacing(1));
fprintf(fileID,'parm tstop %.4e \n', ChangeTimeAt(end));




fclose(fileID);





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

%Given the radii and density, will calculate the mass difference between
%zones and return the mean of the square of this value.
function [MaxDiff, ZoneMass, MassDiff] = MassDifference(Radii, Density)
Vol = (4/3) * pi() * (Radii.^3);
ZoneMass = diff(Vol).*Density;
MassDiff = 100*(diff(ZoneMass)./ZoneMass(2:end));
%Good results for the mean of the quadrature. This function takes the
%mean of quadrature of only those values with a mass difference above
%1.5 (bear in mind the vapour layer is not optimised, so early high
%values are constant).
MaxDiff = mean((abs(MassDiff).*(abs(MassDiff)>2)).^2);
end

%Optimisation function. Given the Ratios, Radii and Layers, will construct
%the mesh and calculate the mass difference using above functions.
function [MaxDiff, Radii, Density] = BoundarySolver(RatioVapour, RatioSplit, RatioIce, RatioCH, LayersVapour, LayersSplit, LayersIce, LayersCH, RadiusVapour,RadiusSplit, RadiusIce, RadiusCH, DensityVapour, DensityIce, DensityCH)

BoundaryVapour = LayersVapour+1;
BoundarySplit = BoundaryVapour + LayersSplit;
BoundaryIce = BoundarySplit + LayersIce;
BoundaryCH = BoundaryIce + LayersCH;

VapourRadii = RatioIncrement(1, BoundaryVapour, 0, RadiusVapour, RatioVapour);
SplitRadii = RatioIncrement(BoundaryVapour, BoundarySplit, RadiusVapour, RadiusSplit, RatioSplit);
IceRadii = RatioIncrement(BoundarySplit, BoundaryIce, RadiusSplit, RadiusIce, RatioIce);
CHRadii = RatioIncrement(BoundaryIce, BoundaryCH, RadiusIce, RadiusCH, RatioCH);

Density = [DensityVapour*ones(1, LayersVapour), DensityIce*ones(1, (LayersIce+LayersSplit)), DensityCH*ones(1, LayersCH)];

Radii = [0, VapourRadii, SplitRadii, IceRadii, CHRadii];
MaxDiff = MassDifference(Radii, Density);

end


