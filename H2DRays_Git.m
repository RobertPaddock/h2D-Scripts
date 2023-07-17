Mode = 5
NumberOfRays = 300;
FocusPosition = 0;
FocalSpotRadius = 0.02;
LensRadius = 50; %For mode 3
LensPosition = -100;
LensHeight = 45; %For mode 4
LensLength = 11; %For mode 4
PowerIndex = 1; %1 for linear increase with radius, 2 for squared etc.

%Mode 1 - Uniform Rays
%Produces a laser beam of a given number of rays, all travelling parallel
%to long axis of cylinder (no dephasing)
if Mode==1 
RayRadiusAtFocus = [1:NumberOfRays]*FocalSpotRadius/NumberOfRays;
RayAngle = repmat(180,1,NumberOfRays);
RayPower = RayRadius.^PowerIndex;

%Mode 2 - Random rays
%Produces a beam of constant radius from lens to focus, but with random
%mapping of rays.
elseif Mode==2  
RayRadiusAtFocus = [1:NumberOfRays]*FocalSpotRadius/NumberOfRays;
RayRadiusAtLens = (2*FocalSpotRadius*rand(NumberOfRays,1) - FocalSpotRadius).';
RayAngle = atand((RayRadiusAtFocus - RayRadiusAtLens)./(FocusPosition-LensPosition));
Rays = [RayRadiusAtLens; RayRadiusAtFocus];
% figure 
% plot([LensPosition FocusPosition], Rays)
RayRadius = RayRadiusAtFocus;
RayStartingZ = FocusPosition;
RayAngle = RayAngle+180;
figure
plot([LensPosition FocusPosition FocusPosition-LensPosition], [RayRadius-(tand(RayAngle)*(FocusPosition-LensPosition)); RayRadius; RayRadius+(tand(RayAngle)*(FocusPosition-LensPosition))])
RayPower = RayRadius.^PowerIndex;

%Mode 3 - Random rays, larger lens
%Produces a beam with random mapping, but with a much larger radius at lens
%than at beam
elseif Mode==3
RayRadiusAtFocus = [1:NumberOfRays]*FocalSpotRadius/NumberOfRays;
RayRadiusAtLens = (2*LensRadius*rand(NumberOfRays,1) - LensRadius).';
RayAngle = atand((RayRadiusAtFocus - RayRadiusAtLens)./(FocusPosition-LensPosition));
Rays = [RayRadiusAtLens; RayRadiusAtFocus];
% figure
% plot([LensPosition FocusPosition], Rays)
RayRadius = RayRadiusAtFocus;
RayStartingZ = FocusPosition;
RayAngle = RayAngle+180;
figure
plot([LensPosition FocusPosition FocusPosition-LensPosition], [RayRadius-(tand(RayAngle)*(FocusPosition-LensPosition)); RayRadius; RayRadius+(tand(RayAngle)*(FocusPosition-LensPosition))])
RayPower = RayRadius.^PowerIndex;

%Mode 4 - Beam at 1m distance with a 26 degree angle, and 50cm diameter
elseif Mode == 4
%For each beam, calculate a random series of x and y positions along the
%lens.
NumberOfRaysPerBeam = NumberOfRays/2;
RayRadiusAtFocus = [1:NumberOfRaysPerBeam]*FocalSpotRadius/NumberOfRaysPerBeam;
Lens1Angle = 25;
Lens1XCenter = -108*cosd(Lens1Angle);
Lens1XHeight = 108*sind(Lens1Angle);
PosnAlongLens = (LensLength*rand(NumberOfRaysPerBeam,1) - LensLength/2).';
RayLens1YPosn = Lens1XHeight + (PosnAlongLens * cosd(Lens1Angle));
RayLens1XPosn = Lens1XCenter + (PosnAlongLens * sind(Lens1Angle));
%And for beam 2
Lens2Angle = -25;
Lens2XCenter = -108*cosd(Lens2Angle);
Lens2XHeight = 108*sind(Lens2Angle);
PosnAlongLens = (LensLength*rand(NumberOfRaysPerBeam,1) - LensLength/2).';
RayLens2YPosn = Lens2XHeight + (PosnAlongLens * cosd(Lens2Angle));
RayLens2XPosn = Lens2XCenter + (PosnAlongLens * sind(Lens2Angle));
%Calculate evenly spaced points at focus. Assign each one of the randomly
%generated lens positions. Combine the lens positions into long lists, and
%randomise order to assign each grid point at focus a corresponding lens
%point. Find the angle required for the corresponding ray.
RayRadiusAtFocus = [1:NumberOfRays]*FocalSpotRadius/NumberOfRays;
CombinedRayLensXPosition = [RayLens1XPosn RayLens2XPosn];
CombinedRayLensYPosition = [RayLens1YPosn RayLens2YPosn];
Randomiser = randperm(NumberOfRays);
CombinedRayLensXPosition = CombinedRayLensXPosition(Randomiser);
CombinedRayLensYPosition = CombinedRayLensYPosition(Randomiser);
RayAngles = atand((RayRadiusAtFocus - CombinedRayLensYPosition)./(CombinedRayLensXPosition));
RayAngles = RayAngles+180;
figure
plot([CombinedRayLensXPosition; repmat(FocusPosition,1, NumberOfRays)], [RayRadiusAtFocus-(tand(RayAngles).*(CombinedRayLensXPosition)); RayRadiusAtFocus])
%And save key data
RayRadius = RayRadiusAtFocus;
RayAngle = RayAngles;
RayStartingZ = FocusPosition;
RayPower = RayRadius.^PowerIndex;

%Mode 5 - Beam at 1m distance with a 26 degree angle, and 11cm lens
%diameter, plus a central beam.
elseif Mode == 5
NumberOfRaysPerBeam = ceil(NumberOfRays/3);
NumberOfRays = NumberOfRaysPerBeam*3;
RayRadiusAtFocus = [1:NumberOfRaysPerBeam]*FocalSpotRadius/NumberOfRaysPerBeam;
Lens1Angle = 25;
Lens1XCenter = -108*cosd(Lens1Angle);
Lens1XHeight = 108*sind(Lens1Angle);
PosnAlongLens = (LensLength*rand(NumberOfRaysPerBeam,1) - LensLength/2).';
RayLens1YPosn = Lens1XHeight + (PosnAlongLens * cosd(Lens1Angle));
RayLens1XPosn = Lens1XCenter + (PosnAlongLens * sind(Lens1Angle));
%And for beam 2
Lens2Angle = -25;
Lens2XCenter = -108*cosd(Lens2Angle);
Lens2XHeight = 108*sind(Lens2Angle);
PosnAlongLens = (LensLength*rand(NumberOfRaysPerBeam,1) - LensLength/2).';
RayLens2YPosn = Lens2XHeight + (PosnAlongLens * cosd(Lens2Angle));
RayLens2XPosn = Lens2XCenter + (PosnAlongLens * sind(Lens2Angle));
%And for beam 3
Lens3Angle = 0;
Lens3XCenter = -108*cosd(Lens3Angle);
Lens3XHeight = 108*sind(Lens3Angle);
PosnAlongLens = (LensLength*rand(NumberOfRaysPerBeam,1) - LensLength/2).';
RayLens3YPosn = Lens3XHeight + (PosnAlongLens * cosd(Lens3Angle));
RayLens3XPosn = Lens3XCenter + (PosnAlongLens * sind(Lens3Angle));
%Calculate evenly spaced points at focus. Assign each one of the randomly
%generated lens positions. Combine the lens positions into long lists, and
%randomise order to assign each grid point at focus a corresponding lens
%point. Find the angle required for the corresponding ray.
RayRadiusAtFocus = [1:NumberOfRays]*FocalSpotRadius/NumberOfRays;
CombinedRayLensXPosition = [RayLens1XPosn RayLens2XPosn RayLens3XPosn];
CombinedRayLensYPosition = [RayLens1YPosn RayLens2YPosn RayLens3YPosn];
Randomiser = randperm(NumberOfRays);
CombinedRayLensXPosition = CombinedRayLensXPosition(Randomiser);
CombinedRayLensYPosition = CombinedRayLensYPosition(Randomiser);
RayAngles = atand((RayRadiusAtFocus - CombinedRayLensYPosition)./(CombinedRayLensXPosition));
RayAngles = RayAngles+180;
figure
plot([CombinedRayLensXPosition; repmat(FocusPosition,1, NumberOfRays)], [RayRadiusAtFocus-(tand(RayAngles).*(CombinedRayLensXPosition)); RayRadiusAtFocus])
%And save key data
RayRadius = RayRadiusAtFocus;
RayAngle = RayAngles;
RayStartingZ = FocusPosition;
RayPower = RayRadius.^PowerIndex;

end


filename = sprintf('Rays');
filepath = sprintf(['MultiFile/', filename, '.inf']);
fileID = fopen(filepath,'w');
fprintf(fileID, 'source laser 0.527 \r\n')
  for i=1:NumberOfRays
    %Ray StartingRadius StartingZ StartingAngle
  fprintf(fileID,'RAY %.6e %.6e %.6e %.6e \r\n', RayRadius(i), RayStartingZ, RayAngle(i), RayPower(i));
end
fclose(fileID);