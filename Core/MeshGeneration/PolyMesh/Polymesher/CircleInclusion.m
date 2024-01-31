function [x] = CircleInclusion(Demand,Arg)
  global Dati
  BdBox = Dati.domain;
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
   global Dati
   d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
   d2 = dCircle(P,Dati.circle.Center(1),Dati.circle.Center(2),Dati.circle.Radius);  
   Dist = dDiff(d1,d2);

%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  x = cell(2,1); %No boundary conditions specified for this problem
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%
