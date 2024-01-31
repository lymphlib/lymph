function [x] = DoubleCircle(Demand,Arg)
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
  d1 = dCircle(P,Dati.Circle1.Center(1),Dati.Circle1.Center(2),Dati.Circle1.Radius);  
  d2 = dCircle(P,Dati.Circle2.Center(1),Dati.Circle2.Center(2),Dati.Circle2.Radius);  
%   d1 = dCircle(P,-0.5,0,0.75);  
%   d2 = dCircle(P,0.5,0,0.75);  
  d3 = dLine(P,0,-1.5,0,1.5);
  d4 = dIntersect(d1,d3); 
  d5 = dIntersect(d2,-d3); 
  Dist = dUnion(d4,d5);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  x = cell(2,1); %No boundary conditions specified for this problem
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%