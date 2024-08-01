%> @file  Rectangle.m
%> @author Ilario Mazzieri
%> @date 8 March 2023
%> @brief Generation of rectangle based on PolyMesher - version 1.1.
%>
%> The original version is taken from http://paulino.princeton.edu/software.html
%>------------------ PolyMesher  version: 1.1 (Aug13) ---------------------%
%> Ref1: C Talischi, GH Paulino, A Pereira, IFM Menezes,                   %
%>      "PolyMesher: A general-purpose mesh generator for polygonal        %
%>      elements written in Matlab", Struct Multidisc Optim, 2012,         %
%>      DOI 10.1007/s00158-011-0706-z                                      %
%>                                                                         %
%> Ref2: A Pereira, C Talischi, GH Paulino, IFM Menezes, MS Carvalho,      %
%>      "Implementation of fluid flow topology optimization in PolyTop",   %
%>      Struct Multidisc Optim, 2013, DOI XX.XXXX/XXXXXX-XXX-XXX-X         %
%>-------------------------------------------------------------------------%

function [x] = Rectangle(Demand,Arg)  
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
   d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
   Dist = d1;

%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  x = cell(2,1); %No boundary conditions specified for this problem
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%