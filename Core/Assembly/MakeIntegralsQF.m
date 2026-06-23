%> @file  MakeIntegralsQF.m
%> @author Mattia Corti
%> @date 16 September 2025
%> @brief Assembly of the matrices
%>
%==========================================================================
%> @section classMakeIntegralsQF Class description
%==========================================================================
%> @brief            Construction of the integral quantities for matrices
%> assembly with quadrature-free implementation on a specific polygon
%>
%> @param Coeff      Coefficients of the monomial expansion of the bases
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param ie         Polygonal element
%>
%> @retval Integral  Struct containing the integral quantities  
%>
%==========================================================================

function [Integral] = MakeIntegralsQF(Coeff, femregion, ie)
            
        % Construction of the polygon in reference coordinates
        [coords_ie, Jdet, BJ, ~] = GetJacobianPolygon(femregion.bbox(ie,:), femregion.coords_element{ie}');

        % Computation of the monomial integrals in the polygon v
        [I, ~] = IntegralOverPolygon(1, 2*femregion.degree(ie), 2*femregion.degree(ie), coords_ie);

        % Computation of the bases integrals of bilinear forms
        Integral.phiphiC     = Jdet*Coeff.phiphiC*I;
        Integral.gradxgradxC = Jdet/(BJ(1,1))^2*Coeff.gradxgradxC*I;
        Integral.gradygradyC = Jdet/(BJ(2,2))^2*Coeff.gradygradyC*I;
        Integral.gradxgradyC = Jdet/(BJ(1,1)*BJ(2,2))*Coeff.gradxgradyC*I;
        Integral.gradygradxC = Jdet/(BJ(1,1)*BJ(2,2))*Coeff.gradygradxC*I;
        Integral.gradxphiC   = Jdet/BJ(1,1)*Coeff.gradxphiC*I;
        Integral.gradyphiC   = Jdet/BJ(2,2)*Coeff.gradyphiC*I;
        Integral.phigradxC   = Jdet/BJ(1,1)*Coeff.phigradxC*I;
        Integral.phigradyC   = Jdet/BJ(2,2)*Coeff.phigradyC*I;

        % Computation of the monomial integrals in the polygon v
        [I3, ~] = IntegralOverPolygon(1, 3*femregion.degree(ie), 3*femregion.degree(ie), coords_ie);

        % Computation of the bases integrals of trilinear forms
        Integral.phiphiphiC  = Jdet*Coeff.phiphiphiC*I3;

end
