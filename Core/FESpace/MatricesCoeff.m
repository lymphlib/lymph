%> @file  MatricesCoeff.m
%> @author Mattia Corti
%> @date 11 November 2023 
%> @brief The function constructs the coefficients of the basis terms for
%> matrices construction.
%> 
%> The function constructs the 2D-monomial coefficients for the terms in 
%> the matrices.
%>
%======================================================================
%> @section classMatricesCoeff Class description
%======================================================================
%> @brief The function constructs the coefficients of the basis terms for
%> matrices construction.
%>
%> @param femregion   Finite Element struct (see CreateDOF.m). 
%> @param Lx          Coefficients of the Legendre polynomials in x-direction.
%> @param Ly          Coefficients of the Legendre polynomials in y-direction.
%> @param dLx         Coefficients of the gradient of Legendre polynomials in x-direction.
%> @param dLy         Coefficients of the gradient of Legendre polynomials in y-direction.
%>
%> @retval Coeff      Coefficients of the terms of the matrices.
%======================================================================

function [Coeff] = MatricesCoeff(femregion, Lx, Ly, dLx, dLy)

    % Initialization of Legendre polynomials products in 1D
    LxLx       = zeros(2*(femregion.degree+1)-1,femregion.nbases^2);
    LyLy       = zeros(2*(femregion.degree+1)-1,femregion.nbases^2);
    dLxdLx     = zeros(2*(femregion.degree+1)-1,femregion.nbases^2);
    dLydLy     = zeros(2*(femregion.degree+1)-1,femregion.nbases^2);
    LxdLx      = zeros(2*(femregion.degree+1)-1,femregion.nbases^2);
    LydLy      = zeros(2*(femregion.degree+1)-1,femregion.nbases^2);
    dLxLx      = zeros(2*(femregion.degree+1)-1,femregion.nbases^2);
    dLyLy      = zeros(2*(femregion.degree+1)-1,femregion.nbases^2);
        
    % Legendre polynomials products in 1D repeated for double combinations
    Lx2 = repmat(Lx,1,femregion.nbases);
    Ly2 = repmat(Ly,1,femregion.nbases);
    dLx2 = repmat(dLx,1,femregion.nbases);
    dLy2 = repmat(dLy,1,femregion.nbases);
    
    Lx3 = reshape(repmat(Lx,femregion.nbases,1),[femregion.degree+1,femregion.nbases^2]);
    Ly3 = reshape(repmat(Ly,femregion.nbases,1),[femregion.degree+1,femregion.nbases^2]);
    dLx3 = reshape(repmat(dLx,femregion.nbases,1),[femregion.degree+1,femregion.nbases^2]);
    dLy3 = reshape(repmat(dLy,femregion.nbases,1),[femregion.degree+1,femregion.nbases^2]);

    %% Cycle over the polynomial degrees
    for ii = 1:femregion.degree+1    

            % Construction of the index values for the coefficients
            ii_index = (ii:ii+femregion.degree); 

            LxLx(ii_index,:)    = LxLx(ii_index,:)   + Lx3(ii,:).*Lx2;
            LyLy(ii_index,:)    = LyLy(ii_index,:)   + Ly3(ii,:).*Ly2;

            dLxdLx(ii_index,:)  = dLxdLx(ii_index,:) + dLx3(ii,:).*dLx2;
            dLydLy(ii_index,:)  = dLydLy(ii_index,:) + dLy3(ii,:).*dLy2;

            LxdLx(ii_index,:)   = LxdLx(ii_index,:)  + Lx3(ii,:).*dLx2;
            LydLy(ii_index,:)   = LydLy(ii_index,:)  + Ly3(ii,:).*dLy2;
            dLxLx(ii_index,:)   = dLxLx(ii_index,:)  + dLx3(ii,:).*Lx2;
            dLyLy(ii_index,:)   = dLyLy(ii_index,:)  + dLy3(ii,:).*Ly2;
    end

    % Reshape structures for the output creation
    LxLx   = reshape(LxLx,   [size(LxLx,1) 1 size(LxLx,2)]);
    LyLy   = reshape(LyLy,   [1 size(LyLy,1) size(LyLy,2)]);
    dLxdLx = reshape(dLxdLx, [size(dLxdLx,1) 1 size(dLxdLx,2)]);
    dLydLy = reshape(dLydLy, [1 size(dLydLy,1) size(dLydLy,2)]);
    dLxLx  = reshape(dLxLx,  [size(dLxLx,1) 1 size(dLxLx,2)]);
    dLyLy  = reshape(dLyLy,  [1 size(dLyLy,1) size(dLyLy,2)]);
    LxdLx  = reshape(LxdLx,  [size(LxdLx,1) 1 size(LxdLx,2)]);
    LydLy  = reshape(LydLy,  [1 size(LydLy,1) size(LydLy,2)]);

    %% Creation of the outuput structures
    phiphiC     = pagemtimes(LxLx,LyLy);
    
    gradxgradxC = pagemtimes(dLxdLx, LyLy);
    gradygradyC = pagemtimes(LxLx, dLydLy);
    gradxgradyC = pagemtimes(dLxLx, LydLy);
    gradygradxC = pagemtimes(LxdLx, dLyLy);

    gradxphiC   = pagemtimes(dLxLx, LyLy);
    gradyphiC   = pagemtimes(LxLx, dLyLy);
    
    phigradxC   = pagemtimes(LxdLx, LyLy);
    phigradyC   = pagemtimes(LxLx, LydLy);

    % Reshape of the final structures for optimization purposes
    Coeff.phiphiC     = reshape(phiphiC,[(2*(femregion.degree+1)-1)^2, femregion.nbases^2])';
    Coeff.gradxgradxC = reshape(gradxgradxC,[(2*(femregion.degree+1)-1)^2, femregion.nbases^2])';
    Coeff.gradygradyC = reshape(gradygradyC,[(2*(femregion.degree+1)-1)^2, femregion.nbases^2])';
    Coeff.gradxgradyC = reshape(gradxgradyC,[(2*(femregion.degree+1)-1)^2, femregion.nbases^2])';
    Coeff.gradygradxC = reshape(gradygradxC,[(2*(femregion.degree+1)-1)^2, femregion.nbases^2])';
    Coeff.gradxphiC   = reshape(gradxphiC,[(2*(femregion.degree+1)-1)^2, femregion.nbases^2])';
    Coeff.gradyphiC   = reshape(gradyphiC,[(2*(femregion.degree+1)-1)^2, femregion.nbases^2])';
    Coeff.phigradxC   = reshape(phigradxC,[(2*(femregion.degree+1)-1)^2, femregion.nbases^2])';
    Coeff.phigradyC   = reshape(phigradyC,[(2*(femregion.degree+1)-1)^2, femregion.nbases^2])';

end
