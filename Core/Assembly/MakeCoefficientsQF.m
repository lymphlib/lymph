%> @file  MakeCoefficientsQF.m
%> @author Mattia Corti
%> @date 16 September 2025
%> @brief Assembly of the matrices
%>
%==========================================================================
%> @section classMakeCoefficientsQF Class description
%==========================================================================
%> @brief            Construction of the coefficients of the monomial
%> expansions of the bases for the quadrature-free implementations
%>
%> @param deg        Polynomial degree
%> @param nbases     Number of basis functions
%>
%> @retval Coeff     Struct containing the coefficients expansions  
%>
%==========================================================================

function [Coeff] = MakeCoefficientsQF(deg, nbases)

    %% Construction of the basis functions
    [Lx, Ly, dLx, dLy] = Evalshape2DCoeff(deg);
    
    %% Initialization of Legendre polynomials products in 1D for bilinear forms
    LxLx       = zeros(2*(deg+1)-1,nbases^2);
    LyLy       = zeros(2*(deg+1)-1,nbases^2);
    dLxdLx     = zeros(2*(deg+1)-1,nbases^2);
    dLydLy     = zeros(2*(deg+1)-1,nbases^2);
    LxdLx      = zeros(2*(deg+1)-1,nbases^2);
    LydLy      = zeros(2*(deg+1)-1,nbases^2);
    dLxLx      = zeros(2*(deg+1)-1,nbases^2);
    dLyLy      = zeros(2*(deg+1)-1,nbases^2);
        
    % Legendre polynomials products in 1D repeated for double combinations
    Lx2 = repmat(Lx,1,nbases);
    Ly2 = repmat(Ly,1,nbases);
    dLx2 = repmat(dLx,1,nbases);
    dLy2 = repmat(dLy,1,nbases);
    
    Lx3 = reshape(repmat(Lx,nbases,1),[deg+1,nbases^2]);
    Ly3 = reshape(repmat(Ly,nbases,1),[deg+1,nbases^2]);
    dLx3 = reshape(repmat(dLx,nbases,1),[deg+1,nbases^2]);
    dLy3 = reshape(repmat(dLy,nbases,1),[deg+1,nbases^2]);

    %% Cycle over the polynomial degrees
    for ii = 1:deg+1    

            % Construction of the index values for the coefficients
            ii_index = (ii:ii+deg); 

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

    %% Creation of the outuput structures for bilinear forms
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
    Coeff.phiphiC     = sparse(reshape(phiphiC,[(2*(deg+1)-1)^2, nbases^2])');
    Coeff.gradxgradxC = sparse(reshape(gradxgradxC,[(2*(deg+1)-1)^2, nbases^2])');
    Coeff.gradygradyC = sparse(reshape(gradygradyC,[(2*(deg+1)-1)^2, nbases^2])');
    Coeff.gradxgradyC = sparse(reshape(gradxgradyC,[(2*(deg+1)-1)^2, nbases^2])');
    Coeff.gradygradxC = sparse(reshape(gradygradxC,[(2*(deg+1)-1)^2, nbases^2])');
    Coeff.gradxphiC   = sparse(reshape(gradxphiC,[(2*(deg+1)-1)^2, nbases^2])');
    Coeff.gradyphiC   = sparse(reshape(gradyphiC,[(2*(deg+1)-1)^2, nbases^2])');
    Coeff.phigradxC   = sparse(reshape(phigradxC,[(2*(deg+1)-1)^2, nbases^2])');
    Coeff.phigradyC   = sparse(reshape(phigradyC,[(2*(deg+1)-1)^2, nbases^2])');

    %% Initialization of Legendre polynomials products in 1D for trilinear forms
    LxLx       = zeros(2*deg+1,nbases^3);
    LyLy       = zeros(2*deg+1,nbases^3);
    LxLxLx     = zeros(3*deg+1,nbases^3);
    LyLyLy     = zeros(3*deg+1,nbases^3);
        
    % Legendre polynomials products in 1D repeated for triple combinations
    Lx2 = repmat(Lx,1,nbases^2);
    Ly2 = repmat(Ly,1,nbases^2);

    Lx3 = reshape(repmat(Lx,nbases^2,1),[deg+1,nbases^3]);
    Ly3 = reshape(repmat(Ly,nbases^2,1),[deg+1,nbases^3]);

    Lx4 = repmat(reshape(repmat(Lx,nbases,1),[deg+1,nbases^2]),1,nbases);
    Ly4 = repmat(reshape(repmat(Ly,nbases,1),[deg+1,nbases^2]),1,nbases);

    %% 1st-cycle over the polynomial degrees
    for ii = 1:deg+1    

            % Construction of the index values for the coefficients
            ii_index = (ii:ii+deg); 

            LxLx(ii_index,:)    = LxLx(ii_index,:)   + Lx3(ii,:).*Lx2;
            LyLy(ii_index,:)    = LyLy(ii_index,:)   + Ly3(ii,:).*Ly2;

    end

    %% 2nd-cycle over the polynomial degrees
    for ii = 1:deg+1    

            % Construction of the index values for the coefficients
            ii_index = (ii:ii+2*deg); 

            LxLxLx(ii_index,:)    = LxLxLx(ii_index,:)   + Lx4(ii,:).*LxLx;
            LyLyLy(ii_index,:)    = LyLyLy(ii_index,:)   + Ly4(ii,:).*LyLy;

    end
    
    % Reshape structures for the output creation
    LxLxLx  = reshape(LxLxLx, [size(LxLxLx,1) 1 size(LxLxLx,2)]);
    LyLyLy  = reshape(LyLyLy, [1 size(LyLyLy,1) size(LyLyLy,2)]);
  
    %% Creation of the outuput structures for trilinear forms
    phiphiphiC = pagemtimes(LxLxLx,LyLyLy);

    % Reshape of the final structures for optimization purposes
    Coeff.phiphiphiC     = sparse(reshape(phiphiphiC,[(3*deg+1)^2, nbases^3])');

end
