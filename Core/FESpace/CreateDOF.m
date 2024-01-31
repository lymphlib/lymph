%> @file  CreateDOF.m
%> @author Mattia Corti, Paola F. Antonietti, Ilario Mazzieri
%> @date 26 March 2023 
%> @brief Construction of the the femregion structure.
%> 
%> The function creates the femregion stucture, which is at the basis of
%> the code. Indeed the structure contains all the information about the
%> PolyDG discretization. 
%> 
%==========================================================================
%> @section classCreateDOF Class description
%==========================================================================
%> @brief  Construction of the the femregion structure.
%>
%> @param Data           Structure containing the information about the
%> discretization and the physical parameters of the problem.
%> @param region         Structure containing the mesh information.
%>
%> @retval femregion     Structure containing all the information 
%> about the finite element approximation. We report the specific
%> structure:
%> -# femregion.degree: Polynomial degree of the PolyDG approximation.
%> -# femregion.nedges: Number of edges for each cell of the mesh.
%> -# femregion.nbases: Number of basis for each element of the mesh.
%> -# femregion.nel: Number of elements of the mesh.
%> -# femregion.nel_p: Number of porous elements of the mesh.
%> -# femregion.nel_a: Number of acoustic elements of the mesh.
%> -# femregion.nel_e: Number of elastic elements of the mesh.
%> -# femregion.nel_f: Number of fluid elements of the mesh.
%> -# femregion.ndof: Total number of degrees of freedom.
%> -# femregion.ndof_p: Total number of porous degrees of freedom.
%> -# femregion.ndof_a: Total number of acoustic degrees of freedom (scalar).
%> -# femregion.ndof_e: Total number of elastic degrees of freedom (scalar).
%> -# femregion.ndof_f: Total number of fluid degrees of freedom (scalar).
%> -# femregion.nqn: Number of quadrature nodes for each mesh element.
%> -# femregion.coord: Coordinates of the mesh element vertices.
%> -# femregion.bbox: Coordinates of the bounding box for eache mesh element.
%> -# femregion.area: Computed area of each mesh element.
%> -# femregion.coords_element: Coordinates of the vertices for each mesh element.
%> -# femregion.max_kb: Max kb parameter for each mesh element \cite cangiani2017.
%> -# femregion.connectivity: Connectivity of the edges for each mesh element. 
%> -# femregion.connectivity_bc: Connectivity of the boundary edges for each mesh element. 
%> -# femregion.bc_tag: Tag for boundary edges. 
%> -# femregion.id: id for mesh elements (integer).
%> -# femregion.tag: Tag for mesh elements (character). 
%==========================================================================

function [femregion] = CreateDOF(Data, region)

    % Computation of the number of required basis
    nln = 0.5.*(Data.degree+1).*(Data.degree+2);
    
    % Computation of the number of quadrature nodes
    nqn = 2*Data.degree + 1;

    % Struct construction
    femregion = struct('degree',    Data.degree, ...
        'nedges',                   region.nedges', ...
        'nbases',                   nln, ...
        'nel',                      region.ne,...
        'nel_p',                    (size(find(region.tag=='P'),1)),...
        'nel_a',                    (size(find(region.tag=='A'),1)),...
        'nel_e',                    (size(find(region.tag=='E'),1)),...
        'nel_f',                    (size(find(region.tag=='F'),1)),...
        'ndof',                     nln*region.ne,...
        'ndof_p',                   nln*(size(find(region.tag=='P'),1)),...
        'ndof_a',                   nln*(size(find(region.tag=='A'),1)),...
        'ndof_e',                   nln*(size(find(region.tag=='E'),1)),...
        'ndof_f',                   nln*(size(find(region.tag=='F'),1)),...
        'nqn',                      nqn,...
        'coord',                    region.coord,...
        'bbox',                     region.BBox, ...
        'area',                     region.area, ...
        'id',                       region.id, ...
        'tag',                      region.tag, ...
        'bc_tag',                   region.bc_tag, ...
        'connectivity_bc',          region.connectivity_bc);

    % Assignment of cell elements to the struct
    femregion.coords_element = region.coords_element;
    femregion.max_kb         = region.max_kb;
    femregion.connectivity   = region.connectivity;
    
    

end



