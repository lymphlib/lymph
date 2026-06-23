%> @file  CreateDOF.m
%> @author Mattia Corti, Paola F. Antonietti, Ilario Mazzieri, Caterina
%> Leimer Saglio
%> @date 29 May 2026
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
%> -# femregion.degree: Vector of polynomial degrees of the PolyDG 
%> approximation for each element.
%> -# femregion.nedges: Number of edges for each cell of the mesh.
%> -# femregion.nbases: Vector of number of basis for each element of the mesh.
%> -# femregion.nel: Number of elements of the mesh.
%> -# femregion.nel_phys: vector of number of elements for each physics
%> -# femregion.ndof: Total number of degrees of freedom.
%> -# femregion.ndof_phys: vector of number of dofs for each physics.
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
%> -# femregion.id_phys: rescaled element IDs for each physics so that their local numbering starts from 1 (shifted version of the global volume element IDs in femregion.id).
%> -# femregion.label: Tag for mesh elements (character). 
%> -# femregion.phys: vector of physics labels (character).

%==========================================================================

function [femregion] = CreateDOF(Data, region)

    % Computation of degree vector
    if isscalar(Data.degree)
        degree_vector = Data.degree*ones(region.ne,1);
    else
        degree_vector = Data.degree;
    end

    % Computation of the number of required basis
    nln = 0.5.*(degree_vector+1).*(degree_vector+2);
   
    % Computation of the number of quadrature nodes
    nqn = max(degree_vector) + 1;

    % Struct construction
    femregion = struct( ...
        'degree',                   degree_vector, ...
        'nedges',                   region.nedges', ...
        'nbases',                   nln, ...
        'nel',                      region.ne, ...
        'ndof',                     sum(nln), ...
        'nqn',                      nqn, ...
        'coord',                    region.coord, ...
        'bbox',                     region.BBox, ...
        'area',                     region.area, ...
        'id',                       region.id, ...
        'label',                    region.label, ...
        'bc_tag',                   region.bc_tag, ...
        'connectivity_bc',          region.connectivity_bc);

    % Multiphysics: number of elements and DOFs for each physics
    femregion.nel_phys  = zeros(length(Data.LabEl),1);
    femregion.ndof_phys = zeros(length(Data.LabEl),1);
    femregion.id_phys   = femregion.id;

    for ii = 1:length(Data.LabEl)
	% Volume label associated with the specific physics
        femregion.phys(ii,1)    = Data.LabEl{ii};
	% Number of elements associated with the specific physics
        femregion.nel_phys(ii)  = length(find(region.label==Data.LabEl{ii}));
        % Number of dofs associated with the specific physics
	femregion.ndof_phys(ii) = sum(nln(region.label==Data.LabEl{ii}));
        % Rescale element IDs for physics ii so that their local numbering starts from 1 (shifted version of the global volume element IDs in femregion.id).
	femregion.id_phys(region.label==Data.LabEl{ii}) = femregion.id_phys(region.label==Data.LabEl{ii})-min(Data.TagEl{ii})+1;
    end

    % Assignment of cell elements to the struct
    femregion.coords_element = region.coords_element;
    femregion.max_kb         = region.max_kb;
    femregion.connectivity   = region.connectivity;
    
end



