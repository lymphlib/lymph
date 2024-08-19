%> @file  SaveSolutionPoroAcuEla.m
%> @author Ilario Mazzieri
%> @date 28 June 2024
%> @brief  Save solution vector in a struct
%>
%==========================================================================
%> @section classSaveSolutionPoroAcuEla Class description
%> @brief  Save solution vector in a struct
%
%> @param Uh     Solution vector
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval Solutions  Struct with solution vector
%>
%==========================================================================

function  [Solutions]  = SaveSolutionPoroAcuEla(Uh,femregion)

np = femregion.ndof_p;
ne = femregion.ndof_e;
na = femregion.ndof_a;

up_h  = [];
wp_h  = [];
phi_h = [];
ue_h  = [];
dot_up_h  = [];
dot_wp_h  = [];
dot_phi_h = [];
dot_ue_h  = [];

    % poro - acoustic - elastic
if (np > 0 && na > 0 && ne > 0)

    up_h   = Uh(1         : 2*np);
    wp_h   = Uh(2*np+1    : 4*np);
    phi_h  = Uh(4*np+1    : 4*np+na);
    ue_h   = Uh(4*np+na+1 : 4*np+na+2*ne);

    ndof_vel = 4*np+na+2*ne;

    dot_up_h   = Uh(ndof_vel + 1        : ndof_vel+2*np);
    dot_wp_h   = Uh(ndof_vel+2*np+1     : ndof_vel+4*np);
    dot_phi_h  = Uh(ndof_vel+4*np+1     : ndof_vel+4*np+na);
    dot_ue_h   = Uh(ndof_vel+4*np+na+1  : ndof_vel+4*np+na+2*ne);


    % poro - acoustic
elseif (np > 0 && na > 0 && ne == 0)

    up_h   = Uh(1         : 2*np);
    wp_h   = Uh(2*np+1    : 4*np);
    phi_h  = Uh(4*np+1    : 4*np+na);

    ndof_vel = 4*np+na;

    dot_up_h   = Uh(ndof_vel + 1        : ndof_vel+2*np);
    dot_wp_h   = Uh(ndof_vel+2*np+1     : ndof_vel+4*np);
    dot_phi_h  = Uh(ndof_vel+4*np+1     : ndof_vel+4*np+na);

    % poro  - elastic
elseif (np > 0 && na == 0 && ne > 0)

    up_h   = Uh(1         : 2*np);
    wp_h   = Uh(2*np+1    : 4*np);
    ue_h   = Uh(4*np+1    : 4*np+2*ne);

    ndof_vel = 4*np+2*ne;

    dot_up_h   = Uh(ndof_vel + 1        : ndof_vel+2*np);
    dot_wp_h   = Uh(ndof_vel+2*np+1     : ndof_vel+4*np);
    dot_ue_h   = Uh(ndof_vel+4*np+1     : ndof_vel+4*np+2*ne);



    % acoustic - elastic
elseif (np == 0 && na > 0 && ne > 0)

    phi_h  = Uh(1    : na);
    ue_h   = Uh(na+1 : na+2*ne);

    ndof_vel = na+2*ne;

    dot_phi_h  = Uh(ndof_vel+1     : ndof_vel+na);
    dot_ue_h   = Uh(ndof_vel+na+1  : ndof_vel+na+2*ne);


    % poro
elseif (np > 0 && na == 0 && ne == 0)

    up_h   = Uh(1         : 2*np);
    wp_h   = Uh(2*np+1    : 4*np);

    ndof_vel = 4*np;

    dot_up_h   = Uh(ndof_vel + 1        : ndof_vel+2*np);
    dot_wp_h   = Uh(ndof_vel+2*np+1     : ndof_vel+4*np);

    % acoustic
elseif (np == 0 && na > 0 && ne == 0)

    phi_h  = Uh(1    : na);

    ndof_vel = na;

    dot_phi_h  = Uh(ndof_vel+1     : ndof_vel+na);


    % elastic
elseif (np == 0 && na == 0 && ne > 0)

    ue_h   = Uh(1 : 2*ne);

    ndof_vel = 2*ne;

    dot_ue_h   = Uh(ndof_vel+1  : ndof_vel+2*ne);


end


Solutions.up_h      = up_h;
Solutions.wp_h      = wp_h;
Solutions.phi_h     = phi_h;
Solutions.ue_h      = ue_h;
Solutions.dot_up_h  = dot_up_h;
Solutions.dot_wp_h  = dot_wp_h;
Solutions.dot_phi_h = dot_phi_h;
Solutions.dot_ue_h  = dot_ue_h;

