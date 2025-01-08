clear all;
close all;
clc;
addpath('Pre\','Post\','Solver\');
% Set up the problem (para & boundary etc.)
E = 1e9;
nu = 0.3;

problem_type = 0; % 0 for plane stress / 1 for plane strain

n_int_xi  = 3; % quadrature points
n_int_eta = 3;
n_int_l  = 3; % quadrature points for liner edge

analy_u = @(x,y) [-x^2; % analytical u solution
                  0;];  % analytical v solution

%f = @(x, y) [2*E/(nu^2-1);  % fx
%             0;];           % fy
f = @(x, y) [0;     % fx
             0;];   % fy

% Init
D  = [1,  nu, 0;
      nu, 1,  0;
      0,  0,  (1-nu)/2;];
D  = E/(1 - nu^2) * D;
if problem_type == 1
    D  = [1-nu, nu,   0;
          nu,   1-nu, 0;
          0,    0,    (1-2*nu)/2;];
    D  = E/(1 + nu)/(1 - 2*nu) * D;
end

n_int = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);
[xl, weight_l] = Gauss(n_int_l, -1, 1);


% Read the mesh (IEN)

    % Complicated to read .msh
    % Replace it with manually generated mesh instead
    % mesh generation
    n_en   = 4;    % number of nodes in an element
    n_el_x = 60;               % number of elements in x-dir
    n_el_y = 60;               % number of elements in y-dir
    n_el   = n_el_x * n_el_y;  % total number of elements
    
    n_np_x = n_el_x + 1;      % number of nodal points in x-dir
    n_np_y = n_el_y + 1;      % number of nodal points in y-dir
    n_np   = n_np_x * n_np_y; % total number of nodal points
    
    x_coor = zeros(n_np, 1);
    y_coor = x_coor;
    
    hx = 1.0 / n_el_x;        % mesh size in x-dir
    hy = 1.0 / n_el_y;        % mesh size in y-dir
    
    % generate the nodal coordinates
    for ny = 1 : n_np_y
        for nx = 1 : n_np_x
            index = (ny-1)*n_np_x + nx; % nodal index
            x_coor(index) = (nx-1) * hx;
            y_coor(index) = (ny-1) * hy;
        end
    end
    
    % IEN array
    IEN = zeros(n_el, n_en);
    for ex = 1 : n_el_x
        for ey = 1 : n_el_y
            ee = (ey-1) * n_el_x + ex; % element index
            IEN(ee, 1) = (ey-1) * n_np_x + ex;
            IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
            IEN(ee, 3) =  ey    * n_np_x + ex + 1;
            IEN(ee, 4) =  ey    * n_np_x + ex;
        end
    end

% Set the boundary (ID)

    % Complicated to read .msh
    % Replace it with manually generated mesh instead
    % ID array
    ID = zeros(n_np, 2);
    g = zeros(n_np, 2);
    h = zeros(n_np, 2);
    node_type = zeros(n_np, 1); % 0 for interior nodes / 1 for Dirichlet / 2 for Neumann

    counter = 0;
    for ny = 2 : n_np_y - 1
        for nx = 2 : n_np_x - 1
            index = (ny-1)*n_np_x + nx;
            counter = counter + 2;
            ID(index,1) = counter - 1;
            ID(index,2) = counter;
        end
    end
    for nx = 2 : n_np_x - 1
        %counter = counter + 2;
        %ID(nx,1) = counter - 1;
        %ID((n_np_y-1)*n_np_x + nx,1) = counter;

        %g(nx, 1) = x_coor(nx) * (x_coor(nx) - 1);
        %g((n_np_y-1)*n_np_x + nx, 1) = -x_coor((n_np_y-1)*n_np_x + nx) * (x_coor(nx) - 1);

        node_type(nx) = 2;
        node_type((n_np_y-1)*n_np_x + nx) = 2;
        h(nx, 2) = 10;
    end
    n_eq = counter;
    
    LM = zeros(2, n_el, n_en);
    for i = 1 : n_el
        for j  = 1 : n_en
            LM(1, i, j) = ID(IEN(i, j), 1);
            LM(2, i, j) = ID(IEN(i, j), 2);
        end
    end


% Assemble K & F
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );
  
    k_ele = zeros(n_en*2, n_en*2); % element stiffness matrix
    f_ele = zeros(n_en*2, 1);    % element load vector
  
    for ll = 1 : n_int
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end
    
        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        
        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
          
            Ba = [Na_x, 0, Na_y; 0, Na_y, Na_x;];
            for ii = 1 : 2
                pp = 2*(aa-1)+ii;
                tmp = f(x_l, y_l);
                f_ele(pp) = f_ele(pp) + weight(ll) * tmp(ii) * Na * detJ;
    
                for bb = 1 : n_en
                    Nb = Quad(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                    Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
                
                    Bb = [Nb_x, 0; 0, Nb_y; Nb_y, Nb_x;];
                
                    ed = [1,0; 0,1;];
                    for jj = 1 : 2
                        qq = 2*(bb-1)+jj;
                        k_ele(pp, qq) = k_ele(pp, qq) + weight(ll) * ed(ii, :) * Ba * D * Bb * ed(:, jj) * detJ;
                    end
             
                end % end of bb loop
            end
        end % end of aa loop
        
    end % end of quadrature loop

    for ii = 1 : n_en
        if (node_type(IEN(ee, ii)) == 2) && (node_type(IEN(ee, mod(ii, n_en)+1)) == 2)
            switch ii
                case 1
                    xi_l = xl;
                    eta_l = ones(1, n_int_l);
                case 2
                    xi_l = -1 * ones(1, n_int_l);
                    eta_l = xl;
                case 3
                    xi_l = xl;
                    eta_l = -1 * ones(1, n_int_l);
                case 4
                    xi_l = ones(1, n_int_l);
                    eta_l = xl;
            end
            for ll = 1 : n_int_l
                x_l = 0.0; y_l = 0.0;
                dx_dxi = 0.0; dx_deta = 0.0;
                dy_dxi = 0.0; dy_deta = 0.0;
                for aa = 1 : n_en
                    x_l = x_l + x_ele(aa) * Quad(aa, xi_l(ll), eta_l(ll));
                    y_l = y_l + y_ele(aa) * Quad(aa, xi_l(ll), eta_l(ll));    
                    [Na_xi, Na_eta] = Quad_grad(aa, xi_l(ll), eta_l(ll));
                    dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                    dx_deta = dx_deta + x_ele(aa) * Na_eta;
                    dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                    dy_deta = dy_deta + y_ele(aa) * Na_eta;
                end
            
                detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
                
                for aa = 1 : n_en
                    Na = Quad(aa, xi_l(ll), eta_l(ll));
                    [Na_xi, Na_eta] = Quad_grad(aa, xi_l(ll), eta_l(ll));
                    Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                    Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                  
                    for jj = 1 : 2
                        pp = 2*(aa-1)+jj;
                        f_ele(pp) = f_ele(pp) + weight_l(ll) * h(IEN(ee, ii), jj) * Na * detJ;
                    end
                end % end of aa loop
                
            end % end of quadrature loop
        end
    end

    for aa = 1 : n_en
        for ii = 1 : 2
            pp = 2*(aa-1)+ii;
            PP = LM(ii, ee, aa);
            if PP == 0
                continue;
            end
            F(PP) = F(PP) + f_ele(pp);

            for bb = 1 : n_en
                for jj = 1 : 2
                    qq = 2*(bb-1)+jj;
                    QQ = LM(jj, ee, bb);
                    if QQ > 0
                        K(PP, QQ) = K(PP, QQ) + k_ele(pp, qq);
                    else
                        F(PP) = F(PP) - g(IEN(ee, bb), jj) * k_ele(pp, qq);
                    end
                end
            end
        end
    end
end

% Get disp
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 2);

for ii = 1 : n_np
    for jj = 1 : 2
        index = ID(ii, jj);
        if index > 0
            disp(ii, jj) = dn(index);
        else
            disp(ii, jj) = g(ii, jj);
        end
    end
end

save("HEAT", "disp", "n_el_x", "n_el_y");

% Process stress & strain

% Plots

plotVector(n_el_x, n_el_y, disp);

% Error terms