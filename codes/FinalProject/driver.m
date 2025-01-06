clear all; clc;

E  = 1e9;
nu = 0.3;
D  = [1,  nu, 0;
      nu, 1,  0;
      0,  0,  (1-nu)/2;];
D  = E/(1 - nu^2) * D;

% exact solution
exact_u = @(x,y) -x^2;
exact_v = @(x,y) 0;

f_x = @(x,y) 2*E/(nu^2-1); % source term
f_y = @(x,y) 0; 

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

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

% ID array
ID = zeros(n_np,2);
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
    counter = counter + 2;
    ID(nx,1) = counter - 1;
    ID((n_np_y-1)*n_np_x + nx,1) = counter;
end
n_eq = counter;

LM = zeros(2, n_el, n_en);
for i = 1 : n_el
  for j  = 1 : n_en
    LM(1, i, j) = ID(IEN(i, j), 1);
    LM(2, i, j) = ID(IEN(i, j), 2);
  end
end

% allocate the stiffness matrix and load vector
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(n_en, n_en); % element stiffness matrix
  f_ele1 = zeros(n_en, 1);    % element load vector
  f_ele2 = zeros(n_en, 1);    % element load vector
  
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
      
      f_ele1(aa) = f_ele1(aa) + weight(ll) * f_x(x_l, y_l) * Na;
      f_ele2(aa) = f_ele2(aa) + weight(ll) * f_y(x_l, y_l) * Na;
      
      Ba = [Na_x, 0, Na_y; 0, Na_y, Na_x;];

      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        
        Bb = [Nb_x, 0; 0, Nb_y; Nb_y, Nb_x;];
        
        ed = [1,0; 0,1;];
        for ii = 1 : 2
           for jj = 1 : 2
               %pp = 2*(aa-1)+ii;
               %qq = 2*(bb-1)+jj;
               k_ele(aa, bb) = k_ele(aa, bb) + weight(ll) * ed(ii,:) * Ba * D * Bb * ed(:,jj);
           end
        end
        
      end % end of bb loop
    end % end of aa loop
    
  end % end of quadrature loop
 
  for aa = 1 : n_en
    PP = LM(1, ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele1(aa);
      
      for bb = 1 : n_en
        QQ = LM(1, ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
      end  
    end
    PP = LM(2, ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele2(aa);
      
      for bb = 1 : n_en
        QQ = LM(2, ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
      end  
    end
  end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp_x = zeros(n_np, 1);
disp_y = zeros(n_np, 1);

for ii = 1 : n_np
  index = ID(ii, 1);
  if index > 0
    disp_x(ii) = dn(index);
  else
    % modify disp with the g data. Here it does nothing because g is zero
  end

  index = ID(ii, 2);
  if index > 0
    disp_y(ii) = dn(index);
  else
    % modify disp with the g data. Here it does nothing because g is zero
  end
end

% save the solution vector and number of elements to disp with name
% HEAT.mat
save("HEAT", "disp_x", "n_el_x", "n_el_y");

% EOF