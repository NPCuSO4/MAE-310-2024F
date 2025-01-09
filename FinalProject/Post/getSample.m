function [sam_x, sam_y, sam_u, sam_strain, sam_stress] = getSample(n_en, n_el, D, x_coor, y_coor, IEN, disp, n_sam, sam_xi, sam_eta)

sam_x = zeros(n_sam*n_el, 1);
sam_y = zeros(n_sam*n_el, 1);
sam_u = zeros(n_sam*n_el, 2);
sam_strain = zeros(n_sam*n_el, 3);
sam_stress = zeros(n_sam*n_el, 3);

for ee = 1 : n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );
    u_ele =   disp( IEN(ee, 1:n_en), :);
  
    for ii = 1 : n_sam
        jj = n_sam * (ee-1) + ii;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        du_dxi = 0.0; du_deta = 0.0;
        for aa = 1 : n_en
            sam_x(jj) = sam_x(jj) + x_ele(aa) * Quad(aa, sam_xi(ii), sam_eta(ii));
            sam_y(jj) = sam_y(jj) + y_ele(aa) * Quad(aa, sam_xi(ii), sam_eta(ii));
            sam_u(jj, :) = sam_u(jj, :) + u_ele(aa, :) * Quad(aa, sam_xi(ii), sam_eta(ii));
            [Na_xi, Na_eta] = Quad_grad(aa, sam_xi(ii), sam_eta(ii));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
            du_dxi  = du_dxi  + u_ele(aa, :) * Na_xi;
            du_deta = du_deta + u_ele(aa, :) * Na_eta;
        end
        %detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        
        tmp_strain = zeros(3, 1);
        tmp_strain(1) = du_dxi(1)/dx_dxi + du_deta(1)/dx_deta;
        tmp_strain(2) = du_dxi(2)/dy_dxi + du_deta(2)/dy_deta;
        tmp_strain(3) = du_dxi(1)/dy_dxi + du_deta(1)/dy_deta ...
                            + du_dxi(2)/dx_dxi + du_deta(2)/dx_deta;
        tmp_stress = D * tmp_strain;
        sam_strain(jj, :) = tmp_strain(:);
        sam_stress(jj, :) = tmp_stress(:);
    end

end