function [E_L2, E_H1] = calcError(n_en, n_el, x_coor, y_coor, IEN, analy_u, analy_du, disp, n_E_xi, n_E_eta)

[xi, eta, weight] = Gauss2D(n_E_xi, n_E_eta);
n_int_E = n_E_xi * n_E_eta;

L2_top = 0.0; L2_bot = 0.0; H1_top = 0.0; H1_bot = 0.0;

for ee = 1 : n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );
    u_ele =   disp( IEN(ee, 1:n_en), :);
  
    for ll = 1 : n_int_E
        x_l = 0.0; y_l = 0.0; u_l = zeros(2, 1);
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        du_dxi = 0.0; du_deta = 0.0;
        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            u_l(1) = u_l(1) + u_ele(aa, 1) * Quad(aa, xi(ll), eta(ll));
            u_l(2) = u_l(2) + u_ele(aa, 2) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
            du_dxi  = du_dxi  + u_ele(aa, :) * Na_xi;
            du_deta = du_deta + u_ele(aa, :) * Na_eta;
        end
        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        analy_ut = analy_u(x_l, y_l);
        analy_dut = analy_du(x_l, y_l);
        L2_top = L2_top + weight(ll) * (u_l(1) - analy_ut(1))^2 * detJ;
        L2_top = L2_top + weight(ll) * (u_l(2) - analy_ut(2))^2 * detJ;
        L2_bot = L2_bot + weight(ll) * (analy_ut(1)^2 + analy_ut(2)^2) * detJ;

        H1_top = H1_top + weight(ll) * (du_dxi(1)/dx_dxi - analy_dut(1))^2 * detJ;
        H1_top = H1_top + weight(ll) * (du_deta(1)/dy_deta - analy_dut(2))^2 * detJ;
        H1_top = H1_top + weight(ll) * (du_dxi(2)/dx_dxi - analy_dut(3))^2 * detJ;
        H1_top = H1_top + weight(ll) * (du_deta(2)/dy_deta - analy_dut(4))^2 * detJ;
        H1_bot = H1_bot + weight(ll) * (analy_dut(1)^2 + analy_dut(2)^2) * detJ;
        H1_bot = H1_bot + weight(ll) * (analy_dut(3)^2 + analy_dut(4)^2) * detJ;
    end
end

E_L2 = (L2_top/L2_bot)^0.5;
E_H1 = (H1_top/H1_bot)^0.5;