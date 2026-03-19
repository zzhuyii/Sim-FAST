function [pass, modeStr, Pn, phi, phiPn, DCR] = Check_Truss_LRFD(Pu, Ag, An, E, KL, r, ...
                                                                 Fy, Fu, Rp )

%CHECK_TRUSS_LRFD  AASHTO LRFD axial member check 
%
% =========================================================================
% AASHTO LRFD REFERENCES:
%   Resistance factors      : Article 6.5.4.2
%   Tension resistance      : Article 6.8.2.1 (Eqs. 6.8.2.1-1 and 6.8.2.1-2)
%   Tension slenderness     : Article 6.8.4
%   Compression resistance  : Article 6.9.4.1 (Eqs. 6.9.4.1.1-1 and 6.9.4.1.1-2)
%   Compression slenderness : Article 6.9.3
% =========================================================================
%
% Variables:
%   Pu          : axial force (N). Tension positive, compression negative.
%   A           : gross cross-sectional area Ag (m^2)
%   E           : elastic modulus (Pa)
%   KL          : effective length K*L (m)
%   r           : governing radius of gyration (m) — use min(rx,ry)
%   Fy          : specified minimum yield strength (Pa)
%
% Variables: — resistance factors:
%   phi_ty      : phi for tension yielding in gross section
%                 AASHTO 6.5.4.2: phi_y = 0.95
%   phi_cb      : phi for axial compression, steel only
%                 AASHTO 6.5.4.2: phi_c = 0.95
%
% Variables: — net section fracture per 6.8.2.1-2:
%   Fu          : tensile strength (Pa). Default: 427e6 (A500 Gr.C)
%   An          : net area (m^2). Default: An = Ag (welded, no holes)
%   Rp          : hole factor. 1.0 = drilled/reamed; 0.90 = punched full size
%                 Default: 1.0
%   U_shlag     : shear lag reduction factor (Article 6.8.2.2).
%                 1.0 if force transmitted to all elements.
%                 Default: 1.0
%   phi_uf      : phi for tension fracture in net section
%                 AASHTO 6.5.4.2: phi_u = 0.80
%
% Member classification for slenderness limits:
%   memberType  : string — "TopChord", "BottomChord" (primary)
%                           or "Web" (secondary)
%                 Assume all members are primary
%
% OUTPUTS:
%   pass        : true if DCR <= 1.0
%   modeStr     : governing limit state description string
%   Pn          : nominal resistance (N)
%   phi         : resistance factor used for governing limit state
%   phiPn       : factored resistance phi*Pn (N)
%   DCR         : demand-to-capacity ratio = |Pu| / |phi*Pn|
%
% SLENDERNESS WARNINGS (printed but do not alter DCR):
%   Compression primary   : KL/r <= 120  (AASHTO 6.9.3)
%   Compression secondary : KL/r <= 140  (AASHTO 6.9.3)
%   Tension primary       : l/r  <= 200  (AASHTO 6.8.4)
%   Tension secondary     : l/r  <= 240  (AASHTO 6.8.4)



    % Resistance factors — AASHTO LRFD 6.5.4.2
    phi_ty   = 0.95;  % 6.5.4.2: tension yielding
    phi_cb   = 0.95;  % 6.5.4.2: axial compression, steel only
    phi_uf   = 0.80;  % 6.5.4.2: fracture in net section

    % Net section fracture parameters — AASHTO 6.8.2.1-2
    U_shlag  = 1.0;   % shear lag factor for connections, here is 1.0
    
    KLr = KL / r;   % KL/r (dimensionless); used for both compression and tension checks

    if Pu >= 0
        % =====================================================================
        % TENSION MEMBER
        % AASHTO LRFD Article 6.8.2.1
        % Pr = lesser of:
        %   Eq. 6.8.2.1-1: phi_y * Fy * Ag        (yielding on gross section)
        %   Eq. 6.8.2.1-2: phi_u * Fu * An * Rp * U (fracture on net section)
        % =====================================================================

        % Slenderness check — AASHTO 6.8.4
        % NOTE: l/r uses unbraced length only (no K factor for tension members)
        Lr_limit = 200; 
        if KLr > Lr_limit
            fprintf('[WARN] Bar l/r = %.1f > AASHTO 6.8.4 limit = %d\n', ...
                    KLr, Lr_limit);
        end

        % Eq. 6.8.2.1-1: Gross section yielding
        Pny = Fy * Ag;
        Pr_yield    = phi_ty * Pny;

        % Eq. 6.8.2.1-2: Net section fracture
        Pnu = Fu * An * Rp * U_shlag;
        Pr_fracture = phi_uf * Pnu;

        % Governing tensile resistance = lesser of the two (6.8.2.1)
        if Pr_yield <= Pr_fracture
            Pn      = Pny;
            phi     = phi_ty;
            phiPn   = Pr_yield;
            modeStr = sprintf('Tension-Yield      [6.8.2.1-1] phi_y=%.2f | Pr_yield=%.1fN <= Pr_frac=%.1fN', ...
                              phi_ty, Pr_yield, Pr_fracture);
        else
            Pn      = Pnu;
            phi     = phi_uf;
            phiPn   = Pr_fracture;
            modeStr = sprintf('Tension-Fracture   [6.8.2.1-2] phi_u=%.2f | Pr_frac=%.1fN < Pr_yield=%.1fN', ...
                              phi_uf, Pr_fracture, Pr_yield);
        end

    else
        % =====================================================================
        % COMPRESSION MEMBER
        % AASHTO LRFD Article 6.9.4.1
        %
        % Po = Fy * Ag                       (nominal yield resistance)
        % Pe = pi^2 * E * Ag / (KL/r)^2     (elastic critical buckling resistance)
        %
        % If Po/Pe <= 2.25  [Eq. 6.9.4.1.1-1]:  Pn = [0.658^(Po/Pe)] * Po
        % Otherwise         [Eq. 6.9.4.1.1-2]:  Pn = 0.877 * Pe
        %
        % NOTE: Po/Pe = (KL/r)^2 * Fy / (pi^2 * E) = lambda_c^2
        %       The Po/Pe notation is used here to match AASHTO 6.9.4.1 directly.
        % =====================================================================

        % Slenderness check — AASHTO 6.9.3
        KLr_limit = 120;
        if KLr > KLr_limit
            fprintf('[WARN] Bar KL/r = %.1f > AASHTO 6.9.3 limit = %d\n', ...
                    KLr, KLr_limit);
        end

        % Nominal yield and elastic buckling resistances
        Po = Fy * Ag;                             % nominal yield resistance (N)
        Pe = (pi^2 * E * Ag) / (KLr^2);          % elastic critical buckling resistance (N)

        ratio = Po / Pe;   % = (KL/r)^2 * Fy/(pi^2*E)

        if ratio <= 2.25
            % Eq. 6.9.4.1.1-1: Inelastic buckling
            Pn      = (0.658^ratio) * Po;
            modeStr = sprintf('Compression-Inelastic [6.9.4.1.1-1] phi_c=%.2f | Po/Pe=%.3f<=2.25', ...
                              phi_cb, ratio);
        else
            % Eq. 6.9.4.1.1-2: Elastic buckling
            Pn      = 0.877 * Pe;
            modeStr = sprintf('Compression-Elastic   [6.9.4.1.1-2] phi_c=%.2f | Po/Pe=%.3f>2.25', ...
                              phi_cb, ratio);
        end

        phi = phi_cb;
        phiPn = phi * Pn;

    end

    % LRFD DEMAND-TO-CAPACITY RATIO
    DCR  = abs(Pu) / abs(phiPn);
    pass = (DCR <= 1.0);

end