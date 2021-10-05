%
% Tangent stiffness of multiplicative plasticity with linear hardening
%
function [Dtan] = mulPlastTan(mp, D, L, b, alpha, ep)
    %
    EPS = 1E-12;
    Iden = [1 1 1]'; two3 = 2/3; stwo3 = sqrt(two3); %constants
    mu = mp(2); beta = mp(3); H = mp(4); Y0 = mp(5); %material properties
    ftol = Y0 * 1E-6; %tolerance for yield
    R = inv(eye(3) - L); %inc. deformation gradient
    bm = [b(1) b(4) b(6); b(4) b(2) b(5); b(6) b(5) b(3)];
    bm = R * bm * R'; %trial elastic left C-G
    b = [bm(1, 1) bm(2, 2) bm(3, 3) bm(1, 2) bm(2, 3) bm(1, 3)]';
    [~, P] = eig(bm); %eigenvalues
    eigen = sort(real([P(1, 1) P(2, 2) P(3, 3)]))'; %principal stretch
    %
    % Duplicated eigenvalues
    TMP = -1;

    for I = 1:2

        if abs(eigen(1) - eigen(3)) < EPS
            eigen(I) = eigen(I) + TMP * EPS;
            TMP = -TMP;
        end

    end

    if abs(eigen(1) - eigen(2)) < EPS; eigen(2) = eigen(2) + EPS; end;

        if abs(eigen(2) - eigen(3)) < EPS; eigen(2) = eigen(2) + EPS; end;
            %
            % EIGENVECTOR MATRIX N*N = M(6,*)
            M = zeros(6, 3); %eigenvector matrices

            for K = 1:3
                KB = 1 + mod(K, 3);
                KC = 1 + mod(KB, 3);
                EA = eigen(K);
                EB = eigen(KB);
                EC = eigen(KC);
                D1 = EB - EA;
                D2 = EC - EA;
                DA = 1 / (D1 * D2);
                M(1, K) = ((b(1) - EB) * (b(1) - EC) + b(4) * b(4) + b(6) * b(6)) * DA;
                M(2, K) = ((b(2) - EB) * (b(2) - EC) + b(4) * b(4) + b(5) * b(5)) * DA;
                M(3, K) = ((b(3) - EB) * (b(3) - EC) + b(5) * b(5) + b(6) * b(6)) * DA;
                M(4, K) = (b(4) * (b(1) - EB + b(2) - EC) + b(5) * b(6)) * DA;
                M(5, K) = (b(5) * (b(2) - EB + b(3) - EC) + b(4) * b(6)) * DA;
                M(6, K) = (b(6) * (b(3) - EB + b(1) - EC) + b(4) * b(5)) * DA;
            end

            %
            deps = 0.5 * log(eigen); %logarithmic
            sigtr = D * deps; %trial principal stress
            sig = sigtr;
            eta = sigtr - alpha - sum(sigtr) * Iden / 3; %shifted stress
            etat = norm(eta); %norm of eta
            fyld = etat - stwo3 * (Y0 + (1 - beta) * H * ep); %trial yield function

            if fyld >= ftol %yield test
                gamma = fyld / (2 * mu + two3 * H); %plastic consistency param
                N = eta / etat; %unit vector normal to f
                sig = sigtr - 2 * mu * gamma * N; %updated stress
                var1 = 4 * mu^2 / (2 * mu + two3 * H);
                var2 = 4 * mu^2 * gamma / etat; %coefficients
                D = D - (var1 - var2) * (N * N') + var2 * (Iden * Iden') / 3; %tangent stiffness
                D(1, 1) = D(1, 1) - var2; %contr. from 4th-order I
                D(2, 2) = D(2, 2) - var2;
                D(3, 3) = D(3, 3) - var2;
            end

            J1 = sum(eigen);
            J3 = eigen(1) * eigen(2) * eigen(3);
            I2 = [1 1 1 0 0 0]';
            I4 = eye(6); I4(4, 4) = .5; I4(5, 5) = .5; I4(6, 6) = .5;
            Ibb = [0, b(4)^2 - b(1) * b(2), b(6)^2 - b(1) * b(3), 0, b(4) * b(6) - b(1) * b(5), 0;

                b(4) * b(4) - b(1) * b(2), 0, b(5)^2 - b(2) * b(3), 0, 0, b(4) * b(5) - b(2) * b(6);

                b(6)^2 - b(1) * b(3), b(5)^2 - b(2) * b(3), 0, b(5) * b(6) - b(3) * b(4), 0, 0;

                0, 0, b(5) * b(6) - b(3) * b(4), (b(1) * b(2) - b(4)^2) / 2, (b(2) * b(6) - b(4) * b(5)) / 2, (b(1) * b(5) - b(4) * b(6)) / 2;

                b(4) * b(6) - b(1) * b(5), 0, 0, (b(2) * b(6) - b(4) * b(5)) / 2, (b(2) * b(3) - b(5)^2) / 2, (b(3) * b(4) - b(5) * b(6)) / 2;

                0, b(4) * b(5) - b(2) * b(6), 0, (b(1) * b(5) - b(4) * b(6)) / 2, (b(3) * b(4) - b(5) * b(6)) / 2, (b(1) * b(3) - b(6)^2) / 2];
            %
            d1 = 1 / ((eigen(2) - eigen(1)) * (eigen(3) - eigen(1)));
            d2 = 1 / ((eigen(3) - eigen(2)) * (eigen(1) - eigen(2)));
            d3 = 1 / ((eigen(1) - eigen(3)) * (eigen(2) - eigen(3)));
            t11 = -J3 * d1 / eigen(1); t12 = -J3 * d2 / eigen(2); t13 = -J3 * d3 / eigen(3);
            t21 = d1 * eigen(1); t22 = d2 * eigen(2); t23 = d3 * eigen(3);
            t31 = t21 * (J1 - 4 * eigen(1)); t32 = t22 * (J1 - 4 * eigen(2)); t33 = t23 * (J1 - 4 * eigen(3));
            %
            CT1 = d1 * Ibb + t11 * (I4 - (I2 - M(:, 1)) * (I2 - M(:, 1))') + t21 * (b * M(:, 1)' + M(:, 1) * b') + t31 * M(:, 1) * M(:, 1)';
            CT2 = d2 * Ibb + t12 * (I4 - (I2 - M(:, 2)) * (I2 - M(:, 2))') + t22 * (b * M(:, 2)' + M(:, 2) * b') + t32 * M(:, 2) * M(:, 2)';
            CT3 = d3 * Ibb + t13 * (I4 - (I2 - M(:, 3)) * (I2 - M(:, 3))') + t23 * (b * M(:, 3)' + M(:, 3) * b') + t33 * M(:, 3) * M(:, 3)';
            %
            Dtan = M * D * M' + 2 * (CT1 * sig(1) + CT2 * sig(2) + CT3 * sig(3));
