function [q, J] = assemble_newton_schemeApproxChi(var, z, Phi, Chi, dChidz, d2Chidz2, varphi)
%Assembles the Newton scheme for the approximated system, using the approximations for Chi. In

    %Initialises needed variables
    dz = z(2) - z(1);
    N  = numel(var);
    M  = N + 2;

    %Helper functions
    A = @(v_i) Phi*v_i./((1-v_i)*(1-Phi));
    Ap = @(v_i) Phi./((1-v_i).^2*(1-Phi));
    App= @(v_i) 2*Phi./((1-v_i).^3*(1-Phi));

    B = @(v_i) 3*v_i.^3./(1-v_i);
    Bp= @(v_i) 3*v_i.^2.*(3-2*v_i)./(1-v_i).^2;
    Bpp = @(v_i) 6*v_i.*(3-3*v_i+v_i.^2)./(1-v_i).^3;


    %--- Build neighbor arrays ---
    vm = [varphi;      var(1:end-1)];   % var_{i-1}
    v  = var(:);                         % var_i
    vp = [var(2:end);  Phi];            % var_{i+1}

    %--- Pull off derivatives at full-index i+1 ---
    Chi_trunc = Chi(2:end-1);
    DChi_trunc = dChidz(2:M-1);
    DDChi_trunc = d2Chidz2(2:M-1);



    %------------- Residual -------------

    q = A(var).*(1+Chi_trunc).*(vp-2*var+vm)/dz^2 +Ap(var).*(1+Chi_trunc).*( (vp-vm)/(2*dz) ).^2 ...
        +(4*A(var)-Bp(var).*(1+Chi_trunc)).*DChi_trunc.*(vp-vm)/(2*dz) - B(var).*(1+Chi_trunc).*DDChi_trunc-4*B(var).*DChi_trunc.^2;


    %------------- Interior Jacobian diagonals -------------

    jm1 = A(var).*(1+Chi_trunc)/dz^2-2*Ap(var).*(1+Chi_trunc).*(vp-vm)./(4*dz^2) ...
           -(4*A(var)-Bp(var).*(1+Chi_trunc)).*DChi_trunc/(2*dz);    % sub-diagonal J(i,i-1)

    jD  = Ap(var).*(1+Chi_trunc).*(vp-2*var+vm)/dz^2 -2*A(var).*(1+Chi_trunc)/dz^2 ...
           +App(var).*(1+Chi_trunc).*( (vp-vm)/(2*dz) ).^2 + (4*Ap(var)-Bpp(var).*(1+Chi_trunc) ).*DChi_trunc.*(vp-vm)/(2*dz) - Bp(var).*(1+Chi_trunc).*DDChi_trunc-4*Bp(var).*DChi_trunc.^2;  % main diagonal J(i,i)

 
    jp1 = ...
  A(var).*(1+Chi_trunc)/dz^2                          ... % A(v_i)*(1+F_i)/Δz^2
+ 2*Ap(var).*(1+Chi_trunc).*(vp-vm)/(4*dz^2)         ... % 2*A'(v_i)*(1+F_i)*(v_{i+1}-v_{i-1})/(4 Δz^2)
+ (4*A(var) - Bp(var).*(1+Chi_trunc)).*DChi_trunc/(2*dz);      % (2A - B'(1+F))*F'/(2 Δz)


%Efficient sparse indexing
    J = spdiags([[jm1(1:end-1);0],jD,[0;jp1(2:end)]],[-1,0,1], N, N);


end
