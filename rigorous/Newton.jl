function matrix_elements_old(L, phi, phi_inv, basis_size) 
    stack([phi(L(phi_inv(unit_vector(k,basis_size)))) for k=0:basis_size-1])
end

matrix_rows(L, phi, phi_inv, basis_size) = k->phi(L(phi_inv(unit_vector(k, basis_size))))

function matrix_elements(L, phi, phi_inv, basis_size)
    stack(pmap(matrix_rows(L, phi, phi_inv, basis_size), 0:basis_size-1))
end

eye(n) = [(j==k) ? RR(1) : RR(0) for j=1:n,k=1:n]

function newton_one_step(p, F_and_DF, phi, phi_inv, basis_size)
    Fp, DFp = F_and_DF(p)
    M = matrix_elements(DFp, phi, phi_inv, basis_size)
    p .- phi_inv(M \ phi(Fp))
end
    
function newton_F(p0, F_and_DF, phi, phi_inv, basis_size)
    p = p0
    Fp, DFp = F_and_DF(p)
    err = sum(abs.(Fp))
    println(0, " ", err)
    for k=1:20
        p_new = newton_one_step(p, F_and_DF, phi, phi_inv, basis_size)
        Fp_new, DFp_new = F_and_DF(p_new)
        err_new = sum(abs.(Fp_new))
        println(k, " ", err_new)
        if err_new >= err
            break 
        else
            p, err = p_new, err_new
        end
    end
    println("process ended")
    p
end
    