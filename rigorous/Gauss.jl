function gauss(A_::Matrix)
    A = copy(A_)      
    m, n = size(A)
    @assert m <= n

    for k in 1:m
        
        # seek largest nonzero A[j,k] across j
        best, abs_best = k, abs(A[k,k])
        for i in k+1:m
            if abs(A[i,k]) > abs_best
                best, abs_best = i, abs(A[i,k])
            end
        end

        # swap rows to position pivot
        A[k,:], A[best,:] = A[best,:], A[k,:]

        # divide row by pivot
        pivot = A[k,k]
        for j in k:n
            A[k,j] /= pivot
        end

        # clear zeros above and below
        for k_that in 1:m
            if k_that != k
                ratio = A[k_that,k]/A[k,k]
                for j in k:n
                    A[k_that,j] -= ratio*A[k,j]
                end
            end
        end
    end
    A
end
