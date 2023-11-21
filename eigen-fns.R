cijk = function(i, j, k) {
    if(i > j) {
        stop("must have i <= j")
    }
    if(k >= 2^(j - i - 1)) {
        stop("must have k < 2^(j-i-1)")
    }
    if(i == j) {
        cijk = rep(1, 2^j)
        return(cijk)
    }
    cijk = rep(0, 2^j)
    start = 2^(i+1) * k + 1
    end = start + 2^(i+1) - 1
    replacement = c(rep(1, 2^i), rep(-1, 2^i))
    cijk[start:end] = replacement
    return(cijk)
}

get_evec_order = function(i, j, k) {
    if(i == j) return(1)
    if(i == (j - 1)) return(2)
    return(2^(j - i - 1) + k + 1)
}

Dij = function(i, j, value = 1) {
    block_size = 2^i
    matrix_size = 2^j
    m = matrix(value, nrow = block_size, ncol = block_size)
    n_blocks = matrix_size / block_size
    matrix_list = rep(list(m), n_blocks)
    return(as.matrix(Matrix::bdiag(matrix_list)))
}

get_evals = function(r, evals) {
    if(r < 1 && r > 0) {
        return((rep(r^(-1), length(evals)) + (1-r)^(-1) * evals^(-1))^(-1))
    } else if(r == 1) {
        return(evals)
    } else if(r == 0){
        return(rep(1, length(evals)))
    }
}
