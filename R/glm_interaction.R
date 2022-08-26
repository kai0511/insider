#' @export
glm_interaction <- function(residual, train_indicator, interaction_indicator, column_factor, tol = 1e-10, n_cores =10){

    unique_ita <- unique(interaction_indicator)
    coeff_matrix <- matrix(0, nrow = length(unique_ita), ncol = nrow(column_factor))
    pval_matrix <- matrix(0, nrow = length(unique_ita), ncol = nrow(column_factor))

    for(i in unique_ita) {

        ids <- which(interaction_indicator == i);

        st_idx <- 1; ed_idx <- 1
        nonzero_num <- length(ids) * ncol(column_factor);
        outcomes = rep(0,nonzero_num);
        features = matrix(0, nrow = nonzero_num, ncol = nrow(column_factor))

        for(k in ids){
            ed_idx = st_idx + ncol(column_factor) - 1;
            features[st_idx:ed_idx, ] = t(column_factor);
            outcomes[st_idx:ed_idx] = residual[k,];
            st_idx = ed_idx + 1;
        }

        data <- data.frame(response = outcomes, features)
        fit <- glm(response ~ . - 1, family = gaussian(), data = data)
        coeff_matrix[i,] <- unname(coefficients(fit))
        pval_matrix[i,] <- coef(summary(fit))[,4]
    }
    return(list(coeff_matrix, pval_matrix))
}


        
