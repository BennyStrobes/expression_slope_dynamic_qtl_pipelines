args = commandArgs(trailingOnly=TRUE)
library(lme4)



run_linear_model <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[4,4]
        coef <- summary(fit)$coefficients[4,1]
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[4,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_all") {
        fit_full <- lmer(expr ~ genotype + time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines) + (1|time_steps), REML=FALSE)
        fit_null <- lmer(expr ~ genotype + time_steps + (1|cell_lines) + (1|time_steps), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[4,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_cell_line_slope_correlated") {
        fit_full <- lmer(expr ~ genotype + time_steps + genotype_interaction:time_steps_interaction + (time_steps | cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~ genotype + time_steps + (time_steps | cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[4,1]
        #tvalue <- coefs[4,3]
        #pvalue2 <- 2*pt(-abs(tvalue),df=num_samp-1)
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_cell_line_slope_uncorrelated") {
        fit <- lmer(expr ~ genotype + time_steps + genotype_interaction:time_steps_interaction + (time_steps || cell_lines))
        coefs <- data.frame(coef(summary(fit)))
        coef <- coefs[4,1]
        tvalue <- coefs[4,3]
        pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
    }

    return(list(coef=coef, pvalue=pvalue))
}


run_linear_model_one_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1:time_steps_interaction + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[5,4]
        coef <- summary(fit)$coefficients[5,1]
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1:time_steps_interaction + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1:time_steps_interaction + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[5,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_all") {
        fit_full <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + genotype_interaction:time_steps_interaction + (1|cell_lines) + (1|time_steps), REML=FALSE)
        fit_null <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + (1|cell_lines) + (1|time_steps), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[5,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_cell_line_slope_correlated") {
        fit_full <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + genotype_interaction:time_steps_interaction + (time_steps | cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + (time_steps | cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[5,1]
        #tvalue <- coefs[4,3]
        #pvalue2 <- 2*pt(-abs(tvalue),df=num_samp-1)
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 

    return(list(coef=coef, pvalue=pvalue))
}

run_linear_model_two_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[6,4]
        coef <- summary(fit)$coefficients[6,1]
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[6,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_all") {
        fit_full <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + genotype_interaction:time_steps_interaction + (1|cell_lines) + (1|time_steps), REML=FALSE)
        fit_null <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + (1|cell_lines) + (1|time_steps), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[6,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_cell_line_slope_correlated") {
        fit_full <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + genotype_interaction:time_steps_interaction + (time_steps | cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + (time_steps | cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[6,1]
        #tvalue <- coefs[4,3]
        #pvalue2 <- 2*pt(-abs(tvalue),df=num_samp-1)
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 

    return(list(coef=coef, pvalue=pvalue))
}

run_linear_model_three_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2, cov3) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + cov3:time_steps_interaction + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[7,4]
        coef <- summary(fit)$coefficients[7,1]
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + cov3:time_steps_interaction + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + cov3:time_steps_interaction + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[7,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_all") {
        fit_full <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + cov3:time_steps_interaction + genotype_interaction:time_steps_interaction + (1|cell_lines) + (1|time_steps), REML=FALSE)
        fit_null <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + cov3:time_steps_interaction + (1|cell_lines) + (1|time_steps), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[7,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_cell_line_slope_correlated") {
        fit_full <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + cov3:time_steps_interaction + genotype_interaction:time_steps_interaction + (time_steps | cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~ genotype + time_steps + cov1:time_steps_interaction + cov2:time_steps_interaction + cov3:time_steps_interaction + (time_steps | cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[7,1]
        #tvalue <- coefs[4,3]
        #pvalue2 <- 2*pt(-abs(tvalue),df=num_samp-1)
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 

    return(list(coef=coef, pvalue=pvalue))
}


null_response <- function(){
    return(list(coef=0,pvalue=1))
}

# Compute a slope for each cell line that represents the direction of expression for this gene across time
compute_cell_line_slopes <- function(time_steps, expr, permute) {
    # Initialize output vector
    slopes <- c()
    num_lines <- length(time_steps)
    # Compute cell line slope for each line (loop through lines)
    for (line_num in 1:num_lines) {
        cell_line_time_steps <- as.numeric(strsplit(time_steps[line_num], ',')[[1]])
        cell_line_expr <- as.numeric(strsplit(expr[line_num], ',')[[1]])
        # Standardize Expression Data
        cell_line_expr <- (cell_line_expr - mean(cell_line_expr))/sd(cell_line_expr)
        # if permute == True, shuffle time steps
        if (permute == "True") {
            cell_line_time_steps <- sample(cell_line_time_steps)
        }
        # Fit linear model
        fit <- lm(cell_line_expr ~ cell_line_time_steps)
        # Extract slope of time step variable
        cell_line_slope <- summary(fit)$coefficients[2,1]
        # Add slope of time step variable to vector keep track of slopes from all cell lines
        slopes <- c(slopes, cell_line_slope)
    }
    return(as.numeric(slopes))
}


fit_genotype_to_cell_line_slopes_no_cov <- function(cell_line_slopes,genotype) {
    # Standardize cell_line_slopes
    cell_line_slopes <- (cell_line_slopes - mean(cell_line_slopes))/sd(cell_line_slopes)
    # Fit linear model
    fit <- lm(cell_line_slopes ~ genotype)
    # Extract relevent data from LM
    coef <- summary(fit)$coefficients[2,1]
    pvalue <- summary(fit)$coefficients[2,4]
    return(list(coef=coef,pvalue=pvalue))
}

fit_genotype_to_cell_line_slopes_one_cov <- function(cell_line_slopes,genotype, cov1) {
    # Standardize cell_line_slopes
    cell_line_slopes <- (cell_line_slopes - mean(cell_line_slopes))/sd(cell_line_slopes)
    # Fit linear model
    fit <- lm(cell_line_slopes ~ genotype + cov1)

    # Extract relevent data from LM
    coef <- summary(fit)$coefficients[2,1]
    pvalue <- summary(fit)$coefficients[2,4]
    return(list(coef=coef,pvalue=pvalue))
}

fit_genotype_to_cell_line_slopes_two_cov <- function(cell_line_slopes,genotype, cov1, cov2) {
    # Standardize cell_line_slopes
    cell_line_slopes <- (cell_line_slopes - mean(cell_line_slopes))/sd(cell_line_slopes)
    # Fit linear model
    fit <- lm(cell_line_slopes ~ genotype + cov1 + cov2)

    # Extract relevent data from LM
    coef <- summary(fit)$coefficients[2,1]
    pvalue <- summary(fit)$coefficients[2,4]
    return(list(coef=coef,pvalue=pvalue))
}

fit_genotype_to_cell_line_slopes_three_cov <- function(cell_line_slopes,genotype, cov1, cov2, cov3) {
    # Standardize cell_line_slopes
    cell_line_slopes <- (cell_line_slopes - mean(cell_line_slopes))/sd(cell_line_slopes)
    # Fit linear model
    fit <- lm(cell_line_slopes ~ genotype + cov1 + cov2 + cov3)

    # Extract relevent data from LM
    coef <- summary(fit)$coefficients[2,1]
    pvalue <- summary(fit)$coefficients[2,4]
    return(list(coef=coef,pvalue=pvalue))
}


#####################
# Command line args
#####################


input_data_file = args[1]
covariate_method = args[2]
output_file = args[3]
permute = args[4]



stop = FALSE
count = 0
f = file(input_data_file, "r")

sink(output_file)

while(!stop) {
    # Read current line
    next_line = readLines(f, n = 1)
    # Read current line


    # Parse the line
    data = strsplit(next_line,'\t')[[1]]    
    rs_id = data[1]
    ensamble_id = data[2]
    cell_lines = as.factor(strsplit(data[3],';')[[1]])
    pc1 = as.numeric(strsplit(data[4],';')[[1]])
    pc2 = as.numeric(strsplit(data[5],';')[[1]])
    pc3 = as.numeric(strsplit(data[6],';')[[1]])
    genotype = as.numeric(strsplit(data[7],';')[[1]])

    time_steps = strsplit(data[8],';')[[1]]
    expr = strsplit(data[9],';')[[1]]



    tryCatch(
    {
        # Compute a slope for each cell line that represents the direction of expression for this gene across time
        cell_line_slopes <- compute_cell_line_slopes(time_steps, expr, permute)

        # Fit linear model of genotype onto cell_line_slopes
        if (covariate_method == "none") {
            lm_results <- fit_genotype_to_cell_line_slopes_no_cov(cell_line_slopes,genotype)
        } else if (covariate_method == "pc1") {
            lm_results <- fit_genotype_to_cell_line_slopes_one_cov(cell_line_slopes,genotype, pc1)
        } else if (covariate_method == "pc1_2") {
            lm_results <- fit_genotype_to_cell_line_slopes_two_cov(cell_line_slopes,genotype, pc1, pc2)
        } else if (covariate_method == "pc1_3") {
            lm_results <- fit_genotype_to_cell_line_slopes_three_cov(cell_line_slopes,genotype, pc1, pc2, pc3)
        }

        # print result to output file!!
        new_line <- paste0(next_line,"\t",lm_results$coef,"\t", lm_results$pvalue,"\n")
        cat(new_line)
    },
    error = function(e){
        new_line <- paste0(next_line,"\t",0.0,"\t", 1.0,"\n")
        cat(new_line)
    }
    )

    count = count + 1

    if(length(next_line) == 0) {
        stop = TRUE
        close(f)
    }

}
# close output file handle
sink()

