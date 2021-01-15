#' Selection Operator for Jointly analyzing multiple variants (SOJO)
#' 
#' This function computes penalized Selection Operator for JOintly analyzing multiple variants (SOJO) within a mapped locus,
#'  based on LASSO regression derived from GWAS summary statistics. 
#' 
#' @param sum.stat.discovery A data frame including GWAS summary statistics of genetic variants within a mapped locus. 
#' The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS; N, sample size.
#' @param sum.stat.validation A data frame including GWAS summary statistics from a validation dataset. It should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' Freq1, the allele frequency of Allele1; b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS; N, sample size.
#' @param LD_ref The reference LD correlation matrix including SNPs at the locus. The row names and column names of the matrix should be SNP names in reference sample.
#' @param snp_ref The reference alleles of SNPs in the reference LD correlation matrix. The names of the vector should be SNP names in reference sample.
#' @param v.y The phenotypic variance of the trait. Default is 1.
#' @param lambda.vec The tuning parameter sequence given by user. If not specified, the function will compute its own tuning parameter sequence
#' ,which is recommended.
#' @param standardize Logical value for genotypic data standardization, prior to starting the algorithm. 
#' The coefficients in output are always transformed back to the original scale. Default is \code{standardize = TRUE}.
#' @param nvar The number of variants aiming to be selected in the model. If \code{sum.stat.validation} is provided, \code{nvar} is the maximum number
#' of variants in the model.
#' For example, if \code{nvar = 5}, then the algorithm will stop before the sixth variant is selected. Default is 50.
#' 
#' @note Users can download reference LD correlation matrices from https://www.dropbox.com/home/sojo\%20reference\%20ld\%20matrix. 
#' These LD matrices are based on 612,513 chip markers in Swedish Twin Registry. If chip markers are only a small subset of the analysis, LD matrix from the 1000 Genomes Project
#' can be used (see the GitHub tutorial). The function will then take overlapping SNPs between summary statistics and reference LD matrix. 
#' 
#' The function returns results along the whole LASSO path when tuning parameter changes. Users can specify several tunining parameters or how many variants
#' should be selected.    
#' 
#' The optimal tuning parameter can be suggested by validation. If the GWAS summary statistics from a validation dataset are provided in \code{sum.stat.validation}, then the out of sample R^2 for each tuning parameter in \code{lambda.v} 
#' will be computed. The tuning parameter gives the largest out of sample R^2 will be considered as optimal. The optimal tuning parameter and the variants and their effect sizes
#' at this tuning parameter will be reported in \code{beta.opt} and \code{lambda.opt}.
#' 
#' When a tiny \code{lambda.vec} is specified, the LASSO solution is similar to the standard multiple regression, 
#' which may cause error due to complete LD between variants. 
#' 
#' Note the length of lambda.v in result may be longer than \code{nvar}. Because a lambda will be recorded when a variant is 
#' added into or removed from the model. 
#' 
#' @return A list is returned with:
#' \itemize{
#' \item{beta.opt }{The optimal variants and their effect sizes in terms of out of sample R^2. Only available when \code{sum.stat.validation} is provided.}
#' \item{lambda.opt }{The optimal tuning parameter in terms of out of sample R^2. Only available when \code{sum.stat.validation} is provided.}
#' \item{R2 }{The out of sample R^2 for each tuning parameter in \code{lambda.v}. Only available when \code{sum.stat.validation} is provided.}
#' \item{lambda.v }{The tuning parameter sequence actually used.}
#' \item{beta.mat }{The LASSO estimates at the tuning parameters in \code{lambda.v} stored in sparse matrix format. The reference alleles in results are same as those in the discovery gwas results.}
#' \item{selected.markers }{The vector of selected variants. The variants being ahead are selected earlier in LASSO path.}
#' }
#' 
#' @author Zheng Ning
#' 
#' @references 
#' Ning Z, Lee Y, Joshi PK, Wilson JF, Pawitan Y, Shen X (2017). A selection operator for summary association statistics 
#' reveals locus-specific allelic heterogeneity of complex traits. \emph{Submitted}.
#' 
#' @seealso 
#' sojo tutorial: https://github.com/zhenin/sojo
#' @import Matrix
#' 
#' @examples 
#'\dontrun{
#' ## The GWAS summary statistics of SNPs in 1 MB window centred at rs11090631 
#' data(sum.stat.discovery)
#' head(sum.stat.discovery)
#' 
#' ## The reference matrix and corresponding reference alleles 
#' download.file("https://www.dropbox.com/s/ty1udfhx5ohauh8/LD_chr22.rda?raw=1", destfile = paste0(find.package('sojo'), "example.rda"))
#' load(file = paste0(find.package('sojo'), "example.rda"))
#' 
#' res <- sojo(sum.stat.discovery, LD_ref = LD_mat, snp_ref = snp_ref, nvar = 20)
#' 
#' ## LASSO path plot
#' matplot(log(res$lambda.v), t(as.matrix(res$beta.mat)), lty = 1, type = "l", xlab = expression(paste(log, " ",lambda)), 
#' ylab = "Coefficients", main = "Summary-level LASSO")
#' 
#' ## LASSO solution for user supplied tuning parameters
#' res2 <- sojo(sum.stat.discovery = sum.stat.discovery, LD_ref = LD_mat, snp_ref = snp_ref, lambda.vec = c(0.004,0.002))
#' 
#' 
#' ## LASSO solution and the optimal tuning parameter when validation dataset is available
#' data(sum.stat.validation)
#' head(sum.stat.validation)
#' 
#' res.valid <- sojo(sum.stat.discovery, sum.stat.validation = sum.stat.validation, LD_ref = LD_mat, snp_ref = snp_ref, nvar = 20)
#' res.valid$beta.opt  # the optimal variants and their effect sizes
#' res.valid$lambda.opt  # the optimal tuning parameter
#' res.valid$R2  # out of sample R^2
#' }
#' @export
#' 
sojo <- 
  function(sum.stat.discovery, sum.stat.validation = NULL, LD_ref, snp_ref, v.y=1, lambda.vec=NA, standardize = T, nvar = 50){
    
    colnames_input <- c("SNP", "A1", "A2","b", "se", "N")
    colnames_lack <- setdiff(colnames_input, intersect(colnames(sum.stat.discovery), colnames_input))
    
    
    if(length(colnames_lack) > 0){
      colnames_lack <- paste(colnames_lack, collapse = ', ')
      stop(paste("The following columns are missing:", colnames_lack))
    }
    
    if(ncol(LD_ref) != length(snp_ref)){
      stop("The SNPs in reference LD matrix and its reference allele vector does't match! Please check.")
    }
    
    
    rownames(LD_ref) <- colnames(LD_ref) <- names(snp_ref)
    
    if(is.null(sum.stat.validation)){
      snps.overlap <- intersect(sum.stat.discovery$SNP,names(snp_ref))
      if(length(snps.overlap) == 0){
        stop("There is no overlapping SNPs between summary statistics and reference LD matrix! Please check.")
      }
      rownames(sum.stat.discovery) <- sum.stat.discovery$SNP
      sum.stat <- sum.stat.discovery[snps.overlap,]
    } else{
      snps.overlap <- intersect(intersect(sum.stat.discovery$SNP,names(snp_ref)),sum.stat.validation$SNP)
      if(length(snps.overlap) == 0){
        stop("There is no overlapping SNPs between discovery sample, validation sample and reference sample! Please check.")
      }
      rownames(sum.stat.discovery) <- sum.stat.discovery$SNP
      rownames(sum.stat.validation) <- sum.stat.validation$SNP
      sum.stat <- sum.stat.discovery[snps.overlap,]
      sum.stat.valid <- sum.stat.validation[snps.overlap,]
    }
    
    
    
    ##### Map reference allele (follow LD_ref at here, change back to ref in discovery sample at the end)#####
    
    LD_mat_save <- LD_ref[snps.overlap,snps.overlap]
    LD_mat_save[lower.tri(LD_mat_save,diag = T)] <- 0
    LD_use <- LD_mat_save + t(LD_mat_save)
    diag(LD_use) <- 1
    rownames(LD_use) <- colnames(LD_use) <- snps.overlap
    snp_ref_use <- snp_ref[snps.overlap]
    
    
    index <- sum.stat$A2 != snp_ref_use
    tmp <- sum.stat$A1[index]
    sum.stat$A1[index] <- sum.stat$A2[index]
    sum.stat$A2[index] <- tmp
    sum.stat$b[index] <- -sum.stat$b[index]
    
    betas_meta <-  sum.stat$b
    betas_se <-  sum.stat$se
    n.vec <- sum.stat$N
    
    p <- length(betas_meta)
    var.X <- v.y / n.vec / betas_se^2
    if(standardize == T){
      B <- LD_use
      #Xy <- betas_meta * v.y / betas_se^2 / sqrt(var.X) / n.vec
      Xy <- betas_meta * sqrt(v.y) / betas_se / sqrt(n.vec)
    }else{
      B <- diag(sqrt(var.X)) %*% LD_use %*% diag(sqrt(var.X))
      Xy <- betas_meta * v.y / betas_se^2 / n.vec
    }
    
    ## X'X
    lambda.v <- c()
    beta <- numeric(p)
    
    lambda <- max(abs(Xy))
    beta.mat <- sA.mat <- matrix(0,p,0)
    lambda.v <- c(lambda.v,lambda)
    j1 <- which.max(abs(Xy))
    A <- c(j1)
    nA <- (1:p)[-A]
    sA <- sign(Xy[A])
    sA.v <- numeric(p)
    sA.v[A] <- sA
    sA.mat <- cbind(sA.mat, sA.v)
    XaXa_inv <- solve(B[A,A])
    beta[A] <- XaXa_inv %*% (Xy[A] - lambda * sA)
    beta.mat <- cbind(beta.mat, beta)
    
    # Next hitting time
    XjXa <- B[nA,A]
    ##### Compute the hitting times #####
    
    
    while(lambda > 0 & length(A) < (nvar + 1)){
      
      # LASSO solution as lamb decreases from lamb_k
      
      temp <- XjXa %*% XaXa_inv
      
      posi1 <- (Xy[nA] - temp %*% Xy[A]) / (1 - temp %*% sA)
      nega1 <- (Xy[nA] - temp %*% Xy[A]) / (-1 - temp %*% sA)
      both <- cbind(posi1,nega1)
      rownames(both) <- nA
      both[which(is.na(both))] <- -Inf
      hit <- max(both[both < lambda - 1e-10])
      
      sign_j <- (-1)^(which(both == hit, arr.ind = TRUE)[2] - 1)  # 1 for 1st col, -1 for 2nd col
      ind <- nA[which(both == hit, arr.ind = TRUE)[1]]
      
      # Next crossing time
      
      cross_all <- (XaXa_inv %*% Xy[A]) / (XaXa_inv %*% sA)
      
      ind_cross <- which(cross_all < lambda - 1e-10)
      if(length(ind_cross) == 0){
        cross <- -Inf
      } else{
        cross <- max(cross_all[ind_cross])
      }
      ind2 <- which(cross_all == cross)
      
      # Update lambda
      lambda <- max(hit, cross)
      
      if (cross < hit){
        A <- c(A, ind)
        sA <- c(sA, sign_j)
        sA.v <- numeric(p)
        sA.v[A] <- sA
        sA.mat <- cbind(sA.mat, sA.v)
        nA <- (1:p)[-A]
      } else{
        beta[A[ind2]] <- 0
        A <- A[-ind2]
        sA <- sA[-ind2]
        sA.v <- numeric(p)
        sA.v[A] <- sA
        sA.mat <- cbind(sA.mat, sA.v)
        nA <- (1:p)[-A]
      }
      
      if(length(A) == p){
        lambda.v <- c(lambda.v,lambda)
        XaXa_inv <- solve(B[A,A])
        beta[A] <- XaXa_inv %*% (Xy[A] -lambda * sA)
        beta.mat <- cbind(beta.mat, beta)
        break
      }
      
      XaXa_inv <- try(solve(B[A,A]), silent = T)
      if(class(XaXa_inv) == "try-error"){
        A <- A[-length(A)]
        break
      }
      beta[A] <- XaXa_inv %*% (Xy[A] - lambda * sA)
      beta.mat <- cbind(beta.mat, beta)
      
      # Next hitting time
      XjXa <- B[nA,A]
      lambda.v <- c(lambda.v,lambda)
      
    }
    
    if(standardize == T){
      beta.mat <- diag(1 / sqrt(var.X)) %*% beta.mat 
    }
    
    ##### Use validation set to select optimum lambda #####
    
    if(!is.null(sum.stat.validation)){
      r2_sum <- function(beta_est, b_uni, se_uni, LD_ref, var.X, N){
        cov.y_hat.y <- crossprod(var.X * b_uni, beta_est)  ## cov(Xb_hat, y)
        var.y_hat <- t(sqrt(var.X) * beta_est) %*% LD_ref %*% (sqrt(var.X) * beta_est)  ## var(y_hat)
        var.y <- median(N*var.X*se_uni^2)  # var(y)
        return(cov.y_hat.y^2 / var.y_hat / var.y)
      }
      
      index2 <- sum.stat.valid$A2 != snp_ref_use
      sum.stat.valid$b[index2] <- -sum.stat.valid$b[index2]
      var.X.valid <- 2*sum.stat.valid$Freq1*(1-sum.stat.valid$Freq1)
      
      R2 <- numeric(ncol(beta.mat))
      for(i in 2:ncol(beta.mat)){
        R2[i] <- r2_sum(beta_est=beta.mat[,i], b_uni=sum.stat.valid$b, se_uni=sum.stat.valid$se, LD_ref=LD_use, var.X=var.X.valid, N=sum.stat.valid$N)
      }
      
    }
    
    if(is.na(lambda.vec)){
      beta.mat <- diag(1 - index*2, nrow = length(index)) %*% beta.mat
      rownames(beta.mat) <- rownames(LD_use)
      selected.markers <- rownames(beta.mat)[A]
      if(is.null(sum.stat.validation)){
        return(list(lambda.v = lambda.v, beta.mat = Matrix(beta.mat,sparse = TRUE), selected.markers = selected.markers[1:nvar]))
      } else{
        lambda.opt <- lambda.v[which.max(R2)]
        snps.opt <- which(abs(beta.mat[,which.max(R2)]) > 1e-10)
        beta.opt <- beta.mat[snps.opt,which.max(R2)]
        return(list(beta.opt = beta.opt, lambda.opt = lambda.opt, R2 = R2, lambda.v = lambda.v, beta.mat = Matrix(beta.mat,sparse = TRUE), selected.markers = selected.markers[1:nvar]))
      }
    }
    
    
    
    
    ##### Give a tuning parameter, get result #####
    
    lap <- function(lambda){
      if(lambda > max(lambda.v))
        return(numeric(p))
      if(lambda < lambda.v[length(lambda.v)] || lambda == lambda.v[length(lambda.v)])
        return(XaXa_inv %*% (Xy[A] - lambda * sA))
      k <- length(which(lambda < lambda.v))
      beta <- (beta.mat[,k+1]-beta.mat[,k])/(lambda.v[k+1]-lambda.v[k])*(lambda - lambda.v[k]) + beta.mat[,k]
      return(beta)
    }
    if(min(lambda.vec)< min(lambda.v))
      stop(paste("Too many variants will be selected. Please set a larger nvar or a larger lambda."))
    bm <- matrix(0,p,length(lambda.vec))
    for(i in 1:length(lambda.vec)){
      bm[,i] <- lap(lambda.vec[i])
    }
    bm <- diag(1 - index*2, nrow = length(index)) %*% bm
    rownames(bm) <- rownames(LD_use)
    return(list(lambda.v = lambda.vec, beta.mat = Matrix(bm,sparse = TRUE)))
  }