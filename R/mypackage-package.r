#' mypackage
#'
#' @name mypackage
#' @docType package
#' @export cross.cor
#' @export corstarsl
#' @export multilevel_parameters
#' @import Hmisc pander

cross.cor <- function(x, y, type = "pearson") {
  Rnew <- corstarsl(data.frame(x, y), type)
  # x - number of columns
  # y - number of rows
  Rnew.ncol <- 1 + length(x)
  
  Rnew.nrow <- 1 + length(y)
  
  Rnew <- Rnew[c(Rnew.ncol:nrow(Rnew)),c(1:Rnew.nrow)]
  return(Rnew)
}

corstarsl <- function(x, type = "pearson", pvalue = "stars"){ 
  #library(Hmisc) 
  x <- as.matrix(x) 
  R <- Hmisc::rcorr(x, type = type)$r 
  p <- Hmisc::rcorr(x, type = type)$P 
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1] 
  p <- format(round(cbind(rep(-1.11, ncol(x)), p), 2))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  if (pvalue == "stars") {
    Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  }
  if (pvalue == "values") {
    Rnew <- matrix(paste(R, ", (_p_=", p,")", sep=""), ncol=ncol(x)) 
  }
  
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  row.names(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew) 
  
  ## remove last column and return the matrix (which is now a data frame)
  #Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew) 
} 

multilevel_parameters=function(modlist) {
  library(lme4)
  # May be working with (2 group varibles with one varience):
  # nullmodel <- lmer(score ~ cohort90 + female + sclass + (1 + cohort90 | schoolid) + schtype + schurban, data = mydata, REML = FALSE)
  # Didn't work, need to be adopted
  # nullmodel <- lmer(score ~ -1 + (1 | schoolid) + (1 | cohort90) + caseid, data = mydata, REML = FALSE)
  if (class(modlist) != "list") {
    modlist <- list(modlist)
  }
  
  modclass=unlist(lapply(modlist,class))
  ifelse(
    all(modclass[1]==modclass)==FALSE,return("Error: Objects in list need to all be of the same class"),modclass )
  if(modclass[1]=="mer") { #For models fit using lmer
    # Get formulas
    VarFormula_short.list=lapply(modlist,function(i) {
      formula_array <- deparse(formula(i), width = 30L)
      VarFormula_short = formula_array[[1]]
    } )
    
    #Get variance of fixed effects by multiplying coefficients by design matrix
    VarF.list=lapply(modlist,function(i) varF=var(as.vector(lme4::fixef(i) %*% t(i@X))) )
    
    #Get variance of random effects by extracting variance components
    # In our term this means Intercept variance ( between-school (level 2))
    VarRand.list=lapply(modlist,function(i) do.call(rbind,lapply(lme4::VarCorr(i),function(j) j[1])) )
    
    #Get residual variance
    # In our term this means residual Varience (within-school between-student (level 1))
    VarResid.list=lapply(modlist,function(i) attr(lme4::VarCorr(i), "sc")^2 )
    
    # To calcaulate VPC we first need total varience (as sum of variences), then devide Intercept varience on Total varience
    # 1. Total Varience
    TotalVarience = lapply(seq_along(modlist),function(i) (colSums(VarRand.list[[i]])+VarResid.list[[i]]) )
    # 2. VPC
    VPC = lapply(seq_along(modlist),function(i) (VarRand.list[[i]]/TotalVarience[[i]]) )
    #Calculate marginal R-squared (fixed effects/total variance)
    Rm.list=lapply(seq_along(modlist),function(i) VarF.list[[i]]/(VarF.list[[i]]+colSums(VarRand.list[[i]])+VarResid.list[[i]]) )
    #Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc.list=lapply(seq_along(modlist),function(i) (VarF.list[[i]]+colSums(VarRand.list[[i]]))/(VarF.list[[i]]+colSums(VarRand.list[[i]])+VarResid.list[[i]]) )
    #Bind R^2s into a matrix and return
    Rsquared.mat=cbind(VarFormula_short.list, Rm.list,Rc.list, VPC);
    colnames(Rsquared.mat)=c("Formula","Marginal","Conditional", "VPC or ICC1 or t00 ")
    return(Rsquared.mat)
  }
}

vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

rsquared.lme=function(modlist) {  
  modclass=unlist(lapply(modlist,class))
  ifelse(
    all(modclass[1]==modclass)==FALSE,return("Error: Objects in list need to all be of the same class"),modclass )
  if(modclass[1]=="mer" | modclass[1] == "lmerMod") { #For models fit using lmer
    #Get variance of fixed effects by multiplying coefficients by design matrix
    VarF.list=lapply(modlist,function(i) varF=var(as.vector(fixef(i) %*% t(i@X))) )
    #Get variance of random effects by extracting variance components
    VarRand.list=lapply(modlist,function(i) do.call(rbind,lapply(VarCorr(i),function(j) j[1])) )
    #Get residual variance
    VarResid.list=lapply(modlist,function(i) attr(VarCorr(i), "sc")^2 )
    #Calculate marginal R-squared (fixed effects/total variance)
    Rm.list=lapply(seq_along(modlist),function(i) VarF.list[[i]]/(VarF.list[[i]]+colSums(VarRand.list[[i]])+VarResid.list[[i]]) )
    #Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc.list=lapply(seq_along(modlist),function(i) (VarF.list[[i]]+colSums(VarRand.list[[i]]))/(VarF.list[[i]]+colSums(VarRand.list[[i]])+VarResid.list[[i]]) )
    #Bind R^2s into a matrix and return
    Rsquared.mat=do.call(cbind,list(Rm.list,Rc.list)); colnames(Rsquared.mat)=c("Marginal","Conditional")
    return(Rsquared.mat)
  }    
  if(modclass[1]=="lme") {#For models fit using lme
    #Get design matrix of fixed effects from model
    Fmat.list=lapply(modlist,function(i) model.matrix(eval(i$call$fixed)[-2],i$data) )
    #Get variance of fixed effects by multiplying coefficients by design matrix
    VarF.list=lapply(seq_along(modlist),function(i) var(as.vector(fixef(modlist[[i]]) %*% t(Fmat.list[[i]]))) )
    #Get variance of random effects by extracting variance components
    VarComp.list=lapply(modlist,function(i) VarCorr(i) )
    VarRand.list=lapply(VarComp.list,function(i) as.numeric(i[rownames(i)!="Residual" & rownames(i)=="(Intercept)","Variance"]) )
    #Get residual variance
    VarResid.list=lapply(VarComp.list,function(i) as.numeric(i[rownames(i)=="Residual","Variance"]) )
    #Calculate marginal R-squared (fixed effects/total variance)
    Rm.list=lapply(seq_along(modlist),function(i) VarF.list[[i]]/(VarF.list[[i]]+sum(VarRand.list[[i]])+VarResid.list[[i]]) )
    Rc.list=lapply(seq_along(modlist),function(i) (VarF.list[[i]]+ifelse(length(VarRand.list[[i]][1])==1,VarRand.list[[i]][1],colSums(VarRand.list[[i]])))/(VarF.list[[i]]+ifelse(length(VarRand.list[[i]][1])==1,VarRand.list[[i]][1],colSums(VarRand.list[[i]]))+VarResid.list[[i]]) )
    #Bind R^2s into a matrix and return
    Rsquared.mat=do.call(cbind,list(Rm.list,Rc.list)); colnames(Rsquared.mat)=c("Marginal","Conditional")
    return(Rsquared.mat)
    #Or return error if objects are not linear mixed effects models
  } else { print("Error: Function requires list of objects of class 'mer' or 'lme'") }
} 

#' @export xtables_html
xtables_html <- function (table, pander = TRUE, ...) {
  if (pander == TRUE) {    
    result <- pandoc.table(table, style = "rmarkdown")
    return (result)
  }
  else {
    result <- print.xtable(xtable(table, ...), type="html")
  }
}

plot.cohenD <- function (var1, var2) {
  # Standardized Mean Difference (Cohen's d)
  ES <- 0.8
  # get mean2 depending on value of ES from d = (u1 - u2)/sd
  mean1 <- ES*1 + 1
  # create x sequence
  x <- seq(1 - 3*1, mean1 + 3*1, .01)
  # generate normal dist #1
  y1 <- dnorm(x, 1, 1)
  # put in data frame
  df1 <- data.frame("x" = x, "y" = y1)
  # generate normal dist #2
  y2 <- dnorm(x, mean1, 1)
  # put in data frame
  df2 <- data.frame("x" = x, "y" = y2)
  # get y values under overlap
  y.poly <- pmin(y1,y2)
  # put in data frame
  poly <- data.frame("x" = x, "y" = y.poly)
  
  # Cohen's U3, proportion of control > 50th perc. treatment
  u3 <- 1 - pnorm(1, mean1,1)
  u3 <- round(u3,3)
  
  # plot with ggplot2
  ggplot(df1, aes(x,y, color="treatment")) +
    # add line for treatment group
    geom_line(size=1) + 
    # add line for control group
    geom_line(data=df2, aes(color="control"),size=1) +
    # shade overlap
    geom_polygon(aes(color=NULL), data=poly, fill="red", alpha=I(4/10),
                 show_guide=F) +
    # add vlines for group means
    geom_vline(xintercept = 1, linetype="dotted") + 
    geom_vline(xintercept = mean1, linetype="dotted") + 
    # add plot title
    opts(title=paste("Visualizing Effect Sizes 
                     (Cohen's d = ",ES,"; U3 = ",u3,")", sep="")) +
    # change colors and legend annotation
    scale_color_manual("Group", 
                       values= c("treatment" = "black","control" = "red")) +
    # remove axis labels
    ylab(NULL) + xlab(NULL)
}

#' @export z0
z0 <- function(x) {
  return (round(x,0))
}

#' @export p1
p1 <- function(x) {
  return (round(x,0))
}

#' @export lmp
lmp <- function (modelobject) {
  if (class(modelobject) == "lm") {
    f <- summary(modelobject)$fstatistic
  }
  else if (class(modelobject) == "summary.lm") {
    f <- modelobject$fstatistic
  } else {
    stop("Not an object of class 'lm' ")    
  }
  p <- unname(pf(f[1],f[2],f[3],lower.tail=F))  #attributes(p) <- NULL
  return(p)
}

#' @export pv
pv <- function (x, digits=2, eps = 0.001, prefix = TRUE) {
  if (x <= eps) {
    formatted = format.pval(x, digits=digits, eps = eps)
  } else {
    formatted = paste("=", round(x, digits), sep="")
  }
  if (prefix == TRUE) {
    return (formatted)    
  } else {
    return (substr(formatted,2,8))
  }
  
} 

#' @export x_pvalue
x_pvalue <- function (x, pv, format="pvalues") {
  #x <- trunc(x, 4)
  #pv <- trunc(pv, 4)
  if (format == "stars") {
    mystars <- ifelse(pv <= .001, "***", ifelse(pv <= .01, "**", ifelse(pv <= .05, "*", "")))
    result <- paste(x, mystars, sep="")
  }
  if (format == "pvalues") {
    result <- paste(x, ", (p=", pv(pv),")", sep="")
  }
  if (format == "values") {
    result <- paste(x, ", (p=", pv,")", sep="")
  }
  
  if (format == "significant_values") {
    if (pv < 0.05) {
      result <- paste(x, ", (p=", pv,")", sep="")
    } else {
      result <- x
    }
  }
  return(result)
}

#' @export regression_table
regression_table <- function(scales, predictors, data) {
  
Rresult <- data.frame()
for (output in scales) {
  
  model <- lm(as.formula(paste(output," ~ ", paste(predictors, collapse = " + "))), data = data)
  coeff <- summary(model)$coefficients
  model_coefficients <- as.data.frame(t(summary(model)$coefficients))
  model_beta_coefficients <- as.data.frame(t(QuantPsyc::lm.beta(model)))
  table_header <- c(names(model_coefficients), "R2", "F")
  
  Rnew <- data.frame(rbind(table_header))
  
  colnames(Rnew) <- table_header
  for (i in 1:length(model_coefficients)) {
    value <- x_pvalue(round(model_coefficients[i][,1][1],2), model_coefficients[i][,1][4])
    Rnew[table_header[i]] <- value
  }
  Rnew["R2"] <- round(summary(model)$r.squared, 2)
  Rnew["F"] <- paste0("F(",summary(model)$fstatistic["numdf"], ",", summary(model)$fstatistic["dendf"], ")=", x_pvalue(round(summary(model)$fstatistic["value"],2),pv(lmp(model))))

  
  rownames(Rnew) <- output
  Rresult <- rbind(Rresult, Rnew)
  }
return (Rresult)
}

#' @export regression_table2
regression_table2 <- function(scales, predictors, data) {
  
  Rresult <- data.frame()
  for (output in scales) {
    
    model <- lm(as.formula(paste(output," ~ ", paste(predictors, collapse = " + "))), data = data)
    coeff <- as.data.frame(summary(model)$coefficients)
    model_beta_coefficients <- as.data.frame(t(QuantPsyc::lm.beta(model)))
    
    #Standart <- t(cbind(NA,Intercept = model_beta_coefficients))
    coeff <- cbind(coeff, Standart = NA)
    for (row in rownames(coeff)) {
      if (row == "(Intercept)") {
        coeff["(Intercept)","Standart"] <- NA
      } else {
        coeff[row,"Standart"] <- model_beta_coefficients[row]
      }
      
    }
    coeff <- as.data.frame(t(coeff))
    
    Rnew <- data.frame()  
    for (index in 1:length(coeff)) {
      i <- coeff[index]
      cell <- data.frame()
      b <-  x_pvalue(round(t(i)[1,"Estimate"],2), t(i)[1,"Pr(>|t|)"])
      #if (t(i)[1,"Standart"] != NA
      b_std <- x_pvalue(round(t(i)[1,"Standart"],2), t(i)[1,"Pr(>|t|)"])
      SE <- round(t(i)[1,"Std. Error"],2)
      t <- round(t(i)[1,"t value"],2)
      Rnew[names(coeff)[index],output] <- paste0("b=",b," _b_=",b_std," SE=", SE," _t_=",t)
      
    }
    Rnew <- as.data.frame(t(Rnew))
    Rnew["R2"] <- round(summary(model)$r.squared, 2)
    Rnew["F"] <- paste0(
      x_pvalue(round(summary(model)$fstatistic["value"],2),pv(lmp(model))),
      " _df_=(",summary(model)$fstatistic["numdf"], ",", summary(model)$fstatistic["dendf"], ")")
    
    
    #rownames(Rnew) <- output
    #Rresult <- Rnew
    Rresult <- rbind(Rresult, Rnew)
  }
  return (Rresult)
}

#' @export regression_table3
regression_table3 <- function(output, predictors, data) {
  
  model <- lm(as.formula(paste(output," ~ ", paste(predictors, collapse = " + "))), data = data)
  coeff <- as.data.frame(summary(model)$coefficients)
  model_beta_coefficients <- as.data.frame(t(QuantPsyc::lm.beta(model)))
  
  #Standart <- t(cbind(NA,Intercept = model_beta_coefficients))
  coeff <- cbind(coeff, Standart = NA)
  for (row in rownames(coeff)) {
    if (row == "(Intercept)") {
      coeff["(Intercept)","Standart"] <- NA
    } else {
      coeff[row,"Standart"] <- model_beta_coefficients[row]
    }
    
  }
  coeff <- as.data.frame(t(coeff))
  
  Rnew <- data.frame()  
  for (index in 1:length(coeff)) {
    #Rnew[names(coeff)[index],output] <- 1
    i <- coeff[index]
    cell <- data.frame()
    
    cell["b","Value"] <-  x_pvalue(round(t(i)[1,"Estimate"],2), t(i)[1,"Pr(>|t|)"])
    #if (t(i)[1,"Standart"] != NA
    cell["b_std","Value"] <- x_pvalue(round(t(i)[1,"Standart"],2), t(i)[1,"Pr(>|t|)"])
    cell["SE","Value"] <- round(t(i)[1,"Std. Error"],2)
    cell["t","Value"] <- round(t(i)[1,"t value"],2)
    cell <- cbind(Predictor=c(names(coeff)[index], "", "", ""), Variable = c("b","$\beta$","_SE_","_t_"), cell)
    Rnew <- rbind(Rnew, cell)
  }
  #Rnew <- as.data.frame(t(Rnew))
  #rbind(Rnew, R2 = c(0))
  #Rnew["R2"] <- 2
  Rnew <- rbind(Rnew, R2 = cbind(Predictor="", Variable = "$R^2$", Value = round(summary(model)$r.squared, 2) ))
  
  Rnew <- rbind(Rnew, F = cbind(Predictor="", Variable = "$F$", Value = 
                                  paste0(x_pvalue(round(summary(model)$fstatistic["value"],2),pv(lmp(model))),
                                         " _df_=(",summary(model)$fstatistic["numdf"], ",", summary(model)$fstatistic["dendf"], ")")))
  return (Rnew)
}

#' @export regression_multi_table
regression_multi_table <- function(output, predictors, data) {
  Rresult <- regression_table3(output[1], predictors, data)[, c("Predictor","Variable")]
  for (scale in output) {
    Rtemp <- subset(regression_table3(scale, predictors, data), select = "Value")
    colnames(Rtemp) <- scale
    Rresult <- cbind(Rresult,Rtemp)
  }
  rownames(Rresult) <- NULL
  return (Rresult)
}

#' @export rev
rev <- function(input, length) {
  fun <- function(x) {
    return (abs(x - (length+1)))
  }
  sapply(input, fun)
}
