#My functions




# returns df with all dupes of specified variable
dup<-
  function(df,id){
    df[unlist(df[,id]) %in% unlist(df[duplicated(df[,id]),id]) ,]
  }


# takes lm object and outputs a OR/CI/p value table
cilm <-
  function(mod){
    coef <- summary(mod)$coefficients
    ci <- confint.default(mod)
    coef[,2:3] <- ci
    colnames(coef) <- c('Estimate', '95% LB', '95% UB','p value')
    
    return(round(coef,3))
    
  }


# takes glm object (logistic regression model) and outputs a OR/CI/p value table
cilog <-
  function(mod){
    sm <- summary(mod)
    coef <- sm$coefficients[,1]
    p <- sm$coefficients[,4]
    ci <- confint.default(mod)
    
    ctab <- cbind(coef,ci,p)
    ctab[,1:3] <- exp(ctab[,1:3])
    
    colnames(ctab)[1] <- 'OR'
    
    return(round(ctab,3))
    
  }





# takes a multinom model object (pacckage nnet)
# --  creates OR table with CIs and p values

cimn <- function(mod){
  
  #coefficients
  coef <- summary(mod)$coefficients %>% t %>% as.data.frame
  coef$param <- rownames(coef)
  coefb <- coef%>% pivot_longer(cols=1:(ncol(coef)-1), names_to = 'level', values_to = 'OR')
  
  #p values
  z <- summary(mod)$coefficients/summary(mod)$standard.errors
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  pa <- p %>% t %>% as.data.frame
  pa$param <- rownames(pa)
  pb <- pa%>% pivot_longer(cols=1:(ncol(pa)-1), names_to = 'level', values_to = 'p.value')
  
  #CIs
  lb <- summary(mod)$coefficients - summary(mod)$standard.errors *1.959964
  ub <- summary(mod)$coefficients + summary(mod)$standard.errors *1.959964
  
  lba <- lb %>% t %>% as.data.frame
  uba <- ub %>% t %>% as.data.frame
  
  lba$param <- row.names(lba)
  uba$param <- row.names(uba)
  
  lbb <- lba%>% pivot_longer(cols=1:(ncol(lba)-1), names_to = 'level', values_to = '2.5%')
  ubb <- uba%>% pivot_longer(cols=1:(ncol(uba)-1), names_to = 'level', values_to = '97.5%')
  
  ci <- merge(lbb,ubb, by=c('param','level'))
  
  
  #all together
  
  coci <- merge(coefb, ci, by=c('param','level'))
  tab <- merge(coci, pb, by=c('param','level'))
  
  tab[,3:5] <- exp(tab[,3:5])
  
  tab[,3:6] <- round(tab[,3:6],3)
  
  return(tab)
}



# You can compute the high leverage observation by looking at the ratio of number
# of parameters estimated in model and sample size. If an observation has a ratio
# greater than 2 -3 times the average ratio, then the observation considers as 
# high-leverage points. I personally like to use this simple function to identify 
# high-leverage observations.
# From https://towardsdatascience.com/how-to-detect-unusual-observations-on-your-regression-model-with-r-de0eaa38bc5b

highleverage <- 
  function(fit){
    p <- length(coefficients(fit))
    n <- length(fitted(fit))
    ratio <-p/n
    plot(hatvalues(fit), main="Index Plot of Ratio")
    abline(h=c(2,3)*ratio, col="red", lty=2)
    text(1:n, hatvalues(fit), rownames(fit$model),adj=1.2)
  }


############################# modcsv ###########################################
#currently only works for lm and glm binomial, and multinom models
modcsv<-
  function(list, file){
    
    for(i in 1:length(list)){
      
      #glm models have class object that is length 2
      if(length(class(list[[i]])) > 1){
        if(class(list[[i]])[1]=='glm'){
          
          #use function cilog to get ORs for logistic model
          if(list[[i]]$family$family =='binomial'){
            
            write.table(data.frame(t(rep('----------------------------------------',5))), file, sep=',', col.names = F, row.names=F, append=T)
            write.table(data.frame(t(c('Outcome:',names(list)[[i]]))), file,row.names=F,col.names=F,append=T)
            write.table(data.frame(t(c('',colnames(cilog(list[[i]]))))), file, sep=',',col.names=F,row.names=F, append=T)
            write.table(cilog(list[[i]])%>%round(3), file,col.names=F, sep=',', append=T)
          }
        }
        
        
        #use function cimn to get ORs for multinomial models
        if(class(list[[i]])[1]=='multinom'){
          
          write.table(data.frame(t(rep('------------------------------------------',5))), file, sep=',', col.names = F, row.names=F, append=T)
          write.table(data.frame(t(c('Outcome:',names(list)[[i]]))), file,row.names=F,col.names=F,append=T)
          write.table(cimn(list[[i]]), file, row.names=F, sep=',', append=T)
          
        }
        
      }
      
      #lm models have class object that is length 1
      if(length(class(list[[i]])) == 1){
        if(class(list[[i]])=='lm'){
          
          write.table(data.frame(t(rep('------------------------------------------',5))), file, sep=',', col.names = F, row.names=F, append=T)
          write.table(data.frame(t(c('Outcome:',names(list)[[i]]))), file,row.names=F,col.names=F,append=T)
          write.table(data.frame(t(c('',colnames(cilm(list[[i]]))))), file, sep=',',col.names=F,row.names=F, append=T)
          write.table(cilm(list[[i]])%>%round(3), file,col.names=F, sep=',', append=T)
        }
      }
      
    }
    
    
  }




# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut] %>%round(3),
    p = pmat[ut] %>% round(3)
  )
}