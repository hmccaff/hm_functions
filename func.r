#My functions




# returns df with all dupes of specified variable
dup<-
  function(df,id){
    df[unlist(df[,id]) %in% unlist(df[duplicated(df[,id]),id]) ,]
  }


# takes lm object and outputs a OR/CI/p value table
# ******* deprecated *************
cilm <-
  function(mod){
    coef <- summary(mod)$coefficients
    ci <- confint.default(mod)
    coef[,2:3] <- ci
    colnames(coef) <- c('Estimate', '95% LB', '95% UB','p value')
    
    return(round(coef,3))
    
  }


# takes glm object (logistic regression model) and outputs a OR/CI/p value table
# ******* deprecated *************
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



#format p values - used in coeftab and tabsum
fmtp <- 
  function(pvec){
    pvec <- as.numeric(pvec)
    sapply(pvec,
           function(x){
             if(is.na(x)){
               return('')
             }else if(round(x,3)==0){
               return(paste('**',"<0.001",'**',sep=''))
             }
             else if(x<.05){
               return(paste('**',sprintf(x,fmt="%.3f"),'**',sep=''))
             }else{
               return(sprintf(x,fmt="%.3f"))
             }
           }
    )
  }


#for arsenal table summary - generate table summary and format p values
#for an arsenal table objects
tabsum <-
  function(tab){
   st<- summary(tab,
            text=T,
            labelTranslations=labs
            )%>%
      as.data.frame()
   
   
    st$`p value` <- fmtp(st$`p value`)
    
    st
  }
  

#takes a table of joint tests from emmeans and cleans it up
jttab <-
  function(tab){
    newtab<-
      tab%>%
      mutate(
        `model term` = paste('**',`model term`,'**',sep=''),
        F.ratio = sprintf(F.ratio,fmt="%.3f"),
        Chisq = sprintf(Chisq,fmt="%.3f"),
        p.value = fmtp(p.value)
        
      )
    
    rownames(newtab) <- newtab$`model term`
    newtab[,-1]
  }



#takes a model object and outputs a coefficient table with basic information
coeftab <-
  function(fit,cimethod='profile',type='response'){
    
    #cox models
    if('coxph' %in% class(fit)){
      est <- summary(fit)$conf.int[,1] %>% sprintf(fmt="%.2f")
      lb <- summary(fit)$conf.int[,3] %>% sprintf(fmt="%.2f")
      ub <- summary(fit)$conf.int[,4] %>% sprintf(fmt="%.2f")
      pval <- summary(fit)$coefficients[,5] %>% fmtp
      
      df <- 
        data.frame(
          `HR`=est,
          `confint`= paste('(',lb,',',ub,')',sep =''),
          `pval`=pval
        )
      rownames(df) <- names(fit$coefficients)
      
      bot<-
        data.frame(
          `HR`= c('',
                        round(fit$n,0)),
          `confint`='',
          `pval`=''
        )
      rownames(bot) <- c('---',
                         'N')
      
      df <- rbind(df,bot)
    }
    
    #lmer objects, linear mixed models
    if('lmerModLmerTest' %in% class(fit)){
      
      #this class is found with lmer objects
      if(type=='lp'){
        #point estimates -remove intercept, keep var components
        est <- 
          c(fit@beta)%>%
          sprintf(fmt="%.3f")
        
        ci <-
          confint(fit,method=cimethod,parm=rownames(summary(fit)$coefficients))%>%
          as.data.frame()%>%
          mutate(across(c(1,2),function(x){sprintf(fmt="%.3f",x)}))
        
        #pval
        pval <- 
          summary(fit)$coefficients[,5]%>%
          fmtp
        
        #on response scale with fixed effects only
      }else if(type=='response'){
        #point estimates -remove intercept, keep var components
        est <- 
          fit@beta[-1]%>%
          sprintf(fmt="%.2f")
        
        ci <-
          confint(fit,method=cimethod,
                  parm=rownames(summary(fit)$coefficients)[-1])%>%
          as.data.frame()%>%
          mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
        
        #pval -remove intercept
        pval <- 
          summary(fit)$coefficients[-1,5]%>%
          fmtp
      }
      
      if(is.null(dim(ci))){
        df <- 
          data.frame(
            `Estimate`=est,
            `confint`= paste('(',ci[1],',',ci[2],')',sep =''),
            `pval`=pval
          )
      }
      else{
        df <- 
          data.frame(
            `Estimate`=est,
            `confint`= paste('(',ci[,1],',',ci[,2],')',sep =''),
            `pval`=pval
          )
      }
      
      #add random effects to bottom of table
      vc <- as.data.frame(VarCorr(fit))
      #random effects std deviations
      ranef <- vc[is.na(vc$var2) & vc$grp != 'Residual',]
      mid <- 
        data.frame(Estimate=c('','',sprintf(ranef$sdcor,fmt="%.3f")),
                   confint='',
                   pval=''
        )
      rn<-paste(ranef$grp,'|',ranef$var1,' SD',sep='')
      rn <- c('------','variance components:',rn)
      rownames(mid) <- rn
      
      #random effects correlations
      ranef <- vc[!is.na(vc$var2),]
      if(nrow(ranef)>0){
        mid2<-
          data.frame(Estimate=sprintf(ranef$sdcor,fmt="%.3f"),
                     confint='',
                     pval=''
          )
        rn <- paste(ranef$grp,'|',ranef$var1,'*',ranef$var2,' Corr',sep='')
        rownames(mid2) <- rn
        
        mid <- rbind(mid,mid2)
      }
      #residual SD last
      resid <- vc[vc$grp == 'Residual',]
      mid3 <- 
        data.frame(Estimate=sprintf(resid$sdcor,fmt="%.3f"),
                   confint='',
                   pval=''
        )
      rn <- 'Residual SD'
      rownames(mid3) <- rn
      mid <- rbind(mid,mid3)
      
      #add N and fit statistics to bottom of table
      bot<-
        data.frame(
          `Estimate`= c('',round(nrow(fit@frame),0),
                        sprintf(summary(fit)$AICtab[1],fmt="%.2f")),
          `confint`='',
          `pval`=''
        )
      rownames(bot) <- c('---',
                         'N','AIC')
      
      df <- rbind(df,mid,bot)
      
      
      
    }
    
    if('glmerMod' %in% class(fit)){
      
      #on linear predictor scale with var comps
      if(type=='lp'){
        #point estimates -remove intercept, keep var components
        est <- 
          c(fit@beta)%>%
          sprintf(fmt="%.3f")
        
        ci <-
          confint(fit,method=cimethod,parm=rownames(summary(fit)$coefficients))%>%
          as.data.frame()%>%
          mutate(across(c(1,2),function(x){sprintf(fmt="%.3f",x)}))
        
        #pval
        pval <- 
          summary(fit)$coefficients[,4]%>%
          fmtp
        
        #on response scale with fixed effects only
      }else if(type=='response'){
        #point estimates -remove intercept, keep var components
        est <- 
          exp(fit@beta[-1])%>%
          sprintf(fmt="%.2f")
        
        ci <-
          confint(fit,method=cimethod,
                  parm=rownames(summary(fit)$coefficients)[-1])%>%
          as.data.frame()%>%
          exp%>%
          mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
        
        #pval -remove intercept
        pval <- 
          summary(fit)$coefficients[-1,4]%>%
          fmtp
      }
      
      if(is.null(dim(ci))){
        df <- 
          data.frame(
            `Estimate`=est,
            `confint`= paste('(',ci[1],',',ci[2],')',sep =''),
            `pval`=pval
          )
      }
      else{
        df <- 
          data.frame(
            `Estimate`=est,
            `confint`= paste('(',ci[,1],',',ci[,2],')',sep =''),
            `pval`=pval
          )
        rownames(df) <- rownames(ci)
      }
      
      #add random effects to bottom of table
      vc <- as.data.frame(VarCorr(fit))
      #random effects std deviations
      ranef <- vc[is.na(vc$var2),]
      mid <- 
        data.frame(Estimate=c('','',sprintf(ranef$sdcor,fmt="%.3f")),
                   confint='',
                   pval=''
        )
      rn<-paste(ranef$grp,'|',ranef$var1,' SD',sep='')
      rn <- c('------','variance components:',rn)
      rownames(mid) <- rn
      
      #random effects correlations
      ranef <- vc[!is.na(vc$var2),]
      if(nrow(ranef)>0){
        mid2<-
          data.frame(Estimate=sprintf(ranef$sdcor,fmt="%.3f"),
                     confint='',
                     pval=''
          )
        rn <- paste(ranef$grp,'|',ranef$var1,'*',ranef$var2,' Corr',sep='')
        rownames(mid2) <- rn
        
        mid <- rbind(mid,mid2)
      }
      
      #add N and fit statistics to bottom of table
      bot<-
        data.frame(
          `Estimate`= c('',round(nrow(fit@frame),0),
                        sprintf(summary(fit)$AICtab[1],fmt="%.2f"),
                        sprintf(summary(fit)$AICtab[2],fmt="%.2f")),
          `confint`='',
          `pval`=''
        )
      rownames(bot) <- c('---',
                         'N','AIC','BIC')
      
      df <- rbind(df,mid,bot)
      
      if(type=='response'){
        names(df)[1] <- 'OR'
      }
      
      
      
    }
    
    if('gee' %in% class(fit)){
      
      #on linear predictor scale with var comps
      if(type=='lp'){
        #point estimates -remove intercept
        est <- 
          summary(fit)$coefficients[-1,1]%>%
          sprintf(fmt="%.3f") 
        pt <- summary(fit)$coefficients[-1,1]
        robse <- summary(fit)$coefficients[-1,4]
        lb <- (pt - 1.959964*robse) %>% sprintf(fmt="%.3f")
        ub <- (pt + 1.959964*robse)%>% sprintf(fmt="%.3f")
        ci <- data.frame(lb=lb,ub=ub)
        
        #pval -remove intercept
        z<- summary(fit)$coefficients[-1,5]
        p <- pnorm(-abs(z))*2
        
        pval <-  fmtp(p)
        
        
        #on response scale 
      }else if(type=='exp'){
        #point estimates -remove intercept
        
          
          est <- 
            summary(fit)$coefficients[-1,1]%>%
            exp%>%
            sprintf(fmt="%.2f") 
          pt <- summary(fit)$coefficients[-1,1]
          robse <- summary(fit)$coefficients[-1,4]
          lb <- (pt - 1.959964*robse) %>% exp %>% sprintf(fmt="%.2f")
          ub <- (pt + 1.959964*robse) %>% exp %>% sprintf(fmt="%.2f")
          ci <- data.frame(lb=lb,ub=ub)
          
          #pval -remove intercept
          z<- summary(fit)$coefficients[-1,5]
          p <- pnorm(-abs(z))*2
          
          pval <-  fmtp(p)
          
      }else if(type=='response'){
        #point estimates -remove intercept
        if(fit$family[[1]]=='gaussian'){
          # do not exponentiate if response is gaussian
          est <- 
            summary(fit)$coefficients[-1,1]%>%
            sprintf(fmt="%.2f") 
          pt <- summary(fit)$coefficients[-1,1]
          robse <- summary(fit)$coefficients[-1,4]
          lb <- (pt - 1.959964*robse) %>% sprintf(fmt="%.2f")
          ub <- (pt + 1.959964*robse) %>% sprintf(fmt="%.2f")
          ci <- data.frame(lb=lb,ub=ub)
          
          #pval -remove intercept
          z<- summary(fit)$coefficients[-1,5]
          p <- pnorm(-abs(z))*2
          
          pval <-  fmtp(p)
        }else{
          
          est <- 
            summary(fit)$coefficients[-1,1]%>%
            exp%>%
            sprintf(fmt="%.2f") 
          pt <- summary(fit)$coefficients[-1,1]
          robse <- summary(fit)$coefficients[-1,4]
          lb <- (pt - 1.959964*robse) %>% exp %>% sprintf(fmt="%.2f")
          ub <- (pt + 1.959964*robse) %>% exp %>% sprintf(fmt="%.2f")
          ci <- data.frame(lb=lb,ub=ub)
          
          #pval -remove intercept
          z<- summary(fit)$coefficients[-1,5]
          p <- pnorm(-abs(z))*2
          
          pval <-  fmtp(p)
        }
      }
      
      if(is.null(dim(ci))){
        df <- 
          data.frame(
            `Estimate`=est,
            `confint`= paste('(',ci[1],',',ci[2],')',sep =''),
            `pval`=pval
          )
      }
      else{
        df <- 
          data.frame(
            `Estimate`=est,
            `confint`= paste('(',ci[,1],',',ci[,2],')',sep =''),
            `pval`=pval
          )
        rownames(df) <- names(p)
      }
      #add N and fit statistics to bottom of table
      
      #working corr NA if max clust size =1
      if(nrow(fit$working.correlation)>1){wc <- sprintf(fit$working.correlation[2,1],fmt="%.3f")}
      if(nrow(fit$working.correlation)==1){wc <- as.character(NA)}
      
      bot<-
        data.frame(
          `Estimate`= c('',
                        round(fit$nobs,0),length(unique(fit$id)),
                        as.character(fit$call[3]),
                        fit$family[[1]],
                        summary(fit)$model$corstr,
                        wc),
          `confint`='',
          `pval`=''
        )
      rownames(bot) <- c('---',
                         'N','N groups','grouping var.','dist.',
                         'corstr.','working correlation')
      
      df <- rbind(df,bot)
      
      if(type=='response'){
        if(fit$family[[1]]=='gaussian'){
          names(df)[1] <- 'Estimate'
        }
        if(fit$family[[1]]=='binomial'){
          names(df)[1] <- 'OR'
        }
        if(fit$family[[1]]=='poisson'){
          names(df)[1] <- 'IRR'
        }
        #assuming log link is used for gamma
        if(fit$family[[1]]=='Gamma'){
          names(df)[1] <- 'exp(beta)'
        }
      }
      
      
    }
    
    else if('glm' %in% class(fit)){
      
      #on linear predictor scale 
      if(type=='lp'){
        #point estimates 
        est <- summary(fit)$coefficients[,1]%>% sprintf(fmt="%.3f")
        
        ci <- 
          confint(fit,method=cimethod)%>%
          as.data.frame()%>%
          mutate(across(c(1,2),function(x){sprintf(fmt="%.3f",x)}))
        
        #format p values
        pval <- summary(fit)$coefficients[,4]%>%
          fmtp
        
        #exponentiated
      }else if(type=='exp'){
          #remove intercept
          est <- 
            summary(fit)$coefficients[-1,1]%>%
            exp%>% 
            sprintf(fmt="%.2f")
          
          ci <- 
            confint(fit,method=cimethod)[-1,] %>% 
            as.data.frame()%>%
            exp%>%
            mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
          
          #pval
          pval <-  
            summary(fit)$coefficients[-1,4]%>% 
            fmtp
          
        }else if(type=='response'){
          
          if(fit$family[1]=='binomial'){
            #remove intercept
            est <- 
              summary(fit)$coefficients[-1,1]%>%
              exp%>% 
              sprintf(fmt="%.2f")
            
            ci <- 
              confint(fit,method=cimethod)[-1,] %>% 
              as.data.frame()%>%
              exp%>%
              mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
            
            #pval
            pval <-  
              summary(fit)$coefficients[-1,4]%>% 
              fmtp
          }else if(fit$family[1]=='gaussian'){
            #remove intercept
            est <- 
              summary(fit)$coefficients[-1,1]%>%
              sprintf(fmt="%.2f")
            
            ci <- 
              confint(fit,method=cimethod)[-1,] %>% 
              as.data.frame()%>%
              mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
            
            #pval
            pval <-  
              summary(fit)$coefficients[-1,4]%>% 
              fmtp
          }
        
        
      }
      
      if(is.null(dim(ci))){
        df <- 
          data.frame(
            `Estimate`=est,
            `confint`= paste('(',ci[1],',',ci[2],')',sep =''),
            `pval`=pval
          )
      }
      else{
        df <- 
          data.frame(
            `Estimate`=est,
            `confint`= paste('(',ci[,1],',',ci[,2],')',sep =''),
            `pval`=pval
          )
      }
      #add N and fit statistics to bottom of table
      bot<-
        data.frame(
          `Estimate`= c('',nrow(fit$model),summary(fit)$family[[1]]),
          `confint`='',
          `pval`=''
        )
      rownames(bot) <- c('---','N','dist.')
      
      df <- rbind(df,bot)
      
      if(type=='response'){
        if(fit$family[[1]]=='gaussian'){
          names(df)[1] <- 'Estimate'
        }
        if(fit$family[[1]]=='binomial'){
          names(df)[1] <- 'OR'
        }
        if(fit$family[[1]]=='poisson'){
          names(df)[1] <- 'IRR'
        }
      }
      
      
      
    }
    
    df
    
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