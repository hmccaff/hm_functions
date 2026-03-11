#My functions

# Utility function to time-stamp filenames
timestamp <- function(prefix, ext = "csv") {
  paste0(prefix, "_", format(Sys.Date(), "%y%m%d"), ".", ext)
}


# returns df with all dupes of specified variable
dup <-
  function(df, id) {
    df[unlist(df[, id]) %in% unlist(df[duplicated(df[, id]), id]), ]
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



# Function to add a flextable with a header to the Word document
add_table_to_doc <- function(doc, df, table_title,width, fontsize=12, orient='none') {
  # Convert the data frame to a flextable
  ft <- flextable(df)
  
  ft <- width(ft, j=1:length(width), width=width)
  ft<- hrule(ft, rule='exact')
  ft <- height(ft, height=.25)
  ft <- fontsize(ft, size = fontsize, part = "all")
  
  # make header and ** cells bold
  ft <- ft %>% bold(part = "header")
  
  for(i in 1:nrow(df)){
    for(j in 1:ncol(df)){
      if(grepl("\\*",df[i,j])){
        
        clean_text <- gsub("\\*", "", df[i,j])
        
        # Apply the cleaned text back to the cell
        ft <- compose(ft, i = i, j = j, part = "body", value = as_paragraph(as_chunk(clean_text)))
        ft <-
          ft %>%
          bold(i=i,j=j)
      }
    }
  }
  
  # Add the table title
  doc <- doc %>% body_add_fpar(
    fpar(ftext(table_title, fp_text(bold = TRUE, font.size = 12)))
  ) 
  
  
  # Add the flextable to the document
  doc <- doc %>% body_add_flextable(ft)
  
  
  #use given orientation
  if(orient=='landscape'){
  doc <- body_end_section_landscape(doc)
  }
  if(orient=='portrait'){
    doc <- body_end_section_portrait(doc)
  }
  
  return(doc)
}

#format p values - used in coeftab and tabsum
fmtp <- 
  function(pvec){
    pvec <- as.numeric(pvec)
    sapply(pvec,
           function(x){
             if(is.na(x)){
               return('')
             }else if(x < .001){
               return(paste("<0.001",'***',sep=''))
             }else if(x<.01){
                return(paste(sprintf(x,fmt="%.3f"),'**',sep=''))
             }else if(x<.05){
               return(paste(sprintf(x,fmt="%.3f"),'*',sep=''))
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
    
    if(tab$control$test==T){
    st$`p value` <- fmtp(st$`p value`)
    }
    
    names(st)[which(names(st)=='')] <- 'Variable'
    
    st
  }


#takes a table of joint tests from emmeans and cleans it up
jttab <-
  function(tab){
    newtab<-
      tab%>%
      mutate(
        `model term` = `model term`,
        F.ratio = sprintf(F.ratio,fmt="%.3f"),
        Chisq = sprintf(Chisq,fmt="%.3f"),
        p.value = fmtp(p.value)
        
      )
    
    rownames(newtab) <- newtab$`model term`
    newtab[,-1]%>%
      rownames_to_column(var = 'Joint Tests')
  }

#binds two data frames containing model output of different row count, filling blanks with whitespace
cbind_fill <-
  function(x,y){
    #separate x and y into chunks divided by '---'
    xdiv <- which(x[,1]=='---')
    xchunks<- list(x[1:(xdiv[1]-1),]%>%list)
    if(length(xdiv)>1){
      for(i in 2:length(xdiv)){
        xchunks <- append(xchunks,
                          x[(xdiv[i-1]):(xdiv[i]-1),]%>%list
        )
      }
    }
    xchunks <- append(xchunks,x[c((xdiv[length(xdiv)]):nrow(x)),]%>%list)
    
    ydiv <- which(y[,1]=='---')
    ychunks<- list(y[1:(ydiv[1]-1),]%>%list)
    if(length(ydiv)>1){
      for(i in 2:length(ydiv)){
        ychunks <- append(ychunks,
                          y[(ydiv[i-1]):(ydiv[i]-1),]%>%list
        )
      }
    }
    ychunks <- append(ychunks,y[c((ydiv[length(ydiv)]):nrow(y)),]%>%list)
    
    if(length(xchunks) != length(ychunks)){
      break('Model result tables have different number of chunks')
    }else{
      chunklist <- list()
      for(i in 1:length(xchunks)){
        xc <- xchunks[[i]] %>% as.data.frame
        yc <- ychunks[[i]] %>% as.data.frame
        #if x and y have equal rows
        if((nrow(xc) - nrow(yc)) == 0){
          
          #divider
          div <- data.frame('.'=rep('',nrow(xc)))
          
          dfc <- cbind(xc,div,yc)
          
        }
        #if x is longer
        else if((nrow(xc) - nrow(yc)) > 0){
          
          #add a chunk of blank rows to match length of y to length of x
          drow <- nrow(xc)-nrow(yc)
          bot <- 
            as.data.frame(
              lapply(names(yc), function(col) rep("", drow)),
              stringsAsFactors = FALSE
            )
          names(bot) <- names(yc)
          
          ycnew <- rbind(yc,bot)
          
          #divider
          div <- data.frame('.'=rep('',nrow(xc)))
          
          dfc <- cbind(xc,div,ycnew)
          
        }
        #if y is longer
        else if((nrow(xc) - nrow(yc)) < 0){
          
          #add a chunk of blank rows to match length of x to length of y
          drow <- nrow(yc)-nrow(xc)
          bot <- 
            as.data.frame(
              lapply(names(xc), function(col) rep("", drow)),
              stringsAsFactors = FALSE
            )
          names(bot) <- names(xc)
          
          xcnew <- rbind(xc,bot)
          
          #divider
          div <- data.frame('.'=rep('',nrow(yc)))
          
          dfc <- cbind(xcnew,div,yc)
        }
        
        chunklist <- append(chunklist,dfc%>%list())
        
        
      }
      
      #recombine chunks into a single df
      df <- chunklist[[1]]
      names(df) <- names(chunklist[[2]])
      for(i in 2:length(chunklist)){
        df <- rbind(df,chunklist[[i]])
      }
      df
      
    }
  }


#uses cbind_fill to make a list of mods into a wide "stack" in order
modstack<-
  function(mods){
    stack <- mods[[1]]
    for(i in 2:length(mods)){
      stack <- cbind_fill(stack,mods[[i]])
    }
    stack
  }

# #binds two data frames of different row count, filling blanks with whitespace
# cbind_fill <-
#   function(x,y){
#     
#     #if x and y have equal rows
#     if((nrow(x) - nrow(y)) == 0){
#       
#       #divider
#       div <- data.frame('.'=rep('',nrow(x)))
#       
#       df <- cbind(x,div,y)
#       
#       df
#     }
#     
#     #if x has more rows: big model -> small model 
#     else if((nrow(x) - nrow(y)) > 0){
#       extrarows <- nrow(x) - nrow(y)
#       
#       #top of y - estimates
#       top <- y[1:which(y[,1]=='---')-1,]
#       
#       #bottom of y - summary stats
#       bot <- y[which(y[,1]=='---'):nrow(y),]
#       
#       #empty cells between estimates and summary stats
#       midrow <- rep('',ncol(y))
#       mid <- data.frame(t(replicate(extrarows, midrow)))
#       names(mid) <- names(y)
#     
#       newy <- rbind(top,mid,bot)
#       
#       #divider
#       div <- data.frame('.'=rep('',nrow(x)))
#       
#       df <- cbind(x,div,newy)
#       
#       df
#     }
#     
#     #if y has more rows: small model -> big model
#     else if((nrow(x) - nrow(y)) < 0){
#       extrarows <- nrow(y) - nrow(x)
#       
#       #top of x - estimates
#       top <- x[1:which(x[,1]=='---')-1,]
#       
#       #bottom of y - summary stats
#       bot <- x[which(x[,1]=='---'):nrow(x),]
#       
#       #empty cells between estimates and summary stats
#       midrow <- rep('',ncol(x))
#       mid <- data.frame(t(replicate(extrarows, midrow)))
#       names(mid) <- names(x)
#       
#       newx <- rbind(top,mid,bot)
#       
#       #divider
#       div <- data.frame('.'=rep('',nrow(y)))
#       
#       df <- cbind(newx,div,y)
#       
#       df
#     }
#     
#     
#   }
# 
# #uses cbind_fill to make a list of mods into a wide "stack" in order
# modstack<-
#   function(mods){
#     stack <- mods[[1]]
#     for(i in 2:length(mods)){
#       stack <- cbind_fill(stack,mods[[i]])
#     }
#     stack
#   }



#takes a posterior object and creates esetimate, CI and pseudo p
conttab<-
  function(cont,name,type='or'){
    if(type=='or'){
      median = median(cont)%>%exp%>% sprintf(fmt="%.2f")
      lb = posterior_interval(cont,prob=.95)[[1]]%>% exp%>% sprintf(fmt="%.2f")
      ub = posterior_interval(cont,prob=.95)[[2]]%>% exp%>% sprintf(fmt="%.2f")
    }
    
    
    
    pseudo_p <- function(x){
                        prob <- mean(x>0)
                        p <- 2* min(prob,1-prob)
                        fmtp(p)
    }
    
    p <- pseudo_p(cont)
    
    
    newtab<-
      data.frame(
        contrast = name,
        median = median,
        credint = paste('(',lb,', ',ub,')',sep=''),
        pseudo_p = p
        
      )
    newtab
  }


#takes a table of pairwise tests from emmeans and cleans it up with confints
pwtab <-
  function(tab,type='lp'){
    
    contrast <- tab %>% as.data.frame %>% select('contrast')%>%unlist()
    
    if(length(grep('\\|',contrast))>0){
    lvl <- tab %>% as.data.frame %>% select(2)%>%unlist()
    contrast <- paste0(contrast,' | ',lvl)
    }
    
    if(type=='hr'){
      lb <- tab %>% confint %>% as.data.frame %>% select('asymp.LCL') %>% unlist %>% sprintf(fmt="%.2f")
      ub <- tab %>% confint %>% as.data.frame %>% select('asymp.UCL') %>% unlist%>% sprintf(fmt="%.2f")
      est <- tab %>% as.data.frame %>% select('ratio')%>% unlist%>% sprintf(fmt="%.2f")
      pval <- tab %>% as.data.frame %>% select('p.value')%>% unlist %>% fmtp
      
    }
    
    if(type=='or'){
      lb <- tab %>% confint %>% as.data.frame %>% select('asymp.LCL') %>% unlist %>% sprintf(fmt="%.2f")
      ub <- tab %>% confint %>% as.data.frame %>% select('asymp.UCL') %>% unlist%>% sprintf(fmt="%.2f")
      est <- tab %>% as.data.frame %>% select('odds.ratio')%>% unlist%>% sprintf(fmt="%.2f")
      pval <- tab %>% as.data.frame %>% select('p.value')%>% unlist %>% fmtp
      
    }
    if(type=='lp'){
      lb <- tab %>% confint %>% as.data.frame %>% select('lower.CL') %>% unlist %>% sprintf(fmt="%.2f")
      ub <- tab %>% confint %>% as.data.frame %>% select('upper.CL') %>% unlist%>% sprintf(fmt="%.2f")
      est <- tab %>% as.data.frame %>% select('estimate')%>% unlist%>% sprintf(fmt="%.2f")
      pval <- tab %>% as.data.frame %>% select('p.value')%>% unlist %>% fmtp
      
    }
    else if(type == 'exp'){
    lb <- tab %>% confint %>% as.data.frame %>% select('lower.CL') %>% unlist %>%exp %>% sprintf(fmt="%.2f")
    ub <- tab %>% confint %>% as.data.frame %>% select('upper.CL') %>% unlist %>%exp %>% sprintf(fmt="%.2f")
    est <- tab %>% as.data.frame %>% select('estimate')%>% unlist %>%exp %>% sprintf(fmt="%.2f")
    pval <- tab %>% as.data.frame %>% select('p.value')%>% unlist %>% fmtp
    }
    
    
    newtab<-
      data.frame(
        contrast = contrast,
        estimate = est,
        confint = paste('(',lb,', ',ub,')',sep=''),
        p.value = pval
        
      )
    
  }


#add contrasts (from pwtab object) to model output table
add_cont <-
  function(mod,cont){
    colnames(cont) <- colnames(mod)
    cont[,1] <- as.character(cont[,1])
    
    top <- mod[-c(which(mod[,1]=='---'):nrow(mod)),]
    mid <- 
      rbind(
        c('---','','',''),
        c('Model Contrasts','','',''),
        cont
      )
    bot <- mod[c(which(mod[,1]=='---'):nrow(mod)),]
    
    rbind(top,mid,bot)
    
  }



#takes a model object and outputs a coefficient table with basic information
coeftab <-
  function(fit,cimethod='profile',type='response',labs=T,zph=T,pool=F,vif=F){
   
    if(pool==T){
      pfit<-pool(fit)
      sfit <- summary(pfit,conf.int=T,exponentiate=F)
      trim<-sfit[-1,c(1,2,7,8,6)]
      
      
      
      if('gee' %in% class(fit[[1]])){
        df <- 
          trim%>%
          mutate(
            Estimate = sprintf(fmt="%.2f",estimate),
            confint = paste('(',sprintf(fmt="%.2f",`2.5 %`),',',sprintf(fmt="%.2f",`97.5 %`),')',sep=''),
            pval = fmtp(`p.value`)
          )%>%
          select(term,Estimate,confint,pval)
        
        rownames(df) <- df$term
        df <- df %>% select(-term)
        
        #working corr NA if max clust size =1
        if(nrow(fit[[1]]$working.correlation)>1){wc <- sprintf(fit[[1]]$working.correlation[2,1],fmt="%.3f")}
        if(nrow(fit[[1]]$working.correlation)==1){wc <- as.character(NA)}
        
        bot<-
          data.frame(
            `Estimate`= c('',
                          as.character(fit[[1]]$call$formula[2]),
                          round(fit[[1]]$nobs,0),
                          length(unique(fit[[1]]$id)),
                          as.character(fit[[1]]$call[3]),
                          fit[[1]]$family[[1]],
                          summary(fit[[1]])$model$corstr,
                          wc,
                          pfit$m),
            `confint`='',
            `pval`=''
          )
        rownames(bot) <- c('---','outcome',
                           'N','N groups','grouping var.','dist.',
                           'corstr.','working correlation','imputed datasets')
        df <- rbind(df,bot)
      }else{
        df <- 
          trim%>%
          mutate(
            Estimate = sprintf(fmt="%.2f",estimate),
            confint = paste('(',sprintf(fmt="%.2f",`2.5 %`),',',sprintf(fmt="%.2f",`97.5 %`),')',sep=''),
            pval = fmtp(`p.value`)
          )%>%
          select(term,Estimate,confint,pval)
        
        rownames(df) <- df$term
        df <- df %>% select(-term)
        
        
        bot<-
          data.frame(
            `Estimate`= c('',
                          as.character(fit[[1]]$call$formula[2]),
                          fit[[1]]$family[[1]],
                          round(nrow(fit[[1]]$model),0),
                          pfit$m),
            `confint`='',
            `pval`=''
          )
        rownames(bot) <- c('---','outcome','dist.',
                           'N','imputed datasets')
        
        
        df <- rbind(df,bot)
        
      }
      
    }
    
    #stanreg mixed models
    if(('stanreg' %in% class(fit)) & ('glm' %in% class(fit)) & ('lmerMod' %in% class(fit))){
      #select posterior columns of interest
      post <- as.matrix(fit)
      post <- post[,which(colnames(post)%in% names(modlabs)),drop=F] 
      postdf <- post %>% as.data.frame
      
      if(type=='response'){
        
        if(fit$family[1]=='binomial'){
          est <- 
            lapply(postdf,median)%>%
            unlist%>%
            exp%>% 
            sprintf(fmt="%.2f")
          
          ci <-
            posterior_interval(post,prob=.95)%>%
            exp%>%
            as.data.frame()%>%
            mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
          
          
          #pseudo p value - probability that beta is as extreme as 0 in either direction
          pseudo_pval <- lapply(postdf,
                                function(x){
                                  prob <- mean(x>0)
                                  p <- 2* min(prob,1-prob)
                                  fmtp(p)
                                }
          )%>% unlist()
          
        }else if(fit$family[1]=='gaussian'){
          ############################ needs code
        }
        
        
      }
      
      if(is.null(dim(ci))){
        df <- 
          data.frame(
            `Estimate`=est,
            `cred_int`= paste('(',ci[1],',',ci[2],')',sep =''),
            `pseudo_pval`=pseudo_pval
          )
      }
      else{
        df <- 
          data.frame(
            `Estimate`=est,
            `cred_int`= paste('(',ci[,1],',',ci[,2],')',sep =''),
            `pseudo_pval`=pseudo_pval
          )
        rownames(df) <- rownames(ci)
      }
      #add N and fit statistics to bottom of table
      vc<- VarCorr(fit)
      
      bot<-
        data.frame(
          `Estimate`= c('', nrow(fit$data),as.character(fit$formula[[2]]),as.character(fit$formula)[[3]],fit$family[[1]],fit$prior.info$prior$dist),
          `cred_int`='',
          `pseudo_pval`=''
        )
      rownames(bot) <- c('---','N','Outcome','formula RHS.','dist.','prior dist.')
      
      # if(!is.null(fit$call$offset)){
      #   bot <- rbind(bot, c('as.character(c(fit$call$offset))','',''))
      #   bot<-
      #     data.frame(
      #       `Estimate`= c('',nrow(fit$model),as.character(fit$formula[[2]]),summary(fit)$family[[1]],
      #                     as.character(c(fit$call$offset))),
      #       `cred_int`='',
      #       `pseudo_pval`=''
      #     )
      #   rownames(bot) <- c('---','N','Outcome','dist.','Offset')
      # }
      
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
    
    #stanreg glm models
    else if('stanreg' %in% class(fit)){
      if(type=='response'){
        
        if(fit$family[1]=='binomial'){
          #remove intercept
          est <- 
            fit$coefficients[-1]%>%
            exp%>% 
            sprintf(fmt="%.2f")
          
          ci <-
            posterior_interval(fit,prob=.95)[-1,,drop=F]%>%
            exp%>%
            as.data.frame()%>%
            mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
          
          
          #pseudo p value - probability that beta is as extreme as 0 in either direction
          posterior <- as.data.frame(fit)%>% select(-'(Intercept)')
          pseudo_pval <- lapply(posterior,
                      function(x){
                        prob <- mean(x>0)
                        p <- 2* min(prob,1-prob)
                        fmtp(p)
                      }
          )%>% unlist()
          
        }else if(fit$family[1]=='gaussian'){
          ############################ needs code
        }
        
        
      }
      
      if(is.null(dim(ci))){
        df <- 
          data.frame(
            `Estimate`=est,
            `cred_int`= paste('(',ci[1],',',ci[2],')',sep =''),
            `pseudo_pval`=pseudo_pval
          )
      }
      else{
        df <- 
          data.frame(
            `Estimate`=est,
            `cred_int`= paste('(',ci[,1],',',ci[,2],')',sep =''),
            `pseudo_pval`=pseudo_pval
          )
        rownames(df) <- rownames(ci)
      }
      #add N and fit statistics to bottom of table
      bot<-
        data.frame(
          `Estimate`= c('',nrow(fit$model),as.character(fit$formula[[2]]),fit$family[[1]],fit$prior.info$prior$dist),
          `cred_int`='',
          `pseudo_pval`=''
        )
      rownames(bot) <- c('---','N','Outcome','dist.','prior dist.')
      
      if(!is.null(fit$call$offset)){
        bot <- rbind(bot, c('as.character(c(fit$call$offset))','',''))
        bot<-
          data.frame(
            `Estimate`= c('',nrow(fit$model),as.character(fit$formula[[2]]),summary(fit)$family[[1]],
                          as.character(c(fit$call$offset))),
            `cred_int`='',
            `pseudo_pval`=''
          )
        rownames(bot) <- c('---','N','Outcome','dist.','Offset')
      }
      
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
    
    #cox models
    else if('coxph' %in% class(fit)){
      est <- summary(fit)$conf.int[,1] %>% sprintf(fmt="%.2f")
      lb <- summary(fit)$conf.int[,3] %>% sprintf(fmt="%.2f")
      ub <- summary(fit)$conf.int[,4] %>% sprintf(fmt="%.2f")
      
      #if model uses clusters, or is svycoxph, use robust pvals (extra column)
      if('id' %in% names(fit$call)| 'svycoxph' %in% class(fit)){
        pval <- summary(fit)$coefficients[,6] %>% fmtp
      }else if('cluster' %in% names(fit$call)){
        pval <- summary(fit)$coefficients[,6] %>% fmtp
      }else{
        pval <- summary(fit)$coefficients[,5] %>% fmtp
      }
      
      df <- 
        data.frame(
          `HR`=est,
          `confint`= paste('(',lb,',',ub,')',sep =''),
          `pval`=pval
        )
      rownames(df) <- names(fit$coefficients)
      
      #add n clusters to bottom if model uses clusters
      if('id' %in% names(fit$call)){
        bot<-
          data.frame(
            `HR`= c('',
                    round(fit$n,0),
                    as.character(fit$call$formula[[2]])[[2]],
                    'Cox PH',
                    round(fit$n.id,0),
                    as.character(fit$call$id)),
            `confint`='',
            `pval`=''
          )
        rownames(bot) <- c('---',
                           'N',
                           'Outcome',
                           'Model',
                           'N groups',
                           'grouping var')
      }else if('cluster' %in% names(fit$call)){
        bot<-
          data.frame(
            `HR`= c('',
                    round(fit$n,0),
                    as.character(fit$call$formula[[2]])[[2]],
                    'Cox PH',
                    as.character(fit$call$cluster)),
            `confint`='',
            `pval`=''
          )
        rownames(bot) <- c('---',
                           'N',
                           'Outcome',
                           'Model',
                           'cluster var')
         
      }else{
      bot<-
        data.frame(
          `HR`= c('',
                  round(fit$n,0),
                  as.character(fit$call$formula[[2]])[[2]],
                  'Cox PH'
                  ),
          `confint`='',
          `pval`=''
        )
      rownames(bot) <- c('---',
                         'N',
                         'Outcome',
                         'Model')
      
      }
      
      #tests of proportionality
      if(zph){
        zph_df <- cox.zph(fit)[[1]] %>% as.data.frame
        
        zphbot <- 
          data.frame(
            `HR`=c('',fmtp(zph_df$p)),
            `confint`='',
            `pval`=''
          )
        rownames(zphbot) <- c('ZPH pvals:',rownames(zph_df))
      }
      
      if(zph){
        df <- rbind(df,bot,zphbot)
      }
      if(!zph){
      df <- rbind(df,bot)
      }
    }
    
    
    #cox models
    else if('coxme' %in% class(fit)){
      est <- summary(fit)$coefficients[,2] %>% sprintf(fmt="%.2f")
      lb <- confint(fit)[,1] %>% exp() %>% sprintf(fmt="%.2f")
      ub <- confint(fit)[,2] %>% exp() %>% sprintf(fmt="%.2f")
      
      pval <- summary(fit)$coefficients[,5] %>% fmtp
      
      df <- 
        data.frame(
          `HR`=est,
          `confint`= paste('(',lb,',',ub,')',sep =''),
          `pval`=pval
        )
      rownames(df) <- names(fit$coefficients)
      
      
      #add random effects to bottom of table
      vc <- unlist(fit$vcoef)
      
      mid <- 
        data.frame(HR=c('','',sprintf(vc,fmt="%.3f")),
                   confint='',
                   pval=''
        )
      rn<-names(vc)
      rn <- c('------','variance components:',rn)
      rownames(mid) <- rn
      
      
      bot<-
        data.frame(
          `HR`= c('',
                  fit$n[1],
                  fit$n[2],
                  as.character(fit$call$formula[[2]])[[2]],
                  'Cox PH'
                  ),
          `confint`='',
          `pval`=''
        )
      rownames(bot) <- c('---',
                         'N events',
                         'N',
                         'Outcome',
                         'Model')
      
      #tests of proportionality
      if(zph){
        zph_df <- cox.zph(fit)[[1]] %>% as.data.frame
        
        zphbot <- 
          data.frame(
            `HR`=c('',fmtp(zph_df$p)),
            `confint`='',
            `pval`=''
          )
        rownames(zphbot) <- c('ZPH pvals:',rownames(zph_df))
      }
      
      if(zph){
        df <- rbind(df,mid,bot,zphbot)
      }
      if(!zph){
      df <- rbind(df,mid,bot)
      }
    }
    
    #lmer objects, linear mixed models
    else if('lmerModLmerTest' %in% class(fit)){
      
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
        rownames(df) <- rownames(ci)
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
    
    else if('glmerMod' %in% class(fit)){
      
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
    
    else if('gee' %in% class(fit)){
      
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
      }
      #add N and fit statistics to bottom of table
      
      #working corr NA if max clust size =1
      if(nrow(fit$working.correlation)>1){wc <- sprintf(fit$working.correlation[2,1],fmt="%.3f")}
      if(nrow(fit$working.correlation)==1){wc <- as.character(NA)}
      
      bot<-
        data.frame(
          `Estimate`= c('',
                        as.character(fit$call$formula[2]),
                        round(fit$nobs,0),length(unique(fit$id)),
                        as.character(fit$call[3]),
                        fit$family[[1]],
                        summary(fit)$model$corstr,
                        wc),
          `confint`='',
          `pval`=''
        )
      rownames(bot) <- c('---','outcome',
                         'N','N groups','grouping var.','dist.',
                         'corstr.','working correlation')
      
      
      rownames(df) <- rownames(summary(fit)$coefficients)[-1]
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
      
      if(type=='exp'){
        if(fit$family[[1]]=='gaussian'){
          names(df)[1] <- 'Exp(beta)'
        }
      }
      
    }
    
    
    #Firth bias-corrected logistic regression
    else if('logistf' %in% class(fit)){
      
      #on linear predictor scale 
      if(type=='lp'){
        #point estimates 
        est <- fit$coefficients%>% sprintf(fmt="%.3f")
        
        ci <- 
          confint(fit,method=cimethod)%>%
          as.data.frame()%>%
          mutate(across(c(1,2),function(x){sprintf(fmt="%.3f",x)}))
        
        
        #format p values
        pval <- fit$prob%>%fmtp
        
        #exponentiated
      }else if(type=='response'){
        #point estimates 
        est <- fit$coefficients%>% exp %>% sprintf(fmt="%.3f")
        est <- est[-1]
        
        ci <- 
          confint(fit,method=cimethod)%>%
          as.data.frame()%>% exp %>% 
          mutate(across(c(1,2),function(x){sprintf(fmt="%.3f",x)}))
        ci <- ci[-1,]
        
        #format p values
        pval <- fit$prob%>%fmtp
        pval <- pval[-1]
        
        
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
      #add N and fit statistics to bottom of table
      bot<-
        data.frame(
          `Estimate`= c('',nrow(fit$model),as.character(fit$formula[[2]]),"Firths bias-reduced logistic regression"),
          `confint`='',
          `pval`=''
        )
      rownames(bot) <- c('---','N','Outcome','dist.')
      
      df <- rbind(df,bot)
      
      if(type=='response'){
          names(df)[1] <- 'OR'
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
        
        
        #for simple models (intercept + coef), confint df has different dimensions
        if(length(fit$coefficients)==2){
          ci<-
            confint(fit,method=cimethod)[-1,] %>% 
            exp%>%
            t()%>%
            as.data.frame()%>%
            mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
          rownames(ci) <- names(fit$coefficients)[2]
        }else{
          ci <- 
            confint(fit,method=cimethod)[-1,] %>% 
            as.data.frame()%>%
            exp%>%
            mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
        }
        
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
          
          
          #for simple models (intercept + coef), confint df has different dimensions
          if(length(fit$coefficients)==2){
            ci<-
            confint(fit,method=cimethod)[-1,] %>% 
              exp%>%
              t()%>%
              as.data.frame()%>%
              mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
            rownames(ci) <- names(fit$coefficients)[2]
          }else{
          ci <- 
            confint(fit,method=cimethod)[-1,] %>% 
            as.data.frame()%>%
            exp%>%
            mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
          }
          
          #pval
          pval <-  
            summary(fit)$coefficients[-1,4]%>% 
            fmtp
        }else if(fit$family[1]=='gaussian'){
          #remove intercept
          est <- 
            summary(fit)$coefficients[-1,1]%>%
            sprintf(fmt="%.2f")
          
          if(length(fit$coefficients)==2){
            ci<-
              confint(fit,method=cimethod)[-1,] %>% 
              t()%>%
              as.data.frame()%>%
              mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
            rownames(ci) <- names(fit$coefficients)[2]
          }else{
          ci <- 
            confint(fit,method=cimethod)[-1,] %>% 
            as.data.frame()%>%
            mutate(across(c(1,2),function(x){sprintf(fmt="%.2f",x)}))
          }
          
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
        rownames(df) <- rownames(ci)
      }
      #add N and fit statistics to bottom of table
      bot<-
        data.frame(
          `Estimate`= c('',nrow(fit$model),as.character(fit$formula[[2]]),summary(fit)$family[[1]]),
          `confint`='',
          `pval`=''
        )
      rownames(bot) <- c('---','N','Outcome','dist.')
      
      if(!is.null(fit$call$offset)){
        bot <- rbind(bot, c('as.character(c(fit$call$offset))','',''))
        bot<-
          data.frame(
            `Estimate`= c('',nrow(fit$model),as.character(fit$formula[[2]]),summary(fit)$family[[1]],
                          as.character(c(fit$call$offset))),
            `confint`='',
            `pval`=''
          )
        rownames(bot) <- c('---','N','Outcome','dist.','Offset')
      }
      
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
      
      if(type=='exp'){
        if(fit$family[[1]]=='gaussian'){
          names(df)[1] <- 'Exp(beta)'
        }
      }
      
      
    }
    
    
    #ad vif to bottom
    if(vif==T){
      if('lmerModLmerTest' %in% class(fit) | 'glmerMod' %in% class(fit)){
        v <- car::vif(lm(fit@call$formula, data=fit@frame))
      }else if(pool==T){
        v <- car::vif(fit[[1]])
        if(is.matrix(v)){
          v <- v[,1]
        }
      }else{
        v <- car::vif(fit)
      }
      
      vifdf <- data.frame(Estimate= c('',sprintf(v,fmt="%.2f")),
                          confint='',
                          pval=''
      )
      rownames(vifdf) <- c('Diagnostics:',paste('VIF', names(v), sep=' '))
      
      
      df <- rbind(df, vifdf)
    }
    
    
    #replace row names with labels from modlabs vector
    if(labs==T){
      if(exists("modlabs", envir = .GlobalEnv)){
          rownames(df) <- ifelse(rownames(df) %in% names(modlabs), 
                           modlabs[rownames(df)], 
                           rownames(df))  # Replace with corresponding new name
    }else{stop('error: no modlabs provided in environment')}
    }
    
    df%>%
      rownames_to_column(var = 'Model Results')
    
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