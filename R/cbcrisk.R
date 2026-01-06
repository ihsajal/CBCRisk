cbcrisk <- function(mastectomy = 0, race = 0, profile, start.age, pred.year=5, print.output=T){
  cbcrisk_nonBlack <- function(profile, start.age, pred.year=5, print.output=T){
    
    ####### start age check #########
    if(start.age <18|start.age >89)
    {stop("current age has to be within 18-89")}
    
    ###### profile check ##########
    
    if (length(profile)<8|length(profile)>8)
    {stop("Please enter a valid profile (8 variables)")}
    
    if(profile[1]> 3 | profile[1]<1)
    {stop("Please enter a valid code for 'Age at first BC diagnosis'")}
    
    if(profile[2]> 3 | profile[2]<1)
    {stop("Please enter a valid code for 'Anti-estrogen therapy'")}
    
    if(profile[3]> 3 | profile[3]<1)
    {stop("Please enter a valid code for 'First degree family history of BC'")}
    
    if(profile[4]> 2 | profile[4]<1)
    {stop("Please enter a valid code for 'High risk pre-neoplasia'")}
    
    if(profile[5]> 5 | profile[5]<1)
    {stop("Please enter a valid code for 'Breast density'")}
    
    if(profile[6]> 3 | profile[6]<1)
    {stop("Please enter a valid code for 'ER status'")}
    
    if(profile[7]> 3 | profile[7]<1)
    {stop("Please enter a valid code for 'First BC type'")}
    
    if(profile[8]> 4 | profile[8]<1)
    {stop("Please enter a valid code for 'Age at first child birth'")}
    
    
    ########### consistency of the age variables check #########
    
    if (profile[1] ==2 & start.age <30)
    {stop("Current age has to be equal to or greater than the age at first breast cancer diagnosis")}
    
    
    if (profile[1] ==3 & start.age <40)
    {stop("Current age has to be equal to or greater than the age at first breast cancer diagnosis")}
    
    
    if (profile[8] ==2 & start.age <30)
    {stop("Age at first child birth has to be equal to or less than the current age, If there is no birth before current age please choose Nulliparous")}
    
    
    if (profile[8] ==3 & start.age <40)
    {stop("Age at first child birth has to be equal to or less than the current age, if there is no birth before current age please choose Nulliparous")}
    
    
    #######################################
    
    ########### To avoid installing a relatively large number of dependencies, this package uses squeezeBlanks and recode functions from car package by copying the functions directly from the package. 
    squeezeBlanks <- function(text){
      gsub(" *", "",  text)
    }
    
    recode <- function (var, recodes, as.factor, as.numeric = TRUE, levels) 
    {
      lo <- -Inf
      hi <- Inf
      recodes <- gsub("\n|\t", " ", recodes)
      recode.list <- rev(base::strsplit(recodes, ";")[[1]])
      is.fac <- is.factor(var)
      if (missing(as.factor)) 
        as.factor <- is.fac
      if (is.fac) 
        var <- as.character(var)
      result <- var
      for (term in recode.list) {
        if (0 < length(grep(":", term))) {
          range <- base::strsplit(base::strsplit(term, "=")[[1]][1], ":")
          low <- try(eval(parse(text = range[[1]][1])), silent = TRUE)
          if (class(low) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 low)
          }
          high <- try(eval(parse(text = range[[1]][2])), silent = TRUE)
          if (class(high) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 high)
          }
          target <- try(eval(parse(text = base::strsplit(term, "=")[[1]][2])), 
                        silent = TRUE)
          if (class(target) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 target)
          }
          result[(var >= low) & (var <= high)] <- target
        }
        else if (0 < length(grep("^else=", squeezeBlanks(term)))) {
          target <- try(eval(parse(text = base::strsplit(term, "=")[[1]][2])), 
                        silent = TRUE)
          if (class(target) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 target)
          }
          result[1:length(var)] <- target
        }
        else {
          set <- try(eval(parse(text = base::strsplit(term, "=")[[1]][1])), 
                     silent = TRUE)
          if (class(set) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 set)
          }
          target <- try(eval(parse(text = base::strsplit(term, "=")[[1]][2])), 
                        silent = TRUE)
          if (class(target) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 target)
          }
          for (val in set) {
            if (is.na(val)) 
              result[is.na(var)] <- target
            else result[var == val] <- target
          }
        }
      }
      if (as.factor) {
        result <- if (!missing(levels)) 
          factor(result, levels = levels)
        else as.factor(result)
      }
      else if (as.numeric && (!is.numeric(result))) {
        result.valid <- na.omit(result)
        opt <- options(warn = -1)
        result.valid <- as.numeric(result.valid)
        options(opt)
        if (!any(is.na(result.valid))) 
          result <- as.numeric(result)
      }
      result
    }
    
    
    ###################################################
    
    
    abs_risk_cbc=function(profile, start.age, pred.year, h1star.seer, h2.seer)
    {
      require(survival) 
      
      ### Preparing the data:
      data=data.frame(t(profile))
      ## Changing the variable names:
      colnames(data) <- c("age_1std","hormf","famhx_bc","lcis_atyph","brd_density","ER","bc_typ","cat_age1stb")
      
      ## Recoding the inputs:
      data$age_1std=recode(data$age_1std,"1='<30';2='30-40';3='40+'",as.factor=T)
      data$hormf=recode(data$hormf,"1='Yes';2='No';3='Unk'",as.factor=T)
      data$famhx_bc=recode(data$famhx_bc,"1='Yes';2='No';3='Unk'",as.factor=T)
      data$lcis_atyph=recode(data$lcis_atyph,"1='Yes';2='Unk'",as.factor=T)
      data$brd_density=recode(data$brd_density, "4='all_fat';3='scattered';2='heterog_dense';1='extrm_dense';5='Unk'",as.factor=T)
      data$ER=recode(data$ER,"1='Neg';2='Pos';3='Unk'",as.factor=T)
      data$bc_typ=recode(data$bc_typ,"1='Pure DCIS';2='Invasive_DCIS';3='Pure Invasive'",as.factor=T)
      data$cat_age1stb=recode(data$cat_age1stb,"1='<30/nulli';2='30-39';3='40+';4='Unk'",as.factor=T)
      
      data$grp=1010
      
      
      
      
      
      cbc_coeff <- c( 0.78496379 ,0.27152492,-0.05094703,-0.2511127, 0.18167327,0.44951127,0.45036138,0.70396452,0.53252559,0.42360416,0.39372522,0.12897412,-0.17737098,0.30621685, 0.50522831,0.26240513 ,
                      1.31195801,-0.02886171)
      
      names(cbc_coeff) <- c("age_1std<30", "age_1std30-40" ,"hormfUnk","hormfYes","famhx_bcUnk", "famhx_bcYes" ,"lcis_atyphYes","brd_densityextrm_dense","brd_densityheterog_dense","brd_densityscattered","brd_densityUnk","ERNeg",                   
                            "ERUnk","bc_typInvasive_DCIS","bc_typPure DCIS","cat_age1stb30-39", "cat_age1stb40+","cat_age1stbUnk" )
      
      
      x1=model.matrix(fmod11.bcsc,data)
      RR=exp(x1%*%cbc_coeff)       ### Getting the RR#####  
      
      aa=c(18,30,35,40,45,50,55,60,65,70,75,80,85)
      bb=c(29,34,39,44,49,54,59,64,69,74,79,84,89)+1
      
      
      
      ##### Choosing the hazard rates:    
      
      if(h1star.seer==T){
        h1star=c(0.002642137,0.003639010,0.003480063,0.003355744,0.003450679,0.003212342,0.003505955,0.003559329,0.004015681,0.003960725,0.004364455,0.003938050,0.003553980)}
      else {h1star=c(0.005938242,0.004942339,0.004235975,0.003648188,0.003079530,0.003297890,0.003486713,0.003210002,0.004035436,0.00431467,0.003866516,
                     0.004041648,0.003070979)}
      
      if (h2.seer==T)
      {other.hazard=c(0.03465443,0.03191816,0.02610106,0.02000662,0.01701598,0.01677001,0.01791297,0.01960874,0.02232529,0.02948460,0.04226024,0.06412088,0.10045459)}
      else {other.hazard=c(0.01672640,0.02119907,0.01616915,0.01288945,0.01173593,0.01086192,0.01095524,0.01292321,0.01443535,0.02118252,0.02893834,0.04330877,0.05676192)}
      
      ar=c(0.7073222,0.7031257,0.6135916,0.5476586,0.4961982,0.4727390,0.4566321,0.4531747,0.4097558,0.4398585,0.4400653,0.4078442,0.4004540)
      # 
      cbc.hazard=h1star*(1-ar)
      
      
      
      st=findInterval(start.age,aa)
      
      if((start.age+pred.year) %in% aa){ 
        
        fn= (findInterval(start.age+pred.year,aa))-1
      } else {fn=findInterval(start.age+pred.year,aa)}
      
      
      surv.cbc=rep(0,length=fn+1)     ##### Survival function vector for cbc
      surv.other=rep(0,length=fn+1)     #####  Survival function vector for non-cbc
      surv.cbc[1]=1
      surv.other[1]=1
      surv.cbc_a=0
      surv.other_a=0
      
      abs_prob=rep(0,length=fn)
      
      ### Calculating the survival function and absolute risk###    
      for (i in 1:fn)
      {
        
        if (start.age>=aa[i] & start.age<=bb[i] & start.age+pred.year<=bb[i]){   #### this condition applies when the prediction ends in the same interval as the current age
          abs_prob[i] = ((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(start.age+pred.year-start.age)))
          
        }
        else if (start.age>=aa[i] & start.age<=bb[i] & start.age+pred.year>bb[i]){ ### This allows the prediction to go further
          abs_prob[i] = ((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(bb[i]-start.age)))
          
          surv.cbc[i+1]=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(bb[i]-aa[i]))   ### Survival function till the previous interval
          surv.other[i+1]= surv.other[i]*exp(-other.hazard[i]*(bb[i]-aa[i]))
          
          surv.cbc_a=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(start.age-aa[i])) #### Survival function till start.age
          surv.other_a=surv.other[i]*exp(-other.hazard[i]*(start.age-aa[i]))     #### Survival function till start.age
        } 
        else if (start.age<aa[i] & start.age+pred.year>aa[i] & start.age+pred.year<=bb[i]){   ### This is for the probability calculated in the last interval of the prediction
          abs_prob[i]=((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(surv.cbc[i]/surv.cbc_a)*(surv.other[i]/surv.other_a)*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(start.age+pred.year-aa[i]))) 
          surv.cbc[i+1]=0
          surv.other[i+1]=0
        }
        else if (start.age<aa[i] & start.age+pred.year>bb[i]){      #### Probability for an interval within the prediction years
          abs_prob[i]=((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(surv.cbc[i]/surv.cbc_a)*(surv.other[i]/surv.other_a)*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(bb[i]-aa[i])))
          
          surv.cbc[i+1]=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(bb[i]-aa[i]))
          surv.other[i+1]= surv.other[i]*exp(-other.hazard[i]*(bb[i]-aa[i]))  
        }
        else 
        {abs_prob[i]=0
        surv.cbc[i+1]=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(bb[i]-aa[i]))
        surv.other[i+1]= surv.other[i]*exp(-other.hazard[i]*(bb[i]-aa[i])) 
        }
        
      } 
      
      result <- c(start.age+pred.year, round(100*sum(abs_prob),2))
      #names(result) <- c("current.age", "nyears", "abs.risk(%)")
      return(result)
    }
    pred_year=seq(pred.year,89,pred.year)
    output=data.frame(matrix(  ,nrow = 0, ncol = 2))
    colnames(output)=c("by age", "CBC risk(%)")
    
    if(start.age + pred.year[1] >89)
    {stop("current age + nyears.pred < 89")}
    
    for (i in 1:length(pred_year))
    {  
      if(start.age+pred_year[i] <=89)
      {output[i,]=abs_risk_cbc(profile,start.age,pred.year=pred_year[i],h1star.seer=T,h2.seer=T)
      }
      else {break}
      
    }
    prof=data.frame(t(profile))
    
    colnames(prof) <- c("age_1std","hormf","famhx_bc","lcis_atyph","brd_density","ER","bc_typ","cat_age1stb")
    
    ## Recoding the inputs:
    prof$age_1std=recode(prof$age_1std,"1='<30';2='30-40';else='40+'",as.factor=T)
    prof$hormf=recode(prof$hormf,"1='Yes';2='No';else='Unk'",as.factor=T)
    prof$famhx_bc=recode(prof$famhx_bc,"1='Yes';2='No';else='Unk'",as.factor=T)
    prof$lcis_atyph=recode(prof$lcis_atyph,"1='Yes';else='Unk'",as.factor=T)
    prof$brd_density=recode(prof$brd_density, "4='Almost ent. fat';3='Scattered';2='Heterog dense';1='Extreme dense';else='Unk'",as.factor=T)
    prof$ER=recode(prof$ER,"1='Neg';2='Pos';else='Unk'",as.factor=T)
    prof$bc_typ=recode(prof$bc_typ,"1='Pure DCIS';2='Invasive_DCIS';else='Pure Invasive'",as.factor=T)
    prof$cat_age1stb=recode(prof$cat_age1stb,"1='<30/nulli';2='30-39';3='40+';else='Unk'",as.factor=T)
    
    colnames(prof)=c("Age at first BC diagnosis", "Anti-estrogen therapy","Family history of BC","High Risk Preneoplasia","Breast density","ER","First BC type","Age at first birth")
    
    row.names(prof) <- NULL
    row.names(output) <- NULL
    
    out.list=list("profile" = prof, "current_age" = start.age, "risk" = output)
    
    if (print.output==T)
    {return(out.list)}
    else {invisible(out.list)}
  }
  cbcriskBlack <- function (profile, start.age, pred.year, print.output = T){  
    
    #profile = c(1,2,1,2); start.age = 55; pred.year= 3
    ####### start age check #########
    if (start.age < 18 | start.age > 89) {
      stop("current age has to be within 18-89")
    }
    
    ###### profile check ##########
    
    if (length(profile) < 4 | length(profile) > 4) {
      stop("Please enter a valid profile (8 variables)")
    }
    if (profile[1] > 3 | profile[1] < 1) {
      stop("Please enter a valid code for 'Breast Density'")
    }
    if (profile[2] > 3 | profile[2] < 1) {
      stop("Please enter a valid code for 'First degree family history of BC'")
    }
    if (profile[3] > 4 | profile[3] < 1) {
      stop("Please enter a valid code for 'Tumor Size'")
    }
    if (profile[4] > 2 | profile[4] < 1) {
      stop("Please enter a valid code for 'Age at first diagnosis'")
    }
    
    ########### consistency of the age variables check #########
    
    if (profile[4] == 2 & start.age < 40) {
      stop("Current age has to be equal to or greater than the age at first breast cancer diagnosis")
    }
    ########### To avoid installing a relatively large number of dependencies, 
    ########### this package uses squeezeBlanks and recode functions from car package
    ########### by copying the functions directly from the package. 
    
    squeezeBlanks <- function(text) {
      gsub(" *", "", text)
    }
    recode <- function(var, recodes, as.factor, as.numeric = TRUE, 
                       levels) {
      lo <- -Inf
      hi <- Inf
      recodes <- gsub("\n|\t", " ", recodes)
      recode.list <- rev(base::strsplit(recodes, ";")[[1]])
      is.fac <- is.factor(var)
      if (missing(as.factor)) 
        as.factor <- is.fac
      if (is.fac) 
        var <- as.character(var)
      result <- var
      for (term in recode.list) {
        if (0 < length(grep(":", term))) {
          range <- base::strsplit(base::strsplit(term, "=")[[1]][1], ":")
          low <- try(eval(parse(text = range[[1]][1])), 
                     silent = TRUE)
          if (class(low) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", low)
          }
          high <- try(eval(parse(text = range[[1]][2])), silent = TRUE)
          if (class(high) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", high)
          }
          target <- try(eval(parse(text = base::strsplit(term, "=")[[1]][2])), silent = TRUE)
          if (class(target) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", target)
          }
          result[(var >= low) & (var <= high)] <- target
        }
        else if (0 < length(grep("^else=", squeezeBlanks(term)))) {
          target <- try(eval(parse(text = base::strsplit(term, "=")[[1]][2])), silent = TRUE)
          if (class(target) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", target)
          }
          result[1:length(var)] <- target
        }
        else {
          set <- try(eval(parse(text = base::strsplit(term, "=")[[1]][1])), silent = TRUE)
          if (class(set) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", set)
          }
          target <- try(eval(parse(text = base::strsplit(term, "=")[[1]][2])), silent = TRUE)
          if (class(target) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", target)
          }
          for (val in set) {
            if (is.na(val)) 
              result[is.na(var)] <- target
            else result[var == val] <- target
          }
        }
      }
      if (as.factor) {
        result <- if (!missing(levels)) 
          factor(result, levels = levels)
        else as.factor(result)
      }
      else if (as.numeric && (!is.numeric(result))) {
        result.valid <- na.omit(result)
        opt <- options(warn = -1)
        result.valid <- as.numeric(result.valid)
        options(opt)
        if (!any(is.na(result.valid))) 
          result <- as.numeric(result)
      }
      result
    }
    
    
    #fmod.10.black <- readRDS("finalmodel.black.rds")
    
    abs_risk_cbc = function(profile, start.age, pred.year){
      require(survival)
      prof.data = data.frame(t(profile))
      colnames(prof.data) <- c("brd_density.m", "famhx_bc", "tm_siz.v2", "age_dx.m")
      prof.data$brd_density.m = recode(prof.data$brd_density.m, "1='scattered/all_fat';2='heter/ext_dense';3='Unk'", as.factor = TRUE)
      prof.data$famhx_bc = recode(prof.data$famhx_bc, "1='Yes';2='No';3='Unk'", as.factor = TRUE)
      prof.data$tm_siz.v2 = recode(prof.data$tm_siz.v2, "1='T0/T1/T2';2='T3/T4';3='TIS';4='Unk'", as.factor = TRUE)
      prof.data$age_dx.m = recode(prof.data$age_dx.m, "1='<40';2='40+'", as.factor = TRUE)
      prof.data$grp = 186
      
      cbc_coeff <- c(0.75796954, 0.14835361, -0.03749922, 0.82421254, 0.76240886, 0.44189823, 0.00000000, 0.34543988)
      names(cbc_coeff) <- c("brd_density.mheter/ext_dense", "brd_density.mUnk", "famhx_bcUnk", 
                            "famhx_bcYes", "tm_siz.v2T3/T4", "tm_siz.v2TIS", "tm_siz.v2Unk", "age_dx.m<40")
      x1 = model.matrix(fmod.10.black, prof.data)
      RR = exp(x1 %*% cbc_coeff)
      aa = c(18, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85)
      bb = c(29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, 89) + 1
      
      h1star = c(0.001984127, 0.004812479, 0.004800216,  
                 0.005273998, 0.004413701, 0.004120600, 0.004997154,  
                 0.004771886, 0.004669028, 0.004312749, 0.005003362, 
                 0.004077472, 0.004118386)
      
      other.hazard = c(0.05978261, 0.05466557, 0.04661755, 
                       0.03654126, 0.03306543, 0.03335578, 0.03399081,  
                       0.03482456, 0.03849821, 0.04764718, 0.06135798, 
                       0.08719388, 0.12253724)
      
      ar = c(0.5082915, 0.5082915, 0.6041760, 0.6162957, 0.3894131, 0.4605934, 0.4106439, 
             0.3021153, 0.3339779, 0.3987406, 0.4696751, 0.2429663, 0.2059945)
      
      cbc.hazard = h1star * (1 - ar)
      st = findInterval(start.age, aa)
      
      fn = ifelse((start.age + pred.year) %in% aa, (findInterval(start.age + pred.year, aa)) - 1, 
                  findInterval(start.age + pred.year, aa))
      surv.cbc = rep(0, length = fn + 1)
      surv.other = rep(0, length = fn + 1)
      surv.cbc[1] = 1
      surv.other[1] = 1
      surv.cbc_a = 0
      surv.other_a = 0
      abs_prob = rep(0, length = fn)
      for (i in 1:fn) {
        if (start.age >= aa[i] & start.age <= bb[i] & start.age + 
            pred.year <= bb[i]) {
          abs_prob[i] = ((cbc.hazard[i] * RR)/(cbc.hazard[i]*RR + other.hazard[i])) * (1 - exp(-(cbc.hazard[i] * 
                                                                                                   RR + other.hazard[i]) * (start.age + pred.year -  start.age)))
        }
        else if (start.age >= aa[i] & start.age <= bb[i] & 
                 start.age + pred.year > bb[i]) {
          abs_prob[i] = ((cbc.hazard[i] * RR)/(cbc.hazard[i] *  RR + 
                                                 other.hazard[i])) * (1 - exp(-(cbc.hazard[i] * RR + other.hazard[i]) * (bb[i] - start.age)))
          surv.cbc[i + 1] = surv.cbc[i] * exp(-cbc.hazard[i] * RR * (bb[i] - aa[i]))
          surv.other[i + 1] = surv.other[i] * exp(-other.hazard[i] * (bb[i] - aa[i]))
          surv.cbc_a = surv.cbc[i] * exp(-cbc.hazard[i] * RR * (start.age - aa[i]))
          surv.other_a = surv.other[i] * exp(-other.hazard[i] * (start.age - aa[i]))
        }
        else if (start.age < aa[i] & start.age + pred.year > 
                 aa[i] & start.age + pred.year <= bb[i]) {
          abs_prob[i] = ((cbc.hazard[i] * RR)/(cbc.hazard[i] *  RR + other.hazard[i])) * (surv.cbc[i]/surv.cbc_a) * 
            (surv.other[i]/surv.other_a) * (1 - exp(-(cbc.hazard[i] * RR + other.hazard[i]) * (start.age + pred.year - aa[i])))
          surv.cbc[i + 1] = 0
          surv.other[i + 1] = 0
        }
        else if (start.age < aa[i] & start.age + pred.year > bb[i]) {
          abs_prob[i] = ((cbc.hazard[i] * RR)/(cbc.hazard[i] *  RR + other.hazard[i])) * (surv.cbc[i]/surv.cbc_a) * 
            (surv.other[i]/surv.other_a) * (1 - exp(-(cbc.hazard[i] * RR + other.hazard[i]) * (bb[i] - aa[i])))
          surv.cbc[i + 1] = surv.cbc[i] * exp(-cbc.hazard[i] *   RR * (bb[i] - aa[i]))
          surv.other[i + 1] = surv.other[i] * exp(-other.hazard[i] * (bb[i] - aa[i]))
        }
        else {
          abs_prob[i] = 0
          surv.cbc[i + 1] = surv.cbc[i] * exp(-cbc.hazard[i] *  RR * (bb[i] - aa[i]))
          surv.other[i + 1] = surv.other[i] * exp(-other.hazard[i] * (bb[i] - aa[i]))
        }
      }
      result <- c(start.age + pred.year, round(100 * sum(abs_prob), 2))
      return(result)
    }
    
    pred_year = seq(pred.year, 89, pred.year)
    output = data.frame(matrix(, nrow = 0, ncol = 2))
    colnames(output) = c("by age", "CBC risk(%)")
    if (start.age + pred.year[1] > 89) {
      stop("current age + nyears.pred < 89")
    }
    for (i in 1:length(pred_year)) {
      if (start.age + pred_year[i] <= 89) {
        output[i, ] = abs_risk_cbc(profile, start.age, pred.year = pred_year[i])
      }
      else {
        break
      }
    }
    
    prof.data = data.frame(t(profile))
    colnames(prof.data) <- c("brd_density.m", "famhx_bc", "tm_siz.v2", "age_dx.m")
    prof.data$brd_density.m = recode(prof.data$brd_density.m, "1='scattered/all_fat';2='heter/ext_dense';3='Unk'", as.factor = TRUE)
    prof.data$famhx_bc = recode(prof.data$famhx_bc, "1='Yes';2='No';3='Unk'", as.factor = TRUE)
    prof.data$tm_siz.v2 = recode(prof.data$tm_siz.v2, "1='T0/T1/T2';2='T3/T4';3='TIS';4='Unk'", as.factor = TRUE)
    prof.data$age_dx.m = recode(prof.data$age_dx.m, "1='<40';2='40+'", as.factor = TRUE)
    prof.data$grp = 186
    
    colnames(prof.data) = c("Breast Density", "Family history of BC", "Tumor Size", "Age at first diagnosis")
    row.names(prof.data) <- NULL
    row.names(output) <- NULL
    out.list = list(profile = prof.data[,-5], current_age = start.age, risk = output)
    if (print.output == T) {  return(out.list)}
    else {  invisible(out.list)}
  }
  cbcrisk_mastectomy <- function(profile, start.age, pred.year=5, print.output=T){
    
    ####### start age check #########
    if(start.age <18|start.age >89)
    {stop("current age has to be within 18-89")}
    
    ###### profile check ##########

    if (length(profile)<9|length(profile)>9)
    {stop("Please enter a valid profile (8 variables)")}
    
    if(profile[1]> 2 | profile[1]<1)
    {stop("Please enter a valid code for 'Age at first child birth'")}
    
    if(profile[2]> 3 | profile[2]<1)
    {stop("Please enter a valid code for 'Age at first BC diagnosis'")}
    
    if(profile[3]> 3 | profile[2]<1)
    {stop("Please enter a valid code for 'BMI'")}
    
    if(profile[4]> 4 | profile[4]<1)
    {stop("Please enter a valid code for 'Breast density'")}
    
    if(profile[5]> 2 | profile[5]<1)
    {stop("Please enter a valid code for 'First degree family history of BC'")}
    
    if(profile[6]> 2 | profile[6]<1)
    {stop("Please enter a valid code for 'ER status'")}
    
    if(profile[7]> 2 | profile[7]<1)
    {stop("Please enter a valid code for 'First BC type'")}
    
    if(profile[8]> 2 | profile[8]<1)
    {stop("Please enter a valid code for 'Lobular carcinoma in situ'")}
    
    if(profile[9]> 2 | profile[9]<1)
    {stop("Please enter a valid code for 'Tumor stage'")}
    
    ########### consistency of the age variables check #########
    
    if (profile[2] == 2 & start.age <30)
    {stop("Current age has to be equal to or greater than the age at first breast cancer diagnosis")}
    
    if (profile[2] == 3 & start.age <40)
    {stop("Current age has to be equal to or greater than the age at first breast cancer diagnosis")}
    
    if (profile[1] == 2 & start.age <40)
    {stop("Age at first child birth has to be equal to or less than the current age, If there is no birth before current age please choose Nulliparous")}
    
    
    #######################################
    
    ########### To avoid installing a relatively large number of dependencies, this package uses squeezeBlanks and recode functions from car package by copying the functions directly from the package. 
    squeezeBlanks <- function(text){
      gsub(" *", "",  text)
    }
    
    recode <- function (var, recodes, as.factor, as.numeric = TRUE, levels) 
    {
      lo <- -Inf
      hi <- Inf
      recodes <- gsub("\n|\t", " ", recodes)
      recode.list <- rev(base::strsplit(recodes, ";")[[1]])
      is.fac <- is.factor(var)
      if (missing(as.factor)) 
        as.factor <- is.fac
      if (is.fac) 
        var <- as.character(var)
      result <- var
      for (term in recode.list) {
        if (0 < length(grep(":", term))) {
          range <- base::strsplit(base::strsplit(term, "=")[[1]][1], ":")
          low <- try(eval(parse(text = range[[1]][1])), silent = TRUE)
          if (class(low) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 low)
          }
          high <- try(eval(parse(text = range[[1]][2])), silent = TRUE)
          if (class(high) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 high)
          }
          target <- try(eval(parse(text = base::strsplit(term, "=")[[1]][2])), 
                        silent = TRUE)
          if (class(target) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 target)
          }
          result[(var >= low) & (var <= high)] <- target
        }
        else if (0 < length(grep("^else=", squeezeBlanks(term)))) {
          target <- try(eval(parse(text = base::strsplit(term, "=")[[1]][2])), 
                        silent = TRUE)
          if (class(target) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 target)
          }
          result[1:length(var)] <- target
        }
        else {
          set <- try(eval(parse(text = base::strsplit(term, "=")[[1]][1])), 
                     silent = TRUE)
          if (class(set) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 set)
          }
          target <- try(eval(parse(text = base::strsplit(term, "=")[[1]][2])), 
                        silent = TRUE)
          if (class(target) == "try-error") {
            stop("\n  in recode term: ", term, "\n  message: ", 
                 target)
          }
          for (val in set) {
            if (is.na(val)) 
              result[is.na(var)] <- target
            else result[var == val] <- target
          }
        }
      }
      if (as.factor) {
        result <- if (!missing(levels)) 
          factor(result, levels = levels)
        else as.factor(result)
      }
      else if (as.numeric && (!is.numeric(result))) {
        result.valid <- na.omit(result)
        opt <- options(warn = -1)
        result.valid <- as.numeric(result.valid)
        options(opt)
        if (!any(is.na(result.valid))) 
          result <- as.numeric(result)
      }
      result
    }
    
    
    ###################################################
    
    
    abs_risk_cbc=function(profile, start.age, pred.year, h1star.seer, h2.seer)
    {
      require(survival) 
      
      ### Preparing the data:
      data=data.frame(t(profile))
      ## Changing the variable names:

      colnames(data) <- c("cat_age1stb_new", "age_1std", "bmi_new", "brd_density", "famhx_bc", "ER", "FTYPNew", "lcis", "stage_new")
      
      ## Recoding the inputs:
      data$cat_age1stb_new=recode(data$cat_age1stb_new,"1='<40/nulli';2='40+'",as.factor=T)
      data$age_1std=recode(data$age_1std,"1='<30';2='30-40';3='40+'",as.factor=T)
      data$bmi_new=recode(data$bmi_new, "1='normal/Undw';2='Ovw';3='Obs'",as.factor=T)
      data$brd_density=recode(data$brd_density, "4='all_fat';3='scattered';2='heterog_dense';1='extrm_dense'",as.factor=T)
      data$famhx_bc=recode(data$famhx_bc,"1='Yes';2='No'",as.factor=T)
      data$ER=recode(data$ER,"1='Neg';2='Pos'",as.factor=T)
      data$FTYPNew=recode(data$FTYPNew,"1='DCIS';2='No DCIS'",as.factor=T)
      data$lcis=recode(data$lcis,"1='No';2='Yes'",as.factor=T)
      data$stage_new=recode(data$stage_new,"1='early_bc';2='adv_primry/metatastatic'",as.factor=T)
      
      data$grp=101
      
      cbc_coeff <- c(0.7275486, 1.1216776, 0.3293037, 0.7178398, 0.4121097, 1.1693814, 1.0885620, 
                     0.6729445, 0.4252677, 0.3001046, 0.3293037, 0.9202828, 0.3784364)
      
      names(cbc_coeff) <- c("cat_age1stb_new40+","age_1std<30","age_1std30-40","bmi_newObs","bmi_newOvw","brd_densityextrm_dense",
                            "brd_densityheterog_dense","brd_densityscattered","famhx_bcYes","ERNeg","FTYPNewDCIS","lcisYes","stage_newadv_primry/metatastatic")
      
      
      x1=model.matrix(fmod_mast,data)
      RR=exp(x1%*%cbc_coeff)       ### Getting the RR#####  
      
      aa=c(18,30,35,40,45,50,55,60,65,70,75,80,85)
      bb=c(29,34,39,44,49,54,59,64,69,74,79,84,89)+1
      
      
      
      ##### Choosing the hazard rates:    
      

      h1star=c(0.002578981, 0.003436104, 0.003289667, 0.003466177, 
               0.003708858, 0.003650760, 0.003734786, 0.004205306, 
               0.004662258, 0.004533459, 0.004893812, 0.004472440, 0.003767342)
      
      other.hazard=c(0.02573451, 0.02597403, 0.02313724, 0.01820899, 0.01642355, 
                     0.01661314, 0.01841396, 0.02081018, 0.02455919, 0.03218279, 
                     0.04547102, 0.06995559, 0.10790407)
      
      
      ar=c(0.9349496, 0.9107625, 0.8927819, 0.8572174, 0.8358139, 0.8071595, 0.7968251, 
           0.8002398, 0.7733977, 0.7796478, 0.7876976, 0.7333668, 0.7872015)
      # 
      cbc.hazard=h1star*(1-ar)
      
      st=findInterval(start.age,aa)
      
      if((start.age+pred.year) %in% aa){ 
        
        fn= (findInterval(start.age+pred.year,aa))-1
      } else {fn=findInterval(start.age+pred.year,aa)}
      
      
      surv.cbc=rep(0,length=fn+1)     ##### Survival function vector for cbc
      surv.other=rep(0,length=fn+1)     #####  Survival function vector for non-cbc
      surv.cbc[1]=1
      surv.other[1]=1
      surv.cbc_a=0
      surv.other_a=0
      
      abs_prob=rep(0,length=fn)
      
      ### Calculating the survival function and absolute risk###    
      for (i in 1:fn)
      {
        
        if (start.age>=aa[i] & start.age<=bb[i] & start.age+pred.year<=bb[i]){   #### this condition applies when the prediction ends in the same interval as the current age
          abs_prob[i] = ((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(start.age+pred.year-start.age)))
          
        }
        else if (start.age>=aa[i] & start.age<=bb[i] & start.age+pred.year>bb[i]){ ### This allows the prediction to go further
          abs_prob[i] = ((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(bb[i]-start.age)))
          
          surv.cbc[i+1]=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(bb[i]-aa[i]))   ### Survival function till the previous interval
          surv.other[i+1]= surv.other[i]*exp(-other.hazard[i]*(bb[i]-aa[i]))
          
          surv.cbc_a=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(start.age-aa[i])) #### Survival function till start.age
          surv.other_a=surv.other[i]*exp(-other.hazard[i]*(start.age-aa[i]))     #### Survival function till start.age
        } 
        else if (start.age<aa[i] & start.age+pred.year>aa[i] & start.age+pred.year<=bb[i]){   ### This is for the probability calculated in the last interval of the prediction
          abs_prob[i]=((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(surv.cbc[i]/surv.cbc_a)*(surv.other[i]/surv.other_a)*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(start.age+pred.year-aa[i]))) 
          surv.cbc[i+1]=0
          surv.other[i+1]=0
        }
        else if (start.age<aa[i] & start.age+pred.year>bb[i]){      #### Probability for an interval within the prediction years
          abs_prob[i]=((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(surv.cbc[i]/surv.cbc_a)*(surv.other[i]/surv.other_a)*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(bb[i]-aa[i])))
          
          surv.cbc[i+1]=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(bb[i]-aa[i]))
          surv.other[i+1]= surv.other[i]*exp(-other.hazard[i]*(bb[i]-aa[i]))  
        }
        else 
        {abs_prob[i]=0
        surv.cbc[i+1]=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(bb[i]-aa[i]))
        surv.other[i+1]= surv.other[i]*exp(-other.hazard[i]*(bb[i]-aa[i])) 
        }
        
      } 
      
      result <- c(start.age+pred.year, round(100*sum(abs_prob),2))
      #names(result) <- c("current.age", "nyears", "abs.risk(%)")
      return(result)
    }
    pred_year=seq(pred.year,89,pred.year)
    output=data.frame(matrix(  ,nrow = 0, ncol = 2))
    colnames(output)=c("by age", "CBC risk(%)")
    
    if(start.age + pred.year[1] >89)
    {stop("current age + nyears.pred < 89")}
    
    for (i in 1:length(pred_year))
    {  
      if(start.age+pred_year[i] <=89)
      {output[i,]=abs_risk_cbc(profile,start.age,pred.year=pred_year[i],h1star.seer=T,h2.seer=T)
      }
      else {break}
      
    }
    prof=data.frame(t(profile))
    
    colnames(prof) <- c("cat_age1stb_new", "age_1std", "bmi_new", "brd_density", "famhx_bc", "ER", "FTYPNew", "lcis", "stage_new")
    
    ## Recoding the inputs:
    prof$cat_age1stb_new=recode(prof$cat_age1stb_new,"1='<40/nulli';2='40+'",as.factor=T)
    prof$age_1std=recode(prof$age_1std,"1='<30';2='30-40';else='40+'",as.factor=T)
    prof$bmi_new=recode(prof$bmi_new, "1='normal/underweight';2='overweight';3='obese'",as.factor=T)
    prof$brd_density=recode(prof$brd_density, "4='Almost ent. fat';3='Scattered';2='Heterog dense';1='Extreme dense'",as.factor=T)
    prof$famhx_bc=recode(prof$famhx_bc,"1='Yes';2='No'",as.factor=T)
    prof$ER=recode(prof$ER,"1='Neg';2='Pos'",as.factor=T)
    prof$FTYPNew=recode(prof$FTYPNew,"1='DCIS';2='Invasive'",as.factor=T)
    prof$lcis=recode(prof$lcis,"1='Yes';2='Yes'",as.factor=T)
    prof$stage_new=recode(prof$stage_new,"1='early BC';2='adv. primary/metastatic'",as.factor=T)
    
    
    colnames(prof)=c("Age at first birth", "Age at first BC diagnosis", "Body Mass Index", 
                     "Breast density", "Family history of BC", "ER", "First BC type","LCIS", "Tumor Stage")
    
    row.names(prof) <- NULL
    row.names(output) <- NULL
    
    out.list=list("profile" = prof, "current_age" = start.age, "risk" = output)
    
    if (print.output==T)
    {return(out.list)}
    else {invisible(out.list)}
  }
  
  if(mastectomy == 1 && length(profile)==9){return(cbcrisk_mastectomy(profile, start.age, pred.year, print.output = T))}
  if(mastectomy == 1 && length(profile)!=9){print("Mastectomy provided is 1 (Patient is scheduled to undergo mastectomy and considering undergoing CPM). Please insert a correct patient profile for this patient")}
  if(mastectomy == 0 && race==1 && length(profile)==4){return(cbcriskBlack(profile, start.age, pred.year, print.output = T))}
  if(mastectomy == 0 && race==1 && length(profile)!=4){print("Race provided is 1 (non-Hispanic Black). Please insert a correct patient profile for this race")}
  if(mastectomy == 0 && race==0 && length(profile)==8){return(cbcrisk_nonBlack(profile, start.age, pred.year, print.output = T))}
  if(mastectomy == 0 && race==0 && length(profile)!=8){print("Race provided is 0 (unspecified). Please insert a correct patient profile for this race")}
}

