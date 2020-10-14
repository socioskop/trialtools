
### load libraries
library(foreach)
library(data.table)
library(magrittr)
library(lubridate)

### raw data pull
grab.it <- function(path, reg, years, pop, vars=NULL, s.id=NULL, show.vars=F){
  # reformat args
  if (substr(path, nchar(path), nchar(path))!="/"){path <- paste0(path, "/")}
  if (is.null(s.id)){s.id <- "pnr"}
  
  # get file list
  files <- list.files(path, pattern=reg)
  files <- files[grepl(paste0(as.character(years), collapse="|"), files)]
  
  # pull to list
  ff <- foreach(file = files) %do% {
    f <- haven::read_sas(paste0(path, file))
    colnames(f) <- tolower(colnames(f))
    if (show.vars==T & file==files[1]){print(colnames(f))}
    f$year <- paste0(stringr::str_extract_all(file, "[0-9]")[[1]][1:4], collapse="")
    f <- f[f[[s.id]] %in% pop,]
    if (!is.null(vars)){f <- f[,vars]}
    f
  }
  f <- dplyr::bind_rows(ff)
  return(f)
}


### formatting: fixed string format with padding for numericals
form.it <- function(x, digits=3, perc=FALSE){
  x <- format(round(x, digits), nsmall=digits)
  if (perc==T){x <- paste0(x, "%")}
  return(x)
}

### reduce: takes last observation before time
last.obs <- function(data, id, time) {
  d <- data.table::data.table(data[!is.na(data[[id]]),])
  d <- d[order(get(id), get(time)),]
  d[ , N :=   .N, by=get(id)]
  d[ , n := 1:.N, by=get(id)]
  d <- d[n==N]
  return(as.data.frame(d)[,!colnames(d) %in% c("n", "N")])
}

### reduces data by some function within a certain time fram (from index looking forward 'interval')
reduce.forward <- function(data, id="pnr", index="index", time="date", vars, interval, func="sum", cens.date=NULL){
  d <- data[!is.na(data[[id]]),]
  
  if (!is.null(cens.date)){
    d[[cens.date]] <- lubridate::as_date(ifelse(is.na(d[[cens.date]]), lubridate::as_date(round(d[[index]]+interval)), d[[cens.date]]))
    d <- d[d[[time]]>d[[index]] & d[[time]]<(lubridate::as_date(round(d[[cens.date]]))) & !is.na(d[[index]]) & !is.na(d[[time]]),]
  } else {
    d <- d[d[[time]]>d[[index]] & d[[time]]<(lubridate::as_date(round(d[[index]]+interval))) & !is.na(d[[index]]) & !is.na(d[[time]]),]
  }
  
  if (length(vars)==1) {is.na(d[,vars]) <- 0} else {
    d[,vars] <- lapply(d[,vars], function(x) replace(x, is.na(x), 0))
  }
  d <- aggregate(as.formula(paste0(".~", id, "+", index)), data=d[,c(id, index, vars)], FUN=func, na.action = na.omit)
  return(d)
}

### Twoways: generates twoway statistics for categorical and numerical variables (+ n-count descriptives)
tt.tabulate <- function(data, xs, treat, weight=NULL, num=NA, cat=NA, bin=NA, num.test="auto", na.is="collapsed", no_wgt=NULL){
  t <- data.frame()
  for (x in xs){
    ## verbose
    print(x)
    
    ## dealing with missings
    if (na.is=="collapsed" & !x %in% num){
      
      ## testing if missing cells are based on fewer than 5 observations
      if (x %in% bin) {
        test.tab <- table(data[[x]]>=1, data[[treat]], useNA="ifany")
      } else {
        test.tab <- table(data[[x]], data[[treat]], useNA="ifany")
      }
      issue <- (min(test.tab[test.tab>0])<5)
      
      ## if this is the case collapsing missing with first category
      if (issue==T){
        na.convert <- names(table(data[[x]]))[1]
        data[[x]][is.na(data[[x]])] <- ifelse(class(data[[x]]) %in% c("numeric", "integer"), as.numeric(na.convert), na.convert)
      }
    }
    
    #if (x %in% ttests){test.type <- "ttest"} else {test.type <- "wilcox"}
    
    ## adding row to summary table depending on variable type
    wgt <- switch((x %in% no_wgt | is.null(weight))+1, weight, NULL)
    if (x %in% c("n", "n_wgt")){try(t <- dplyr::bind_rows(t, twoway.n  (data=data, x, treat, weight=wgt)))} 
    if (x %in% num){try(t <- dplyr::bind_rows(t, twoway.num(data=data, x, treat, weight=wgt, digit.m = 2, digit.sd = 2, test = num.test)))} 
    if (x %in% cat){try(t <- dplyr::bind_rows(t, twoway.chi(data=data, x, treat, weight=wgt)))} 
    if (x %in% bin){try(t <- dplyr::bind_rows(t, twoway.chi(data=data, x, treat, weight=wgt, bin=T)))}
  }
  t <- t%>% dplyr::mutate_if(is.ok, function(x) as.numeric(as.character(x)))
  t <- t%>% dplyr::mutate_if(is.factor, function(x) as.character(x))
  return(t)
}

twoway.n <- function(data, x, group, weight=NULL){
  if (is.null(weight)){data$weight <- 1} else {data$weight <- data[[weight]]; data[[weight]] <- NULL}
  
  data <- data[!is.na(data[[group]]),]
  #tab <- table(data[[x]], data[[group]], useNA="ifany")
  tab <- round(questionr::wtd.table(data[[x]], data[[group]], weights = data$weight, na.show = F))
  
  groups <- colnames(tab)
  if (nrow(tab)>1){stop("two many levels")}
  tab <- c(rbind(tab, form.it(prop.table(tab, 2)*100, 1)))
  tab <- c(x, NA, tab, NA)
  tab <- as.data.frame(t(tab))
  colnames(tab) <- c("var", "level", rbind(paste0(groups, ".mean/n"), paste0(groups, ".sd/%")), "p")
  return(tab)
}

twoway.chi <- function(data, x, group, weight=NULL, bin=F){
  if (is.null(weight)){data$weight <- 1} else {data$weight <- data[[weight]]; data[[weight]] <- NULL}
  
  data <- data[!is.na(data[[group]]),c(x, group, "weight")]
  if (bin==T) {data[[x]] <- as.numeric(data[[x]]>=1)}
  #tab <- table(data[[x]], data[[group]], useNA="ifany")
  show.na <- (sum(is.na(data[[x]]))>0)
  tab <- round(questionr::wtd.table(data[[x]], data[[group]], weights = data$weight, na.show = show.na))
  tab.noround <- questionr::wtd.table(data[[x]], data[[group]], weights = data$weight, na.show = F)
  
  tab         <- tab[rowSums(tab)>0,]
  tab.noround <- tab.noround[rowSums(tab.noround)>0,]
  
  groups <- colnames(tab)
  if (nrow(tab)<2){warning("two few levels")}
  data[[x]][is.na(data[[x]])] <- "N/A"
  p <- form.it(chisq.test(tab.noround[rowSums(tab.noround)>0,])$p.value, 3)
  levels <- rownames(tab)
  if (is.na(levels[length(levels)]) & sum(is.na(levels))==1){levels[length(levels)] <- "N/A"}
  tab <- matrix(rbind(tab, form.it(prop.table(tab, 2)*100, 1)), nrow=nrow(tab)) #[,c(seq(1, ncol(tab)*2, 2), seq(2, ncol(tab)*2, 2))]
  if (nrow(tab)>=2){filling <- c(rep(NA, nrow(tab)-1), p)} else {filling <- p}
  if (nrow(tab)>=2 & bin==F){varname   <- c(x, rep(NA, nrow(tab)-1))} else {varname <- rep(x, nrow(tab))}
  tab <- as.data.frame(cbind(varname, levels, tab, filling))
  colnames(tab) <- c("var", "level", rbind(paste0(groups, ".mean/n"), paste0(groups, ".sd/%")), "p")
  if (bin==T){tab <- tab[nrow(tab),]}
  tab[nrow(tab), "test"] <- "chi^2"
  return(tab)
}

twoway.num <- function(data, x, group, weight=NULL, digit.m=1, digit.sd=1, inf=FALSE, test=NULL, shapiro.p=.01){
  if (is.null(weight)){data$weight <- 1; weight <- "weight"} else {data$weight <- data[[weight]]; data[[weight]] <- NULL; weight <- "weight"}
  
  data <- data[!is.na(data[[group]]),]
  groups <- unique(na.omit(data[[group]]))
  tab <- rep(NA, 2*length(groups)+2)
  tab[1] <- x
  tab[2] <- NA
  k <- 3
  for (i in 1:length(groups)){
    tab[k+0] <- form.it(matrixStats::weightedMean  (data[[x]][data[[group]]==groups[i]], w=data$weight[data[[group]]==groups[i]], na.rm=T), digit.m)
    tab[k+1] <- form.it(matrixStats::weightedSd    (data[[x]][data[[group]]==groups[i]], w=data$weight[data[[group]]==groups[i]], na.rm=T), digit.sd)
    tab[k+2] <- form.it(matrixStats::weightedMedian(data[[x]][data[[group]]==groups[i]], w=data$weight[data[[group]]==groups[i]], na.rm=T), digit.sd)
    k <- k+3
  }
  
  if (test=="auto"){
    if(shapiro.test(resid(glm(paste0(x, "~", group), data=data)))$p.value<shapiro.p){
      test <- "rank"
    } else {test <- "ttest"}
  }
  
  if (test=="rank"){
    tab[length(tab)+1] <- tryCatch({form.it(
      ordinal:::anova.clm(ordinal::clm(paste0("as.factor(", x, ")~as.factor(", group, ")"), weights = data[["weight"]], data=data, link = "logit"))$`Pr(>Chisq)`, 3)}, error=function(err) NA)
  } else if (test=="ttest"){
    yy <- data[,x]
    xx <- data[,group]
    weights <- data[,weight]
    m <- glm(yy ~ xx, weights = weights)
    tab[length(tab)+1] <- tryCatch({form.it(as.numeric(na.omit(lmtest::lrtest(m)$`Pr(>Chisq)`)), 3)}, error=function(err) NA)
  }
  
  tab <- as.data.frame(t(tab))
  colnames(tab) <- c("var", "level", rbind(paste0(groups, ".mean/n"), paste0(groups, ".sd/%"), paste0(groups, ".median")), "p")
  tab$test <- test
  return(tab)
}

is.ok <- function(x) {
  (sum(!is.na(x))==(sum(!is.na(as.numeric(as.character(x))))) & length(x[grepl(paste0(letters, collapse="|"), x)])==0)
}

### Survival analysis and plotting
ggsplot <- function(model, titles, table=F){
  p <- survminer::ggsurvplot(
    fit=model, conf.int=T, palette=c(pals::okabe()[2:7]),
    title=titles, legend=c(.1, .25), legend.labs=c("Controls", "MeHAS"),
    risk.table=table, break.time.by=25, #palette=rep("black", 2),
    size=.66, linetype=c("solid", "dashed"), censor=F,
    xlab="Days", ylab="Proportion unemployed", dots=T, xlim=c(0, 371))
  
  #p$plot <- p$plot+labs(x=NULL)+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  return(p)
}

splot <- function(model, title){
  rms::survplot(fit=model, conf=c("none", "bands", "bars")[2],
           xlab="Days", ylab="Proportion of unemployed", xlim=c(0, 371),
           label.curves=list(keys="lines"),
           levels.only=T, abbrev.label=F, 
           loglog=F, logt=F, time.inc=25, dots=T,
           nrisk=F, cex.n.risk=.5, adj.subtitle=F, main=title)
}

