#' IDEA - NOT FOR USE -
#' 
#' Should baggr include power calculations?
#' Danny Toomey, April 2023
#' 
#' One of the main reasons researchers choose to perform meta-analyses is the potential 
#' to detect effects individual studies cannot by increaseing statistical power. While  
#' this is true of fixed-effects models, it is not necessarily true of random effects 
#' models, which are far more common. 
#'
#' There are multiple ways to evaluate the possible effects of insufficient sample size 
#' or excessive heterogeniety in a random effects model, which are the key variables that 
#' would effect its statistical power. An argument against including a power measure would
#' be that heterogeniety is already given in baggr models and evaluating datasets with low 
#' sample sizes is one of the key benefits of Bayesian analyses. It is possible that a low
#' statistical power is why a user would be interesting in building a Bayesian model in the 
#' first place. 
#'
#' However, it is also possible that a user finds a power analysis informative in explaining 
#' unexpected results. Consider a user is expecting to find significant effects in a Bayesian 
#' model of Yusuf et al 1985, and they are surprised when they do not. However, a power 
#' analysis shows the meta-analysis has a power of 0.26. The user may then infer that although 
#' a Bayesian model may be less sensitive to small sample size, there may still be significant
#' limitations in its representation of findings in the field.
#'
#' The [power_calc_example] function illustrates an example of this use case and proposes
#' an implementation for using power analyses to contextualize baggr models. 
#'
#' In general, a power metric may not be useful to all models. However, in cases where users observe 
#' an unexpected effect, the power of the meta-analysis may provide useful context.
#'
#' References: 
#' 1. 	Quintana, D. How to calculate statistical power for your meta-analysis. Medium. 2017 Jul. 
#' 		Available online: https://towardsdatascience.com/how-to-calculate-statistical-power-for-your-meta-analysis-e108ee586ae8
#' 2.	Jackson D, Turner R. Power analysis for random‐effects meta‐analysis. 
#'		Research synthesis methods. 2017 Sep;8(3):290-302. DOI: 10.1002/jrsm.1240. 
#'
#' @import metafor

power_calc_example <- function(){
	bg <- baggr(yusuf,effect='logOR',silent=TRUE)

	total_es <- (sum(yusuf$ai)/sum(yusuf$n1i-yusuf$ai))/(sum(yusuf$ci)/sum(yusuf$n2i-yusuf$ci))
	es <- abs(log(total_es)*sqrt(3)/pi)
	as <- mean(c(yusuf$n1i,yusuf$n2i))
	mk <- nrow(yusuf)-1
	es_y <- metafor::escalc(measure = "OR",ai=yusuf$ai,ci=yusuf$ci,n1i=yusuf$n1i,n2i = yusuf$n2i,data=yusuf)
	rma_yusuf <- metafor::rma(data=yusuf,yi=es_y$yi,vi=es_y$vi)
	hg <- sqrt(rma_yusuf$H2)
		
	message("\n",
	"Following the article referenced in the documentation, here is an example showing  \n",
	"how a power analysis could be used along side a baggr model for Yusuf et al 1985.\n\n",
	"A baggr analysis of Yusuf 1985 reveals no significant effects: \n\n",
	" --- --- --- \n") 
	print(bg)
	message("\n --- --- --- \n\n",
	"Which may surprise a researcher in the field who was expecting to observe an effect.\n",
	"However, additional context is given by the power of the meta analysis: \n")
	cat(paste0(" --- Power = ",format(round(power(es,as,mk,hg),3),nsmall=3)," --- \n"))
	message("\n",
	"Which informs the researcher that the meta-analysis does not meet the conventional 80%\n",
	"benchmark for power, and so may not be indicative of other findings broadly.\n\n",
	"While this type of analysis may not be available for all datasets passed to baggr, \n",
	"(e.g. datasets like Rubin's 8 Schools which precalculates tau and SE), it may be \n",
	"useful to provide a message to the user if a primary dataset is available, such as\n",
	"the Yusuf et al dataset.\n\n",
	"An example implementation is given by the function `prep_ma_with_power`, which performs\n",
	"a power analysis and prepares the dataset for baggr via `prepare_ma`. If power is below 0.8, \n",
	"a warning message is given to inform the user. See the warning messages below for an exmaple\n",
	"using Yusuf et al.\n"
	)
	prep <- prep_ma_with_power(yusuf)

}

power <- function(es,as,mk,hg){
	eq1 <- ((as+as)/((as)*(as))) + ((es^2)/(2*(as+as)))
	eq2 <- hg*(eq1)
	eq3 <- eq2+eq1
	eq4 <- eq3/mk
	eq5 <- (es/sqrt(eq4))
	Power <- (1-pnorm(1.96-eq5))
	return(Power)
}

prep_ma_with_power <- function(data){
	total_es <- (sum(data$ai)/sum(data$n1i-data$ai))/(sum(data$ci)/sum(data$n2i-data$ci))
	es <- abs(log(total_es)*sqrt(3)/pi)
	as <- mean(c(data$n1i,data$n2i))
	mk <- nrow(data)-1
	es_y <- metafor::escalc(measure = "OR",ai=data$ai,ci=data$ci,n1i=data$n1i,n2i = data$n2i,data=data)
	rma_data <- metafor::rma(data=data,yi=es_y$yi,vi=es_y$vi)
	hg <- sqrt(rma_data$H2)

	power <- power(es,as,mk,hg)
	prep <- prepare_ma(data,effect='logOR')
	
	if(power<0.8){
		warning(paste0("The power of your meta anlysis is below 80%, (P = ",format(round(power,3),nsmall=3),"), \n",
		"which may limit the generalizability of your model. \n"))
	} else {
		cat(paste0("Power = ",format(round(power,3),nsmall=3),"\n"))
	}
	return(prep)
	
}


