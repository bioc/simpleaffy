library("methods")

#holds the results of a pairwise comparison
setClass("QCStats",representation(scale.factors="numeric",target="numeric",percent.present="numeric",average.background="numeric",minimum.background="numeric",maximum.background="numeric",bioB="character",bioC="character",bioD="character",creX="character",gapdh3="numeric",gapdhM="numeric",gapdh5="numeric",actin3="numeric",actinM="numeric",actin5="numeric"));

#accessor methods
setGeneric("sfs", function(object) standardGeneric("sfs"))
setMethod("sfs","QCStats",function(object) object@scale.factors)

setGeneric("target", function(object) standardGeneric("target"))
setMethod("target","QCStats",function(object) object@target)

setGeneric("percent.present", function(object) standardGeneric("percent.present"))
setMethod("percent.present","QCStats",function(object) object@percent.present)

setGeneric("avbg", function(object) standardGeneric("avbg"))
setMethod("avbg","QCStats",function(object) object@average.background)

setGeneric("minbg", function(object) standardGeneric("minbg"))
setMethod("minbg","QCStats",function(object) object@minimum.background)

setGeneric("maxbg", function(object) standardGeneric("maxbg"))
setMethod("maxbg","QCStats",function(object) object@maximum.background)

setGeneric("bioB", function(object) standardGeneric("bioB"))
setMethod("bioB","QCStats",function(object) object@bioB)

setGeneric("bioC", function(object) standardGeneric("bioC"))
setMethod("bioC","QCStats",function(object) object@bioC)

setGeneric("bioD", function(object) standardGeneric("bioD"))
setMethod("bioD","QCStats",function(object) object@bioD)

setGeneric("creX", function(object) standardGeneric("creX"))
setMethod("creX","QCStats",function(object) object@creX)

setGeneric("gapdh3", function(object) standardGeneric("gapdh3"))
setMethod("gapdh3","QCStats",function(object) object@gapdh3)

setGeneric("gapdhM", function(object) standardGeneric("gapdhM"))
setMethod("gapdhM","QCStats",function(object) object@gapdhM)

setGeneric("gapdh5", function(object) standardGeneric("gapdh5"))
setMethod("gapdh5","QCStats",function(object) object@gapdh5)

setGeneric("actin5", function(object) standardGeneric("actin5"))
setMethod("actin5","QCStats",function(object) object@actin5)

setGeneric("actinM", function(object) standardGeneric("actinM"))
setMethod("actinM","QCStats",function(object) object@actinM)

setGeneric("actin3", function(object) standardGeneric("actin3"))
setMethod("actin3","QCStats",function(object) object@actin3)

setGeneric("gapdh35", function(object) standardGeneric("gapdh35"))
setMethod("gapdh35","QCStats",function(object) object@gapdh3 - object@gapdh5);

setGeneric("actin35", function(object) standardGeneric("actin35"))
setMethod("actin35","QCStats",function(object) object@actin3 - object@actin5);

setGeneric("gapdh3M", function(object) standardGeneric("gapdh3M"))
setMethod("gapdh3M","QCStats",function(object) object@gapdh3 - object@gapdhM);

setGeneric("actin3M", function(object) standardGeneric("actin3M"))
setMethod("actin3M","QCStats",function(object) object@actin3 - object@actinM);

.nms <- c("atgenomecdf",	
"ath1121501cdf",	
"drosgenome1cdf",	
"hcg110cdf",		
"hgu133acdf",		
"hgu133a2cdf",	
"hgu133atagcdf",	
"hgu133bcdf",		
"hgu133plus2cdf",	
"hgu95acdf",		
"hgu95av2cdf",	
"hgu95bcdf",		
"hgu95ccdf",		
"hgu95dcdf",		
"hgu95ecdf",		
"mgu74acdf",		
"mgu74av2cdf",	
"mgu74bcdf",		
"mgu74bv2cdf",	
"mgu74ccdf",		
"mgu74cv2cdf",	
"mouse4302cdf",	
"mouse430acdf",		
"mouse430a2cdf",	
"mouse430b2cdf",		
"mu11ksubacdf",	
"mu11ksubbcdf",	
"rae230acdf",		
"rae230bcdf",		
"rgu34acdf",		
"rgu34bcdf",		
"rgu34ccdf",
"ygs98cdf")		

.alpha1 <- c(0.04,	
	     0.05,	
	     0.04,	
	     0.04,	
	     0.05,	
	     0.05,	
	     0.05,	
	     0.05,	
	     0.05,	
	     0.04,	
	     0.04,	
	     0.04,	
	     0.04,	
	     0.04,	
	     0.04,	
	     0.04,	
	     0.04,	
	     0.04,	
	     0.04,	
	     0.04,	
	     0.04,	
	     0.05,	
	     0.05,	
	     0.05,	
	     0.05,	
	     0.04,	
	     0.04,	
	     0.05,	
	     0.05,	
	     0.04,	
	     0.04,	
	     0.04,
	     0.04)	

.alpha2 <- c(0.06,	
	     0.065,	
	     0.06,	
	     0.06,	
	     0.065,	
	     0.065,	
	     0.065,	
	     0.065,	
	     0.065,	
	     0.06,	
	     0.06,	
	     0.06,	
	     0.06,	
	     0.06,	
	     0.06,	
	     0.06,	
	     0.06,	
	     0.06,	
	     0.06,	
	     0.06,	
	     0.06,	
	     0.065,	
	     0.065,	
	     0.065,	
	     0.065,	
	     0.06,	
	     0.06,	
	     0.065,	
	     0.065,	
	     0.06,	
	     0.06,	
	     0.06,
             0.06)

.actin3 <- c("AFFX-r2-At-Actin-3_s_at",		
	     "AFFX-r2-At-Actin-3_s_at",		
	     "AFFX-Dros-ACTIN_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-HSAC07/X00351_3_at",		
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX-b-ActinMur/M12481_3_at",	
	     "AFFX_Rat_beta-actin_3_at",	
	     "AFFX_Rat_beta-actin_3_at",	
	     "AFFX_Rat_beta-actin_3_at",	
	     "AFFX_Rat_beta-actin_3_at",	
	     "AFFX_Rat_beta-actin_3_at",
             "AFFX-YFL039C3_at")


.actinM <- c("AFFX-r2-At-Actin-M_s_at",		
	     "AFFX-r2-At-Actin-M_s_at",		
	     "AFFX-Dros-ACTIN_M_r_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-HSAC07/X00351_M_at",		
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX-b-ActinMur/M12481_M_at",	
	     "AFFX_Rat_beta-actin_M_at",	
	     "AFFX_Rat_beta-actin_M_at",	
	     "AFFX_Rat_beta-actin_M_at",	
	     "AFFX_Rat_beta-actin_M_at",	
	     "AFFX_Rat_beta-actin_M_at",
             "AFFX-YFL039CM_at")


.actin5 <- c("AFFX-r2-At-Actin-5_s_at",		
	     "AFFX-r2-At-Actin-5_s_at",		
	     "AFFX-Dros-ACTIN_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-HSAC07/X00351_5_at",		
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX-b-ActinMur/M12481_5_at",	
	     "AFFX_Rat_beta-actin_5_at",	
	     "AFFX_Rat_beta-actin_5_at",	
	     "AFFX_Rat_beta-actin_5_at",	
	     "AFFX_Rat_beta-actin_5_at",	
	     "AFFX_Rat_beta-actin_5_at",
             "AFFX-YFL039C5_at")

.gapdh3 <- c("AFFX-Athal-GAPDH_3_s_at",		
	     "AFFX-Athal-GAPDH_3_s_at",		
	     "AFFX-Dros-GAPDH_3_at",		
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-HUMGAPDH/M33197_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX-GapdhMur/M32599_3_at",	
	     "AFFX_Rat_GAPDH_3_at",		
	     "AFFX_Rat_GAPDH_3_at",		
	     "AFFX_Rat_GAPDH_3_at",		
	     "AFFX_Rat_GAPDH_3_at",		
	     "AFFX_Rat_GAPDH_3_at",
             "AFFX-YER022w3_at")		

.gapdhM <- c("AFFX-Athal-GAPDH_M_s_at",		
	     "AFFX-Athal-GAPDH_M_s_at",		
	     "AFFX-Dros-GAPDH_M_at",		
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-HUMGAPDH/M33197_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX-GapdhMur/M32599_M_at",	
	     "AFFX_Rat_GAPDH_M_at",		
	     "AFFX_Rat_GAPDH_M_at",		
	     "AFFX_Rat_GAPDH_M_at",		
	     "AFFX_Rat_GAPDH_M_at",		
	     "AFFX_Rat_GAPDH_M_at",
	     "AFFX-YER022w3_at")		

.gapdh5 <- c("AFFX-Athal-GAPDH_5_s_at",		
	     "AFFX-Athal-GAPDH_5_s_at",		
	     "AFFX-Dros-GAPDH_5_at",		
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-HUMGAPDH/M33197_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX-GapdhMur/M32599_5_at",	
	     "AFFX_Rat_GAPDH_5_at",		
	     "AFFX_Rat_GAPDH_5_at",		
	     "AFFX_Rat_GAPDH_5_at",		
	     "AFFX_Rat_GAPDH_5_at",		
	     "AFFX_Rat_GAPDH_5_at",
	     "AFFX-YER022w3_at")		

.biob <- c("AFFX-BioB-3_at",      
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-r2-Ec-bioB-3_at",
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",	    
	    "AFFX-BioB-3_at",
            "AFFX-BioB-3_at")	  

.bioc <- c("AFFX-BioC-3_at",      
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-r2-Ec-bioC-3_at",
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",	    
	    "AFFX-BioC-3_at",
            "AFFX-BioC-3_at")	    

.biod <- c("AFFX-BioD-3_at",      
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioD-3_at",	    
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-r2-Ec-bioD-3_at",
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",	    
	    "AFFX-BioDn-3_at",
            "AFFX-BioDn-3_at")

.crex <- c("AFFX-CreX-3_at",      
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-r2-P1-cre-3_at",
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",	    
	    "AFFX-CreX-3_at",
	    "AFFX-CreX-3_at")	    

names(.alpha1) <- .nms
names(.alpha2) <- .nms

names(.actin3) <- .nms
names(.actinM) <- .nms
names(.actin5) <- .nms

names(.gapdh3) <- .nms
names(.gapdhM) <- .nms
names(.gapdh5) <- .nms

names(.biob) <- .nms
names(.bioc) <- .nms
names(.biod) <- .nms
names(.crex) <- .nms


getTao <- function(name) {
  0.015;
}

getAlpha1 <- function(name) {
  .alpha1[name]
}

getAlpha2 <- function(name) {
  .alpha2[name]
}

getActin3 <- function(name) {
  .actin3[name]
}

getActinM <- function(name) {
  .actinM[name]
}

getActin5 <- function(name) {
  .actin5[name]
}

getGapdh3 <- function(name) {
  .gapdh3[name]
}

getGapdhM <- function(name) {
  .gapdhM[name]
}

getGapdh5 <- function(name) {
  .gapdh5[name]
}

getBioB <- function(name) {
  .biob[name]
}


getBioC <- function(name) {
  .bioc[name]
}


getBioD <- function(name) {
  .biod[name]
}


getCreX <- function(name) {
  .crex[name]
}




qc.affy <-function(unnormalised,normalised=NULL,logged=TRUE,
                   tau=getTao(cleancdfname(cdfName(unnormalised))),
                   alpha1=getAlpha1(cleancdfname(cdfName(unnormalised))),
                   alpha2=getAlpha2(cleancdfname(cdfName(unnormalised))),
	           bioB=getBioB(cleancdfname(cdfName(unnormalised))),
                   bioC=getBioC(cleancdfname(cdfName(unnormalised))),
                   bioD=getBioD(cleancdfname(cdfName(unnormalised))),
                   creX=getCreX(cleancdfname(cdfName(unnormalised))),
                   gapdh3=getGapdh3(cleancdfname(cdfName(unnormalised))),
                   gapdhM=getGapdhM(cleancdfname(cdfName(unnormalised))),
                   gapdh5=getGapdh5(cleancdfname(cdfName(unnormalised))),
                   actin3=getActin3(cleancdfname(cdfName(unnormalised))),
                   actinM=getActinM(cleancdfname(cdfName(unnormalised))),
                   actin5=getActin5(cleancdfname(cdfName(unnormalised)))) {
  if(is.null(normalised)) {
    normalised <- call.exprs(unnormalised,"mas5");
  }


  if(!is.element(cleancdfname(cdfName(unnormalised)),.nms)) {
	stop(paste("I'm sorry I do not know about chip type:",cleancdfname(cdfName(unnormalised))))
  }

  x <- exprs(normalised);

  det <- detection.p.val(unnormalised,tau=tau,alpha1=alpha1,alpha2=alpha2);

  dpv<-apply(det$call,2,function(x) { 
            x[x!="P"] <- 0;
	    x[x=="P"] <- 1;
            x<-as.numeric(x);			       
            return(100 * sum(x)/length(x));
  });

  sfs    <- normalised@description@preprocessing$sfs;
  target <- normalised@description@preprocessing$tgt;

  if(!logged) { x <- log2(x); }

  bgsts<-.bg.stats(unnormalised)$zonebg

  meanbg <- apply(bgsts,1,mean);

  minbg  <- apply(bgsts,1,min);

  maxbg  <- apply(bgsts,1,max);

  stdvbg <- sqrt(apply(bgsts,1,var));
  if(!gapdh3 == "NA") {
    gapdh3 <- x[gapdh3,];
    gapdhM <- x[gapdhM,];
    gapdh5 <- x[gapdh5,];
  }
  else {
    print("Warning couldn't find gapdh probesets... setting ratio to 1");
    gapdh3 <- 1;
    gapdhM <- 1;
    gapdh5 <- 1;
  }
  if(!gapdh3 == "NA") {
    actin3 <- x[actin3,];
    actinM <- x[actinM,];
    actin5 <- x[actin5,];
  }
  else {
    print("Warning couldn't find actin probesets... setting ratio to 1");
    actin3 <- 1;
    actinM <- 1;
    actin5 <- 1;
  }

  bioB <- det$call[bioB,];
  bioC <- det$call[bioC,];
  bioD <- det$call[bioD,];
  creX <- det$call[creX,];

  return(new("QCStats",scale.factors=sfs,target=target,percent.present=dpv,average.background=meanbg,minimum.background=minbg,maximum.background=maxbg,
              gapdh3=gapdh3,gapdhM=gapdhM,gapdh5=gapdh5,actin3=actin3,actinM=actinM,actin5=actin5,bioB=bioB,bioC=bioC,bioD=bioD,creX=creX));
}


setGeneric("qc", function(unnormalised,...) standardGeneric("qc"))
setMethod("qc","AffyBatch",function(unnormalised,...) qc.affy(unnormalised,...)) 

setGeneric("getQCParams", function(x) standardGeneric("getQCParams") )

setMethod("getQCParams","AffyBatch",function(x) {
  n <- cleancdfname(cdfName(x))

  res <- list(getGapdh3(n),getGapdhM(n),getGapdh5(n),getActin3(n),getActinM(n),getActin5(n),getBioB(n),getBioC(n),getBioD(n),getCreX(n),getAlpha1(n),getAlpha2(n),getTao(n))
  names(res) <- c("Gapdh3","GapdhM","Gapdh5","Actin3","ActinM","Actin5","BioB","BioC","BioD","CreX","Alpha1","Alpha2","Tao")
  return(res)
})


.bg.stats <- function(unnormalised, grid=c(4,4)) {
pms         <- unlist(pmindex(unnormalised))
mms         <- unlist(mmindex(unnormalised))
all         <- c(pms,mms)
intensities <- exprs(unnormalised)
rws <- nrow(unnormalised)
cls <- ncol(unnormalised)
zonebg <- c();
zonesd <- c();
for(no in 1:length(unnormalised)){
  this.array <- intensities[,no];
  result <- .C("bgmas",as.integer(as.vector(all)),as.integer(length(all)),
       as.double(as.vector(this.array)),as.integer(length(this.array)),
       as.integer(rws),
       as.integer(cls),
       as.integer(grid[1]),as.integer(grid[2]),
       zonebg=double(grid[1] * grid[2]),
       zonesd=double(grid[1] * grid[2]),corrected=double(length(this.array)),PACKAGE="simpleaffy");
  zonesd <- rbind(zonesd, result$zonesd);
  zonebg <- rbind(zonebg, result$zonebg);
  }
  colnames(zonesd) <- paste("zone",1:16,sep=".");
  colnames(zonebg) <- paste("zone",1:16,sep=".");
  rownames(zonesd) <- sampleNames(unnormalised);
  rownames(zonebg) <- sampleNames(unnormalised);
  return(list(zonebg=zonebg,zonesd=zonesd))
}



plot.qc.stats<-function(x,fc.line.col="black",sf.ok.region="light blue",chip.label.col="black",sf.thresh = 3.0,gdh.thresh = 3.0,ba.thresh = 3.0,present.thresh=10,bg.thresh=20,label=NULL,...) {

  sfs    <- log2(sfs(x))

  n      <- length(sfs)

  meansf <- mean(sfs)

  dpv <- percent.present(x)
  dpv <- (round(100*dpv))/100;

  abg <- avbg(x)
  abg <- (round(100*abg))/100;
	
  sfs <- sfs + 6.0;
  if(is.null(label)) { label <- 1:n }
  col=c("red","green");
  d1 <- 0.0;
  d2 <- 0.0;
  d3 <- 0.0;

  for(i in 1:n) {
    for(j in 1:n) { 
      d1 <- max(abs(sfs[i] - sfs[j]),d1);
      d2 <- max(abs(dpv[i] - dpv[j]),d1);
      d3 <- max(abs(abg[i] - abg[j]),d3);
    }
  }
  i<-1:101;
  arc <- 2 * pi / 100;
  x1 <- (7.5 +  meansf) * cos(i * arc);
  y1 <- (7.5 +  meansf) * sin(i * arc);
  x2 <- (4.5 +  meansf) * cos(i * arc);
  y2 <- (4.5 +  meansf) * sin(i * arc);

  plot(x1[1],y1[1],pch=".",xlim=range(-12,12),ylim=range(-12,12),xaxt="n",yaxt="n",xlab="",ylab="",col=sf.ok.region,...);

  polygon(c(x1,rev(x2),x1[1]),c(y1,rev(y2),y1[1]),col=sf.ok.region,border=sf.ok.region);

  x1 <- 3 * cos(i * arc);
  y1 <- 3 * sin(i * arc);
  points(x1,y1,pch=".",col=fc.line.col);
  x1 <- 6 * cos(i * arc);
  y1 <- 6 * sin(i * arc);
  points(x1,y1,pch=".",col=fc.line.col);
  x1 <- 9 * cos(i * arc);
  y1 <- 9 * sin(i * arc);
  points(x1,y1,pch=".",col=fc.line.col);
  text(0,3,"-3",col=fc.line.col)
  text(0,6,"0",col=fc.line.col)
  text(0,9,"+3",col=fc.line.col)
  arc <- 2 * pi / n;

  gdh <- gapdh35(x);
  ba  <- actin35(x)
  bb  <- bioB(x)

  for(i in 1:n) {
    if(d1 > sf.thresh) { col = "red" } else {col="blue"}
     x1 <- sfs[i] * cos(i * arc);
     y1 <- sfs[i] * sin(i * arc);
     x2 <- 6 * cos(i * arc);
     y2 <- 6 * sin(i * arc);
     lines(c(x2,x1),c(y2,y1),col=col);
     points(x1,y1,col=col,pch=20);
     text(x1,y1,col=chip.label.col,label=label[i],adj=0.2 * c(cos(i * arc),sin(i * arc)));
     x2 <- (6 + gdh[i]) * cos(i * arc);
     y2 <- (6 + gdh[i]) * sin(i * arc);
     if(gdh[i] > gdh.thresh) { col = "red" } else {col="blue"}	
     points(x2,y2,pch=1,col=col);
     x2 <- (6 + ba[i]) * cos(i * arc);
     y2 <- (6 + ba[i]) * sin(i * arc);
     if(ba[i] > ba.thresh) { col = "red" } else {col="blue"}	
     points(x2,y2,pch=2,col=col);

     if(d2 > present.thresh) { col = "red" } else {col="blue"}
     x2 <- (9 * cos(i * arc));
     y2 <- (9 * sin(i * arc));
     text(x2,y2,label=paste(dpv[i],"%",sep=""),col=col);

     if(d3 > bg.thresh) { col = "red" } else {col="blue"}
     x2 <- (10 * cos(i * arc));
     y2 <- (10 * sin(i * arc));
     text(x2,y2,label=abg[i],col=col);

     if(bb[i]!="P") {
       x2 <- (11 * cos(i * arc));
       y2 <- (11 * sin(i * arc));
       text(x2,y2,label="bioB",col="red");
     }
  }
  legend(-10,10,pch=1:2,c("GAPDH","Beta Actin"))
}

setMethod("plot","QCStats",function(x,y) plot.qc.stats(x,...))
