NOAH.expression <- fread("/Users/why/Desktop/uw/course/2020 fall/biost 544/data/expression_data.csv", header = TRUE, sep = ',')[,-1]
NOAH.clinical <- read.csv("/Users/why/Desktop/uw/course/2020 fall/biost 544/data/clinical_data.csv", header = TRUE,sep = ',')[,-1]
NOAH.exp.keep <- NOAH.expression[,c("centerid", "patid", "TRIO", "CHMP2B")]
NOAH.exp.keepcenterid <- as.numeric(NOAH.exp.keep$centerid)
NOAH.exp.keep$patid <- as.numeric(NOAH.exp.keep$patid)

NOAH <- inner_join(NOAH.exp.keep, NOAH.clinical, by=c("centerid","patid"))

set.seed(1)

NOAH.cv <- NOAH %>% mutate(fold = sample(1:5, nrow(NOAH), replace = TRUE))
evaluate.mse <- function(predictive.model, test.data){
  predictions <- predict(predictive.model, test.data)
  errors <- test.data$invasive_tumor_cells.pct - predictions
  mse <- mean(errors^2)
  return(mse)
}

spans<-seq(1,5,1)
mse.cv <- matrix(0, nrow = length(spans), ncol = 5)
for(i in 1:length(spans)){
  for(k in 1:5){
    pred.mod <- lm(invasive_tumor_cells.pct~poly(TRIO,spans[i],raw=T)+poly(CHMP2B,spans[i],raw=T),
                      data=NOAH.cv %>% filter(fold != k))
                      ## the control arg just allows for extrapolation
    mse.cv[i,k] <- evaluate.mse(pred.mod,
                                      NOAH.cv %>% filter(fold == k))
  }
}
apply(mse.cv,1, mean) ## This line calculates the mean of each row