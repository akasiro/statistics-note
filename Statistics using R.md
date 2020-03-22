
[toc]

# Chapter 1 Linear regression [^ 1]
## 1. Basic Model and Definition
* Population form
  $$ y_i  = {\beta}_0 + {\beta}_i x_i + {\varepsilon}_i $$

* Estimate model
  $$ \hat{y_i} = {\beta}_0 + {\beta}_i x_i $$

* Residual

  $$e_i = y_i - \hat{y_i}$$
## 2. Distributional assumptions underlying linear regression
* **Assumption 1**: the relationship between *y~i~* and *x~i~* is linear
* **Assumption 2**: variance of residuals is constant regardless of *x~i~*  [*homoscedasticity and the homogeneity of error variance assumption*]
* **Assumption 3**: residuals are normally distributed
* **Assumption 4**: *x* measured without error and the error is unrelated to model error
* **Assumption 5**: residuals for any two individuals are assumed to be independent of one another.

## 3. An example regression with R

1. Data source: ***swiss*** the built-in dataset of R
2. code

```R
library(stargazer) #package offer a formated output
library(car) #package using for plot

# 0. a quick view of data
head(swiss)
# 1. descriptive table
stargazer(swiss)
# 2. a quick kde distribution plot of dependent variable
d <- density(swiss$Fertility)
plot(d)
# 3. Build a multiple regression model
Model1 <- lm(Fertility~Examination + Education, data=swiss)
# 4. collinearity check
vif(Model1)
# 5. summary of model
summary(Model1)
# 6. check the entire F test
anova(Model1)
# 7. Checking regression assumption
# 7.1 Residual distribution (linearity assumption and homogeneity of variance) and interaction effect
residualPlots(Model1)
# 7.2 residual normal distribution assumption check
qqPlot(Model1)
# 8 format output
stargazer(Model1, type = 'text')
#stargazer(Model1, type = 'html', out = 'Model1.htm')
```

3. Result interpretation





[^ 1]:most ideas from *Finch, W. H., Bolin, J. E., & Kelley, K. (2019). Multilevel Modeling Using R (2nd Edition ed.): Chapman and Hall/CRC.*


 <div STYLE="page-break-after: always;"></div>
# Chapter 2 Multilevel linear modeling
## Two-level Multilevel linear modeling
### 1. Why and when we use multilevel linear model
* what problem it solve (idea from Hox, Moerbeek and Schoot, 2017 chapter1)
***Multilevel research***: *Research that individuals interact with the social context to which whty belong,  that individuals are influenced by the contexts or groups to which they belong, and that those groups are in turn influenced by the individuals who make up that group*
***Multilevel research*** concerns a population with a hierarchical structure
***Multilevel problem***: *a problem that concerns the relationships between variables that are measured at a number of different hierarchical levels.*
* why use MLM
1. solve intraclass correlation in the hierarchical data structure sample
> the average correlation (intraclass correlation) between variables in individual level from the same group will be higher than that from different group, which violates the assumption of independence of the observations and causes a strong biases (Walsh, 1947) and produces spurious 'significant' results.
> 


* what's the advantage compared to other modeling
* In what context we use it
a hierarchical data set consisting of individuals nested within groups, with one single outcome or response variable that is measured at the lowest level, and explanatory variables at all existing levels.
### 2. What's the assumption
1. level2 residuals are independent between clusters
2. level2 intercepts and coefficients are independent to the level1 residuals
3. level1 residuals are normally distributed and have a constant variance
4. level2 intercept and slopes have a multivariable normal distribution with a constant covariance matrix
### 3. data structure and preclean before model
#### 3.1 data structure

| level2 id | level1 id | level1 dv | level2 iv1 | level2 iv1 | level1 iv1 | level1 iv2 |
| --------- | --------- | --------- | ---------- | ---------- | ---------- | ---------- |
| 1         | 1         |           |            |            |            |            |
| 1| 2||||||

***Depend variable must be level 1***
Example datasets can be found in https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book.git

#### 3.2 centering
variable - mean
* **why**
1.reduce collinearity caused by including an interaction term in a regression model (not just for MLM) (e.g. Iversen, 1991), which is important in MLMs (Wooldridge, 2004)
2. centering provides a potential advantage in terms of interpretation of results

* **how**
1. *grand mean centering*: most commonly used (Bickel, 2007) and in most cases work well (Kreft, de Leeuw, and Aiken, 1995)
2. *group mean centering*: an alternative way
### 4. Theory
#### 4.1 Model
Full model (random slopes model):  $$y = {\gamma}_{00}+ {\gamma}_{10} x_{ij} + U_{0j} + U_{1j}  x_{ij} + {\varepsilon}_{ij}$$
> $$({\gamma}_{00}+ {\gamma}_{10} x_{ij})$$ --- ***fixed effect***
> $$(U_{0j} + U_{1j}  x_{ij})$$ --- ***random effect***

Level 1 model: $$y_{ij} = {\beta}_{0j} + {\beta}_{1j} x + {\varepsilon}_{ij}$$
> $$y_{ij}$$ --- $$y$$ of group $$j$$  individual $$i$$

Level 2 model:
$${\beta}_{0j} = {\gamma}_{00} + U_{0j}$$
$${\beta}_{1j} = {\gamma}_{10} + U_{1j} $$

##### 4.1.2 standard notation of the model (from Hox, Moerbeek and Schoot, 2017 chapter2)
$$Y_{ij} = {\gamma}_{00} + {\gamma}_{p0}X_{pij} + {\gamma}_{0q}Z_{qj} + {\gamma}_{pq}Z_{qj}X_{pij}+{\mu}_{pj}X_{pij}+{\mu}_{0j} + e_{ij}$$

or

$$Y_{ij} = {\gamma}_{00} + \sum_p{\gamma}_{0q}Z_{qj} + \sum_q{\gamma}_{p0}X_{pij} + \sum_p\sum_q{\gamma}_{pq}Z_{qj}X_{pij}+\sum_p{\mu}_{pj}X_{pij}+{\mu}_{0j} + e_{ij}$$

> subscript $$i$$: individual $$j%%
>
> subscript $$j$$: group $$j$$
>
> subscript $$p$$: the iv in the lowest (individual) level
>
> subscript $$q$$: the iv in the highest (group) level

##### Deducing
**Simple linear model**
$$y = {\beta}_0 + {\beta}_1 x + {\varepsilon}$$
> $${\beta}_0$$ --- intercept, conditional mean of $$y$$ when $$x=0$$
> $${\beta}_1$$ --- coefficient
> $${\varepsilon}$$ --- random variation

**Level 1 model**
$$y_{ij} = {\beta}_{0j} + {\beta}_{1j} x + {\varepsilon}_{ij}$$
> $$y_{ij}$$ --- $$y$$ of group $$j$$  individual $$i$$

**Random Intercept**
Level 1 model +
$${\beta}_{0j} = {\gamma}_{00} + U_{0j}$$
=> $$y = {\gamma}_{00} + U_{0j} + {\beta}_1 x + {\varepsilon}_{ij}$$

> $${\gamma}_{00}$$ --- ***fixed effect***, average or general intercept across group
> $$U_{0j}$$ --- ***random effect***, group-specific effect on the intercept
> $${\tau}^2$$ --- variance of $$U_{0j}$$.
> $${\sigma}^2$$ -- variance of $${\varepsilon}$$

**Random intercept: null model**
$$y = {\gamma}_{00} + U_{0j}  + {\varepsilon}_{ij}$$
Usually use for calculate ICC (intraclass correlation)
$$\hat{{\rho}_I} = {{\tau}^2\over {\tau}^2+{\sigma}^2}$$

**Random slopes**
Random intercept model + 
$${\beta}_{1j} = {\gamma}_{10} + U_{1j} $$
=> $$y = {\gamma}_{00}+ {\gamma}_{10} x_{ij} + U_{0j} + U_{1j}  x_{ij} + {\varepsilon}_{ij}$$
> $$({\gamma}_{00}+ {\gamma}_{10} x_{ij})$$ --- ***fixed effect***
> $$(U_{0j} + U_{1j}  x_{ij})$$ --- ***random effect***

#### 4.2 Estimation
1. Maximum likelihood estimation (MLE)
> robust and produces estimates that are asymptotically efficient and consistent (Hox, Moerbeek and Schoot, 2017)

2. Restricted maximum likelihood estimation (REML)
> more accurate with regard to the extimation of variance parameters than is MLE (Kreft & De Leeuw, 1998)
> generally the preferred method, though for testing variance parameters (or any random effect) it is necessary to use MLE (Snijders & Bosker, 1999)

3. Generalized least squares (GLS)
> GLS estimation approximates ML estimates
> faster to compute than MLE
> use when MLE fail to converge
> but is less efficient and the GLS-derived standard errors are inaccurate

4. Generalized estimatiing equations (GEE)
> faster than MLE
> less efficient than MLE

5. Bayesian methods
> deal with multicollinearity (Can et al., 2014)
> deal with non-normality
> deal with smaller sample size on the highest level

6. Bootstrapping
### 5  how to conduct MLM in R
#### 5.0 packages
```R
library(lme4)
library(lmerTest)
```
#### 5.1 Centering
```R
# grand mean centering
var_center <- dataset$var - mean(dataset$var)
```
#### 5.2 Check ICC
```R
# conduct a null model
model0 <- lmer(dv ~ 1 + (1|level2id), data = dataset)
summary(model0)
```
calculate ICC using $$\hat{{\rho}_I} = {{\tau}^2\over {\tau}^2+{\sigma}^2}$$

#### 5.3 modeling
```R
# function lmer(), random effect in parentheses, default estimation REML

# null model
model0 <- lmer(dv ~ 1 + (1|level2id), data = dataset)
# random intercept model 1: level1 iv
model1 <- lmer(dv ~ iv_level1 + (1|level2id), data = dataset)
# random intercept model 2: level2 iv
model2 <- lmer(dv ~ iv_level1 + iv_level2 + (1|level2id), data = dataset)
# random intercept model 3: cross-level interactions
model3 <- lmer(dv ~ iv_level1 + iv_level2 + iv_level1*iv_level2 + (1|level2id), data = dataset)
# random coefficients model 1: basic
## no need to add (1|level2id), it's implicit included, if don't want random intercept use (-1+iv_level1|level2id)
model4 <- lmer(dv ~ iv_level1 + (iv_level1|level2id), data = dataset)
# random coefficients model 2: random slopes are corelated
model5 <- lmer(dv ~ iv1_level1 + iv2_level2 + (iv1_level1+iv2_level2|level2id), data = dataset)
# random coefficients model 3.1: random slopes are uncorelated (estimate intercept seperated)
model6.1 <- lmer(dv ~ iv1_level1 + iv2_level2 + (iv1_level1|level2id) + (iv2_level1|level2id), data = dataset)
# random coefficients model 3.2: random slopes are uncorelated (without separate intercept estimated)
model6.2 <- lmer(dv ~ iv1_level1 + iv2_level2 + (1|level2id) + (-1 + iv1_level1|level2id) + (-1 + iv2_level1|level2id), data = dataset)
```
#### 5.4 model fit comparing
```R
# function anova()
anova(model1, model2)
```
#### 5.5 lme4 and hypothesis testing
As lme4 library doesn't provide p-values, we can test using ***bootstrapped confidence intervals*** to estimate the significance of fixed and random effect in MLM
> How it work
> step 1. 从现在的样本中进行随机抽样，抽取样本数为B的样本（default B = 500)
> step 2. 利用重新得到的样本进行模型估计
> step 3. 重复step1 和 step2得到多组估计值的样本，利用估计值的样本分布估计95%置信区间的估计值
> $$H_0: {\theta} = {\theta}_0$$
> $${\theta}$$: real value of the parameter of interest (e.g. random intercept variance)
> $${\theta}_0$$: the estimate of the parameter in the model(e.g. random intercept variance in summary(model4))
> * For random effect:
> if 0 not in the interval, it means there are random effect or significant
> * For fixed effect:
> if the estimate of the model in the interval, it means significant
```R
# use confint() function and there are three kind of bootstrap options
# 1. percentile bootstrap(perc)
confint(model4, method=c("boot"), boot.type=c("perc"))
# 2. the standard error bootstrap(basic)
confint(model4, method=c("boot"), boot.type=c("basic"))
# 3. the normal bootstrap(norm)
confint(model4, method=c("boot"), boot.type=c("norm"))
# 4
confint(model4, method=c("Wald"))
# 5
confint(model4, method=c("profile"))
```
### 6. how to explain the result
A example using data https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/raw/master/chapter%202/popularity/SPSS/popular2.sav
#### 6.0 load packages
```R
library(foreign) # for read data
library(lme4) # MLM modeling
library(lmerTest) # povide p-value for modeling
```
#### 6.1 read data
```R
popular2 = read.spss("https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/raw/master/chapter%202/popularity/SPSS/popular2.sav", to.data.frame=TRUE)
head(popular2)
```
**output**
| index | pupil | class | extrav | sex | texp | popular | popteach | Zextrav | Zsex | Ztexp | Zpopular | Zpopteach | Cextrav | Ctexp | Csex |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 1 | 5 | girl | 24 | 6.3 | 6 | -0.1703149 | 0.9888125 | 1.486153 | 0.8850133 | 0.66905609 | -0.215 | 9.737 | 0.5 |
| 2 | 2 | 1 | 7 | boy | 24 | 4.9 | 5 | 1.4140098 | -1.0108084 | 1.486153 | -0.1276291 | -0.04308451 | 1.785 | 9.737 | -0.5 |
| 3 | 3 | 1 | 4 | girl | 24 | 5.3 | 6 | -0.9624772 | 0.9888125 | 1.486153 | 0.1616973 | 0.66905609 | -1.215 | 9.737 | 0.5 |
| 4 | 4 | 1 | 3 | girl | 24 | 4.7 | 5 | -1.7546396 | 0.9888125 | 1.486153 | -0.2722923 | -0.04308451 | -2.215 | 9.737 | 0.5 |
| 5 | 5 | 1 | 5 | girl | 24 | 6.0 | 6 | -0.1703149 | 0.9888125 | 1.486153 | 0.6680185 | 0.66905609 | -0.215 | 9.737 | 0.5 |
| 6 | 6 | 1 | 4 | boy | 24 | 4.7 | 5 | -0.9624772 | -1.0108084 | 1.486153 | -0.2722923 | -0.04308451 | -1.215 | 9.737 | -0.5 |

#### 6.2 check ICC
**code**
```
m1 <- lmer(popular~1+(1|class), data=popular2)
summary(m1)
```
**output**
> Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
> Formula: popular ~ 1 + (1 | class)
>   Data: popular2
>
> REML criterion at convergence: 6330.5
>
> Scaled residuals: 
>    Min      1Q  Median      3Q     Max 
> -3.5655 -0.6975  0.0020  0.6758  3.3175 
>
> Random effects:
>  Groups   Name        Variance Std.Dev.
>  class    (Intercept) 0.7021   0.8379  
>  Residual             1.2218   1.1053  
> Number of obs: 2000, groups:  class, 100
> 
> Fixed effects:
>             Estimate Std. Error       df t value Pr(>|t|)    
> (Intercept)  5.07786    0.08739 98.90973    58.1   <2e-16 \*\*\*
> 
> Signif. codes:  0 ‘\*\*\*’ 0.001 ‘\*\*’ 0.01 ‘\*’ 0.05 ‘.’ 0.1 ‘ ’ 1

**explain**

>  $$ICC = {class var\over class var + residual var} = {0.7021\over 0.7021+1.2218} = 0.36$$
>
> The intraclass correlation is 0.36. Thus, 36 percent of the variance of the popularity scores is at the group level, which is very high for social science data

#### 6.3 modeling
**code**
```R
# a cross-level interaction random slope model
m2 <- lmer(popular~sex+extrav+texp+extrav*texp+(extrav|class), data=popular2)
summary(m2)
```
**output**
> Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
> Formula: popular ~ sex + extrav + texp + extrav \* texp + (extrav | class)
>    Data: popular2
> 
> REML criterion at convergence: 4780.5
> 
> Scaled residuals: 
>      Min       1Q   Median       3Q      Max 
> -3.12872 -0.63857 -0.01129  0.67916  3.05006 
> 
> Random effects:
>  Groups   Name        Variance Std.Dev. Corr 
>  class    (Intercept) 0.478639 0.69184       
>           extrav      0.005409 0.07355  -0.64
>  Residual             0.552769 0.74348       
> Number of obs: 2000, groups:  class, 100
> 
> Fixed effects:
>               Estimate Std. Error         df t value Pr(>|t|)    
> (Intercept) -1.210e+00  2.719e-01  1.093e+02  -4.449 2.09e-05 \*\*\*
> sexgirl      1.241e+00  3.623e-02  1.941e+03  34.243  < 2e-16 \*\*\*
> extrav       8.036e-01  4.012e-02  7.207e+01  20.031  < 2e-16 \*\*\*
> texp         2.262e-01  1.681e-02  9.851e+01  13.458  < 2e-16 \*\*\*
> extrav:texp -2.473e-02  2.555e-03  7.199e+01  -9.679 1.15e-14 \*\*\*
> 
> Signif. codes:  0 ‘\*\*\*’ 0.001 ‘\*\*’ 0.01 ‘\*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> 
> Correlation of Fixed Effects:
>             (Intr) sexgrl extrav texp  
> sexgirl      0.002                     
> extrav      -0.867 -0.065              
> texp        -0.916 -0.047  0.801       
> extrav:texp  0.773  0.033 -0.901 -0.859

**explain**

> ***fixed effect***
>
> extraverted pupils are more popular (0.80 and significant)
>
> the regression coefficient for the cross-level interaction is -0.03, which is small but significant. Thus, the difference between extraverted and introverted pupils is smaller with more experienced teachers
>
> ***random effect***
>
> the coefficient of extrav (0.8) should not be interpretted without considering the variance (0.005). the slope of extrav are assumed to follow a normal distribution (mean=0.8, variance=0.005)

#### 6.4 model comparing

**code**

```R
anova(m1,m2)
```

**output**
> Data: popular2
> Models:
> m1: popular ~ 1 + (1 | class)
> m2: popular ~ sex + extrav + texp + extrav \* texp + (extrav | class)
>    Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
> m1  3 6333.5 6350.3 -3163.7   6327.5                             
> m2  9 4765.6 4816.0 -2373.8   4747.6 1579.8      6  < 2.2e-16 \*\*\*
> 
> Signif. codes:  0 ‘\*\*\*’ 0.001 ‘\*\*’ 0.01 ‘\*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#### 6.5 hypothesis testing

**code**

```R
confint(m2, method=c("profile"))
```

**output**
| | 2.5% | 97.5%|
|---|---|---|
| .sig01 | 0.43943525 | 0.91581706|
| .sig02 | -1.00000000 | 1.00000000|
| .sig03 | 0.00000000 | 0.12214512|
| .sigma | 0.71972711 | 0.76832074|
| (Intercept) | -1.74053910 | -0.67551715|
| sexgirl | 1.16942204 | 1.31181392|
| extrav | 0.72501143 | 0.88239025|
| texp | 0.19316323 | 0.25903297|
| extrav:texp | -0.02975121 | -0.01972519|

**explain**

fixed effect: all estimate are in the 95% interval, hence all significant

sig01 test variance of random intercept, since 0 not in interval, there is statistically significant variation in the intercept across class CI95[0.44, 0.92]

sig02 is the correlation between the random intercept and random slope, which is not significant

sig03 is the test of random slope, since 0 in the interval, it's marginal significant


### 7. graphing

### 8. formatting results

### Reference
Kreft, I.G.G. & de Leeuw, J. (1998). Introducing Multilevel Modeling. Thousand Oaks, CA: Sage.
Snijders, T. & Bosker, R. (1999). Multilevel Analysis: An Introduction to Basic andAdvanced Multilevel Modeling, 1st edition. Thousand Oaks, CA: Sage.
Iversen, G. (1991). Contextual Analysis. Newbury Park, CA: Sage.
Wooldridge, J. (2004). Fixed Effects and Related Estimators for Correlated Random Coefficient and Treatment Effect Panel Data Models. East Lansing: Department of Economics, Michigan State University.
Bickel, R. (2007). Multilevel Analysis for Applied Research: It’s Just Regression! New York: Guilford Press.
Kreft, I.G.G., de Leeuw, J., & Aiken, L. (1995). The Effect of Different Forms of Centering in Hierarchical Linear Models. Multivariate Behavioral Research, 30, 1–22.
Hox, J. J., Moerbeek, M., & Schoot, R. v. d. (2017). Multilevel Analysis: Techniques and Applications, Third Edition (3rd Edition ed.). New York: Routledge.
Walsh, J.E. (1947). Concerning the effect of intraclass correlation on certain significance tests. Annals of Mathematical Statistics, 18, 88–96.
Can, S., van de Schoot, R., & Hox, J. (2014). Collinear latent variables in multilevel confirmatory factor analysis: A comparison of maximum likelihood and Bayesian estimations. Educational and Psychological Measurement, 75(3), 406–427.


 <div STYLE="page-break-after: always;"></div>
# R tips
## 1. read SPSS data file
```R
# activate package (this package is installed in base distribute, no need to install)
library(foreign)
# find the file path (if you already know, skip this step)
file.choose()
# read file
dataset = read.spss("C:\\path\\file.sav", to.data.frame=TRUE)
# check input
# 1. check size
dim(dataset)
# 2. check head
head(dataset)
# 3. check structure
str(dataset)
# 4. view the dataset
view(dataset)
# 5. manually edit
fix(dataset)

```