\name{cbcrisk}
\alias{cbcrisk}
\alias{CBCRisk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predicts the risk of contralateral breast cancer of a woman.}
\description{
\code{cbcrisk} is used to predict the chance of developing contralateral breast cancer for a woman whose current age (i.e., age at counseling) is within 18-89. This function gives the predictions starting at the current age (\code{start.age}) and then in increments of a specified number of years (\code{pred.year}) up to age 89.
}
\usage{
cbcrisk(race, profile, start.age, pred.year, print.output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{race}{
  race of the woman: all or unspecified race (1); non-Hispanic Black (2).
  }
  \item{profile}{
    a vector containing the values of 8 variables for all or unspecified race or a vector of 4 variables for non-Hispanic Black patient. The variables and their possible values (categories) are listed below in the Details section.
}
  \item{start.age}{
  current age (i.e., age at couselling) of the woman. It must be equal to or greater than the age at first breast cancer diagnosis.
}
  \item{pred.year}{
increments (in years) at which the predictions will be made, starting at \code{start.age} up to age 89. 
e.g; 2, 5, 10, ... . The default is 5 years. 
}
  \item{print.output}{
a logical indicating whether the output should be printed on R console. Defaults to TRUE. Setting it to FALSE suppresses the output.
}
}
\details{
The \code{cbcrisk} function needs a vector of length 8 as a patient profile for a patient with unspecified race. The elements of the patient profile vector must consist of the numerical codes of the categories (noted in parentheses below) of the following risk factors:

Age at first BC diagnosis: <30 (1); 30-39 (2); 40+ (3).   	
	
Anti-estrogen therapy:	Yes	(1); No (2); Unknown	(3).
		
1st degree family history of BC: 	Yes	(1); No (2); Unknown (3).
		
High Risk Preneoplasia: Yes	(1); No/Unknown	(2).
		
Breast density:	Extreme dense (1); Heterogeneously dense (2); Scattered	(3); Almost Entirely fat (4); Unknown	(5).
		
ER status:	Negative (1); Positive (2); Unknown	(3).
		
First BC type: Pure DCIS (1); Mixed	(2); Pure Invasive (3).
		
Age at first birth:	<30/Nulliparous	(1); 30-39 (2); 40+ (3); Unknown (4).

The \code{cbcrisk} function needs a vector of length 4 as a patient profile for a non-Hispanic Black patient. The elements of the patient profile vector must consist of the numerical codes of the categories (noted in parentheses below) of the following risk factors:

Breast density:	Scattered or Almost Entirely fat (1); Heterogeneous or Extreme dense (2);  Unknown (3).

1st degree family history of BC: 	Yes	(1); No (2); Unknown (3).

Tumor size: T0/T1/T2 (1); T3/T4 (2); TIS (3); Unknown (4)

Age at first BC diagnosis: <40 (1); 40+ (2).   	
	
}
\value{
\code{cbcrisk} returns a list of three objects--- the patient profile provided by the user, the current age of the woman provided by the user, and the contralateral breast cancer risks until age 89.

 \item{profile}{the patient profile containing information about the variables discussed above.}
 
 \item{current.age}{current age of the Black woman.}
 
 \item{risk}{the risks of contralateral breast cancer up to age 89.}
}

\references{
1. Chowdhury, M., Euhus, D., Onega, T., Biswas, S. and Choudhary, P. K. (2017) A model for individualized risk prediction of contralateral breast cancer, Breast Cancer Research and Treatment, 161, 153-160. 

2. Chowdhury, M., Euhus, D., Banu, A., Umbricht, C., Biswas, S. and Choudhary, P. K. (2018) Validation of a Personalized Risk Prediction Model for Contralateral Breast Cancer, Breast Cancer Research and Treatment, 170, 415-423.

3. Chowdhury, M., Euhus, D., O'Donnell, M., Onega, T., Choudhary, P. K. and Biswas, S. (2018) Dose-Dependent Effect of Mammographic Breast Density on the Risk of Contralateral Breast Cancer, Breast Cancer Research and Treatment, 170, 143-14.

4. Sajal, I. H., Chowdhury, M., Wang, T., Euhus, D., Choudhary, P. K. and Biswas, S. (2022) CBCRisk-Black: a personalized contralateral breast cancer risk prediction model for black women, Breast Cancer Research and Treatment, 194, 179-186.   
}
\author{
Marzana Chowdhury, 
Ibrahim Hossain Sajal
}

\examples{
## The profile vector
patient_profile_1= c(3, 2, 1, 2, 3, 2, 1, 2)
patient_profile_2= c(2, 2, 1, 2)

## install packages "survival" (if not already installed)
## install the "cbcrisk" package

library(cbcrisk)

## risks of a patient of any or unspecified race at every 5 years (by default) starting at age 55
cbcrisk(1, patient_profile_1, 55)  

## risks of a black patient at every 5 years (by default) starting at age 55
cbcrisk(2, patient_profile_2, 55)  

## risks of a black patient at every 10 years starting at age 55
cbcrisk(2, patient_profile_2, 55, 10) 
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
