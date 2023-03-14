\name{dataIEA}
\alias{dataIEA}
\docType{data}
\title{
Data for understanding of good citizenship behaviour
}
\description{
Data set from the International Civic and Citizenship Education Study 2016 (Schulz et al., 2018). As part of a comprehensive evaluation of education systems, the IEA conducted surveys in 1999, 2009 and 2016 in school classes of 14-year olds to investigate civic education with the same scientific rigor as the evaluation of more traditional educational skills of language and mathematics. The present study focuses on the third wave of the survey that was conducted in 2016.

Questions regarding citizenship norms in all three waves asked respondents to explain their understanding of what a good adult citizen is or does. The survey then lists a variety of activities for respondents to rate in terms of how important these activities are in order to be considered a good adult citizen. The twelve items range from obeying the law and voting in elections, to protecting the environment and defending human rights.

Covariates included are customary determinants of citizenship norms from the literature at the individual-level of socio-economic measures and country-level measure of gross domestic product (GDP) per capita.
}
\usage{data("dataIEA")}
\format{
  A data frame with 90221 observations on the following 28 variables.
  \describe{
    \item{\code{ICCS_year}}{Year of survey}
    \item{\code{COUNTRY}}{Country}
    \item{\code{IDSTUD}}{Study ID}
    \item{\code{TOTWGTS}}{Study weight}
    \item{\code{obey}}{Always obeying the law}
    \item{\code{rights}}{Taking part in activities promoting human rights}
    \item{\code{local}}{Participating in activities to benefit people in the local community}
    \item{\code{work}}{Working hard}
    \item{\code{envir}}{Taking part in activities to protect the environment}
    \item{\code{vote}}{Voting in every national election}
    \item{\code{history}}{Learning about the country's history}
    \item{\code{respect}}{Showing respect for government representatives}
    \item{\code{news}}{Following political issues in the newspaper, on the radio, on TV, or on the Internet}
    \item{\code{protest}}{Participating in peaceful protests against laws believed to be unjust}
    \item{\code{discuss}}{Engaging in political discussions}
    \item{\code{party}}{Joining a political party}
    \item{\code{female}}{Female}
    \item{\code{books}}{Number of books at home}
    \item{\code{edexp}}{Educational expectations}
    \item{\code{ed_mom}}{Mother education}
    \item{\code{ed_dad}}{Father education}
    \item{\code{nonnat_born}}{Non-native born}
    \item{\code{immigrantfam}}{Immigrant family}
    \item{\code{nonnat_lang}}{Non-native language level}
    \item{\code{gdp_constant}}{GDP}
    \item{\code{log_gdp_constant}}{Log GDP}
    \item{\code{gdp_currentusd}}{GDP in USD}
    \item{\code{log_gdp_currentusd}}{Log GDP in USD}
  }
}
\references{
Schulz, W., Ainley, J., Fraillon, J., Losito, B., Agrusti, G., & Friedman, T. (2018). Becoming
citizens in a changing world: IEA International Civic and Citizenship Education Study 2016
international report. Springer.
}
\keyword{datasets}