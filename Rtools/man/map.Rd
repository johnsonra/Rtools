\name{map}

\alias{map}

\docType{data}

\title{2004 MALD map}

\description{Mapping by Admixture Linkage Disequilibrium (MALD) map from
the Smith et al paper in 2004.}

\usage{data(map)}

\format{
  A data frame with 3011 observations on the following 43 variables.
  \describe{
    \item{snp}{SNP ID (dbSNP ID provided if available)}
    \item{submitter}{Submitter}
    \item{chr}{Chromosome (23 == X)}
    \item{gen.pos}{Interpoated deCODE genetic position}
    \item{phys.pos}{Pysical position in HG16}
    \item{direction}{Direction of strand in HG 16}
    \item{genotyping.site}{Site of revalidation genotyping}
    \item{fasta}{FASTA sequence, length = 201 bp (variant whose
      frequency is quoted is given to the right of the `/')}
    \item{snp.code}{SNP (same direction as FASTA; variant whose
      frequency is quoted is given to the right of the `/')} 
    \item{notes}{Notes about possible mapping problems}
    \item{rank}{rank in map (missing if not in 2154 well-spaced set)}
    \item{f.Euro.WAfri}{Frequency; 21\% European, 79\% West African}
    \item{f.Euro.AmIn}{Frequency; 50\% European, 50\% Amerindian}
    \item{f.WAfri.AmIn}{Frequcncy; 50\% West African, 50\% Amerindian}
    \item{f.Euro.EAsian}{Frequency; 50\% European, 50\% East Asian}
    \item{f.WAfri}{West African frequency (no Botswana)}
    \item{n.WAfri}{West African count (no Botswana)}
    \item{f.EuroAm}{European American frequency}
    \item{n.EuroAm}{European American count}
    \item{f.AmIn}{Amerindian frequency}
    \item{n.AmIn}{Amerindian count}
    \item{f.EAsian}{East Asian frequency}
    \item{n.EAsian}{East Asian count}
    \item{f.nigerian}{Beni Nigerian frequency}
    \item{n.nigerian}{Beni Nigerian count}
    \item{f.ghana}{Ghana frequency}
    \item{n.ghana}{Ghana count}
    \item{f.cameroon}{Cameroon frequency}
    \item{n.cameroon}{Cameroon count}
    \item{f.senegal}{Senegal frequency}
    \item{n.senegal}{Senegal count}
    \item{f.botswana}{Botswana frequency}
    \item{n.botswana}{Botswana count}
    \item{f.utah}{European Americans from Utah (CEPH) frequency}
    \item{n.utah}{European Americans from Utah (CEPH) count}
    \item{f.baltimore}{European Americans from Baltimore frequency}
    \item{b.baltimore}{European Americans from Baltimore count}
    \item{f.chicago}{European Americans from Chicago frequency}
    \item{n.chicago}{European Americans from Chicago count}
    \item{f.mexico}{Zapotec Amerindians from Mexico frequency}
    \item{n.mexico}{Zapotec Amerindians from Mexico count}
    \item{f.cantonese}{Cantonese Chinese frequency}
    \item{n.cantonese}{Cantonese Chinese count}
  }
}

\source{Smith et al; A High-Density Admixture Map for Disease Gene
  Discovery in African Americans; American Journal of Human Genetics,
  2004, 74, 1001-13}

\examples{data(map)}

\keyword{datasets}
