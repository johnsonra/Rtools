# install.all.R
# install all packages I normally use (getting the most up-to-date version)
# Randall Johnson
# CCR Collaborative Bioinformatics Resource at Frederick National Laboratory
# Leidos Biomedical Research, Inc

install.all <- function(repos = "http://watson.nci.nih.gov/cran_mirror/")
{
    # install cran packages
    to.install <- c('gtools', 'haplo.stats', 'hapsim', 'hwde', 'maps', 'Matrix',
                    'plotrix', 'Rserve', 'RODBC', 'RJDBC', 'sm',
                    'statmod', 'gdata', 'Hmisc', 'lme4', 'mapdata',
                    'mapproj', 'maptools', 'SASmixed', 'TeachingDemos', 'vioplot',
                    'genetics', 'gmodels', 'gplots', 'RgoogleMaps', 'tseries',
                    'exactci', 'rmeta', 'inline', 'ncdf', 'RSQLite', 'MSIseq',
                    'Epi', 'longpower', 'gap', 'VGAM', 'wordcloud',
                    'XML', 'LEAPFrOG', 'DirichletReg', 'pwr', 'WriteXLS', 'meta',
                    'sqldf', 'doBy', 'devtools', 'pbkrtest', 'lpSolve', 'diptest',
                    'mclust', 'missForest', 'magrittr')

    install.packages(to.install, repos = repos, dependencies = TRUE)

    # check that they all loaded
    loaded <- sapply(to.install, require, character.only = TRUE)

    if(sum(!loaded) > 0)
        print(paste('The following packages were not loaded:',
              paste(to.install[!loaded], collapse = ', ')))

    # try to install some Github packages
    if(require(devtools))
    {
        devtools::install_github("johnsonra/Rtools", subdir = 'Rtools')
        devtools::install_github("johnsonra/ALDdata", subdir = 'ALDdata')
        devtools::install_github("johnsonra/ALDsuite", subdir = 'ALDsuite')
    }else{
        warning("devtools package required for installation of packages from Github.")
    }

    # install some from Bioconductor
    source("http://bioconductor.org/biocLite.R")
    biocLite()
    biocLite(c('AnnotationDbi', 'hexbin', 'GEOquery', 'simpleaffy', 'SeqVarTools',
               'qvalue', 'EBImage', 'rhdf5'))
}
