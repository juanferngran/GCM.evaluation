# Hypothesis Test for the Difference of Two Population Proportions ---------

# wts.cccma <- clusterGrid(cccma, type = "lamb")
# clusters.era.interim <- clusterGrid(era.interim, type = "lamb")
# cccma.CMIP5.transitionProb <- transitionProb(wts.cccma)
# eraIn.transitionProb <- transitionProb(clusters.era.interim)

cccma.cmip5.test <- prop.test(cbind(as.vector(cccma.CMIP5.transitionProb), as.vector(eraIn.transitionProb)))

str(cccma.cmip5.test)

# List of 9
# $ statistic  : Named num 15.6
# ..- attr(*, "names")= chr "X-squared"
# $ parameter  : Named num 298
# ..- attr(*, "names")= chr "df"
# $ p.value    : num 1
# $ estimate   : Named num [1:299] 0.555 0.799 0.833 0.26 0.449 ...
# ..- attr(*, "names")= chr [1:299] "prop 1" "prop 2" "prop 3" "prop 4" ...
# $ null.value : NULL
# $ conf.int   : NULL
# $ alternative: chr "two.sided"
# $ method     : chr "299-sample test for equality of proportions without continuity correction"
# $ data.name  : chr "cbind(as.vector(cccma.CMIP5.transitionProb), as.vector(eraIn.transitionProb))"
# - attr(*, "class")= chr "htest"