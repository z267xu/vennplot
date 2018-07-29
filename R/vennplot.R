#' Draw Venn and Euler diagram in 2D or 3D
#'
#' @param disjoint.combinations Named numeric vector or data.frame where each column should be factor.  See Details.
#' @param vars Extract specific variables of data.frame as \code{disjoint.combinations}. If \code{vars = NULL}, all the information of data.frame will be extracted.
#' @param Delta The length of step for method "lineSearch" or the initial interval of test points for method "NelderMead".
#' @param ThreeD Draw Venn diagram in 3D. See Examples.
#' @param lambda It can be \code{NULL} or a numeric vector. If \code{lambda = NULL}, the loss function optimize lambda, else, based on the given lambda, loss function will calculate stress respectively then return the minimum one and corresponding lambdas.
#' @param stressWay If data set can be separated into a few groups, there will be two ways to express stress: one is to sum up all the stress (named "sum"; default), the other is to use total TSS divide by total RSS (named "combine").
#' @param delta Closeness between groups.
#' @param weight The weight of \code{disjoint.combinations}. It should have the same length with \code{disjoint.combinations}.
#' @param expand If some balls should not intersect and the code fails to detect it. It is possible to be fixed manually but sacrificing stress.
#' @param twoWayGenerate Boolean factor, if false, any missing intersections are set as zero.
#' @param scaleSearch Provide multiple methods to optimize scale lambda. The default method is "NelderMead". See Details.
#' @param twoWaySearch If two way intersections are missing, multiple methods are available to generate two way intersections. The default method is "lineSearch". See Details.
#' @param scaleSeachTolerance A list with tolerance value and boolean factor " proportional". The loop of NelderMead and lineSearch in scaleSearch will end when the difference or  proportional difference matches the tolerance value.
#' @param distanceTolerance A list with tolerance value and boolean factor " proportional". The Newton method of finding distance will end when the difference or  proportional difference matches the tolerance value.
#' @param lossTolerance A list with ToleranceofLoss, maximumStep, ALPHA, ToleranceofStepsize and boolean factor "proportional". If ALPHA is null, the step size will be searched through Newton method and it will stop when step reaches the maximum step or the difference matches ToleranceofStepsize; else step size will be fixed with ALPHA . The loss will end when the difference or proportional difference or the total loss value matches the "ToleranceofLoss".
#' @param stressBound The loop of method NelderMead will stop when stress is beyond the stressBound.
#' @param maximumStep The maximum searching step for method NelderMead and Newton method of calculating distance.
#' @param planeSize The plane size of calculating disjoint intersections numerically.
#' @param lower The lower bound of the interval to be searched for the "goldenSectionSearch" and "L-BFGS-B". See Details.
#' @param upper The upper bound of the interval to be searched for the "goldenSectionSearch" and "Brent". See Details.
#' @param control A list of control parameters. See Details
#' @param hessian Logical. A numerically differentiated Hessian matrix be returned or not. See Details.
#' @param mar Plot margins.
#' @param cols Color of balls. If \code{NULL}, rainbow color will be set.
#' @param alpha Color darkness.
#' @param smooth For 3D plot, if true, the balls will be much more smoother. However, based on the high resolution, if the number of balls is too much, when rotating, the new window stumbles.
#' @param ... Any further graphical parameters to be passed to the \code{plot} function.
#'
#' @details
#' 1. One way sets must be given in \code{disjoint.combination}. e.g.\code{disjoint.combination = c( B=2, AB=0.5)} is not allowed.  \code{disjoint.combination = c(A = 0, B=2, AB=0.5)} works.
#' 2. Except "NelderMead" and "lineSearch", "goldenSectionSearch" in \code{scaleSearch} and \code{twoWaySearch} is based on \code{\link{optimize}} and the rest methods are based on \code{\link{optim}}.
#' 3. \code{lower}, \code{upper}, \code{control} and \code{hessian} share the same parameters with \code{\link{optim}}, and \code{lower}, \code{upper} can also be used in \code{\link{optimize}}
#'
#' @return An object of the class \code{vennplot} with following components:
#' \describe{
#'   \item{xy}{centres of the balls (columns are (\code{x}, \code{y}) or (\code{x}, \code{y}, \code{z}) coordinates).}
#'   \item{radius}{radii of the balls.}
#'   \item{loss}{total loss of \code{vennplot}.}
#'   \item{stress}{stress value for solution.}
#' }
#' @examples
#' # 3D Venn plot with arbitray sets
#' disjoint.combinations <- c(A=80, B=50,C=100, D = 100,E = 100,
#'                           "A&C"=30, "A&D"= 30,"B&E" = 30, "A&E" = 40, h = 40, "B&h" = 10)
#' ve <- vennplot(disjoint.combinations, ThreeD = TRUE)
#'
#' # data frame
#' vennplot(disjoint.combinations = sharks, vars = c("Au","USA","Fa","Sex"),
#'          scaleSearch = "lineSearch", expand = 1.1)
#'
#' @export
vennplot <- function(disjoint.combinations = NULL, vars = NULL, Delta = 0.1,
                     ThreeD = FALSE, lambda = NULL, stressWay = c("sum","combine"),
                     delta = 0.01, weight = NULL, expand = NULL, twoWayGenerate = FALSE,
                     scaleSearch = c("NelderMead", "lineSearch", "goldenSectionSearch",
                                     "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                     twoWaySearch = c("lineSearch", "NelderMead", "goldenSectionSearch",
                                      "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                     scaleSeachTolerance = list(value = 1e-5,  proportional = FALSE),
                     distanceTolerance = list(value = 1e-5,  proportional = FALSE),
                     lossTolerance = list(ToleranceofLoss = 1e-10, maximumStep = 10, ALPHA = 1e-2,
                                          ToleranceofStepsize = 1e-5,  proportional = FALSE),
                     stressBound = 1e-3, maximumStep = 50, planeSize = 50,
                     lower=  -Inf, upper = Inf, control = list(), hessian = FALSE,
                     mar = rep(1,4), cols = NULL, alpha = 0.3, smooth = FALSE, ...){

  if(is.null(disjoint.combinations)) stop("combinations should not be empty")

  if (is.data.frame(disjoint.combinations)){
    disjoint.combinations <- extractCombinations(disjoint.combinations, vars = vars)
  }
  combinations <- disjoint2combinations(disjoint.combinations)
  # reorder disjoint.combinations with ways
  reOrder <- reorderDisjointCombinations(combinations, disjoint.combinations)
  disjoint.combinations <- reOrder$newDisjointCombinations
  # scale proportionately
  combProp <- combinations/sum(disjoint.combinations)
  disProp <- disjoint.combinations/sum(disjoint.combinations)
  # Check weight vector
  if(is.null(weight)){
    weight <- rep(1,length(disjoint.combinations))
  } else {
    if(length(weight)!=length(disjoint.combinations)){
      stop("weight must be NULL or the same length as the number of disjoint combinations")
    }
  }
  weight <- weight[reOrder$newOrder]
  disjointSetNames <- names(disProp)
  names(weight) <- disjointSetNames

  # Get the number of intersections
  numWays <- str_count(disjointSetNames, pattern = "&") + 1

  if(max(numWays)==2){
    nonEmptyTwoWays <- which(combProp!=0)
    weight <- weight[nonEmptyTwoWays]
    disjointSetNames <- disjointSetNames[nonEmptyTwoWays]
    numWays <- numWays[nonEmptyTwoWays]
    combProp <- combProp[nonEmptyTwoWays]
    disProp <- disProp[nonEmptyTwoWays]
  }

  oneWays <- which(numWays == 1)
  # Check colours
  if (is.null(cols)) {
    cols <- rainbow(length(oneWays), alpha = alpha)
  }
  oneWaySetName <- disjointSetNames[oneWays]
  oneWaySet <- combProp[oneWays]
  m <- length(oneWays)

  # and for those that are larger than one way set
  if (length(disjointSetNames) == length(oneWays)) {
    largerThanOneWaySetName <- NULL
  } else {
    largerThanOneWaySetName <- disjointSetNames[-oneWays]
    # All names appearing in two and higher order ways must appear
    # as individual input sets too
    if (!all(unique(unlist(str_split(largerThanOneWaySetName,"&"))) %in% oneWaySetName)){
      stop("Some intersection sets contain sets which do not appear individually.")
    }
  }
  # Detect disjoint groups of sets and return groups in a list
  groups <- groupDetection(largerThanOneWaySetName = largerThanOneWaySetName,
                           oneWaySetName= oneWaySetName)
  # Order groups from largest to smallest
  groupOrder <- order(sapply(groups, length), decreasing = T)
  groups <- groups[groupOrder]
  ngroups <- length(groups)

  # get the combinations and proportions separated by groups

  combPropGroup <- list()
  disPropGroup <- list()
  for(i in 1:ngroups){
    groupMember <- groups[[i]]
    groupWay <- oneWaySetName[groupMember]
    for(j in 1:length(disjointSetNames)){
      if(any(groupWay %in% str_split(disjointSetNames[j],"&")[[1]])){
        groupWay <- c(groupWay, disjointSetNames[j])}
    }
    groupWay <- unique(groupWay)
    combPropGroup[[i]] <- combProp[which(disjointSetNames %in% groupWay)]
    disPropGroup[[i]] <- disProp[which(disjointSetNames %in% groupWay)]
  }

  # Calculate radius
  radius <- if (ThreeD){
    (3*oneWaySet/(4*pi))^(1/3)
  } else {
    sqrt(oneWaySet/pi)
  }
  # Get radii within each group
  radiusGroup <- lapply(groups,
                        function(grp){radius[grp]})

  weightGroup <- lapply(disPropGroup,
                        function(prop){
                          weight[which(names(weight) %in% names(prop)==TRUE)]})
  # lossFunction
  groupCentre <- list()
  if(length(lower) != length(upper)){
    stop("the length of lower vector must be equal to the length of upper vector")
  }
  if(length(lower) == 1 && length(upper) == 1){
    lower = rep(lower,2)
    upper = rep(upper,2)
  }
  if(length(hessian) == 1){hessian = rep(hessian, 2)}

  scaleSearch <- match.arg(scaleSearch)
  twoWaySearch <- match.arg(twoWaySearch)
  stressWay <- match.arg(stressWay)
  stressGroup <- rep(0, ngroups)
  RSSGroup <- rep(0, ngroups)
  TSSGroup <- rep(0, ngroups)
  loss <- 0
  for (gn in 1:ngroups){
    lossFunctionOutput <- lossFunction(gn = gn, combProp = combPropGroup[[gn]], Delta = Delta,
                                       radius = radiusGroup[[gn]], disProp = disPropGroup[[gn]], lambda = lambda,
                                       weight = weightGroup[[gn]], ThreeD = ThreeD,
                                       method = c(scaleSearch, twoWaySearch), twoWayGenerate = twoWayGenerate,
                                       lower=  lower, upper = upper, control = control, hessian = hessian,
                                       expand = expand, scaleSeachTolerance = scaleSeachTolerance,
                                       distanceTolerance = distanceTolerance, lossTolerance = lossTolerance,
                                       stressBound = stressBound, maximumStep = maximumStep, planeSize = planeSize)
    groupCentre[[gn]] <- lossFunctionOutput$centre
    stressGroup[gn] <- lossFunctionOutput$stress
    RSSGroup[gn] <- lossFunctionOutput$RSS
    TSSGroup[gn] <- lossFunctionOutput$TSS
    loss <- loss + lossFunctionOutput$loss
  }
  if(stressWay == "sum"){
    stress <- sum(stressGroup)
  } else {
    #combine stress
    stress <- sum(RSSGroup)/sum(TSSGroup)
  }
  #combine all groups
  if(ngroups > 1){
    combinedGroups <- combineGroups(groupCentre = groupCentre, radiusGroup = radiusGroup,delta = delta)
    xy <- combinedGroups$xy
    radius <- combinedGroups$radius
    oneWaySetName <- combinedGroups$namexy
  } else {
    xy <- groupCentre[[1]]
    radius <- radiusGroup[[1]]
    oneWaySetName <- names(radius)
  }
  if(ThreeD){
    sphere(xy = xy, radius = radius, cols = cols, alpha = alpha, oneWaySetName = oneWaySetName, smooth = smooth)
  } else {
    plot.new()
    plotCircle(xy, radius, oneWaySetName, cols, mar = mar)
  }
  list(xy = xy, radius = radius, loss = loss, stress = stress)
}
#--- helper functions -----------------------------------------------------

#library(rgl)
#library(stringr)
#library(Rcpp)
#sourceCpp(file = "your path/vennCpp2.cpp")

#lossFunction, the main function of vennplot; return centres of each group
lossFunction <- function(gn, combProp, radius, Delta, lambda,disProp,
                         weight, ThreeD, method, twoWayGenerate,
                         lower, upper, control, hessian, expand, scaleSeachTolerance,
                         distanceTolerance, lossTolerance, stressBound, maximumStep, planeSize){
  disjointSetNames <- names(combProp)
  numWays <- str_count(disjointSetNames,pattern = "&")+1
  #one way set
  oneWaySetName <- disjointSetNames[which(numWays==1)]
  m <- length(oneWaySetName)
  if(m == 1){
    if(ThreeD){centre <- matrix(c(0,0,0),nrow = 1)}else{centre <- matrix(c(0,0),nrow=1)}
    rownames(centre) <- oneWaySetName
    loss <- 0
    oneWayStress <- calculateStress(centre, radius, disProp, weight, ThreeD, planeSize,twoWayGenerate)
    RSS <- oneWayStress$RSS
    TSS <- oneWayStress$TSS
    stress <- oneWayStress$stress
  }else{
    complete <- TRUE
    if(max(numWays)>=3){
      firstGenerate <- twoWayGeneration(combProp = combProp, mu = 1)
      #To determine which example we use: Example 1 or Example 2
      if(firstGenerate$resam != 0){complete <- FALSE}
      newTworWaySet <- firstGenerate$newTworWaySet
    }else{newTworWaySet <- combProp[which(numWays == 2)]}
    EDandInitialLocation <- EuclideanDistance(newTworWaySet = newTworWaySet, oneWaySetName = oneWaySetName,
                                              radius = radius, ThreeD = ThreeD, initial = TRUE,
                                              expand = expand, distanceTolerance = distanceTolerance,
                                              maximumStep = maximumStep)
    ED <- EDandInitialLocation$ED
    xy <- EDandInitialLocation$xy
    if(twoWayGenerate == TRUE &&  complete == FALSE){
      out <- solveWithMu(ED = ED, xy = xy, combProp = combProp, ThreeD = ThreeD,
                         Delta = Delta,  radius = radius, disProp = disProp, weight = weight, method = method,
                         firstGenerate = firstGenerate,
                         lower=  lower, upper = upper, control = control, hessian = hessian,
                         expand = expand, scaleSeachTolerance = scaleSeachTolerance,
                         distanceTolerance = distanceTolerance, lossTolerance = lossTolerance,
                         stressBound = stressBound, maximumStep = maximumStep, planeSize = planeSize,
                         twoWayGenerate = twoWayGenerate)
      ED <- out$ED
      xy <- out$xy
    }
    result <- solveWithLambda(ED = ED, xy = xy, ThreeD = ThreeD, lambda = lambda, Delta = Delta,
                              radius = radius, disProp = disProp, weight = weight, method = method,
                              lower=  lower, upper = upper, control = control, hessian = hessian,
                              scaleSeachTolerance = scaleSeachTolerance, lossTolerance = lossTolerance,
                              stressBound = stressBound, maximumStep = maximumStep, planeSize = planeSize,
                              twoWayGenerate = twoWayGenerate)
    centre <- result$xy
    rownames(centre) <- oneWaySetName
    stress <- result$stress
    loss <- result$loss
    RSS <- result$RSS
    TSS <- result$TSS
  }
  list(centre = centre, loss = loss, stress = stress, RSS = RSS, TSS = TSS)
}
#optimize lambda to get centre
solveWithLambda <- function(ED , xy , ThreeD , lambda, Delta, radius ,
                            disProp, weight, method, lower, upper, control,
                            hessian, scaleSeachTolerance, lossTolerance,
                            stressBound, maximumStep, planeSize, twoWayGenerate) {
  stress_0_Calculation <- calculateStress(xy = xy,radius = radius,
                                          disProp = disProp,
                                          weight = weight,ThreeD = ThreeD, planeSize = planeSize,
                                          twoWayGenerate = twoWayGenerate)
  stress_0 <- stress_0_Calculation$stress
  lambda_1_Stress <- findOptimalStress(lambda = 1, xy = xy,weight = weight,
                                       radius = radius, disProp = disProp, ED = ED,
                                       ThreeD = ThreeD, lossTolerance = lossTolerance, planeSize = planeSize,
                                       twoWayGenerate = twoWayGenerate)
  stress_1 <- lambda_1_Stress$stress
  if(ThreeD){offset <- 3}else{offset <- 2}
  if(min(stress_0,stress_1) == stress_0 && dim(xy)[1]<=offset){
    stress <- stress_0
    loss <- 0
    RSS <- stress_0_Calculation$RSS
    TSS <- stress_0_Calculation$TSS
  }else{
    if(is.null(lambda)){
      if(method[1] == "NelderMead"){
        Centre <- list();
        lambda_2_Stress <- findOptimalStress(lambda = 1+Delta, xy = xy, weight = weight,
                                             radius = radius, disProp = disProp,
                                             ED = ED, ThreeD = ThreeD, lossTolerance = lossTolerance,
                                             planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        stress_2 <- lambda_2_Stress$stress
        lambda_3_Stress <- findOptimalStress(lambda = 1-Delta, xy = xy, weight = weight,
                                             radius = radius, disProp = disProp, ED = ED,
                                             ThreeD = ThreeD, lossTolerance = lossTolerance,
                                             planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        stress_3 <- lambda_3_Stress$stress
        Centre[[1]] <- lambda_1_Stress$xy
        Centre[[2]] <- lambda_2_Stress$xy
        Centre[[3]] <- lambda_3_Stress$xy
        STRESS <- c(stress_1,stress_2,stress_3)
        xy <- Centre[[which(STRESS == min(STRESS))[1]]]
        loss <- c(lambda_1_Stress$loss, lambda_2_Stress$loss, lambda_3_Stress$loss)[which(STRESS == min(STRESS))[1]]
        LAMBDA <- c(1,1+Delta,1-Delta)[order(STRESS)]
        STRESS <- sort(STRESS)
        count <- 0
        while (count< maximumStep && min(STRESS) > stressBound){
          count <- count+1
          lambda_0 <- mean(LAMBDA[1:2])
          lambda_R <- lambda_0 + (lambda_0 - LAMBDA[3])
          lambda_R_stress <- findOptimalStress(lambda = lambda_R, xy = xy, weight = weight,
                                               radius = radius, disProp = disProp, ED = ED,
                                               ThreeD = ThreeD, lossTolerance = lossTolerance,
                                               planeSize = planeSize, twoWayGenerate = twoWayGenerate)
          stress_R <- lambda_R_stress$stress
          #Reflection:
          if(STRESS[1] <= stress_R &&  stress_R< STRESS[2]){
            STRESS <- c(STRESS[1] , stress_R , STRESS[2])
            LAMBDA <- c(LAMBDA[1], lambda_R, LAMBDA[2])
          } else if (stress_R < STRESS[1]){
            #Expansion
            xy <- lambda_R_stress$xy
            loss <- lambda_R_stress$loss
            lambda_E <- lambda_0 + 2*(lambda_R - lambda_0)
            lambda_E_stress <- findOptimalStress(lambda = lambda_E, xy = xy, weight = weight,
                                                 radius = radius, disProp = disProp, ED = ED,
                                                 ThreeD = ThreeD, lossTolerance = lossTolerance,
                                                 planeSize = planeSize, twoWayGenerate = twoWayGenerate)
            stress_E = lambda_E_stress$stress
            if(stress_E<stress_R){
              STRESS <- c(stress_E , stress_R , STRESS[1])
              LAMBDA <- c(lambda_E, lambda_R, LAMBDA[1])
              xy <- lambda_E_stress$xy
              loss <- lambda_E_stress$loss
            }else if(stress_R<=stress_E && stress_E <STRESS[1]){
              STRESS <- c(stress_R, stress_E, STRESS[1])
              LAMBDA <- c(lambda_R, lambda_E, LAMBDA[1])
            }
            else if(STRESS[1]<=stress_E && stress_E <STRESS[2]){
              STRESS <- c(stress_R, STRESS[1],stress_E)
              LAMBDA <- c(lambda_R, LAMBDA[1], lambda_E)
            }else{
              STRESS <- c(stress_R, STRESS[1],STRESS[2])
              LAMBDA <- c(lambda_R, LAMBDA[1],LAMBDA[2])
            }
          } else if(stress_R >= STRESS[2]) {
            #Contraction
            lambda_C <-  lambda_0 + 0.5*(LAMBDA[3] - lambda_0)
            lambda_C_stress <- findOptimalStress(lambda = lambda_C, xy = xy, weight = weight,
                                                 radius = radius, disProp = disProp, ED = ED,
                                                 ThreeD = ThreeD, lossTolerance = lossTolerance,
                                                 planeSize = planeSize, twoWayGenerate = twoWayGenerate)
            stress_C <- lambda_C_stress$stress
            if (stress_C < STRESS[3]){
              STRESS <- c(STRESS[1],STRESS[2],stress_C)
              LAMBDA <- c(LAMBDA[1],LAMBDA[2],lambda_C)[order(STRESS)]
              STRESS <- sort(STRESS)
              if(min(STRESS) == stress_C){
                xy <- lambda_C_stress$xy
                loss <- lambda_C_stress$loss
              }
            } else {
              #shrink
              LAMBDA[2] <- (LAMBDA[1]+LAMBDA[2])/2
              STRESS[2] <- findOptimalStress(lambda = LAMBDA[2], xy = xy, weight = weight,
                                             radius = radius, disProp = disProp, ED = ED,
                                             ThreeD = ThreeD, lossTolerance = lossTolerance,
                                             planeSize = planeSize,twoWayGenerate = twoWayGenerate)$stress
              LAMBDA[3] <- (LAMBDA[1]+LAMBDA[3])/2
              STRESS[3] <- findOptimalStress(lambda = LAMBDA[3], xy = xy, weight = weight,
                                             radius = radius, disProp = disProp, ED = ED,
                                             ThreeD = ThreeD, lossTolerance = lossTolerance,
                                             planeSize = planeSize,twoWayGenerate = twoWayGenerate)$stress
            }
          }
          if(allConnectedCpp(xy = xy, radius = radius, ThreeD = ThreeD) == FALSE){break}
          if(BoolScaleNMCpp( proportional = scaleSeachTolerance$ proportional,
                             value = scaleSeachTolerance$value,
                             LAMBDA = LAMBDA, STRESS = STRESS) == FALSE){break}
        }
        if(allConnectedCpp(xy = xy, radius = radius, ThreeD = ThreeD) == FALSE) {
          lambdaStress <-  findOptimalStress(lambda = LAMBDA[2], xy = xy, weight = weight,
                                             radius = radius, disProp = disProp, ED = ED,
                                             ThreeD = ThreeD, lossTolerance = lossTolerance,
                                             planeSize = planeSize, twoWayGenerate = twoWayGenerate)
          xy <- lambdaStress$xy
          if(allConnectedCpp(xy = xy, radius = radius, ThreeD = ThreeD) == FALSE) {
            lambdaStress <-  findOptimalStress(lambda = LAMBDA[3], xy = xy, weight = weight,
                                               radius = radius, disProp = disProp, ED = ED,
                                               ThreeD = ThreeD, lossTolerance = lossTolerance,
                                               planeSize = planeSize, twoWayGenerate = twoWayGenerate)
            xy <- lambdaStress$xy
          }
          stress <- lambdaStress$stress
          loss <- lambdaStress$loss
          RSS <- lambdaStress$RSS
          TSS <- lambdaStress$TSS
        } else {
          lambdaStress <-  findOptimalStress(lambda = LAMBDA[1], xy = xy, weight = weight,
                                             radius = radius, disProp = disProp, ED = ED,
                                             ThreeD = ThreeD, lossTolerance = lossTolerance,
                                             planeSize = planeSize, twoWayGenerate = twoWayGenerate)
          xy <- lambdaStress$xy
          stress <- lambdaStress$stress
          loss <- lambdaStress$loss
          RSS <- lambdaStress$RSS
          TSS <- lambdaStress$TSS
        }
      } else if (method[1] == "lineSearch"){
        lambda_2_Stress <- findOptimalStress(lambda = 1+Delta, xy = xy, weight = weight,
                                             radius = radius, disProp = disProp, ED = ED,
                                             ThreeD = ThreeD, lossTolerance = lossTolerance,
                                             planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        stress_2 <- lambda_2_Stress$stress
        lambda_3_Stress <- findOptimalStress(lambda = 1-Delta, xy = xy, weight = weight,
                                             radius = radius, disProp = disProp, ED = ED,
                                             ThreeD = ThreeD, lossTolerance = lossTolerance,
                                             planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        stress_3 <- lambda_3_Stress$stress
        if(min(stress_1,stress_2,stress_3) == stress_1){
          xy <- lambda_1_Stress$xy
          loss <- lambda_1_Stress$loss
          stress <- stress_1
          RSS <- lambda_1_Stress$RSS
          TSS <- lambda_1_Stress$TSS
        } else {
          Center <- list()
          Loss <- c()
          if(min(stress_1,stress_2,stress_3) == stress_2){
            stress_n <- stress_2
            stress <- stress_2
            Center[[1]] <- lambda_2_Stress$xy
            Loss[1] <- lambda_2_Stress$loss
            lambda <- 1+Delta
            shrinkage <- TRUE
          } else {
            stress_n <- stress_3
            stress <- stress_3
            Center[[1]] <- lambda_3_Stress$xy
            Loss[1] <- lambda_3_Stress$loss
            lambda <- 1-Delta
            shrinkage <- FALSE
          }
          RSSVec <- c()
          TSSVec <- c()
          n <- 1
          while(stress_n <= stress){
            stress <- stress_n
            if(shrinkage){
              lambda <- lambda + Delta
            } else {
              lambda <- lambda - Delta
            }
            n <- n+1
            lambdaStress <- findOptimalStress(lambda = lambda, xy = Center[[n-1]], weight = weight,
                                              radius = radius, disProp = disProp, ED = ED,
                                              ThreeD = ThreeD, lossTolerance = lossTolerance,
                                              planeSize = planeSize, twoWayGenerate = twoWayGenerate)
            Center[[n]] <- lambdaStress$xy
            Loss[n] <- lambdaStress$loss
            stress_n <- lambdaStress$stress
            RSSVec[n] <- lambdaStress$RSS
            TSSVec[n] <- lambdaStress$TSS
            if(BoolScaleLCpp( proportional = scaleSeachTolerance$ proportional,
                              value = scaleSeachTolerance$value, stress_n = stress_n, stress = stress)){break}
          }
          xy <- Center[[n-1]]
          loss <- Loss[n-1]
          RSS <- RSSVec[n-1]
          TSS <- TSSVec[n-1]
        }
      }else if(method[1] == "goldenSectionSearch"){

        Optimization <- optimize(f = findOptimalStressLambda, lower=  lower[1], upper = upper[1],
                                 xy = lambda_1_Stress$xy, weight = weight,
                                 radius = radius, disProp = disProp, ED = ED,
                                 ThreeD = ThreeD, lossTolerance = lossTolerance,
                                 planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        if(is.null(Optimization$minimum)||is.infinite(Optimization$minimum)){
          stop("This method does not converge")
        }
        lambdaStress <- findOptimalStress(lambda = Optimization$minimum, xy = lambda_1_Stress$xy,weight = weight,
                                          radius = radius, disProp = disProp, ED = ED,
                                          ThreeD = ThreeD, lossTolerance = lossTolerance,
                                          planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        xy <- lambdaStress$xy
        loss <- lambdaStress$loss
        stress <- Optimization$objective
        RSS <- lambdaStress$RSS
        TSS <- lambdaStress$TSS
      } else {
        ## method[1] must be an method appropriate for optim(...)
        Optimization <- optim(par = 1, fn = findOptimalStressLambda,method = method[1],
                              lower = lower[1], upper = upper[1],
                              control = control, hessian = hessian[1],
                              xy = lambda_1_Stress$xy, weight = weight,
                              radius = radius, disProp = disProp, ED = ED,
                              ThreeD = ThreeD, lossTolerance = lossTolerance,
                              planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        if(is.null(Optimization$par)||is.infinite(Optimization$par)){
          stop("This method does not converge")
        }
        lambdaStress <- findOptimalStress(lambda = Optimization$par, xy = lambda_1_Stress$xy,weight = weight,
                                          radius = radius, disProp = disProp, ED = ED,
                                          ThreeD = ThreeD, lossTolerance = lossTolerance,
                                          planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        xy <- lambdaStress$xy
        loss <- lambdaStress$loss
        stress <- Optimization$value
        RSS <- lambdaStress$RSS
        TSS <- lambdaStress$TSS
      }
    } else {
      Centre <- list()
      RSSVec <- c()
      TSSVec <- c()
      STRESS <- rep(0,length(lambda))
      Loss <- rep(1,length(lambda))
      for(i in 1:length(lambda)){
        out <- findOptimalStress(lambda = lambda[i], xy = xy, weight = weight,
                                 radius = radius, disProp = disProp, ED = ED,
                                 ThreeD = ThreeD, lossTolerance = lossTolerance,
                                 planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        Centre[[i]] <- out$xy
        Loss[i] <- out$loss
        STRESS[i] <- out$stress
        RSSVec[i] <- out$RSS
        TSSVec[i] <- out$TSS
      }
      index <- which(STRESS == min(STRESS))[1]
      stress <- STRESS[index]
      xy <- Centre[[index]]
      loss <- Loss[index]
      RSS <- RSSVec[index]
      TSS <- TSSVec[index]
    }
  }
  list(xy = xy, loss = loss, stress = stress, RSS = RSS, TSS = TSS)
}
# if two way intersections missing, optimize mu to get Euclidean ditance
solveWithMu <- function(ED, xy, combProp, ThreeD, radius, disProp,
                        weight, Delta, method, firstGenerate,
                        lower, upper, control, hessian, expand, scaleSeachTolerance,
                        distanceTolerance, lossTolerance,
                        stressBound, maximumStep, planeSize, twoWayGenerate) {
  stress_0 <- calculateStress(xy = xy,radius = radius,disProp = disProp,weight = weight,ThreeD = ThreeD,
                              planeSize = planeSize, twoWayGenerate = twoWayGenerate)$stress
  mu_1_Stress <- findOptimalStress(lambda = 1, xy = xy, weight = weight,
                                   radius = radius, disProp = disProp, ED = ED,
                                   ThreeD = ThreeD, lossTolerance = lossTolerance,
                                   planeSize = planeSize, twoWayGenerate = twoWayGenerate)
  stress_1 <- mu_1_Stress$stress
  if(min(stress_0,stress_1) != stress_0){
    if(method[2] == "NelderMead"){
      EDhat <- list()
      Centre <- list()
      EDplus <- newEuclideanDistance(firstGenerate = firstGenerate, mu = 1+Delta,
                                     combProp = combProp, radius = radius, ThreeD = ThreeD,
                                     expand = expand, distanceTolerance = distanceTolerance,
                                     maximumStep = maximumStep)
      mu_2_Stress <- findOptimalStress(lambda = 1, xy = xy, weight = weight,radius = radius, disProp = disProp,
                                       ED = EDplus, ThreeD = ThreeD, lossTolerance = lossTolerance,
                                       planeSize = planeSize, twoWayGenerate = twoWayGenerate)
      EDplusplus <- newEuclideanDistance(firstGenerate = firstGenerate, mu = 1+2*Delta,
                                         combProp = combProp, radius = radius, ThreeD = ThreeD,
                                         expand = expand, distanceTolerance = distanceTolerance,
                                         maximumStep = maximumStep)
      mu_3_Stress <- findOptimalStress(lambda = 1, xy = xy, weight = weight,radius = radius, disProp = disProp,
                                       ED = EDplusplus, ThreeD = ThreeD, lossTolerance = lossTolerance,
                                       planeSize = planeSize, twoWayGenerate = twoWayGenerate)
      Centre[[1]] <- mu_1_Stress$xy
      Centre[[2]] <- mu_2_Stress$xy
      Centre[[3]] <- mu_3_Stress$xy
      EDhat[[1]] <- ED
      EDhat[[2]] <- EDplus
      EDhat[[3]] <- EDplusplus
      STRESS <- c(stress_1, mu_2_Stress$stress, mu_3_Stress$stress)
      ED <- EDhat[[which(STRESS == min(STRESS))[1]]]
      centre <- Centre[[which(STRESS == min(STRESS))[1]]]
      MU <- c(1,1+Delta,1+2*Delta)[order(STRESS)]
      STRESS <- sort(STRESS)
      count <- 0
      while(count<maximumStep && min(STRESS)> stressBound) {
        #reflection
        count <- count+1
        mu_0 <- mean(MU[1:2])
        mu_R <- mu_0 + (mu_0 - MU[3])
        EDreflection <- newEuclideanDistance(firstGenerate = firstGenerate, mu = mu_R,
                                             combProp = combProp, radius = radius, ThreeD = ThreeD,
                                             expand = expand,distanceTolerance = distanceTolerance,
                                             maximumStep = maximumStep)
        mu_R_stress <- findOptimalStress(lambda = 1, xy = xy, weight = weight, radius = radius, disProp = disProp,
                                         ED = EDreflection,ThreeD = ThreeD, lossTolerance = lossTolerance,
                                         planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        stress_R <- mu_R_stress$stress
        #Reflection:
        if(STRESS[1] <= stress_R &&  stress_R< STRESS[2]){
          STRESS <- c(STRESS[1] , stress_R , STRESS[2])
          MU <- c(MU[1], mu_R, MU[2])
        }else if(stress_R < STRESS[1]){
          #Expansion
          centre <- mu_R_stress$xy
          ED <- EDreflection
          mu_E <- mu_0 + 2*(mu_R - mu_0)
          EDexpansion <- newEuclideanDistance(firstGenerate = firstGenerate, mu = mu_E,
                                              combProp = combProp, radius = radius, ThreeD = ThreeD,
                                              expand = expand, distanceTolerance = distanceTolerance,
                                              maximumStep = maximumStep)
          mu_E_stress <- findOptimalStress(lambda = 1, xy = xy, weight = weight,
                                           radius = radius, disProp = disProp, ED = EDexpansion,
                                           ThreeD = ThreeD, lossTolerance = lossTolerance,
                                           planeSize = planeSize, twoWayGenerate = twoWayGenerate)
          stress_E <- mu_E_stress$stress
          if(stress_E<stress_R){
            STRESS <- c(stress_E , stress_R , STRESS[1])
            MU <- c(mu_E, mu_R, MU[1])
            ED <- EDexpansion
            centre <- mu_E_stress$xy
          }else if(stress_R<=stress_E && stress_E <STRESS[1]){
            STRESS <- c(stress_R, stress_E, STRESS[1])
            MU <- c(mu_R, mu_E, MU[1])
          }
          else if(STRESS[1]<=stress_E && stress_E <STRESS[2]){
            STRESS <- c(stress_R, STRESS[1],stress_E)
            MU <- c(mu_R, MU[1], mu_E)
          }else{
            STRESS <- c(stress_R, STRESS[1],STRESS[2])
            MU <- c(mu_R, MU[1],MU[2])
          }
        }else if(stress_R >= STRESS[2]){
          #Contraction
          mu_C <-  mu_0 + 0.5*(MU[3] - mu_0)
          EDcontraction <- newEuclideanDistance(firstGenerate = firstGenerate, mu = mu_C,
                                                combProp = combProp,
                                                radius = radius, ThreeD = ThreeD,
                                                expand = expand, distanceTolerance = distanceTolerance,
                                                maximumStep = maximumStep)
          mu_C_stress <- findOptimalStress(lambda = 1, xy = xy, weight = weight,radius = radius, disProp = disProp,
                                           ED = EDcontraction, ThreeD = ThreeD, lossTolerance = lossTolerance,
                                           planeSize = planeSize, twoWayGenerate = twoWayGenerate)
          stress_C <- mu_C_stress$stress
          if (stress_C < STRESS[3]){
            STRESS <- c(STRESS[1],STRESS[2],stress_C)
            MU <- c(MU[1],MU[2],mu_E)[order(STRESS)]
            STRESS <- sort(STRESS)
            if(min(STRESS) == stress_C){
              ED <- EDcontraction
              centre <- mu_C_stress$xy}
          }else{
            #shrink
            MU[2] <- (MU[1] + MU[2])/2
            MU[3] <- (MU[1] + MU[3])/2
          }
        }
        if(BoolScaleNMCpp( proportional = scaleSeachTolerance$ proportional,
                           value = scaleSeachTolerance$value,
                           LAMBDA = MU, STRESS = STRESS) == FALSE){break}
      }
      xy <- centre
    } else if(method[2] == "lineSearch") {
      stress_n <- stress_1
      Center <- list()
      Center[[1]] <- mu_1_Stress$xy
      n <- 1
      mu <- 1
      EDhat <- list()
      EDhat[[1]] <- ED
      while (stress_n <= stress_1){
        stress_1 <- stress_n
        mu <- mu + Delta
        ED <- newEuclideanDistance(firstGenerate = firstGenerate, mu = mu,
                                   combProp = combProp, radius = radius, ThreeD = ThreeD,
                                   expand = expand, distanceTolerance = distanceTolerance,
                                   maximumStep = maximumStep)
        L <- findOptimalStress(lambda = 1, xy = Center[[n]], weight = weight,
                               radius = radius, disProp = disProp, ED = ED,
                               ThreeD = ThreeD, lossTolerance = lossTolerance,
                               planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        n <- n+1
        stress_n <- L$stress
        Center[[n]] <- L$xy
        EDhat[[n]] <- ED
      }
      xy <- Center[[n-1]]
      ED <- EDhat[[n-1]]
    } else if(method[2] == "goldenSectionSearch"){
      Optimization <- optimize(f = findOptimalStressMu, lower = lower[2], upper = upper[2],
                               xy = mu_1_Stress$xy, weight = weight,
                               firstGenerate = firstGenerate,combProp = combProp,
                               radius = radius, disProp = disProp, ED = ED,ThreeD = ThreeD,
                               expand = expand, lossTolerance = lossTolerance,
                               distanceTolerance = distanceTolerance,
                               maximumStep = maximumStep,
                               planeSize = planeSize, twoWayGenerate = twoWayGenerate)
      if(is.null(Optimization$minimum) || is.infinite(Optimization$minimum)){
        stop("This method does not converge")
      }else{
        ED <- newEuclideanDistance(firstGenerate = firstGenerate, mu = Optimization$minimum,
                                   combProp = combProp, radius = radius, ThreeD = ThreeD,
                                   expand = expand, distanceTolerance = distanceTolerance,
                                   maximumStep = maximumStep)
        muStress <- findOptimalStress(lambda = 1, xy = mu_1_Stress$xy,weight = weight,
                                      radius = radius, disProp = disProp, ED = ED,
                                      ThreeD = ThreeD, lossTolerance = lossTolerance,
                                      planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        xy <- muStress$xy
      }
    }else{
      Optimization <- optim(par = 1, fn = findOptimalStressMu,method = method[2],
                            lower = lower[2], upper = upper[2],
                            control = control, hessian = hessian[2],
                            xy = mu_1_Stress$xy, weight = weight,
                            firstGenerate = firstGenerate,combProp = combProp,
                            radius = radius, disProp = disProp, ED = ED,ThreeD = ThreeD,
                            expand = expand, lossTolerance = lossTolerance,
                            distanceTolerance = distanceTolerance,
                            maximumStep = maximumStep,
                            planeSize = planeSize, twoWayGenerate = twoWayGenerate)
      if(is.null(Optimization$apr) || is.infinite(Optimization$apr)){
        stop("This method does not converge")
      }else{
        ED <- newEuclideanDistance(firstGenerate = firstGenerate, mu = Optimization$apr,
                                   combProp = combProp, radius = radius, ThreeD = ThreeD,
                                   expand = expand, distanceTolerance = distanceTolerance,
                                   maximumStep = maximumStep)
        L <- findOptimalStress(lambda = 1, xy = mu_1_Stress$xy,weight = weight,
                               radius = radius, disProp = disProp, ED = ED,
                               ThreeD = ThreeD, lossTolerance = lossTolerance,
                               planeSize = planeSize, twoWayGenerate = twoWayGenerate)
        xy <- L$xy
      }
    }
  }
  list(xy = xy, ED = ED)
}
# draw sphere
sphere <- function(xy, radius, cols, alpha, oneWaySetName, smooth){
  open3d()
  if(smooth){
    n <- dim(xy)[1]
    f <- function(s, t) cbind(r * cos(s) * cos(t) + x0,
                              r * sin(s) * cos(t) + y0,
                              r * sin(t) + z0)
    for(i in 1:n){
      x0 <- xy[i,1]
      y0 <- xy[i,2]
      z0 <- xy[i,3]
      r <- radius[i]
      persp3d(f, slim = c(0, pi), tlim = c(0, 2*pi), n = 101, add = T, col = cols[i], alpha = alpha)
    }
  } else {
    spheres3d(xy[,1], xy[,2], xy[,3], radius = radius,
              color = cols, alpha = alpha)
  }
  text3d(xy[,1], xy[,2], xy[,3], oneWaySetName)
}
# if two way intersections missing, given mu to get a new Euclidean distance
newEuclideanDistance <- function(firstGenerate, mu, combProp, radius, ThreeD, expand,
                                 distanceTolerance, maximumStep){
  disjointSetNames <- names(combProp)
  numWays <- str_count(disjointSetNames,pattern = "&")+1
  #one way set
  oneWaySetName <- disjointSetNames[which(numWays==1)]
  oneWaySet <- combProp[which(numWays==1)]
  nextGenerate <- lapply(firstGenerate$New,
                         function(a, oneWaySet, mu = mu){
                           higherWay <- a[1];newGenerate <- a[-1]
                           newGenerateName <- names(newGenerate)
                           for(i in 1:length(newGenerate)){
                             newGenerate[i] <- higherWay*mu^(str_count(names(higherWay),pattern = "&")+1-2)
                             newGenerateNameSeparate <- str_split(newGenerateName,"&")[[i]]
                             minOneWay <- min(oneWaySet[which((names(oneWaySet) %in% newGenerateNameSeparate) == TRUE)])
                             if(newGenerate[i]> minOneWay){
                               newGenerate[i] <- runif(1,min = higherWay,max = minOneWay)
                             }
                           }
                           newGenerate
                         }, oneWaySet = oneWaySet, mu = mu)
  newTworWaySet <- c(combProp[which(numWays == 2)], unlist(nextGenerate))
  EuclideanDistance(newTworWaySet = newTworWaySet,
                    oneWaySetName = oneWaySetName, radius = radius,
                    ThreeD = ThreeD, initial = FALSE, expand = expand,
                    distanceTolerance = distanceTolerance, maximumStep = maximumStep)
}

# given lambda, optimize stress to get centres and corresponding values
findOptimalStress <- function(lambda, xy, weight, radius, disProp, ED, ThreeD, lossTolerance, planeSize, twoWayGenerate){
  if(is.null(lossTolerance$ALPHA)){
    bool <- TRUE
    #Just satisfy ``double'' input in Cpp; ALPHA will be generated through Newton method
    ALPHA <- 0.01
  } else{
    bool <- FALSE
    ALPHA <- lossTolerance$ALPHA
  }
  L <- lossCpp(xy = xy,radius = radius, lambda = lambda, ED = ED,
               ThreeD = ThreeD, ToleranceofLoss = lossTolerance$ToleranceofLoss, maximumStep = lossTolerance$maximumStep,
               ToleranceofStepsize = lossTolerance$ToleranceofStepsize,  proportional = lossTolerance$ proportional,
               ALPHA = ALPHA, Bool = bool)
  stressCalculation <- calculateStress(xy = L$xy,radius = radius,
                                       disProp = disProp, weight = weight, ThreeD = ThreeD,
                                       planeSize = planeSize, twoWayGenerate = twoWayGenerate)
  stress <- stressCalculation$stress
  RSS <- stressCalculation$RSS
  TSS <- stressCalculation$TSS
  ## Return a list
  list(xy = L$xy, stress = stress, loss = L$loss, RSS = RSS, TSS = TSS)
}
# used in optim(...) or optimization(...) function, optimize lambda with minimum stress
findOptimalStressLambda<- function(lambda, xy, weight, radius, disProp, ED, ThreeD, lossTolerance, planeSize, twoWayGenerate){
  stressValues <- findOptimalStress(lambda = lambda, xy = xy,weight = weight,
                                    radius = radius, disProp = disProp,
                                    ED = ED,ThreeD = ThreeD, lossTolerance = lossTolerance,
                                    planeSize = planeSize, twoWayGenerate = twoWayGenerate)
  ## Return the stress only
  stressValues$stress

}
# used in optim(...) or optimization(...) function, optimize mu with minimum stress
findOptimalStressMu<- function(mu, xy, weight, firstGenerate, combProp, radius, disProp, ED,
                               ThreeD, expand, lossTolerance, distanceTolerance, maximumStep, planeSize, twoWayGenerate){
  ED <- newEuclideanDistance(firstGenerate = firstGenerate, mu = mu,
                             combProp = combProp, radius = radius, ThreeD = ThreeD,
                             expand = expand, distanceTolerance = distanceTolerance,
                             maximumStep = maximumStep)
  stressValues <- findOptimalStress(lambda = 1, xy = xy, weight = weight,
                                    radius = radius, disProp = disProp, ED = ED,
                                    ThreeD = ThreeD, lossTolerance = lossTolerance,
                                    planeSize = planeSize, twoWayGenerate = twoWayGenerate)
  ## Return the stress only
  stressValues$stress
}

#transform disjointcombinations to not disjoint combinations
disjoint2combinations <- function(combinations){
  disjointSetNames <- names(combinations)
  #numWays <- str_length(name)
  numWays <- str_count(disjointSetNames,pattern = "&")+1
  combinations <- combinations[order(numWays)]
  numWays <- sort(numWays)
  if(max(numWays)!=1){
    for(i in 1:length(numWays)){
      if(numWays[i] == max(numWays)){break} else {
        Index <- sapply(str_split(names(combinations[-i]),pattern = "&"),
                        function(a){bool <- 0;
                        if(all(str_split(names(combinations[i]),pattern = "&")[[1]]%in%a))
                        {bool <- 1};
                        bool})
        if(length(which(Index == 1))!=0){combinations[i] <- combinations[i] + sum(combinations[-i][which(Index == 1)])}
      }
    }
  }

  combinations
}

#reorder
reorderDisjointCombinations <- function(combinations, disjoint.combinations){
  m <- length(combinations)
  newDisjointCombinations <- rep(0, m)
  combinationsName <- names(combinations)
  newOrder <- rep(0, m)
  disjoint.combinationsName <- names(disjoint.combinations)
  for(i in 1:m){
    wayOrder <- which(disjoint.combinationsName%in%combinationsName[i] == TRUE)
    newDisjointCombinations[i] <- disjoint.combinations[wayOrder]
    newOrder[i] <- wayOrder
  }
  names(newDisjointCombinations) <- combinationsName
  list(newDisjointCombinations = newDisjointCombinations, newOrder = newOrder)
}

# If input is a data list(frame), extract combinations from it
extractCombinations <- function(data, vars) {
  data <- as.data.frame(data)
  #turn character to numeric
  if(is.null(vars)==FALSE){
    allColName <- colnames(data)
    vars <- match.arg(vars , allColName,several.ok = T)
    data <- data[, which(allColName%in%vars) ]
  }
  if(is.null(dim(data))){
    data <- as.data.frame(data)
    colnames(data) <- vars
  }
  if(dim(data)[2] == 1){
    warning("Meaningless factor data frame")
  }
  colReduceForSure <- rep(0,dim(data)[2])
  for(i in 1:dim(data)[2]){
    if(is.numeric(data[,i]) && all(unique(data[,i])%in%c(0,1)) == FALSE){
      colReduceForSure[i] <- i
    } else if(length(unique(data[,i])) == 1){
      colReduceForSure[i] <- i
    }
  }
  # get rid of some numeric columns
  if(all(colReduceForSure==0)==FALSE){
    colReduceForSure <- colReduceForSure[which(colReduceForSure!=0)]
    data <- data[,-colReduceForSure]
    warning(cat(paste("Non-factor column(s)",toString(colReduceForSure),
                      "has(have) been ignored"), "\n"))
  }
  colName <- colnames(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  colAdd <- list()
  colReduce <- rep(0,p)
  for(i in 1:p){
    uniqueName <- unique(data[,i])
    if(all(uniqueName%in%c(0,1))) {next
    } else if (mode(levels(data[,i])) == "character"){
      colReduce[i] <- i
      newDataSet <- matrix(rep(0,n*length(uniqueName)),nrow = n)
      for(j in 1:length(uniqueName)){
        newDataSet[which(data[,i] == uniqueName[j]),j] <- 1
      }
      colnames(newDataSet) <- uniqueName
      colAdd[[i]] <- newDataSet
    }
  }
  colReduce <- colReduce[which(colReduce!=0)]
  if(length(colReduce)!=0) {
    newdata <- as.matrix(data[,-colReduce])
    if(dim(newdata)[2] == 1){
      colnames(newdata) <- colName[-colReduce]
    }
    for(i in 1:length(colAdd)){
      colAddName <- colnames(colAdd[[i]])
      sumAll <- apply(colAdd[[i]], 2, "sum")
      deleteCol <- which(sumAll == min(sumAll))[1]
      newInput <- as.matrix(colAdd[[i]][,-deleteCol])
      if(dim(newInput)[2] == 1){
        colnames(newInput) <- colAddName[-deleteCol]
      }
      newdata <- cbind(newdata, newInput)
    }
  } else {newdata <- data}
  rowSumZero <- which(apply(newdata,1,"sum")==0)
  if(length(rowSumZero)!=0){
    newdata <- as.matrix(newdata[-which(apply(newdata,1,"sum")==0),])
  }
  if(dim(newdata)[2] == 1) {
    disjoint.combinations <- c(sum(newdata))
    names(disjoint.combinations) <- colnames(newdata)
  } else {
    #OUTCOME is not disjoint
    G <- list()
    for(i in 1:dim(newdata)[2]){
      G[[i]] <- t(combn(dim(newdata)[2],i))
    }
    #OUTCOME is disjoint
    newColName <- colnames(newdata)
    nameList <- apply(newdata,1, function(a){paste(newColName[which(a==1)],collapse = "&")})
    disjointOutput <- aggregate(data.frame(count = nameList), list(name = nameList), length)
    disjoint.combinations <- disjointOutput[,2]
    names(disjoint.combinations) <- disjointOutput[,1]
    numWays <- str_count(disjointOutput[,1], pattern = "&")+1
    oneWays <- which(numWays == 1)
    if(length(oneWays) != length(newColName)){
      notIn <- which(newColName %in% disjointOutput[,1][oneWays] == FALSE)
      newIn <- rep(0, length(notIn))
      names(newIn) <- newColName[notIn]
      disjoint.combinations <- c(disjoint.combinations, newIn)
    }
  }
  disjoint.combinations
}

# Calculates distance between two circles
Distance <- function(r1,r2,S,ThreeD, expand, distanceTolerance, maximumStep) {
  if (ThreeD){
    if(S == 0){
      if(is.null(expand)) {
        d <- r1+r2
      } else {
        d <- (r1+r2)*expand
      }
    } else if(abs(S - min(4*pi/3*r1^3,4*pi/3*r2^3)) < distanceTolerance$value) {
      d <- max(r1,r2) - min(r1,r2)
    } else {
      theta <- matrix(c(0,0),nrow = 2)
      thetanew <- theta+1
      f1 <- pi/3*r1^3*(1-cos(theta[1]))^2*(2+cos(theta[1])) +
        pi/3*r2^3*(1-cos(theta[2]))^2*(2+cos(theta[2])) - S
      f2 <- r1*sin(theta[1]) - r2*sin(theta[2])
      k <- 0
      while(k < maximumStep){
        theta <- thetanew
        f1 <- pi/3*r1^3*(1-cos(theta[1]))^2*(2+cos(theta[1])) +
          pi/3*r2^3*(1-cos(theta[2]))^2*(2+cos(theta[2])) - S
        f2 <- r1*sin(theta[1]) - r2*sin(theta[2])
        f <- matrix(c(f1,f2),nrow = 2)
        g <- matrix(c(pi*r1^3*sin(theta[1])^3, pi*r2^3*sin(theta[2])^3,
                      r1*cos(theta[1]),-r2*cos(theta[2])),nrow = 2, byrow = T)
        g <- solve(g)
        thetanew <- theta - g%*%f
        k <- k+1
        if(BoolDistanceCpp( proportional = distanceTolerance$ proportional,
                            value = distanceTolerance$value,
                            f1 = f1, f2 = f2,
                            thetanew = thetanew, theta = theta) == FALSE){break}
      }
      if(any(thetanew>pi) || any(thetanew<0) || k == maximumStep){
        theta1 <- seq(0,pi, length = 200)
        theta2 <- seq(0,pi, length = 200)
        searchMatrix <- distanceCpp(r1, r2,theta1, theta2, S,ThreeD)
        index <- which(searchMatrix== min(searchMatrix),arr.ind=T)[1,]
        d <- r1*cos( theta1[index[1]]) + r2*cos(theta2[index[2]])
      } else {
        d <- r1*cos(thetanew[1]) + r2*cos(thetanew[2])
      }
    }
  }
  else{
    if(S==0){
      if(is.null(expand)){
        d <- r1+r2
      } else {
        d <- (r1+r2)*expand
      }
    }else if( abs(S - min(pi*r1^2,pi*r2^2)) < distanceTolerance$value ){
      d <- max(r1,r2) - min(r1,r2)
    }else{
      theta <- matrix(c(0,0),nrow = 2)
      thetanew <- theta+1
      f1 <- theta[1]*r1^2 - sin(2*theta[1])*r1^2/2 +theta[2]*r2^2 - sin(2*theta[2])*r2^2/2  - S
      f2 <- r1*sin(theta[1]) - r2*sin(theta[2])
      k <- 0
      while(k < maximumStep){
        theta <- thetanew
        f1 <- theta[2]*r2^2 - sin(2*theta[2])*r2^2/2 + theta[1]*r1^2 - sin(2*theta[1])*r1^2/2 - S
        f2 <- r1*sin(theta[1]) - r2*sin(theta[2])
        f <- matrix(c(f1,f2),nrow = 2)
        g <- matrix(c(r1^2-r1^2*cos(2*theta[1]), r2^2-r2^2*cos(2*theta[2]),
                      r1*cos(theta[1]),-r2*cos(theta[2])),nrow = 2, byrow = T)
        g <- solve(g)
        thetanew <- theta - g%*%f
        k <- k+1
        if(BoolDistanceCpp( proportional = distanceTolerance$ proportional,
                            value = distanceTolerance$value,
                            f1 = f1, f2 = f2,
                            thetanew = thetanew, theta = theta) == FALSE){break}
      }
      if(any(thetanew>pi) || any(thetanew<0) || k == maximumStep){
        theta1 <- seq(0,pi, length = 200)
        theta2 <- seq(0,pi, length = 200)
        searchMatrix <- distanceCpp(r1, r2,theta1, theta2, S,ThreeD)
        index <- which(searchMatrix== min(searchMatrix),arr.ind=T)[1,]
        d <- r1*cos( theta1[index[1]]) + r2*cos(theta2[index[2]])
      } else {
        d <- r1*cos(thetanew[1]) + r2*cos(thetanew[2])
      }
    }
  }
  return(d)
}
# Plots circles of Venn diagram
plotCircle <-  function(xy, radius,name, col, mar, ...) {
  par(mar = mar)
  a1 <-  range(c((xy + radius)[,1], (xy - radius)[,1]))
  a2 <- range(c((xy + radius)[,2], (xy - radius)[,2]))
  plot.window(a1, a2, "", asp = 1)
  theta <- seq.int(360)/360*2*pi
  for (i in 1:length(radius)){
    polygon(xy[i,1] +  radius[i]*cos(theta), xy[i,2] + radius[i]*sin(theta), col = col[i],border = col[i])
  }
  text(xy, name)
}

#rotate centres until two groups totally separated
rotateCentres <- function(xy,transxy,radius1,radius2, delta){
  for(i in 1:35){
    theta <- i/18*pi
    if(dim(xy)[2] == 3){
      rotation1 <- matrix(c(1,0,0,0,cos(theta),-sin(theta),0,sin(theta),cos(theta)),nrow = 3,byrow = T)
      rotation2 <- matrix(c(cos(theta),0,sin(theta),0,1,0,-sin(theta),0,cos(theta)),nrow = 3,byrow = T)
      rotation3 <- matrix(c(cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1),nrow = 3,byrow = T)
      rotation <- rotation1%*%rotation2%*%rotation3}
    else{
      rotation <- matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow = 2,byrow = T)
    }
    newxy <- t(rotation%*%t(transxy))
    newxy <- t(transxy[1,] - newxy[1,] + t(newxy))
    Judgement <- allDisjointCpp(xy,newxy,radius1,radius2, delta)
    if (Judgement!=0 ){break}
  }
  list(Judgement =Judgement, newxy = newxy)
}

# after finding centres, combine them with reasonable distance
combineGroups <- function(groupCentre, radiusGroup, delta){
  groupNum <- length(groupCentre)
  radius <- unlist(radiusGroup)
  namexy <- lapply(radiusGroup,
                   function(a){
                     names(a)
                   })
  groupxy <- groupCentre
  radiusxy <- radiusGroup
  for (gn in 1:(groupNum-1)){
    #Randomly find a minimum x, minimum y, maximum x, or maximum y of the second group
    cornerrandom = ceiling(runif(1,0,4))
    if(cornerrandom == 1){
      corner <- which(groupCentre[[gn+1]][,1] == min(groupCentre[[gn+1]][,1]))[1]
    }else if(cornerrandom == 2){
      corner <- which(groupCentre[[gn+1]][,1] == max(groupCentre[[gn+1]][,1]))[1]
    }else if(cornerrandom == 3){
      corner <- which(groupCentre[[gn+1]][,2] == min(groupCentre[[gn+1]][,2]))[1]
    }else{
      corner <- which(groupCentre[[gn+1]][,2] == max(groupCentre[[gn+1]][,2]))[1]
    }
    #corresponding radius
    radiusvec <- radiusGroup[[gn+1]][corner]
    #generate a centre to check which can make this centre is totally disjoint with the former group
    xyvec <- transCpp(xy = groupCentre[[gn]], radius = radiusGroup[[gn]], radiusvec = radiusvec, radiusall = radius)
    #move the second group
    transxy <- t(t(groupCentre[[gn+1]]) + xyvec - groupCentre[[gn+1]][corner,])
    #If two groups are totally seperated
    if(allDisjointCpp(groupCentre[[gn]],transxy,radiusGroup[[gn]],radiusGroup[[gn+1]],delta) == 0){
      out <- rotateCentres(groupCentre[[gn]],transxy,radiusGroup[[gn]],radiusGroup[[gn+1]],delta)
      Judgement <- out$Judgement
      newxy <- out$newxy
      while(Judgement == 0){
        xyvec <- transCpp(xy = groupCentre[[gn]], radius = radiusGroup[[gn]],
                          radiusvec = radiusvec, radiusall = radius)
        transxy <- t(t(groupCentre[[gn+1]]) + xyvec - groupCentre[[gn+1]][corner,])
        out <- rotateCentres(groupCentre[[gn]],transxy,radius[[gn]],radius[[gn+1]],delta)
        Judgement <- out$Judgement
        newxy <- out$newxy
      }
      rownames(newxy) <- namexy[[gn+1]]
    }else{newxy <- transxy}
    #renew the second group centre
    groupxy[[gn+1]] <- newxy
  }
  #Move the last group to the former one with reasonable distance
  for(gn in 1:(groupNum-1)){
    center1 <- rep(0,dim(groupxy[[gn]])[2])
    center2 <- rep(0,dim(groupxy[[gn+1]])[2])
    for(j in 1:length(center1)){
      center1[j] <- (max(groupxy[[gn]][,j]+radiusxy[[gn]])+
                       min(groupxy[[gn]][,j]-radiusxy[[gn]]))/2
      center2[j] <- (max(groupxy[[gn+1]][,j]+radiusxy[[gn+1]])+
                       min(groupxy[[gn+1]][,j]-radiusxy[[gn+1]]))/2
    }
    direction <- center1-center2
    if(delta==0){delta <- (1e-4)}
    # scale it
    direction <- direction/sqrt(sum(direction^2))*delta/10
    newxy <- closeCpp(groupxy[[gn]],groupxy[[gn+1]],radiusxy[[gn]],
                      radiusxy[[gn+1]],delta = delta,direction)$xy
    groupxy[[gn+1]] <- rbind(groupxy[[gn]],newxy)
    radiusxy[[gn+1]] <- c(radiusxy[[gn]],radiusxy[[gn+1]])
    namexy[[gn+1]] <- c(namexy[[gn]],namexy[[gn+1]])
  }
  list(xy = groupxy[[groupNum]], radius = radiusxy[[groupNum]],namexy = namexy[[groupNum]])
}

#given centres, numerically calculate each disjoint part's area and return stress
calculateStress <- function(xy, radius, disProp, weight, ThreeD, planeSize, twoWayGenerate){
  m <- length(radius)
  l <- list()
  oneWaySetName <- names(radius)
  if(ThreeD){
    xuan <- (max(xy[,1]) - min(xy[,1])+2*max(radius))/planeSize
    yuan <- (max(xy[,2]) - min(xy[,2])+2*max(radius))/planeSize
    zuan <- (max(xy[,3]) - min(xy[,3])+2*max(radius))/planeSize
    myList <- list()
    for(k in 1:m){
      for(i in 1:planeSize){
        myList[[i]] <- matrix(rep(0,planeSize^2),nrow =planeSize)
      }
      l[[k]] <- binaryIndexThreeDCpp(myList, xy, radius,k ,yuan, xuan, zuan, planeSize)
    }
    numericArea <- goThroughPixelThreeDCpp(l, m, planeSize)
  }else{
    xuan <- (max(xy[,1]) - min(xy[,1])+2*max(radius))/planeSize
    yuan <- (max(xy[,2]) - min(xy[,2])+2*max(radius))/planeSize
    for (k in 1:m){
      l[[k]] <- binaryIndexCpp(matrix(rep(0,planeSize^2),nrow =planeSize),xy,radius,k,yuan,xuan, planeSize)
    }
    numericArea <- goThroughPixelCpp(myList =  l, m = m, num = planeSize)
  }
  if(m == 1){
    TSS <- (sum(numericArea)/(planeSize*planeSize))^2*weight
    RSS <- 0
    stress <- 0
  } else {
    numericArea <- numericArea[which(getRidofZeroCpp(numericArea)!=0),]
    numericAreaMatrix <- matrix(unlist(unique(as.data.frame(numericArea))),ncol = m)
    numericAreaVector <- countCpp(numericArea,numericAreaMatrix)
    numericAreaVectorName <- apply(numericAreaMatrix,1,
                                   function(a, oneWaySetName)
                                   {paste(oneWaySetName[which(a!=0)],collapse="&")
                                   }, oneWaySetName = oneWaySetName)
    names(numericAreaVector) <- numericAreaVectorName
    lengthDisProp <- length(disProp)
    cstar <- rep(0, lengthDisProp)
    astar <- disProp
    wstar <- weight
    lengthNumericAreaVector <- length(numericAreaVector)
    numericAreaVectorIndex <- rep(0,lengthDisProp)
    numericAreaVectorNameSplit <- str_split(numericAreaVectorName,pattern = "&")
    for(i in 1:lengthDisProp){
      splitDisPropName <- str_split(names(disProp),pattern = "&")[[i]]
      index <- 0
      for (j in 1:lengthNumericAreaVector){
        if(j %in% numericAreaVectorIndex) {next}
        if(all(numericAreaVectorNameSplit[[j]] %in% splitDisPropName,
               length(numericAreaVectorNameSplit[[j]]) == length(splitDisPropName))){
          index <- j
          numericAreaVectorIndex[i] <- index
          break
        }
      }
      if(index != 0){
        cstar[i] <- numericAreaVector[index]
      }
    }
    unnecessaryOverlayIndex <- which(c(1:lengthNumericAreaVector) %in% numericAreaVectorIndex == FALSE)
    lengthofThisIndex <- length(unnecessaryOverlayIndex)
    if(lengthofThisIndex != 0 && twoWayGenerate == FALSE){
      cstar <- c(cstar, numericAreaVector[unnecessaryOverlayIndex])
      astar <- c(disProp, rep(0, lengthofThisIndex))
      wstar <- c(weight, rep(1, lengthofThisIndex))
    }
    fit <- lm(cstar~astar-1, weights = wstar)
    RSS <- sum(fit$residuals^2*wstar)
    TSS <- sum(cstar^2*wstar)
    stress <- RSS/TSS
  }
  list(stress = stress, RSS = RSS, TSS = TSS)
}

#given combinations, detect groups
groupDetection <- function(largerThanOneWaySetName, oneWaySetName){
  vertex <- list()
  if(is.null(largerThanOneWaySetName)){
    for(i in 1:length(oneWaySetName)){
      vertex[[i]] <- i
    }
  } else {
    largerThanOneWaySetLength <- length(largerThanOneWaySetName)
    for(i in 1:largerThanOneWaySetLength){
      Boolean1 <- rep(FALSE, largerThanOneWaySetLength)
      for(j in i:largerThanOneWaySetLength){
        Boolean1[j] <- any((str_split(largerThanOneWaySetName[i],"&")[[1]] %in% str_split(largerThanOneWaySetName[j],"&")[[1]])==TRUE)
      }
      Boolean1True <- which(Boolean1 == TRUE)
      vertex[[i]] <- unique(unlist(str_split(largerThanOneWaySetName[Boolean1True],"&")))
      if(i >1){
        Boolean2 <- rep(FALSE, i-1)
        for(j in 1:(i-1)){
          Boolean2[j] <- any((vertex[[i]]%in%vertex[[j]])==TRUE)
        }
        Boolean2True <- which(Boolean2 == TRUE)
        if(length(Boolean2True) == 1){
          vertex[[Boolean2True]] <- unique(c(vertex[[Boolean2True]],vertex[[i]]))
          vertex[[i]] <- NULL
        }
      }
    }
    if(length(which(sapply(vertex, is.null)))!=0){  vertex <- vertex[-which(sapply(vertex, is.null))]}
    vertexLength <- length(vertex)
    if(length(unlist(vertex)) != length(oneWaySetName)){
      difference <- length(oneWaySetName) - length(unlist(vertex))
      rest <- oneWaySetName[which((oneWaySetName%in%unlist(vertex))== FALSE)]
      for(i in 1:difference){
        vertex[[i+vertexLength]] <- rest[i]
      }
    }
    vertex <- lapply(vertex,
                     function(a,oneWaySetName){
                       which((oneWaySetName%in%a)==TRUE)
                     }, oneWaySetName = oneWaySetName)
  }
  vertex
}
# if two way intersections are missing, generate them
twoWayGeneration <- function(combProp, mu){

  #numWays <- str_length(names(combProp))
  numWays <- str_count(names(combProp),pattern = "&")+1
  oneWaySet <- combProp[which(numWays==1)]
  twoWaySet <- combProp[which(numWays == 2)]
  twoWaySetName <- names(twoWaySet)
  #highWay gives way larger than 2
  highWaySet <- sort(combProp[which(numWays > 2)],decreasing = T)
  highWaySetName <- names(highWaySet)

  highWays <- str_count(highWaySetName,pattern = "&")+1
  resam <- 0
  New <- list()
  newName <- c()
  for (i in 1:length(highWays)){
    com <- combn(highWays[i],2)
    L <- dim(com)[2]
    newTwoWaySet <- rep(0,1e5)
    newTwoWaySetName <- rep(0,1e5)
    for (j in 1:L){
      highWaySetNameSplit <- str_split(highWaySetName[i],"&")[[1]][com[,j]]
      value_a <- oneWaySet[which(names(oneWaySet)==highWaySetNameSplit[1])]
      value_b <- oneWaySet[which(names(oneWaySet)==highWaySetNameSplit[2])]
      valuemin <- min(value_a,value_b)
      highWaySetNamePaste <- paste(highWaySetNameSplit,collapse="&")
      highWaySetNamePasteReorder <- paste(c(highWaySetNameSplit[2], highWaySetNameSplit[1]),collapse="&")
      bool1 <- any( highWaySetNamePaste %in% twoWaySetName, highWaySetNamePasteReorder%in%twoWaySetName)
      if(bool1 == FALSE){
        bool2 <- any( highWaySetNamePaste %in% newName, highWaySetNamePasteReorder %in% newName )
        if(bool2 == FALSE){
          resam <- resam+1
          newName[2*resam - 1] <- highWaySetNamePaste
          newName[2*resam] <- highWaySetNamePasteReorder

          if(highWaySet[i] == 0){
            highWaySet[i] <- min(highWaySet[which(highWaySet!=0)])*0.01
          }
          if((highWaySet[i]*mu^(highWays[i]-2)) < valuemin){
            newTwoWaySet[resam] <- (highWaySet[i]*mu^(highWays[i]-2))
          }else{
            newTwoWaySet[resam] <- runif(1,highWaySet[i],valuemin)
          }
          newTwoWaySetName[resam] <- highWaySetNamePaste
        }
      }
    }
    newTwoWayIndex <- which(newTwoWaySet!=0)
    newTwoWaySet <- newTwoWaySet[newTwoWayIndex]
    if(length(newTwoWaySet) != 0){
      newTwoWaySetName <- newTwoWaySetName[newTwoWayIndex]
      names(newTwoWaySet) <- newTwoWaySetName
      New[[i]] <- c(highWaySet[i], newTwoWaySet)
    }
  }
  if(length(New)== 0){
    newTwoWaySet <- twoWaySet
  } else {
    anyNullList <- which(sapply(New, is.null))
    if(length(anyNullList) !=0 ){
      New <- New[-anyNullList]
    }

    newHighWaySet <- unlist(New)
    if(any(is.na(newHighWaySet))){newHighWaySet <- newHighWaySet[-which(is.na(newHighWaySet))]}
    numWays <- str_count(names(newHighWaySet),pattern = "&")+1
    newTwoWaySet <- c(twoWaySet, newHighWaySet[which(numWays == 2)])
  }
  list(newTworWaySet = newTwoWaySet, resam = resam, New = New)
}

#calculate Euclidean distance matrix and intial location (based on Gram matrix)
EuclideanDistance <- function(newTworWaySet, oneWaySetName, radius, ThreeD, initial, expand, distanceTolerance,
                              maximumStep){
  m <- length(oneWaySetName)
  newTworWaySetName <- names(newTworWaySet)
  ED <- matrix(rep(0, m^2),ncol = m)
  for(i in 1:m){
    for(j in i:m){
      if(j == i){next}else{
        sij1 <- newTworWaySet[which(newTworWaySetName == paste(c(oneWaySetName[i],oneWaySetName[j]),collapse = "&"))]
        sij2 <- newTworWaySet[which(newTworWaySetName == paste(c(oneWaySetName[j],oneWaySetName[i]),collapse = "&"))]
        if(length(sij1) == 0 && length(sij2) == 0 ){
          sij <- 0
        } else { sij <-  c(sij1,sij2) }
        ED[i,j] <- Distance(radius[i], radius[j], sij, ThreeD = ThreeD, expand = expand,
                            distanceTolerance = distanceTolerance, maximumStep = maximumStep)
      }
    }
  }
  ED <- ED + t(ED)
  rownames(ED) <- oneWaySetName
  colnames(ED) <- oneWaySetName
  if(initial == TRUE){
    D2 <- ED^2
    J <- diag(m)-1/m
    G <- -0.5*J%*%D2%*%J
    Em <- svd(G)
    U <- (Em$u) %*% sqrt(diag(Em$d))
    if(ThreeD){
      if (dim(U)[1]>=3){xy <- U[,1:3]
      }else{xy <- cbind(U,rep(0,2))}
    }else{
      xy <- U[,1:2]
    }
    step <- 1
    while((allConnectedCpp(xy,radius,ThreeD) == F)){
      meanvec <- apply(xy,2,"mean")
      deriction <- t(t(xy) - meanvec)
      xy <- xy - 0.1*step*deriction
      step <- step + 1
    }
    out <- list(ED = ED, xy = xy)
  }else{out <- ED}
  out
}

