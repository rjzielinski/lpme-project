# R function ---- "prinSurf" 
#############################################
## Written by Chen Yue on Nov 3rd, 2014
#############################################
## This code is for fitting the principal 
## surface of a data cloud in 3D space
#############################################
#############################################
#############################################
## Main Function "PrinSurf" #################
#############################################
#############################################

# dat        -- N * 3 matrix, each row represents the 3D coordinate of one data point in space
# method     -- choose from 'gb' or 'hs'. 'gb' is our proposed method where 'hs' is the orginal
#               method proposed in Hastie's paper in 1989, then you need to provide hs.neighbor
# sub.size   -- number of data point used to fit the surface, default to all data points. This is
#               used for reducing computing time
# grid.n     -- number of grid point you want to cut in the 2D space, default to 100, which should
#               work in most cases
# max.iter   -- max number of iteration
# err.tol    -- min err.tol
# kn         -- a range for number of knots selection. The output will be length(kn) different surfaces
#               and each one corresponds to one kn number.
# ini        -- choose from 'pca', 'isomap', or 'other'. If you provide 'other', you also need
#               to provide ini.para, if you provide 'isomap', you need to provide 'isomap.n', which
#               defaults to 10
# flag       -- print more information when function is running 

prinSurf <- function(dat, method='gb', sub.size=nrow(dat), grid.n=100, max.iter=10, 
                     err.tol=1e-4, kn=3:10, ini="pca", ini.para=NULL, isomap.n=10, 
                     flag=T, hs.neighbor=0.05){
  gridCurvature = function(grid.fit){
    nGrid <- sqrt(nrow(grid.fit))
    totalCurMat <- 0
    for(highDim in 1:3){
      currMat = matrix(grid.fit[,highDim], nGrid, nGrid)
      partialX2 <- (currMat[3:nGrid,2:(nGrid-1)] - 2 * currMat[2:(nGrid-1), 2:(nGrid-1)] + currMat[1:(nGrid-2),2:(nGrid-1)])
      partialXY <- (currMat[3:nGrid, 3:nGrid] - currMat[3:nGrid,1:(nGrid-2)] - currMat[1:(nGrid-2), 3:nGrid]+currMat[1:(nGrid-2),1:(nGrid-2)])/2
      partialY2 <- (currMat[2:(nGrid-1), 3:nGrid] - 2 * currMat[2:(nGrid-1), 2:(nGrid-1)] + currMat[2:(nGrid-1), 1:(nGrid-2)])
      curvature <- partialX2 + partialXY + partialY2
      totalCurMat <- totalCurMat + curvature^2
    }
    curvatureFinal <- sum(sqrt(totalCurMat))
    return(curvatureFinal)
  }
  if (method == "hs"){
    library(RANN)
    PrinSurf.HS <- function(dat, neighbor.ratio, MaxIter, Error.Tol){
      
      smooth.surface.lpl <- function(x, k, curr.proj, cdat){
        neighbor <- nn2(curr.proj, x, k=k, treetype='kd')
        beta.lpl <- apply(neighbor$nn.idx, 1, function(idx){
          design <- cbind(1,curr.proj[idx,])
          solve(t(design)%*%design)%*%t(design)%*%cdat[idx,]
        })
        fitted.lpl <- apply(cbind(t(beta.lpl),x), 1, function(xx){
          c(1,xx[10:11])%*%matrix(xx[1:9], 3, 3)
        })
        return(list(beta=t(beta.lpl), fitted=t(fitted.lpl)))
      }
      center <- apply(dat, 2, mean)
      cdat  <- cbind( (dat[,1]-center[1]),  (dat[,2]-center[2]),   (dat[,3]-center[3]))
      VarCov <- crossprod(cdat)
      eivec <- eigen(VarCov)$vectors
      score <- cdat%*%eivec[,1:2]
      score <- score * sign(score[1,1])
      prev.proj <- score
      dist.1 <- 1e8
      mse.reduce <- 1
      k <- max(2, floor(dim(dat)[1]*neighbor.ratio))
      Iter <- 0
      while((Iter <= MaxIter)&(mse.reduce > Error.Tol)){
        s.surface <- smooth.surface.lpl(x=prev.proj, k=k, curr.proj=prev.proj, cdat=cdat)
        dist.2 <- sum(apply(dat-s.surface$fitted, 1, function(x){sum(x^2)}))
        mse.reduce <- dist.1-dist.2
        dist.1 <- dist.2
        neighbor <- nn2(s.surface$fitted, cdat, k=1)
        new.proj <- apply(cbind(neighbor$nn.idx, cdat), 1, function(idx){
          A <- t(matrix(s.surface$beta[idx[1], ], 3, 3))
          return(solve(t(A[,2:3])%*%(A[,2:3]))%*%t(A[,2:3])%*%(idx[2:4]-A[,1]))
        })
        new.proj <- t(new.proj)
        prev.proj <- new.proj
        Iter <- Iter + 1
      }
      grid.x <- seq(range(prev.proj[,1])[1],range(prev.proj[,1])[2], length.out=100)
      grid.y <- seq(range(prev.proj[,2])[1],range(prev.proj[,2])[2], length.out=100)
      grid.pred <- expand.grid(grid.x, grid.y)
      grid.surf <- smooth.surface.lpl(x=grid.pred, k=k, curr.proj=prev.proj, cdat=cdat)$fitted
      height <- apply(cdat-s.surface$fitted, 1, function(x){(sum(x^2))^.5})
      return(list(final.proj=prev.proj, Iteration=Iter, proj.surf=s.surface$fitted,
                  grid.surf=grid.surf, height=height))
    } 
    ps.hs <- PrinSurf.HS(dat=dat, neighbor.ratio=hs.neighbor, MaxIter=max.iter, 
                         Error.Tol=err.tol)
    curvature <- gridCurvature(ps.hs$grid.surf)
    mse <- mean(ps.hs$height^2)
    Prinsurf <- list(Proj=ps.hs$final.proj, PS=ps.hs$proj.surf, PS.grid=ps.hs$grid.surf, No.iter=ps.hs$Iteration, 
                     Residual=ps.hs$height, MSE=mse, Curvature=curvature)
  } else if(method == "gb"){
    projection.step <- function(dat, proj, nGrid, kn){
      n.data <- nrow(dat)
      proj <- data.frame(proj)
      names(proj) <- c('X1', 'X2')
      x.gam <- gam(dat[,1]~s(X1,k=kn[1],m=2)+s(X2,k=kn[2],m=2)+s(X1, X2, k=kn[3], m=2), data=proj, method="REML")
      y.gam <- gam(dat[,2]~s(X1,k=kn[1],m=2)+s(X2,k=kn[2],m=2)+s(X1, X2, k=kn[3], m=2), data=proj, method="REML")
      z.gam <- gam(dat[,3]~s(X1,k=kn[1],m=2)+s(X2,k=kn[2],m=2)+s(X1, X2, k=kn[3], m=2), data=proj, method="REML")
      grid.x <- seq(0, 1, length.out=(nGrid))
      pred.grid <- expand.grid(grid.x, grid.x)
      names(pred.grid) <- c("X1", "X2")
      x.pred <- predict(x.gam, pred.grid)
      y.pred <- predict(y.gam, pred.grid)
      z.pred <- predict(z.gam, pred.grid) 
      grid.fit <- cbind(x.pred, y.pred, z.pred)
      newcoor <- apply(dat, 1, function(x){which.min(colSums((t(grid.fit)-x)^2))})
      y.p <- floor((newcoor-1)/nGrid)+1
      x.p <- newcoor-(y.p-1)*nGrid
      projn <- cbind(grid.x[x.p],grid.x[y.p])
      data.fit <- grid.fit[newcoor, ]
      currError <- sum((dat-data.fit)^2)
      projection.step <- list(proj=projn, prin.surf=grid.fit, surfx=x.gam, surfy=y.gam, surfz=z.gam, error=currError)
    }
    # main part of principal surface function
    set.seed(0)
    samp <- sample(1:(dim(dat))[1], sub.size)
    samp.dat <- dat[samp,]
    center <- apply(dat, 2, mean)
    csamp  <- cbind( (samp.dat[,1]-center[1]),  (samp.dat[,2]-center[2]),   (samp.dat[,3]-center[3]))
    cdat  <- cbind( (dat[,1]-center[1]),  (dat[,2]-center[2]),   (dat[,3]-center[3]))
    if(flag==T){
      print("Begin Initialization!")
    }
    if (ini=="pca"){
      VarCov <- crossprod(csamp)
      eivec <- eigen(VarCov)$vectors
      score <- csamp%*%eivec[,1:2]
      score <- score*sign(score[1,1])
      score1 <- score[,1]
      score2 <- score[,2]
      score1 <- (score1-min(score1))/(max(score1)-min(score1))
      score2 <- (score2-min(score2))/(max(score2)-min(score2))
      proj <- cbind(score1, score2)
    } else if (ini=='isomap'){
      library(vegan)
      dis <- vegdist(csamp, method="euclidean")
      ord <- isomap(dis, ndim=2, k=isomap.n)
      ord$points <- ord$points * sign(ord$points[1,1])
      proj <- cbind((ord$points[,1]-min(ord$points[,1]))/diff(range(ord$points[,1])),
                    (ord$points[,2]-min(ord$points[,2]))/diff(range(ord$points[,2])))    
    } else {
      ini.para.1 <- (ini.para[,1]-min(ini.para[,1]))/diff(range(ini.para[,1]))
      ini.para.2 <- (ini.para[,2]-min(ini.para[,2]))/diff(range(ini.para[,2]))
      proj <- cbind(ini.para.1, ini.para.2)
      proj <- proj[samp,]
    }
    if(flag==T){
      print("Initialization Successful!")
    }
    j <- 1
    errDiff <- 100
    errPath <- NULL
    resultList <- NULL
    library(mgcv)
    if((length(kn)==1) & (sum(kn) == -1)){
      while((j<=max.iter)&&(errDiff>err.tol)){
        result <- projection.step(csamp, proj, grid.n+1, rep(kn,3))
        proj.new <- result$proj
        if(j >= 2){
          errDiff <- abs(result$error - errDist)/sub.size
        }
        errDist <- result$error
        errPath <- c(errPath, errDist/sub.size)
        proj <- proj.new
        j <- j+1
      }
      curvaturePen <- gridCurvature(result$prin.surf)
      grid.fit <- result$prin.surf
      newcoor <- apply(cdat, 1, function(x){which.min(colSums((t(grid.fit)-x)^2))})
      n.search <- grid.n+1
      y.p <- floor((newcoor-1)/n.search)+1
      x.p <- newcoor-(y.p-1)*n.search
      grid.x <- seq(0,1, length.out=grid.n+1)
      projection <- cbind(grid.x[x.p],grid.x[y.p])
      projection <- data.frame(projection)
      dat.surf <- grid.fit[newcoor, ]
      height <- apply(cdat-dat.surf, 1, function(x){(sum(x^2))^.5})
      MSE <- sum(height^2)/nrow(dat)
      resultList[[1]] <- list(Proj=projection, PS=dat.surf, PS.grid=result$prin.surf, No.iter=j-1, 
                                 Residual=height, MSE=MSE, errPath=errPath, Curvature=curvaturePen)
    } else {
      for (currKn in kn) {
        j = 1
        errDiff = 100
        errPath = NULL
        projPrev <- proj
        while((j<=max.iter)&&(errDiff>err.tol)){
          result <- projection.step(csamp, projPrev, grid.n+1, kn = c(currKn, currKn, 2*currKn))
          if(j >= 2){
            errDiff <- abs(result$error - errDist)/sub.size
          }
          errDist <- result$error
          errPath <- c(errPath, errDist/sub.size)
          projPrev <- result$proj
          j <- j+1
        }
        curvaturePen <- gridCurvature(result$prin.surf)
        if(flag==T){
          print(paste("Current Df is ", currKn, ", current curvature is ", curvaturePen, sep=""))
        }
        grid.fit <- result$prin.surf
        newcoor <- apply(cdat, 1, function(x){which.min(colSums((t(grid.fit)-x)^2))})
        n.search <- grid.n+1
        y.p <- floor((newcoor-1)/n.search)+1
        x.p <- newcoor-(y.p-1)*n.search
        grid.x <- seq(0,1, length.out=grid.n+1)
        projection <- cbind(grid.x[x.p],grid.x[y.p])
        projection <- data.frame(projection)
        dat.surf <- grid.fit[newcoor, ]
        height <- apply(cdat-dat.surf, 1, function(x){(sum(x^2))^.5})
        MSE <- sum(height^2)/nrow(dat)
        resultList[[currKn]] <- list(Proj=projection, PS=dat.surf, PS.grid=result$prin.surf, No.iter=j-1, 
                                   Residual=height, Curvature=curvaturePen, MSE=MSE, errPath=errPath) 
      } 
    }
    return(resultList)
  } else {
    print("Please choose the method from 'gb' and 'hs'.")
    return(NULL)
  }
}