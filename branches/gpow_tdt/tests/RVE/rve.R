Rve <- function(phengen_recode, lower, upper) {
  pgdata <- as.matrix(read.table(phengen_recode, as.is=T, skip = 1));
  y <- pgdata[,6];
  xdat <- as.matrix(pgdata[,-c(1:6)])
  maf <- colSums(xdat) / length(y);
  xsites <- which(maf < upper & maf > lower);
  xdat <- xdat[,xsites];
  cs <- colSums(xdat[which(y==2),])
  cn <- colSums(xdat[which(y==1),])
  xsites <- c(which(cs==0&cn!=0), which(cn==0&cs!=0))
  xdat <- xdat[,xsites];

  x <- (rowSums(xdat)>0);
  m <- matrix(nrow=2,ncol=2);
  m[1,1] <- sum(x==1 & y==2);
  m[1,2] <- sum(x==0 & y==2);
  m[2,1] <- sum(x==1 & y==1);
  m[2,2] <- sum(x==0 & y==1);
  print(m);
  stat <- fisher.test(m);
  return(c(stat$p.value, stat$estimate));
}
