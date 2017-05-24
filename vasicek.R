r0.1<-5.121 ##r0=given
rinf.1<-6.179 ##  rinfinity=given
lambda.1<-0.229 ##lambda=given
k.1<-0.025  ##kappa = given (mean reversion speed)
sigma.1<-126/12  ##sigma per period
vol.1<-0.0126*(sqrt(1/12))*100  ##volatility per period
theta<-rinf.1+(lambda.1/k.1)  ##theta=rinf+lambda/k

steps1=30    ##number of periods
vasicektree1<-matrix(0,steps1+1,steps1+1)  ##matrix to save the interest rates in
vasicektree1[1,1]<-r0.1  ##r0
vasicektree1[1,2]<-r0.1+k.1*(theta-r0.1)/12+vol.1  ##first period up value
vasicektree1[2,2]<-r0.1+k.1*(theta-r0.1)/12-vol.1  ##first period down value

middleexpec<-r0.1+k.1*(theta-r0.1)/12  ##expected value after first time step

vasicektree1[2,3]<-middleexpec+k.1*(theta-middleexpec)/12  ##expected value after second time step (recombined rate for period 2)

secmiddleup<-vasicektree1[1,2]+k.1*(theta-vasicektree1[1,2])/12   ##expected rate after second step for up value
secmiddledown<-vasicektree1[2,2]+k.1*(theta-vasicektree1[2,2])/12   ##expected rate after second step for down value

##Calulation of probability and r(uu) and r(dd) using quadratic equation factorisation:
##Using factorisation: delta=b^2-4ac && values=(-b+/-Sqrt(delta))/2a
aup<-secmiddleup-vasicektree1[2,3]
bup<-((vasicektree1[2,3]-secmiddleup)^2)-(2*(secmiddleup^2))+(2*secmiddleup*vasicektree1[2,3]-(vol.1^2))
cup<-((secmiddleup-vasicektree1[2,3])*((secmiddleup)^2))-(((vasicektree1[2,3]-secmiddleup)^2)*secmiddleup)+(vasicektree1[2,3]*(vol.1^2))
deltaup<-(bup^2)-(4*aup*cup)
##For Up value taking the maximum
irup<-max((-bup+sqrt(deltaup))/(2*aup),(-bup-sqrt(deltaup))/(2*aup))
vasicektree1[1,3]<-irup

adown<-vasicektree1[2,3]-secmiddledown
bdown<-(vol.1^2)-((vasicektree1[2,3]-secmiddledown)^2)-(2*vasicektree1[2,3]*secmiddledown)+(2*(secmiddledown^2))
cdown<-(secmiddledown*((vasicektree1[2,3]-secmiddledown)^2))+((secmiddledown^2)*(vasicektree1[2,3]-secmiddledown))-(vasicektree1[2,3]*(vol.1^2))
deltadown<-(bdown^2)-(4*adown*cdown)
##For Up value taking the minimum
irdown<-min((-bdown+sqrt(deltadown))/(2*adown),(-bdown-sqrt(deltadown))/(2*adown))
vasicektree1[3,3]<-irdown

for (i in 4:(steps1+1)){ 
  currentrow1=1
  if((i %% 2)==0){ # for even columns
    for(j in 2:(i-2)){ 
      ##rates using 50% probabilities to find the center of tree
      vasicektree1[j,i]=vasicektree1[j,i-1]+k.1*(theta-vasicektree1[j,i-1])/12+(vol.1)   
      vasicektree1[j+1,i]=vasicektree1[j,i-1]+k.1*(theta-vasicektree1[j,i-1])/12-(vol.1) 
    }
    for(m in seq(1,i,(i-1))){
      middleup=vasicektree1[1,(i-1)]+k.1*(theta-vasicektree1[1,(i-1)])/12   ##expected rate after prev step for up value
      ##Solving the quadratic equation:
      aupl = middleup - vasicektree1[2,i]
      bupl = ((vasicektree1[2,i]-middleup)^2)-(2*(middleup^2))+(2*middleup*vasicektree1[2,i]-(vol.1^2))
      cupl = ((middleup-vasicektree1[2,i])*((middleup)^2))-(((vasicektree1[2,i]-middleup)^2)*middleup)+(vasicektree1[2,i]*(vol.1^2))
      deltaupl = (bupl^2)-(4*aupl*cupl)
      irupl = max((-bupl+sqrt(deltaupl))/(2*aupl),(-bupl-sqrt(deltaupl))/(2*aupl))

      middledown=vasicektree1[(i-1),(i-1)]+k.1*(theta-vasicektree1[(i-1),(i-1)])/12   ##expected rate after prev step for down value
      ##Solving the quadratic equation:
      adownl = vasicektree1[(i-1),i]-middledown
      bdownl = (vol.1^2)-((vasicektree1[(i-1),i]-middledown)^2)-(2*vasicektree1[(i-1),i]*middledown)+(2*(middledown^2))
      cdownl = (middledown*((vasicektree1[(i-1),i]-middledown)^2))+((middledown^2)*(vasicektree1[(i-1),i]-middledown))-(vasicektree1[(i-1),i]*(vol.1^2))
      deltadownl = (bdownl^2)-(4*adownl*cdownl)
      irdownl = min((-bdownl+sqrt(deltadownl))/(2*adownl),(-bdownl-sqrt(deltadownl))/(2*adownl))
      ##rates using the new probabilities
      vasicektree1[1,i]=irupl
      vasicektree1[i,i]=irdownl
    }
    
  } else { #odd columns
    for (n in 2:(i-1)){ 
      middleexpec = vasicektree1[(n-1),(i-2)]+k.1*(theta-vasicektree1[(n-1),(i-2)])/12  ##expected value after first time step
      vasicektree1[n,i]=middleexpec+k.1*(theta-middleexpec)/12  ##expected value after second time step (recombined rate for period n+1)
    }
    for(w in seq(1,i,(i-1))){
      
      middleup=vasicektree1[1,(i-1)]+k.1*(theta-vasicektree1[1,(i-1)])/12   ##expected rate after prev step for up value
      ##Solving the quadratic equation:
      aupl = middleup - vasicektree1[2,i]
      bupl = ((vasicektree1[2,i]-middleup)^2)-(2*(middleup^2))+(2*middleup*vasicektree1[2,i]-(vol.1^2))
      cupl = ((middleup-vasicektree1[2,i])*((middleup)^2))-(((vasicektree1[2,i]-middleup)^2)*middleup)+(vasicektree1[2,i]*(vol.1^2))
      deltaupl = (bupl^2)-(4*aupl*cupl)
      irupl = max((-bupl+sqrt(deltaupl))/(2*aupl),(-bupl-sqrt(deltaupl))/(2*aupl))
      
      middledown=vasicektree1[(i-1),(i-1)]+k.1*(theta-vasicektree1[(i-1),(i-1)])/12   ##expected rate after prev step for down value
      ##Solving the quadratic equation:
      adownl = vasicektree1[(i-1),i]-middledown
      bdownl = (vol.1^2)-((vasicektree1[(i-1),i]-middledown)^2)-(2*vasicektree1[(i-1),i]*middledown)+(2*(middledown^2))
      cdownl = (middledown*((vasicektree1[(i-1),i]-middledown)^2))+((middledown^2)*(vasicektree1[(i-1),i]-middledown))-(vasicektree1[(i-1),i]*(vol.1^2))
      deltadownl = (bdownl^2)-(4*adownl*cdownl)
      irdownl = min((-bdownl+sqrt(deltadownl))/(2*adownl),(-bdownl-sqrt(deltadownl))/(2*adownl))
      ##rates using the new probabilities
      vasicektree1[1,i]=irupl
      vasicektree1[i,i]=irdownl
    }
    
  }
  vasicektree1
}

write.csv(vasicektree1,file="Vasicek Tree.csv")  ##saving the matrix as CSV file
