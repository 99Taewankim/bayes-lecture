#정규분포 N(1,1) N(3,1) 반반 합친 target dist?
count=0
init=0
var=5
fn<-function(x){
  ans=1/2*(1/sqrt(2*pi))*exp(-(x-1)^2/2)+1/2*(1/sqrt(2*pi))*exp(-(x-3)^2/2)
  return(ans)
}
save_x<-c()
for(i in 1:10000){
  sampx=rnorm(1,mean=init,var)
  sampunif<-runif(1,0,1)
  if (sampunif<fn(sampx)/fn(init)) {
    init<-sampx
    count<-count+1}
    
  else{
    init=init}
  save_x<-c(save_x,init)
  
}
save_x
acceptrate<-count/10000;acceptrate
hist(save_x)               #proposal pdf
plot(save_x,type='l')    #일정한 잔차를 가지는가?
