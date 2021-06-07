etaclass<-function(patients, beds, quality){
  # calculate the fraction of rejected patients
  # this is 0 when occupied fraction is below quality, 1 when fraction is above
  # 1, and a parable between those 2 values
  # patients: number of patients
  # beds: total number of beds available
  # quality: fraction occupied over which system becomes inefficient
  x<-patients/beds
  a<-1/(2*quality*(quality-1)-(quality^2-1))
  b<-2*a*quality
  c<-a-b
  if(x>1){
    return(1)
  }
  if(x<=quality){
    return(0)
  }
  else{
    return(1-(x*(b-a*x)+c))
  }
}

