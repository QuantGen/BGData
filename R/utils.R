randomString<-function(){
    paste(sample(c(0:9,letters,LETTERS),size=5,replace=TRUE),collapse="")
}
