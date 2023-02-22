for(j in 1:169){
    system2("Rscript", args=c("Src/slimed.R", paste(j)))
}