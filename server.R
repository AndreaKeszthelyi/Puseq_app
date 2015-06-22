# server.R
#Puseq_app

library(shiny)


shinyServer(function(input, output,session) {
  
  ############################################################ sample table #####################################################################   
  
  
  
  sample<-cbind(as.character(c("chr1","chr1", "chr1", "chr1")), as.numeric(c(150, 300, 450, 600)), as.numeric(c(345, 657,148,389)))
  colnames(sample)<-c("chro","pos","count")
  sample<-as.data.frame(sample)
  output$sample <- renderTable({sample})
  
  
  
  ##########################################################Reading in files #####################################################################  
  
  
  dat_df_react<-reactive({
    inFile_df <- input$file_df
    if (is.null(inFile_df))
      return(NULL)
    return(   read.csv(inFile_df$datapath, header= TRUE))})
  
  dat_dr_react<-reactive({
    inFile_dr <- input$file_dr
    if (is.null(inFile_dr))
      return(NULL)
    return(   read.csv(inFile_dr$datapath, header= TRUE))}) 
  
  dat_ef_react<-reactive({
    inFile_ef <- input$file_ef
    if (is.null(inFile_ef))
      return(NULL)
    return(   read.csv(inFile_ef$datapath, header= TRUE))})
  
  dat_er_react<-reactive({
    inFile_er <- input$file_er
    if (is.null(inFile_er))
      return(NULL)
    return(   read.csv(inFile_er$datapath, header= TRUE))})
  
  
  ########################################################## extracting chromosome numbers #####################################################################   
  
  chro_react<-reactive({unique(dat_df_react()$chro)})
    
  
  observe({
    updateSelectInput(
      session,
      "chromo",
    choices=paste(chro_react()))
    
  })
  
  observe({
    updateSelectInput(
      session,
      "chromo_2",
      choices=paste(chro_react()))
    
  })
  
  ########################################################## extracting bin size #####################################################################   
  
  bin_react<-reactive({dat_df_react()$pos[2]-dat_df_react()$pos[1]})
  
  ################################################################## sliders ##########################################################################
  
  output$ui <- renderUI({
   
    dat_df<-dat_df_react()
    dat_df.chr = dat_df[dat_df$chro==input$chromo,]
    if (length(dat_df.chr[,3])>0){maxxval<-((length(dat_df.chr[,3])*(dat_df.chr$pos[2]-dat_df.chr$pos[1])))-input$xrange}else{maxxval<-10000}
  
    switch(input$input_type,
           
           "numeric" =  numericInput("dynamic", "x axis start",
                                     value = 300),
           
           "slider" = sliderInput("dynamic", "x axis start",
                                  min = 1, max = maxxval, value = 1, step = input$xrange))
    
  })

  
  
  
  output$ui_2 <- renderUI({
    dat_df<-dat_df_react()
    dat_df.chr = dat_df[dat_df$chro==input$chromo_2,]
    if (length(dat_df.chr[,3])>0){maxxval_2<-((length(dat_df.chr[,3])*(dat_df.chr$pos[2]-dat_df.chr$pos[1])))-input$xrange}else{maxxval_2<-10000}

    
    switch(input$input_type_2,
           
           
           
           
           "slider" = sliderInput("dynamic_2", "x axis start",
                                  min = 1, max = maxxval_2 , value = 1, step = input$xrange),
            "numeric" =  numericInput("dynamic_2", "x axis start",
                              value = 300))
    
  })
  
  
  
  
  
  
  output$maxy<-renderUI({
    
    dat_df<-dat_df_react()
    dat_dr<-dat_dr_react()
    dat_ef<-dat_ef_react()
    dat_er<-dat_er_react()
    
    if (max(c(dat_df[,3],dat_dr[,3],dat_ef[,3],dat_er[,3]))>0){maxyval<-max(c(dat_df[,3],dat_dr[,3],dat_ef[,3],dat_er[,3]))
    }else{maxyval<-8000}
    
    sliderInput("y.max","y axis max:",
                min = 0,
                max = maxyval,
                step = 50,
                value = 2000)
  })
  
  output$chromoname<-renderUI({
    
    dat_df<-dat_df_react()
    if(is.null(dat_df)){
      
    chromname<-paste("chr")}else{
    chromname<-gsub("[0-9]","",dat_df$chr[1])}
    textInput("chromoname_in", label = ("Enter chromosome name for wig and bedgraph files (match with the chromosome name in your genome browser"),value = paste(chromname))
    
  })
 
  
  
  ###chro
  cur_val.chr <- ""
  observe({
    if (cur_val.chr != input$chromo){
      updateTextInput(session, "chromo_2", NULL, input$chromo)
      cur_val.chr <<- input$chromo
    }
  })
  
  observe({
   if (cur_val.chr != input$chromo_2){
      updateTextInput(session, "chromo", NULL, input$chromo_2)
      cur_val.chr <<- input$chromo_2
    }
  })
  
 
  
  
  ######################################################################### Calculations ##############################################################################################################  
  
  ######################## Counts ####################
  
  
  dat_df.chr.react<-reactive({
    validate(
      need(input$file_df, 'Please check if you have uploaded csv file for delta forward counts if you have please be patient it might take a while'))
    dat_df<-dat_df_react()
    return(dat_df[dat_df$chro==input$chromo,])
  })
  
  dat_dr.chr.react<-reactive({
    validate(
      need(input$file_dr, 'Please check if you have uploaded csv file for delta reverse counts if you have please be patient it might take a while'))
    dat_dr<-dat_dr_react()
    return(dat_dr[dat_dr$chro==input$chromo,])
  })
  dat_ef.chr.react<-reactive({
    validate(
      need(input$file_ef, 'Please check if you have uploaded csv file for epsilon forward counts if you have please be patient it might take a while'))
    dat_ef<-dat_ef_react()
    return(dat_ef[dat_ef$chro==input$chromo,])
  })
  
  dat_er.chr.react<-reactive({
    validate(
      need(input$file_er, 'Please check if you have uploaded csv file for epsilon reverse counts if you have please be patient it might take a while'))
    dat_er<-dat_er_react()
    return(dat_er[dat_er$chro==input$chromo,])
  })
  
  
  
  
  
  
  
  
  ############ ratios table ######################################
  table.ratio.react<-reactive({
    

    
    
    withProgress(message = 'Calculation of ratios in progress',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:300) {
                     incProgress(1/300)
                   }
                 })
    
    chro<-chro_react()
    table.ratio<-data.frame()
    for(chromo in chro[1:length(chro)]){ 
      
      
      dat_df<-dat_df_react()
      dat_dr<-dat_dr_react()
      dat_ef<-dat_ef_react()
      dat_er<-dat_er_react()
      
      if (is.null(dat_df) | is.null(dat_dr) | is.null(dat_ef) | is.null(dat_er))
        return(NULL)
      
      
      dat_df.chr = dat_df[dat_df$chro==chromo,]
      dat_dr.chr = dat_dr[dat_dr$chro==chromo,]
      dat_ef.chr = dat_ef[dat_ef$chro==chromo,]
      dat_er.chr = dat_er[dat_er$chro==chromo,]
      
      
      ### normalization to all counts 
      
      dat_df.chr.norm <- dat_df.chr[,3]/sum(dat_df[,3])
      dat_dr.chr.norm <- dat_dr.chr[,3]/sum(dat_dr[,3])
      dat_ef.chr.norm <- dat_ef.chr[,3]/sum(dat_ef[,3])
      dat_er.chr.norm <- dat_er.chr[,3]/sum(dat_er[,3])
      
      ### Calculating ratios 
      
      dat_df.chr.ratio <- dat_df.chr.norm/(dat_df.chr.norm+dat_ef.chr.norm)
      dat_dr.chr.ratio <- dat_dr.chr.norm/(dat_dr.chr.norm+dat_er.chr.norm)
      dat_ef.chr.ratio <- dat_ef.chr.norm/(dat_df.chr.norm+dat_ef.chr.norm)
      dat_er.chr.ratio <- dat_er.chr.norm/(dat_dr.chr.norm+dat_er.chr.norm)
      
      
      ##### calculating moving average for the ratios
      
      N = input$N.ratio.num   # the parameter for moving ave (2N+1)
      
      dat_df.chr.ratio.ma <- moving.ave.v2(dat_df.chr.ratio, N)
      dat_dr.chr.ratio.ma <- moving.ave.v2(dat_dr.chr.ratio, N)
      dat_ef.chr.ratio.ma <- moving.ave.v2(dat_ef.chr.ratio, N)
      dat_er.chr.ratio.ma <- moving.ave.v2(dat_er.chr.ratio, N)
      
      
      table.ratio.chro<-cbind(as.character(dat_df.chr[,1]),dat_df.chr[,2],dat_df.chr.ratio.ma,dat_dr.chr.ratio.ma,dat_ef.chr.ratio.ma,dat_er.chr.ratio.ma)
      table.ratio<-rbind(table.ratio,table.ratio.chro)
    }
      colnames(table.ratio)<-(c("chromosome","position","delta on forward","delta on reverse","epsilon on forward","epsilon on reverse"))
      return(table.ratio)
    
  })
  
 
 
  

  
  ############################################ origins #################################
  
  
  
  table.ori.react<-reactive({
   
    
    withProgress(message = 'Calculation of origins in progress',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:300) {
                     incProgress(1/300)
                   }
                 })
    chro<-chro_react()
    orilist<-data.frame()
    ratio.table<-table.ratio.react()
    for(chromo in chro[1:length(chro)]){  
      
      ratio.table.chr = ratio.table[ratio.table$chromosome==chromo,]

##calculate differencials
##NN = moving average 


  NN=input$N.ori.num
  
  dat_ef.chr.diff<-moving.ave.v2(diff.sequence(as.numeric(as.vector(ratio.table.chr[,5]))),NN)
  dat_dr.chr.diff<-moving.ave.v2(diff.sequence(as.numeric(as.vector(ratio.table.chr[,4]))),NN)
  
#find local maxima - p= percentile treshold
  p=input$P.ori.num
  
  dat_ef.chr.table<-Findlocalmax(dat_ef.chr.diff,(as.numeric(as.vector(ratio.table.chr[,2]))),p)
  dat_dr.chr.table<-Findlocalmax(dat_dr.chr.diff,(as.numeric(as.vector(ratio.table.chr[,2]))),p)
  
  
#merge peaks within 4 bins
  
  dat_ef.chr.table.merged<-Closeori(as.numeric(as.vector(dat_ef.chr.table[,1])),bin_react())
  dat_dr.chr.table.merged<-Closeori(as.numeric(as.vector(dat_dr.chr.table[,1])),bin_react())
  
#calculate efficiency
  
  dat_ef.chr.orieff<-orieff(dat_ef.chr.table.merged,as.numeric(as.vector(ratio.table.chr[,5])),as.numeric(as.vector(ratio.table.chr[,2])))
  dat_dr.chr.orieff<-orieff(dat_dr.chr.table.merged,as.numeric(as.vector(ratio.table.chr[,4])),as.numeric(as.vector(ratio.table.chr[,2])))
  
#find peaks that are in both datasets within plusminus 4 bins
  
  orilist.chr<-orieff_merge(dat_ef.chr.orieff[,2],dat_ef.chr.orieff[,1],dat_dr.chr.orieff[,2],dat_dr.chr.orieff[,1],chromo,bin_react())
  orilist<-rbind(orilist,orilist.chr)
  
    }
  return(orilist)
  })
  
 output$orilist<-renderTable({table.ori.react()})
  
  ############################################################################# Plotting ##############################################################################################################  
 

######################## Counts ####################

  
  output$countPlot <- renderPlot({
    if(is.null(input$dynamic))
      return()
    
    
    
    
    withProgress(message = 'Plotting in progress',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:150) {
                     incProgress(1/150)
                   }
                 })
    
    
    
    dat_df.chr<-dat_df.chr.react()
    dat_dr.chr<-dat_dr.chr.react()
    dat_ef.chr<-dat_ef.chr.react()
    dat_er.chr<-dat_er.chr.react()
  
    if(is.null(input$y.max)){
      Y_lim<-1000}
    else{Y_lim<-c(0,input$y.max)
      }
    
    

   X_lim<-c((input$dynamic),(input$dynamic+max(input$xrange)))

    

 
  par(mfrow=c(1,1), omd=c(0, 1, 0.2, 0.9), plt=c(0.15, 0.95, 0, 1),lwd=2)
  plot(dat_df.chr[,2],dat_df.chr[,3], xlim = X_lim, col = 1,type="b",ylim = Y_lim, ylab="Counts")
  abline(h=c(input$y.max/4,input$y.max/2,input$y.max/4+input$y.max/2),lwd=1, lty=2)
  par(new=T)
  plot(dat_df.chr[,2],dat_dr.chr[,3], xlim = X_lim, col = 2,type="b",ylab="",ylim = Y_lim,axes = F)
  par(new=T)
  plot(dat_df.chr[,2],dat_ef.chr[,3], xlim = X_lim, col = 7,type="b",ylab="",ylim = Y_lim,axes = F)
  par(new=T)
  plot(dat_df.chr[,2],dat_er.chr[,3], xlim = X_lim, col = 4,type="b",ylab="",ylim = Y_lim,axes = F)
  mtext(paste("Position on ", input$chromo_2), side=1, cex=1.5, line = 3)
  add_legend("topright", legend=c("delta forward", "delta reverse", "epsilon forward", "epsilon reverse"), pch=20, 
             col=c(1,2,7,4),
             horiz=TRUE, bty='n', cex=1.5)
  
  
})

######################### Ratios ########################################

output$ratioPlot <- renderPlot({
  if(is.null(input$dynamic_2))
    return()
 
    withProgress(message = 'Plotting in progress',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                   }
                 })
  ######### plot parameters #########
xi = 50
xi =xi*1000

yl1 = c(0, 1)
yii1 = c(0.25, 0.5, 0.75)


s = 2
p1 = 21
p2 = 21



fs = 1.3
fs2 = 1.5

xlabel= paste("Position on ", input$chromo_2)
ylabel1 = "Number of Reads"
ylabel2 = "Relative amounts to ref."

col1 = "royalblue2"
col2 = "orangered2"
col1.bg = "skyblue2"
col2.bg = "salmon"

col3 = "black"
col3.p = "royalblue2"
col3.m = "orangered2"

x1<-input$dynamic_2
x2<-input$dynamic_2+max(input$xrange_2)
x3<-input$dynamic_2+max(input$xrange_2)/2



#Data to plot

dat_df.chr<-dat_df.chr.react()
dat_dr.chr<-dat_dr.chr.react()
dat_ef.chr<-dat_ef.chr.react()
dat_er.chr<-dat_er.chr.react()

chromo<-input$chromo

ratio.table<-table.ratio.react()
ratio.table.chr = ratio.table[ratio.table$chromosome==chromo,]


dat_df.chr.ratio.mean.ma <-as.numeric(as.vector(ratio.table.chr[,3]))
dat_dr.chr.ratio.mean.ma <- as.numeric(as.vector(ratio.table.chr[,4]))
dat_ef.chr.ratio.mean.ma <-as.numeric(as.vector(ratio.table.chr[,5]))
dat_er.chr.ratio.mean.ma <- as.numeric(as.vector(ratio.table.chr[,6]))


orilist<-table.ori.react()
orilist.chr = orilist[orilist$chromosome==chromo,]

orilist.chr.pos <- as.numeric(as.vector(orilist.chr[,2]))
orilist.chr.eff <- as.numeric(as.vector(orilist.chr[,3]))

##plots

par(mfrow=c(3,1), omd=c(0, 1, 0.2, 0.9), plt=c(0.15, 0.95, 0, 1),lwd=2)
plot(dat_df.chr[,2],dat_df.chr.ratio.mean.ma, cex = s, yaxt = "n", xaxt = "n",ylim = yl1,xlim = c(x1,x2), xlab="", ylab="forward", lwd=3, col=col1, bg=col1.bg, pch=p1, type="l")
par(new=T)
plot(dat_df.chr[,2],dat_ef.chr.ratio.mean.ma, cex = s, yaxt = "n", xaxt = "n",ylim = yl1,xlim = c(x1,x2), xlab="", ylab="", lwd=3, col=col2, bg=col2.bg, pch=p2, type="l",bty = "n")
abline(h=0.5, lwd=1, lty=2)
#axis(3, at=c(x1, x3,x2),labels=c(x1, x3,x2), las=0, cex.axis=fs, lwd=2)
axis(2, at=yii1,labels=yii1, las=1, cex.axis=fs, lwd=2)
legend("topright", legend=c("delta", "epsilon"), pch=20, 
           col=c(col1,col2),
           horiz=TRUE, bty='n', cex=1.3)


plot(dat_df.chr[,2],dat_dr.chr.ratio.mean.ma, cex = s, yaxt = "n", xaxt = "n",ylim = yl1,xlim = c(x1,x2), xlab="", ylab="reverse", lwd=3, col=col1, bg=col1.bg, pch=p1, type="l")
par(new=T)
plot(dat_df.chr[,2],dat_er.chr.ratio.mean.ma, cex = s, yaxt = "n", xaxt = "n",ylim = yl1,xlim = c(x1,x2), xlab="", ylab="", lwd=3, col=col2, bg=col2.bg, pch=p2, type="l",bty = "n")
abline(h=0.5, lwd=1, lty=2)
axis(2, at=yii1,labels=yii1, las=1, cex.axis=fs, lwd=2)
legend("topright", legend=c("delta", "epsilon"), pch=20, 
       col=c(col1,col2),
       horiz=TRUE, bty='n', cex=1.3)




plot(orilist.chr.pos,orilist.chr.eff, cex = s, yaxt = "n", cex.axis = fs, ylim = c(0, 100),xlim = c(x1,x2), xlab="", ylab="origins", lwd=3, col=col1, bg=col1.bg, pch=p1, type="h",lend="square")
axis(2, at=c(0,25, 50, 75, 100),labels=c(0,25, 50, 75, 100), las=1, cex.axis=fs, lwd=2)
par(new=T)
abline(h=c(25,50,75), lwd=1, lty=2)

mtext(xlabel, side=1, cex=fs2, line = 3)


})



########################################################## Files to download ###############################################################################


###Bedgraph file

orieff_bedgraph.react<-reactive({
    chro<-chro_react()
    orieff_bedgraph<-data.frame()
    chromoname_in<-input$chromoname_in
    for(chromo in chro[1:length(chro)]){  
      
      
      table.ori.bedgraph<-table.ori.react()
      
      if (is.null(table.ori.bedgraph))
        return(NULL)
      
      table.ori.bedgraph.chr = table.ori.bedgraph[table.ori.bedgraph$chromosome==chromo,]
      
      chromoname=paste(chromoname_in,gsub("[^0-9]", "",chromo),sep="")
      
  
      pos1<-table.ori.bedgraph.chr$maxpos
      pos2<-as.numeric(as.vector(pos1))+bin_react()
      value<-table.ori.bedgraph.chr$efficiency
      Chromosome<-rep(chromoname, length(table.ori.bedgraph.chr$maxpos))
      orieff_bedgraph.chr<-data.frame(chro=Chromosome, start=pos1, end=pos2, value=value)
      colnames(orieff_bedgraph.chr) <-c("chro", "start", "end", "value")
      orieff_bedgraph<-rbind(orieff_bedgraph,orieff_bedgraph.chr)
    }
  return(orieff_bedgraph)
})



#### wig files 


df_wig.react<-reactive({Wigdata(paste(input$name,"_delta_forward",sep=""), table.ratio.react(),3,bin_react(),"blue",0.5,chro_react(),input$chromoname_in)})
dr_wig.react<-reactive({Wigdata(paste(input$name,"_delta_reverse",sep=""), table.ratio.react(),4,bin_react(),"blue",0.5,chro_react(),input$chromoname_in)})
ef_wig.react<-reactive({Wigdata(paste(input$name,"_epsilon_forward",sep=""), table.ratio.react(),5,bin_react(),"red",0.5,chro_react(),input$chromoname_in)})
er_wig.react<-reactive({Wigdata(paste(input$name,"_epsilon_reverse",sep=""), table.ratio.react(),6,bin_react(),"red",0.5,chro_react(),input$chromoname_in)})

  
  
### file selection  

datasetInput <- reactive({
  switch(input$dataset,
         "all ratios as csv"= table.ratio.react(), 
         "delta forward ratio as wig" = df_wig.react(),
         "delta reverse ratio as wig" = dr_wig.react(),
         "epsilon forward ratio as wig" = ef_wig.react(),
         "epsilon reverse ratio as wig" = er_wig.react(),
         "origin efficiency as csv" = table.ori.react(),
         "origin efficiency as bedgraph" = orieff_bedgraph.react())
})

##### File output

output$downloadData <- downloadHandler(
  filename = function() {
    if (input$dataset=="all ratios as csv" | input$dataset== "origin efficiency as csv"){
    paste(input$name,"_",input$dataset, '.csv', sep='')}
    else if(input$dataset == "origin efficiency as bedgraph"){
      paste(input$name,"_",input$dataset, '.bedgraph', sep='')}
    else{
      paste(input$name,"_",input$dataset, '.wig', sep='')
    }
  },
  content = function(file) {
    
    if (input$dataset=="all ratios as csv" | input$dataset== "origin efficiency as csv"){
    write.csv(datasetInput(), file, row.names = FALSE)}
    else if(input$dataset == "origin efficiency as bedgraph"){
      
    write.table(datasetInput(), file, append = F, quote=F, col.names=F, row.names=F)}
    else{writeLines(unlist(datasetInput()),file,sep="\n")

      }
      
    }
    
  
)





###########################################################functions############################################################
moving.ave.v2 <- function(data, n){ # subroutine to calculate moving averages
  
  dataN <- length(data)
  
  start <- c(rep(1,n), 1:(dataN-n))
  end   <- c((n+1):dataN, rep(dataN,n))
  se <- cbind(start, end) 
  
  average.se <- function(n)mean(data[n[1]:n[2]][!is.nan(data[n[1]:n[2]])])
  r <- apply(se, 1, average.se)
  
  (r)
}




plot.iceberg <- function(x, f, L, col.p, col.m){
  f[is.na(f)]<-L
  pm.area <- plus.minus.area(f-L)
  
  for(i in 1:nrow(pm.area)){
    
    if(pm.area[i,3]==1){col.pm=col.p}
    if(pm.area[i,3]==-1){col.pm=col.m}
    
    polygon(c(x[pm.area[i,1]:pm.area[i,2]],rev(x[pm.area[i,1]:pm.area[i,2]])),
            c(c(rep(L,length(pm.area[i,1]:pm.area[i,2]))),rev(f[pm.area[i,1]:pm.area[i,2]])),
            col=col.pm, border = NA)
  }
}

plus.minus.area <- function(diff){
  
  d=0
  
  s=1
  
  area <- matrix(nrow=0,ncol=3)
  
  
  for(i in 1:length(diff)){
    if(diff[i]*d<0 || i==length(diff)){
      e = i
      if(d>0){area <- rbind(area, c(s,e, 1))}
      if(d<0){area <- rbind(area, c(s,e,-1))}
      #cat(s,e,"\n")
      s = i
      
    }
    d = diff[i]
  }
  
  
  return(area)
}

diff.sequence <- function(vec){
  diff <-c()
  diff[1]=0;
  for(i in 2:length(vec)){
    diff[i] = vec[i]-vec[i-1]
  }
  return(diff)
}





Findlocalmax<-function(diffdata,position,percentile){
  
  max<-c(which(diff(c(TRUE,diff(diffdata)>=0,FALSE))<0 & diffdata>0) )
  tableall<-cbind(position[max], diffdata[max])
  perc<-quantile(tableall[,2],percentile)
  per<-which(tableall[,2]<=perc)
  table<-tableall[-per,]
  return(table)
}



Closeori<-function(pos, bin){
  
  x<-c(1:length(pos)) 
  remove<-c()
  replace<-c()
  for(i in 1:(length(pos)-1)){
    if (abs(pos[i]-pos[i+1])==bin)  
    {remove<-c(remove,x[i])
     next
    }else if
    (abs(pos[i]-pos[i+1])==2*bin) 
    {remove<-c(remove,x[i],x[i+1])
     replace<-c(replace,(pos[x[i]]+bin))
     next
    }else if
    (abs(pos[i]-pos[i+1])==3*bin) 
    {remove<-c(remove,x[i],x[i+1])
     replace<-c(replace,(pos[x[i]]+bin))
     next
    }else if
    (abs(pos[i]-pos[i+1])==4*bin) 
    {remove<-c(remove,x[i],x[i+1])
     replace<-c(replace,(pos[x[i]]+2*bin))
    }
    
  }
  if(length(remove>0)){
  posremove<-c(pos[-remove])}else
  posremove<-pos
  posreplace<-sort(unique(c(posremove,replace)))
  return(posreplace)
  
  
}


#tableclose=peak positions from diff, ratio=pol usage ratio, pos=all position)
orieff<-function(close, ratio, pos){
  
  maxpos<-c(match(close,pos))
  
  ratiomin<-c(tail((which(diff(c(FALSE,diff(ratio[1:maxpos[1]])>0,TRUE))>0)),n=1))
  
  for (i in 1:(length(maxpos)-1)){
    
    ratiomin<-c(ratiomin,(-1+maxpos[i]+tail((which(diff(c(FALSE,diff(ratio[maxpos[i]:maxpos[i+1]])>0,TRUE))>0)),n=1)))
  }
  
  ratiomax<-c()
  
  for (i in 1:(length(maxpos)-1)){
    ratiomax<-c(ratiomax,(-1+maxpos[i]+head((which(diff(c(TRUE,diff(ratio[maxpos[i]:maxpos[i+1]])>=0,FALSE))<0)),n=1)))
  }
  ratiomax<-c(ratiomax, (-1+maxpos[length(maxpos)]+head((which(diff(c(TRUE,diff(ratio[maxpos[length(maxpos)]:length(ratio)])>=0,FALSE))<0)),n=1)))
  #which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
  
  orieff<-c((ratio[ratiomax]-ratio[ratiomin])*100)
  oriefftable<-cbind(close,orieff)
  
  return(oriefftable)
}


orieff_merge<-function(orieff_ef, orieff_ef_pos, orieff_dr, orieff_dr_pos,chro,bin){
  
  value<-c()
  valuepos<-c()
  drpaired<-c()
  efpaired<-c()
  for (i in 1:length(orieff_ef_pos)){
    if 
    (orieff_ef_pos[i] %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[match(orieff_ef_pos[i],orieff_dr_pos)])))
      valuepos<-c(valuepos, orieff_ef_pos[i])
      drpaired<-c(drpaired, orieff_dr_pos[match(orieff_ef_pos[i],orieff_dr_pos)] )
      efpaired<-c(efpaired, orieff_ef_pos[i])
      next
    }else if
    ((orieff_ef_pos[i]+bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]+bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, orieff_ef_pos[i])
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]+bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])   
      next
      
    }else if
    ((orieff_ef_pos[i]-bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]-bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, orieff_ef_pos[i])
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]-bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])   
      next
    }else if
    ((orieff_ef_pos[i]+2*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]+2*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]+bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]+2*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i]) 
      next
    }else if
    ((orieff_ef_pos[i]-2*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]-2*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]-bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]-2*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i]) 
      next
    }else if
    ((orieff_ef_pos[i]+3*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]+3*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]+2*bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]+3*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])
      next
      
    }else if
    ((orieff_ef_pos[i]-3*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]-3*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]-2*bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]-3*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])
      next
    }else if
    ((orieff_ef_pos[i]+4*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]+4*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]+2*bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]+4*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])
      next
      
    }else if
    ((orieff_ef_pos[i]-4*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]-4*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]-2*bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]-4*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])
      
    }
  }
  
  efpairedno<-match(efpaired,orieff_ef_pos)
  drpairedno<-match(drpaired,orieff_dr_pos)
  efunpaired<-orieff_ef_pos[-efpairedno]
  drunpaired<-orieff_dr_pos[-drpairedno]
  
  if((length(efunpaired) != 0) & (length(drunpaired) != 0)){
    
    
    
    for (i in 1:length(drunpaired)){ 
      if 
      (drunpaired[i] %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match(drunpaired[i],orieff_ef_pos)])))
        valuepos<-c(valuepos, drunpaired[i])
        next
      }else if
      ((drunpaired[i]+bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]+bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, drunpaired[i])
        next
        
      }else if
      ((drunpaired[i]-bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]-bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, drunpaired[i])
        next
        
      }else if
      ((drunpaired[i]+2*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]+2*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]+bin))
        next
        
      }else if
      ((drunpaired[i]-2*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]-2*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]-bin))
        next
        
      }else if
      ((drunpaired[i]+3*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]+3*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]+2*bin))
        next
        
      }else if
      ((drunpaired[i]-3*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]-3*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]-2*bin))
        next
      }else if
      ((drunpaired[i]+4*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]+4*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]+2*bin))
        next
        
      }else if
      ((drunpaired[i]-4*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]-4*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]-2*bin))
        
      }} 
      chromosome<-rep(chro,length(valuepos))
      orilist<-as.data.frame(cbind(chromosome,valuepos,value))
      colnames(orilist)<-c("chromosome","maxpos","efficiency")
    
    #### Remove duplicates and zero values ####
    
    orilist_dupl<-which(duplicated(orilist$maxpos))
    if (length(orilist_dupl)>0){
      a<-c()
      for (i in 1:length(orilist_dupl)){
        a<-c(a,(as.numeric(as.vector(orilist$efficiency[orilist_dupl[i]]))+as.numeric(as.vector(orilist$efficiency[orilist_dupl[i]-1])))/2)
        levels(orilist$efficiency)<-c(levels(orilist$efficiency),a)}
      orilist$efficiency[orilist_dupl-1]<-a
      orilist<-orilist[-orilist_dupl,]
      orilist<-orilist[orilist$efficiency !=0,]}
      
    
  }else if(length(efunpaired)==0 | length(drunpaired)==0){
    
    chromosome<-rep(chro,length(valuepos))
    orilist<-as.data.frame(cbind(chromosome,valuepos,value))
    colnames(orilist)<-c("chromosome","maxpos","efficiency")
  
    #### Remove duplicates and zero values ####
  
  orilist_dupl<-which(duplicated(orilist$maxpos))
  if (length(orilist_dupl)>0){
    a<-c()
    for (i in 1:length(orilist_dupl)){
      a<-c(a,(as.numeric(as.vector(orilist$efficiency[orilist_dupl[i]]))+as.numeric(as.vector(orilist$efficiency[orilist_dupl[i]-1])))/2)
      levels(orilist$efficiency)<-c(levels(orilist$efficiency),a)}
      orilist$efficiency[orilist_dupl-1]<-a
      orilist<-orilist[-orilist_dupl,]
      orilist<-orilist[orilist$efficiency !=0,]}}
  
  

  return(orilist) 
}



Wigdata <- function(name, ratio.table, row,  bin, color, h.line,chro,chromoname_in){

  
 header0 = paste('track type=wiggle_0 name="', name,
                  '" description="downloaded from Puseq app Y. Daigaku & A. Keszthelyi, 2015', date(), 
                  '" visibility=full autoScale=off color=', color,
                  ' yLineOnOff=on yLineMark=', h.line, 
                  ' priority=10', sep="")
  wiglist<- list()
  wiglist<-c(wiglist, header0)
  
  for(chromo in chro[1:length(chro)]){ 

    
  header = paste('fixedStep chrom=',paste(chromoname_in,gsub("[^0-9]", "",chromo),sep=""),
                  ' step=', bin, ' span=', bin, sep="")
  
  ratio.table.chr = ratio.table[ratio.table$chromosome==chromo,]
  data.chr<-as.numeric(as.vector(ratio.table.chr[,row]))
  data.chr[is.na(data.chr)]<- 0
  
  wiglist<-c(wiglist,header,list(data.chr))}
  
  return(wiglist)
  
}
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


})
