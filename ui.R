#ui
#Puseq_app
library(shiny)


shinyUI(navbarPage("Puseq data analysis",
  

  tabPanel(("File upload"),
           
           

      sidebarPanel(
        
        
        fileInput('file_df', 'Choose CSV File for delta forward',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
        tags$hr(),
        
        fileInput('file_dr', 'Choose CSV File for delta reverse',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
        tags$hr(),
        
      
       fileInput('file_ef', 'Choose CSV File for epsilon forward',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
       tags$hr(),
       
       
       fileInput('file_er', 'Choose CSV File for epsilon reverse',
                 accept=c('text/csv', 
                          'text/comma-separated-values,text/plain', 
                          '.csv')),
       tags$hr()),
 mainPanel(
 
 tags$blockquote("Please upload a csv file with the counts for each datasets e.g.: delta usage on forward strand, delta usage on reverse strand, epsilon usage on forward strain and epsilon usage on reverse strand"),
       
 tags$blockquote("Sample table for bin size 300:"),
 

 
 tableOutput("sample"))),
 
 
 
 
 
 
tabPanel(("Counts"),
        
         plotOutput("countPlot"),
         hr(),
         
         fluidRow(
           
           column(3,
                  selectInput('chromo', label="Chromosome",""         
                              
                             ) ),
          
                  
        
       
         
      
      column
      (4, offset = 1,
        
        selectInput("input_type", "Input type for x axis start",
                    c( "slider","numeric"
                         )
             ),
             uiOutput("ui"),
             
             sliderInput("xrange",
                  "x axis range:",
                  min = 0,
                  max = 1000000,
                  value = 500000,
                  step = 1000)),
      
      column(4, 
       uiOutput("maxy")
    ))),



tabPanel(("Ratios"),
        
         plotOutput("ratioPlot"),
         hr(),
         
         fluidRow(
           
           column(2,
                  selectInput('chromo_2', label = "Chromosome", "")),
           column(3,offset = 1,
                  numericInput("N.ratio.num", "moving average for ratios", value = 3 ),
                  numericInput("N.ori.num", "moving average for origins", value = 3 ),
                  numericInput("P.ori.num", "percentile treshold for origins", value = 0.3 )),
           column(3, 
              
                selectInput("input_type_2", "Input type for x axis start",
                            c( "slider","numeric")),
                uiOutput("ui_2"),
                
                  sliderInput("xrange_2",
                              "x axis range:",
                              min = 0,
                              max = 1000000,
                              value = 500000,
                             step = 1000)),
         
          column(3, 
                 textInput("name", label = ("Enter name to generate file names"), value = "filename"),
                 uiOutput("chromoname"),
                 
                 selectInput("dataset", "Choose a dataset for csv/wig/bedgraph files:", 
                             choices = c("all ratios as csv",
                                         "delta forward ratio as wig",
                                         "delta reverse ratio as wig",
                                         "epsilon forward ratio as wig",
                                         "epsilon reverse ratio as wig",
                                         "origin efficiency as csv",
                                         "origin efficiency as bedgraph")),
                 downloadButton('downloadData', 'Download'))))
           
          
         
         
)         


)
  

