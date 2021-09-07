## Preliminaries loading data
library(ggpubr)
library(ggrepel)
library(shiny)
library(usmap)
source(file.path("Functions","setup_data.R"))

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("COVID-19 Prevalence and Seroprevalence in U.S. States"),
    helpText("Semi-empirical model utilizing test positivity and reported cases, calibrated using published seroprevalence data."),
    hr(),
    fluidRow(
        column(1),
        column(4,
               sliderInput("daterange",
                           HTML("Select from available Date Range<br/>[Final date in selected range is mapped]:"),
                           min = min(predquant.RE$date),
                           max = max(predquant.RE$date),
                           timeFormat = "%b %d\n%Y",# "%F",
                           value = range(predquant.RE$date))
        ),
        column(1),
        column(4,
               selectizeInput("stateMap", 
                              HTML("Select state(s)<br/>[First state in list is graphed]:"),
                              choices = c("All 50 + DC",statesvec),  
                              selected="All 50 + DC", multiple =TRUE))
    ),
    fluidRow(column(1),
             column(9,
                    radioButtons("Model","Select model:",
                                 choiceNames=list("Random Effects (primary, n fit to each state)",
                                              "Geometric Mean (n=0.5)"),
                                 choiceValues=list("RE","GM"),inline=TRUE))),
    hr(),
    tabsetPanel(
        tabPanel("Seroprevalence Map", plotOutput("SeroMapPlot",height = "600px")),
        tabPanel("Seroprevalence Graph", 
                 fluidRow(column(1),
                          column(9,checkboxInput("AddSeroData","Include Seroprevalence Data?  (individual states only)",
                                                 value=FALSE))),
                 plotOutput("SeroGraph",height = "600px")),
        tabPanel("Total Prevalence Map", plotOutput("TotPrevMapPlot",height = "600px")),
        tabPanel("Total Prevalence Graph", 
                 fluidRow(column(1),
                          column(9,checkboxInput("AddEpiPred","Include Epidemiologic Model Predictions?  (individual states only)",
                                                 value=FALSE))),
                 plotOutput("TotPrevGraph",height = "600px")),
        tabPanel("Undiagnosed Prevalence Map", plotOutput("PrevMapPlot",height = "600px")),
        tabPanel("Undiagnosed Prevalence Graph", plotOutput("PrevGraph",height = "600px")),
        tabPanel("Conceptual Model", imageOutput("Concept",height = "800px"),
                 HTML(paste0("<p>Conceptual model for relationship between test positivity, prevalence of infection, and testing rate. <b>A)</b> Compartmental representation of how the relationships between new infections, undiagnosed and diagnosed prevalence (I<sub>U</sub> and I<sub>D</sub>) and seroprevalence (SP<sub>U</sub> and SP<sub>D</sub>) are modeled for each state, given a bias with power n. All observational inputs are the past τ-day averages of number of positive tests N<sub>+,τ</sub>(t) and number of tests performed N<sub>test,τ</sub>(t), the corresponding test positivity rate P<sub>+,τ</sub>(t) and reported case rate C<sub>+,τ</sub>(t), and the state population size N.  For diagnosed prevalence and seroprevalence, the observational input is the daily reported cases N<sub>+,τ</sub>, and the model parameters are the recovery time after diagnosis T<sub>rec</sub> and the time from infection to seropositivity T<sub>inf</sub>. For undiagnosed prevalence and seroprevalence, our model assumes the test positivity rate is correlated to delayed undiagnosed disease prevalence with a bias parameter b(t) modeled as a negative power function of the testing rate b(t) = [N<sub>test,τ</sub>(t)/N]<sup>–n</sup>. The additional parameters consist of the power parameter n and the initial (missed) seroprevalence SP<sub>o</sub>.  The effective rate parameter 1/T<sub>eff</sub> is time-dependent, and accounts for both T<sub>inf</sub> and ongoing diagnoses so as to not “double count.”  Prevalence and seroprevalence are evaluated with a lag time t<sub>lag</sub>, assumed equal to half the averaging time τ/2. In <b>B)</b>, the diagonal lines represent different values of the bias parameter. In <b>C)</b>, the relationship between testing rate and bias parameter is illustrated. Here the shaded region represents different powers n ranging from 0.1 (lower bound bias) to 0.9 (upper bound bias), the solid line represents n=½.</p>"))),
        tabPanel("Download Excel Calculator",downloadButton("downloadData", "Download Excel Calculator Version"),
                 HTML("<p>This Excel spreadsheet reproduces the calculations given user-supplied data, for fixed (user-modifiable) model parameters.  It is pre-loaded with data for Texas from Covid Tracking Project.</p>"))
    ),
    hr(),
    HTML(paste0("<p>Updated with data from <a href='https://covidtracking.com'>The COVID Tracking Project</a> through ",max(alldat$date),
                              ". Prevalence and seroprevalence estimates calculated through ",
                              max(predquant.RE$date),".</p><b>Note:</b> No further updates due to availability of COVID-19 vaccinations in January 2021 and <a href='https://covidtracking.com/analysis-updates/giving-thanks-and-looking-ahead-our-data-collection-work-is-done'>ending of COVID Tracking Project data collection in March 2021</a>.")),
    hr(),
    HTML("<p><b>Citation: </b>Chiu WA, Ndeffo-Mbah ML (2021) Using test positivity and reported case rates to estimate state-level COVID-19 prevalence and seroprevalence in the United States. PLoS Comput Biol 17(9): e1009374. doi: <a href='https://doi.org/10.1371/journal.pcbi.1009374'>https://doi.org/10.1371/journal.pcbi.1009374</a>.</p>"),
    HTML("<p><b>Source code and data: </b><a href='https://github.com/wachiuphd/COVID-19-US-Semi-Empirical-published'>https://github.com/wachiuphd/COVID-19-US-Semi-Empirical-published</a></p>")
    
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$PrevMapPlot <- renderPlot({
        if (input$Model == "RE") {
            uspredquant <- uspredquant.RE
            predquant <- predquant.RE
        } else {
            uspredquant <- uspredquant.n0.5
            predquant <- predquant.n0.5
        }
        currentdate <-input$daterange[2]
        
        uscurrentpred <- subset(uspredquant,date==currentdate)
        uspastpred <- subset(uspredquant,date==(currentdate-14))
        
        currentpred <- subset(predquant,date==currentdate)
        currentpred <- currentpred[order(currentpred$state),]
        
        pastpred <- subset(predquant,date==(currentdate-14))
        pastpred <- pastpred[order(pastpred$state),]
        
        currentpred$x <- statexy[currentpred$state,"x"]
        currentpred$y <- statexy[currentpred$state,"y"]
        pastpred$x <- statexy[pastpred$state,"x"]
        pastpred$y <- statexy[pastpred$state,"y"]
        
        uscurrentpred$var <- uscurrentpred$I.pct.50.
        uscurrentpred$var.97.5. <- uscurrentpred$I.pct.97.5.
        uscurrentpred$var.2.5. <- uscurrentpred$I.pct.2.5.
        
        uspastpred$var <- uspastpred$I.pct.50.
        currentpred$var <- currentpred$I.pct.50.
        pastpred$var <- pastpred$I.pct.50.
        
        uscurrentpred$trend <- 100*(uscurrentpred$var/uspastpred$var-1)
        uscurrentpred$trend[uscurrentpred$trend>100]<-100
        uscurrentpred$dir<-ifelse(uscurrentpred$trend>0,"Increasing","Decreasing")
        uscurrentpred$dir <- factor(uscurrentpred$dir,levels=c("Increasing","Decreasing"))
        
        currentpred$trend <- 100*(currentpred$var/pastpred$var-1)
        currentpred$trend[currentpred$trend>100]<-100
        currentpred$dir<-ifelse(currentpred$trend>0,"Increasing","Decreasing")
        currentpred$dir <- factor(currentpred$dir,levels=c("Increasing","Decreasing"))
        
        if("All 50 + DC" %in% input$stateMap | length(input$stateMap)==0) {
            mapstates<-statesvec
        } else {
            mapstates<-input$stateMap
        }

        currentpred_states<-subset(currentpred, state %in% mapstates)
        prevmap<-plot_usmap(include=mapstates)+
            geom_point(data=currentpred_states,
                       aes(x = x, y = y, 
                           size = var,
                           shape=dir,
                           color=trend,
                           stroke=2
                       ))+
            geom_label_repel(data=currentpred_states,
                             aes(x=x,y=y,label=paste0(state,"\n",signif(var,2),"%")),
                             point.padding=NA,force=0.1,
                             hjust=0.5,vjust=0.5,size=3,label.size=0.1,alpha=0.8,seed=314)+
            geom_label_repel(data=currentpred_states,
                             aes(x=x,y=y,label=paste0(state,"\n",signif(var,2),"%")),
                             point.padding=NA,force=0.1,
                             hjust=0.5,vjust=0.5,size=3,label.size=0.1,fill=NA,seed=314)+
            # geom_text(data=currentpred_states,
            #           aes(x=x,y=y,label=paste0(state,"\n",signif(var,2),"%")),
            #           hjust=0.5,vjust=0.5,size=3)+
            scale_size(range=c(0.1,24),limits = c(0,ceiling(max(currentpred$var))),
                       name="Semi-empirical\nPrevalence %")+
            scale_color_viridis_c(alpha=0.8,option="plasma",direction=-1,
                                  limits=c(-100,100),
                                  breaks=c(-100,-50,0,50,100),
                                  labels=c("-100","-50","0","+50","+\u2265 100"),
                                  name="2-wk Change (%)")+
            scale_shape_manual(values=c(16,1),name="2-wk Trend",drop=FALSE)+
            guides(size = guide_legend(order = 1),
                   shape = guide_legend(order = 2,
                                        override.aes = list(size = 12))
                   #color = guide_legend(order = 2),
            )+
            ggtitle(paste0("Undiagnosed Infection Prevalence as of ",currentdate,
                           "\nUS overall: ",round(uscurrentpred$var,2),
                           "% [95% CI: ",round(uscurrentpred$var.2.5.,2),
                           "%-",round(uscurrentpred$var.97.5.,2),
                           "%] (2-wk trend: ",uscurrentpred$dir,")"))+
            theme(legend.position = "right",
                  legend.box="vertical",
                  plot.title = element_text(size = 14,face="bold"))
        prevmap
    }
    )
    output$TotPrevMapPlot <- renderPlot({
        if (input$Model == "RE") {
            uspredquant <- uspredquant.RE
            predquant <- predquant.RE
        } else {
            uspredquant <- uspredquant.n0.5
            predquant <- predquant.n0.5
        }
        currentdate <-input$daterange[2]
        
        uscurrentpred <- subset(uspredquant,date==currentdate)
        uspastpred <- subset(uspredquant,date==(currentdate-14))
        
        currentpred <- subset(predquant,date==currentdate)
        currentpred <- currentpred[order(currentpred$state),]
        
        pastpred <- subset(predquant,date==(currentdate-14))
        pastpred <- pastpred[order(pastpred$state),]
        
        currentpred$x <- statexy[currentpred$state,"x"]
        currentpred$y <- statexy[currentpred$state,"y"]
        pastpred$x <- statexy[pastpred$state,"x"]
        pastpred$y <- statexy[pastpred$state,"y"]
        
        uscurrentpred$var <- uscurrentpred$Itot.pct.50.
        uscurrentpred$var.97.5. <- uscurrentpred$Itot.pct.97.5.
        uscurrentpred$var.2.5. <- uscurrentpred$Itot.pct.2.5.
        
        uspastpred$var <- uspastpred$Itot.pct.50.
        currentpred$var <- currentpred$Itot.pct.50.
        pastpred$var <- pastpred$Itot.pct.50.
        
        uscurrentpred$trend <- 100*(uscurrentpred$var/uspastpred$var-1)
        uscurrentpred$trend[uscurrentpred$trend>100]<-100
        uscurrentpred$dir<-ifelse(uscurrentpred$trend>0,"Increasing","Decreasing")
        uscurrentpred$dir <- factor(uscurrentpred$dir,levels=c("Increasing","Decreasing"))
        
        currentpred$trend <- 100*(currentpred$var/pastpred$var-1)
        currentpred$trend[currentpred$trend>100]<-100
        currentpred$dir<-ifelse(currentpred$trend>0,"Increasing","Decreasing")
        currentpred$dir <- factor(currentpred$dir,levels=c("Increasing","Decreasing"))
        
        if("All 50 + DC" %in% input$stateMap | length(input$stateMap)==0) {
            mapstates<-statesvec
        } else {
            mapstates<-input$stateMap
        }
        currentpred_states<-subset(currentpred, state %in% mapstates)
        totprevmap<-plot_usmap(include=mapstates)+
            geom_point(data=currentpred_states,
                       aes(x = x, y = y, 
                           size = var,
                           shape=dir,
                           color=trend,
                           stroke=2
                       ))+
            geom_label_repel(data=currentpred_states,
                             aes(x=x,y=y,label=paste0(state,"\n",signif(var,2),"%")),
                             point.padding=NA,force=0.1,
                             hjust=0.5,vjust=0.5,size=3,label.size=0.1,alpha=0.8,seed=314)+
            geom_label_repel(data=currentpred_states,
                             aes(x=x,y=y,label=paste0(state,"\n",signif(var,2),"%")),
                             point.padding=NA,force=0.1,
                             hjust=0.5,vjust=0.5,size=3,label.size=0.1,fill=NA,seed=314)+
            # geom_text(data=currentpred_states,
            #           aes(x=x,y=y,label=paste0(state,"\n",signif(var,2),"%")),
            #           hjust=0.5,vjust=0.5,size=3)+
            scale_size(range=c(0.1,24),limits = c(0,ceiling(max(currentpred$var))),
                       name="Semi-empirical\nPrevalence %")+
            scale_color_viridis_c(alpha=0.8,option="plasma",direction=-1,
                                  limits=c(-100,100),
                                  breaks=c(-100,-50,0,50,100),
                                  labels=c("-100","-50","0","+50","+\u2265 100"),
                                  name="2-wk Change (%)")+
            scale_shape_manual(values=c(16,1),name="2-wk Trend",drop=FALSE)+
            guides(size = guide_legend(order = 1),
                   shape = guide_legend(order = 2,
                                        override.aes = list(size = 12))
                   #color = guide_legend(order = 2),
            )+
            ggtitle(paste0("Total Infection Prevalence as of ",currentdate,
                           "\nUS overall: ",round(uscurrentpred$var,2),
                           "% [95% CI: ",round(uscurrentpred$var.2.5.,2),
                           "%-",round(uscurrentpred$var.97.5.,2),
                           "%] (2-wk trend: ",uscurrentpred$dir,")"))+
            theme(legend.position = "right",
                  legend.box="vertical",
                  plot.title = element_text(size = 14,face="bold"))
        totprevmap
    }
    )
    output$SeroMapPlot <- renderPlot({
        if (input$Model == "RE") {
            uspredquant <- uspredquant.RE
            predquant <- predquant.RE
        } else {
            uspredquant <- uspredquant.n0.5
            predquant <- predquant.n0.5
        }
        currentdate <-input$daterange[2]
        
        uscurrentpred <- subset(uspredquant,date==currentdate)
        uspastpred <- subset(uspredquant,date==(currentdate-14))
        
        currentpred <- subset(predquant,date==currentdate)
        currentpred <- currentpred[order(currentpred$state),]
        
        pastpred <- subset(predquant,date==(currentdate-14))
        pastpred <- pastpred[order(pastpred$state),]
        
        currentpred$x <- statexy[currentpred$state,"x"]
        currentpred$y <- statexy[currentpred$state,"y"]
        pastpred$x <- statexy[pastpred$state,"x"]
        pastpred$y <- statexy[pastpred$state,"y"]
        
        # uscurrentpred$var <- uscurrentpred$I.pct.50.
        # uscurrentpred$var.97.5. <- uscurrentpred$I.pct.97.5.
        # uscurrentpred$var.2.5. <- uscurrentpred$I.pct.2.5.
        # 
        # uspastpred$var <- uspastpred$I.pct.50.
        # currentpred$var <- currentpred$I.pct.50.
        # pastpred$var <- pastpred$I.pct.50.
        # 
        # uscurrentpred$trend <- 100*(uscurrentpred$var/uspastpred$var-1)
        # uscurrentpred$trend[uscurrentpred$trend>100]<-100
        # uscurrentpred$dir<-ifelse(uscurrentpred$trend>0,"Increasing","Decreasing")
        # uscurrentpred$dir <- factor(uscurrentpred$dir,levels=c("Increasing","Decreasing"))
        # 
        # currentpred$trend <- 100*(currentpred$var/pastpred$var-1)
        # currentpred$trend[currentpred$trend>100]<-100
        # currentpred$dir<-ifelse(currentpred$trend>0,"Increasing","Decreasing")
        # currentpred$dir <- factor(currentpred$dir,levels=c("Increasing","Decreasing"))
        
        if("All 50 + DC" %in% input$stateMap | length(input$stateMap)==0) {
            mapstates<-statesvec
        } else {
            mapstates<-input$stateMap
        }
        
        currentpred_states<-subset(currentpred, state %in% mapstates)
        seromap<-plot_usmap(include=mapstates)+
            geom_point(data=currentpred_states,
                       aes(x = x, y = y, 
                           size = SP.pct.50.,
                           stroke=1,color=SP.pct.50.
                       ),shape=16)+
            geom_label_repel(data=currentpred_states,
                             aes(x=x,y=y,label=paste0(state,"\n",signif(SP.pct.50.,2),"%")),
                             point.padding=NA,force=0.1,
                             hjust=0.5,vjust=0.5,size=3,label.size=0,alpha=0.2,seed=314)+
            geom_label_repel(data=currentpred_states,
                             aes(x=x,y=y,label=paste0(state,"\n",signif(SP.pct.50.,2),"%")),
                             point.padding=NA,force=0.1,
                             hjust=0.5,vjust=0.5,size=3,label.size=0,fill=NA,seed=314)+
            scale_size(range=c(0.1,24),limits = c(0,ceiling(max(currentpred$SP.pct.50.))),
                       name="Semi-empirical\nSeroprevalence %")+
            scale_color_viridis_c(alpha=0.8,option="plasma",begin=1,end=0,name="",limits=c(0,100))+
            guides(size = guide_legend(order = 1),
                   shape = guide_legend(order = 2,
                                        override.aes = list(size = 12))
                   #color = guide_legend(order = 2),
            )+
            ggtitle(paste0("Seroprevalence as of ",currentdate,
                           "\nUS overall: ",round(uscurrentpred$SP.pct.50.,1),
                           "% [95% CI: ",round(uscurrentpred$SP.pct.2.5.,1),
                           "%-",round(uscurrentpred$SP.pct.97.5.,1),
                           "%]"))+
            theme(legend.position = "right",
                  legend.box="vertical",
                  plot.title = element_text(size = 14,face="bold"))
        seromap
    }
    )
    output$PrevGraph <- renderPlot({
        mindate <- input$daterange[1]
        maxdate <- input$daterange[2]
        if (input$Model == "RE") {
            uspredquant <- uspredquant.RE
            predquant <- predquant.RE
        } else {
            uspredquant <- uspredquant.n0.5
            predquant <- predquant.n0.5
        }
        if (length(input$stateMap)==0) {
            statenow <- "US"
            onepred <- subset(uspredquant,date >= mindate & date <= maxdate)
        } else if (input$stateMap[1]=="All 50 + DC") {
            statenow <- "US"
            onepred <- subset(uspredquant,date >= mindate & date <= maxdate)
        } else {
            statenow <- input$stateMap[1]
            onepred <- subset(predquant,state==statenow & date >= mindate & date <= maxdate)
        }
        onepredlast<-tail(onepred,1)
        predtext <- paste0(onepredlast$date,": ",
                           signif(onepredlast$I.pct.50.,3),"% [",
                           signif(onepredlast$I.pct.2.5.,3),"% -",
                           signif(onepredlast$I.pct.97.5.,3),"%]")
        pprev<-ggplot()+
            geom_col(data=onepred,aes(x=date,y=ID.pct.50.,
                                     fill="Diagnosed Cases (smoothed) in\nprevious 10 days"))+
            geom_ribbon(data=onepred,aes(x=date,ymin=I.pct.2.5.,
                                           ymax=I.pct.97.5.,fill=" Posterior 95% CrI"),
                        alpha=0.5)+
            geom_line(data=onepred,aes(x=date,y=I.pct.50.,
                                         color="Posterior median"))+
            scale_fill_viridis_d(begin=0.5,name="")+
            ylim(0,NA)+
            scale_x_date(date_minor_breaks="1 month",date_breaks = "2 months",
                         date_labels = "%b")+  
            theme_bw()+theme(legend.position="bottom")+
            scale_color_viridis_d(begin=0.3,name="Semi-empirical Model")+
            guides(color = guide_legend(order=2),
                   fill = guide_legend(order=3))+
            ggtitle(paste0(statenow," Undiagnosed Infection Prevalence Estimates\n",
                           predtext))+
            xlab("Date")+
            ylab("Undiagnosed Prevalence %")
        pprev
    })
    output$TotPrevGraph <- renderPlot({
        mindate <- input$daterange[1]
        maxdate <- input$daterange[2]
        if (input$Model == "RE") {
            uspredquant <- uspredquant.RE
            predquant <- predquant.RE
        } else {
            uspredquant <- uspredquant.n0.5
            predquant <- predquant.n0.5
        }
        if (length(input$stateMap)==0) {
            statenow <- "US"
            onepred <- subset(uspredquant,date >= mindate & date <= maxdate)
            onestate<-FALSE
        } else if (input$stateMap[1]=="All 50 + DC") {
            statenow <- "US"
            onepred <- subset(uspredquant,date >= mindate & date <= maxdate)
            onestate<-FALSE
        } else {
            statenow <- input$stateMap[1]
            onepred <- subset(predquant,state==statenow & date >= mindate & date <= maxdate)
            onestate<-TRUE
            oneseirpred <- subset(seirpred,state==statenow & date >= mindate & date <= maxdate)
            oneimperialpred <- subset(imperialpred,state==statenow & date >= mindate & date <= maxdate)
        }
        onepredlast<-tail(onepred,1)
        predtext <- paste0(onepredlast$date,": ",
                           signif(onepredlast$Itot.pct.50.,3),"% [",
                           signif(onepredlast$Itot.pct.2.5.,3),"% -",
                           signif(onepredlast$Itot.pct.97.5.,3),"%]")
        ptotprev<-ggplot()+
            geom_col(data=onepred,aes(x=date,y=ID.pct.50.,
                                      fill="Diagnosed Cases (smoothed) in\nprevious 10 days"))+
            geom_ribbon(data=onepred,aes(x=date,ymin=Itot.pct.2.5.,
                                         ymax=Itot.pct.97.5.,fill="  Posterior 95% CrI"),
                        alpha=0.5)+
            geom_line(data=onepred,aes(x=date,y=Itot.pct.50.,
                                       color="Posterior median"))+
            scale_fill_viridis_d(begin=0.5,name="")+
            ylim(0,NA)+
            scale_x_date(date_minor_breaks="1 month",date_breaks = "2 months",
                         date_labels = "%b")+  
            theme_bw()+theme(legend.position="bottom")+
            scale_color_viridis_d(begin=0.3,name="Semi-empirical Model")+
            ggtitle(paste0(statenow," Total Infection Prevalence Estimates\n",
                           predtext))+
            xlab("Date")+
            ylab("Total Prevalence %")
        if (onestate & input$AddEpiPred) {
            ptotprev <- ptotprev + 
                geom_ribbon(data=oneseirpred,
                            aes(x=date,ymin=N_I_tot.2.5./1000,ymax=N_I_tot.97.5./1000,fill=" Extended SEIR CrI"),
                            alpha=0.25)+
                geom_line(data=oneseirpred,
                          aes(x=date,y=N_I_tot.50./1000,linetype="Extended SEIR median"))+
                geom_ribbon(data=oneimperialpred,
                            aes(x=date,ymin=N_I_tot.2.5./1000,ymax=N_I_tot.97.5./1000,fill=" Imperial CrI"),
                            alpha=0.25)+
                geom_line(data=oneimperialpred,
                          aes(x=date,y=N_I_tot.50./1000,linetype="Imperial median"))+
                scale_linetype_manual(name="Epidemiologic Models",values=c(3,2,1))+
                guides(linetype = guide_legend(order=1,nrow=2),
                       color = guide_legend(order=2),
                       fill = guide_legend(order=3,nrow=2))
        } else {
            ptotprev <- ptotprev + 
                guides(color = guide_legend(order=1),
                       fill = guide_legend(order=2))
        }
            
        ptotprev
    })
    output$SeroGraph <- renderPlot({
        mindate <- input$daterange[1]
        maxdate <- input$daterange[2]
        if (input$Model == "RE") {
            uspredquant <- uspredquant.RE
            predquant <- predquant.RE
        } else {
            uspredquant <- uspredquant.n0.5
            predquant <- predquant.n0.5
        }
        if (length(input$stateMap)==0) {
            statenow <- "US"
            onepred <- subset(uspredquant,date >= mindate & date <= maxdate)
            onedat <- subset(usdat,date >= (mindate-14) & date <= (maxdate-14))
            onestate<-FALSE
        } else if (input$stateMap[1]=="All 50 + DC") {
            statenow <- "US"
            onepred <- subset(uspredquant,date >= mindate & date <= maxdate)
            onedat <- subset(usdat,date >= (mindate-14) & date <= (maxdate-14))
            onestate<-FALSE
        } else {
            statenow <- input$stateMap[1]
            onepred <- subset(predquant,state==statenow & date >= mindate & date <= maxdate)
            onedat <- subset(alldat,state==statenow & date >= (mindate-14) & date <= (maxdate-14))
            oneserodat <- subset(serodat,state==statenow & date >= input$daterange[1] & date <= input$daterange[2])
            onecdcseroval <- subset(cdcseroval,state==statenow & date >= input$daterange[1] & date <= input$daterange[2])
            onestate<-TRUE
        }
        onepredlast<-tail(onepred,1)
        predtext <- paste0(onepredlast$date,": ",
                           signif(onepredlast$SP.pct.50.,3),"% [",
                           signif(onepredlast$SP.pct.2.5.,3),"% -",
                           signif(onepredlast$SP.pct.97.5.,3),"%]")
        psero<-ggplot()+
            geom_col(data=onedat,aes(x=(date+14),y=posPct,
                                     fill="Cumulative Reported Cases\n2 weeks Prior"))+
            geom_ribbon(data=onepred,aes(x=date,ymin=SP.pct.2.5.,
                                         ymax=SP.pct.97.5.,fill=" Posterior 95% CrI"),
                        alpha=0.5)+
            geom_line(data=onepred,aes(x=date,y=SP.pct.50.,
                                       color="Posterior median"))+
            scale_fill_viridis_d(begin=0.5,name="")+
            ylim(0,NA)+
            scale_x_date(date_minor_breaks="1 month",date_breaks = "2 months",
                         date_labels = "%b")+  
            theme_bw()+theme(legend.position="bottom")+
            scale_color_viridis_d(begin=0.3,name="Semi-empirical Model")+
            ggtitle(paste0(statenow," Serorevalence Estimates\n",
                           predtext))+
            xlab("Date")+
            ylab("Seroprevalence %")
        if (onestate & input$AddSeroData) {
            psero <- psero+geom_point(data=oneserodat,aes(x=date,y=Data.50.,shape=Source))+
                geom_errorbar(data=oneserodat,aes(x=date,ymin=Data.2.5.,ymax=Data.97.5.),width=0)+
                geom_point(data=onecdcseroval,aes(x=date,y=Data.50.,shape=Source))+
                geom_errorbar(data=onecdcseroval,aes(x=date,ymin=Data.2.5.,ymax=Data.97.5.),width=0)+
                scale_shape_manual(name="Seroprevalence\nData Source",
                                   values=c(18,17,15,16,3,7,1),drop=FALSE)+
                guides(shape = guide_legend(order=1,nrow=4),
                       color = guide_legend(order=2),
                       fill = guide_legend(order=3,nrow=2))
        } else {
            psero <- psero+
                guides(color = guide_legend(order=1),
                       fill = guide_legend(order=2,nrow=2))
        }
        psero
    })
    output$Concept <- renderImage({
        filename <- "ConceptualModel.jpeg"
        list(src = filename,
             contentType="image/jpeg",
             height=800,
             alt = "Conceptual Model")
        
    }, deleteFile = FALSE)
    output$downloadData <- downloadHandler(
        filename = "Semi-Empirical-Model-Spreadsheet-Version.xlsx",
        content = function(con) {
            file.copy("Semi-Empirical-Model-Spreadsheet-Version.xlsx",con)
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
