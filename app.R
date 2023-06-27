library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(ggplot2)
library(plotly)
library(gridExtra)
library(lubridate)
library(reshape2)
library(ggpubr)
library(png)
library(stringr)
library(ggiraph)
library(tidyverse)
library(DT)
library(vcd)

options(spinner.color = "grey", spinner.color.background = "#ffffff", spinner.size = 2, shiny.reactlog=TRUE)

ui <- dashboardPage(skin = "black",
                    dashboardHeader(disable = T),
                    dashboardSidebar(collapsed = TRUE,disable = T),
                    dashboardBody(
                      h2("Vigilância Molecular Baseada em Descarte de Teste Rápido para Covid-19 - UBS-2, Ceilândia"),
                      hr(),
                      column(width = 4,
                       dateRangeInput("date","Data",start = as.Date("2020-01-01"),end = Sys.Date(),format = "dd-mm-yyyy") 
                      ),
                      column(width = 4,
                        pickerInput("virus","Vírus",choices = c("COVID-19 (TR)"="covid_rap",
                                                                "COVID-19 (PCR)"="covid_pcr",
                                                                "fluA"="flu_a","fluB"="flu_b",
                                                                "VSR"="vsr","hPVI"="hpvi","AdV"="adv",
                                                                "HEV"="hev","hMPV"="mpv","α-coronavírus"="acorona",
                                                                "β-coronavírus"="oc43","Não detectado"="ndetect"),
                                    options = list(`actions-box` = TRUE,`deselect-all-text` = "Desmarcar todos",
                                                   `select-all-text` = "Marcar todos",size = 10,
                                                   `selected-text-format` = "count",`count-selected-text` = "{0}/{1} vírus selecionado(s)"),multiple = T,selected = c("covid_rap","covid_pcr","flu_a","flu_b",
                                                                                                                                                                    "vsr","hpvi","adv","hev","mpv","acorona",
                                                                                                                                                                    "oc43","ndetect"))
                      ),
                      column(width = 4,
                        pickerInput("idade","Faixa etária",choices = c("<4","4-10","11-20","21-30","31-40","41-70","70+","Ignorado"),
                                    options = list(`actions-box` = TRUE,`deselect-all-text` = "Desmarcar todos",
                                                   `select-all-text` = "Marcar todos",size = 10,
                                                   `selected-text-format` = "count",`count-selected-text` = "{0}/{1} faixa(s) etária(s) selecionada(s)"),multiple = T,selected = c("<4","4-10","11-20","21-30","31-40","41-70","70+","Ignorado"))
                      ),
                      hr(),
                      fluidPage(
                        box(width=12,height = "auto",
                            tabsetPanel(id="tab",
                              tabPanel("COVID-19",
                                h1(htmlOutput("title1")),
                                valueBoxOutput("top1.1", width = 3),
                                valueBoxOutput("top1.2", width = 3),
                                valueBoxOutput("top1.3", width = 3),
                                valueBoxOutput("top1.4", width = 3),
                                hr(),
                                withSpinner(ggiraphOutput("plot1",width = 'auto',height='auto'),type = 2)
                              ),
                              tabPanel("Outros vírus",
                                h1(htmlOutput("title2")),
                                #valueBoxOutput("top2.1", width = 3),
                                valueBoxOutput("top2.2", width = 4),
                                valueBoxOutput("top2.3", width = 4),
                                valueBoxOutput("top2.4", width = 4),
                                valueBoxOutput("top2.5", width = 4),
                                valueBoxOutput("top2.6", width = 4),
                                valueBoxOutput("top2.7", width = 4),
                                valueBoxOutput("top2.8", width = 4),
                                valueBoxOutput("top2.9", width = 4),
                                valueBoxOutput("top2.10", width = 4),
                                hr(),
                                withSpinner(ggiraphOutput("plot2",width = 'auto',height='auto'),type = 2)
                              ),
                              tabPanel("Concordância",
                                h1(htmlOutput("title3")),
                                tags$hr(),
                                dataTableOutput("table", height = "auto", width = "100%")
                              )
                            )
                            )
                          ),
                      hr(),
                      fluidPage(
                        column(width = 3, align = "center",plotOutput("logo1", height = '260', width = '283')),
                        column(width = 3, align = "center",plotOutput("logo2", height = '92', width = '268')),
                        column(width = 3, align = "center",plotOutput("logo3", height = '128', width = '128')),
                        column(width = 3, align = "center",plotOutput("logo4", height = '180', width = '240'))
                      )
                    )
)


server <- function(input, output, session) {
  
  date=data.frame(dateweek=seq.Date(from=as.Date("2019-12-01"),to=Sys.Date(),by="day"))
  date$ew=paste0(epiweek(date$dateweek),".",epiyear(date$dateweek))
  date$semana=weekdays(date$dateweek)
  date=date[which(date$semana=="domingo"),c("dateweek","ew")]
  
  database <- reactive({
    
    data=data.frame(read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRKSmOKILxFUwcgD4sO7TUAltKsfvTBVV816_H53X3yw2ItPEi4RXZIoFrdmQ5xpZeQNzG4Ybk48Sml/pub?gid=0&single=true&output=csv")) %>% 
      janitor::clean_names() %>% 
      mutate(date=as.Date(data_teste_coleta,format="%d/%m/%Y"),
             ew=paste0(epiweek(date),".",epiyear(date)),
             covid_rap=ifelse(str_detect(toupper(result_covid_teste_rapido),"POS"),1,ifelse(str_detect(toupper(result_covid_teste_rapido),"NEG"),0,9)),
             covid_pcr=ifelse(sars_co_v_2!=1 & sars_co_v_2!=0,9,as.numeric(sars_co_v_2)),
             vsr=ifelse(rsv!=0 & rsv!=1,9,as.numeric(rsv)),
             flu_a=ifelse(flu_a!=0 & flu_a!=1,9,as.numeric(flu_a)),
             flu_b=ifelse(flub_b!=0 & flub_b!=1,9,as.numeric(flub_b)),
             hpvi=ifelse(p2_piv_1==1 | p2_piv_2==1 | p2_piv_3==1 | p2_piv_4==1,1,0),
             acorona=ifelse(p3_229e==1 | p3_nl63==1,1,0),
             adv=ifelse(p2_ad_v==1,1,0),
             hev=ifelse(p2_hev==1,1,0),
             mpv=ifelse(p2_mpv==1,1,0),
             oc43=ifelse(p3_oc43==1,1,0),
             covid_rap=ifelse(is.na(covid_rap),9,covid_rap),
             covid_pcr=ifelse(is.na(covid_pcr),9,covid_pcr),
             vsr=ifelse(is.na(vsr),9,vsr),
             flu_a=ifelse(is.na(flu_a),9,flu_a),
             flu_b=ifelse(is.na(flu_b),9,flu_b),
             hpvi=ifelse(is.na(hpvi),9,hpvi),
             acorona=ifelse(is.na(acorona),9,acorona),
             adv=ifelse(is.na(adv),9,adv),
             hev=ifelse(is.na(hev),9,hev),
             mpv=ifelse(is.na(mpv),9,mpv),
             oc43=ifelse(is.na(oc43),9,oc43),
             resultado=ifelse((covid_rap==1 | covid_pcr==1) & vsr!=1 & flu_a!=1 & flu_b!=1 & hpvi!=1 & acorona!=1 & adv!=1 & hev!=1 & mpv!=1 & oc43!=1,1,
                              ifelse((covid_rap==1 | covid_pcr==1) & (vsr==1 | flu_a==1 | flu_b==1 | hpvi==1 | acorona==1 | adv==1 | hev==1 | mpv==1 | oc43==1),2,
                                     ifelse((covid_rap!=1 & covid_pcr!=1) & (vsr==1 | flu_a==1 | flu_b==1 | hpvi==1 | acorona==1 | adv==1 | hev==1 | mpv==1 | oc43==1),3,9)))) %>% 
      mutate(resultado=factor(resultado,levels=c(1:3,9),labels=c("COVID-19","COVID-19 + Outro vírus respiratórios",
                                                                "Outros vírus respiratórios","Não detectado"))) %>% 
      mutate(ndetect=ifelse(resultado=="Não detectado",1,0)) %>% 
      mutate(idade=as.integer(lubridate::time_length(date-as.Date(data_nascimento,format="%d/%m/%Y"),"years"))) %>% 
      mutate(fxetar=case_when(idade>=0 & idade<4 ~ 1,
                              idade>=4 & idade<11 ~ 2,
                              idade>=11 & idade<21 ~ 3,
                              idade>=21 & idade<31 ~ 4,
                              idade>=31 & idade<41 ~ 5,
                              idade>=41 & idade<71 ~ 6,
                              idade>=71 ~ 7,
                              TRUE ~ 9)) %>% 
      mutate(fxetar=factor(fxetar,levels=c(1:7,9),labels=c("<4","4-10","11-20","21-30","31-40","41-70","70+","Ignorado"))) %>% 
      inner_join(date,by="ew") %>%
      mutate(count=as.numeric(+(apply(across(c(input$virus))==1,1,any)))) %>% 
      filter(date>=min(input$date) & date<=max(input$date)) %>% 
      filter(fxetar %in% c(input$idade))
    
  })
  
  observeEvent(input$tab=="Outros vírus", {
    updatePickerInput(
      session,
      "virus",
      label="Vírus",
      choices = c("fluA"="flu_a","fluB"="flu_b","VSR"="vsr","hPVI"="hpvi","AdV"="adv","HEV"="hev",
                  "hMPV"="mpv","α-coronavírus"="acorona","β-coronavírus"="oc43"),
                  options = list(`actions-box` = TRUE,`deselect-all-text` = "Desmarcar todos",
                                 `select-all-text` = "Marcar todos",size = 10,
                                 `selected-text-format` = "count",`count-selected-text` = "{0}/{1} vírus selecionados"),selected = c("flu_a","flu_b","vsr","hpvi","adv",
                                                                                                                                                  "hev","mpv","acorona","oc43"))
  }, ignoreInit = T)
  
  output$title1 <- renderUI({
    db=database()
    db$count=1
    
    htmlText = paste0("<div style='margin-left: 0.2px'!important;>","Vírus detectados"," (n=",sum(db[db$ndetect==0,]$count),")","</div>")
    HTML(htmlText)
  })
  
  output$title2 <- renderUI({
    db=database()
    db$count=1
    
    htmlText = paste0("<div style='margin-left: 0.2px'!important;>","Outros vírus respiratórios"," (n=",sum(db[db$ndetect!=1 & db$covid_pcr!=1 & db$covid_rap!=1,]$count,na.rm = T),")","</div>")
    HTML(htmlText)
  })
  
  output$title3 <- renderUI({
    
    htmlText = paste0("<div style='margin-left: 0.2px'!important;>","Avaliação de concordância dos testes de Covid-19 (PCR vs Teste rápido)","</div>")
    HTML(htmlText)
  })
  
  output$logo1 <- renderPlot({
    ggplot() +
      background_image(readPNG("/home/silvano/Alex/Logo/logo1.png")) +
      theme(panel.background = element_rect(fill = "#ECF0F5", colour = "#ECF0F5"),
            plot.margin = margin(0, 0, 0, 0, "cm"))
  })
  
  output$logo2 <- renderPlot({
    ggplot() +
      background_image(readPNG("/home/silvano/Alex/Logo/logo2.png")) +
      theme(panel.background = element_rect(fill = "#ECF0F5", colour = "#ECF0F5"),
            plot.margin = margin(0, 0, 0, 0, "cm"))
  })
  
  output$logo3 <- renderPlot({
    ggplot() +
      background_image(readPNG("/home/silvano/Alex/Logo/logo3.png")) +
      theme(panel.background = element_rect(fill = "#ECF0F5", colour = "#ECF0F5"),
            plot.margin = margin(0, 0, 0, 0, "cm"))
  })
  
  output$logo4 <- renderPlot({
    ggplot() +
      background_image(readPNG("/home/silvano/Alex/Logo/logo4.png")) +
      theme(panel.background = element_rect(fill = "#ECF0F5", colour = "#ECF0F5"),
            plot.margin = margin(0, 0, 0, 0, "cm"))
  })
  
  output$top1.1 <- renderValueBox({
    db=database()
    db$res=ifelse(db$resultado!="Não detectado",1,0)
    valueBox(
      paste0(prettyNum(sum(db[db$resultado=="COVID-19",]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[db$resultado=="COVID-19",]$count,na.rm = T)/sum(db$res,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "COVID-19",color = "red",icon = icon("virus"))
  })
  
  output$top1.2 <- renderValueBox({
    db=database()
    db$res=ifelse(db$resultado!="Não detectado",1,0)
    valueBox(
      paste0(prettyNum(sum(db[db$resultado=="COVID-19 + Outro vírus respiratórios",]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[db$resultado=="COVID-19 + Outro vírus respiratórios",]$count,na.rm = T)/sum(db$res,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "COVID-19 + Outro vírus respiratórios",color = "orange",icon = icon("virus"))
  })
  
  output$top1.3 <- renderValueBox({
    db=database()
    db$res=ifelse(db$resultado!="Não detectado",1,0)
    valueBox(
      paste0(prettyNum(sum(db[db$resultado=="Outros vírus respiratórios",]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[db$resultado=="Outros vírus respiratórios",]$count,na.rm = T)/sum(db$res,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "Outros vírus respiratórios",color = "purple",icon = icon("virus"))
  })
  
  output$top1.4 <- renderValueBox({
    db=database()
    valueBox(
      prettyNum(sum(db[db$resultado=="Não detectado",]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
      "Nenhum vírus \n detectado",color = "aqua",icon = icon("virus"))
  })
  
  output$top2.1 <- renderValueBox({
    db=database()
    
    valueBox(
      paste0(prettyNum(sum(db[db$covid_pcr==1 | db$covid_rap==1,]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[db$covid_pcr==1 | db$covid_rap==1,]$count,na.rm = T)/sum(db[db$ndetect!=1,]$count,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "COVID-19",color = "light-blue",icon = icon("virus"))
  })
  
  output$top2.2 <- renderValueBox({
    db=database()
    
    valueBox(
      paste0(prettyNum(sum(db[which(db$flu_a==1),]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[which(db$flu_a==1),]$count,na.rm = T)/sum(db[db$ndetect!=1 & db$covid_pcr!=1 & db$covid_rap!=1,]$count,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "fluA",color = "blue",icon = icon("virus"))
  })
  
  output$top2.3 <- renderValueBox({
    db=database()
    
    valueBox(
      paste0(prettyNum(sum(db[which(db$flu_b==1),]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[which(db$flu_b==1),]$count,na.rm = T)/sum(db[db$ndetect!=1 & db$covid_pcr!=1 & db$covid_rap!=1,]$count,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "fluB",color = "red",icon = icon("virus"))
  })
  
  output$top2.4 <- renderValueBox({
    db=database()
    
    valueBox(
      paste0(prettyNum(sum(db[which(db$vsr==1),]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[which(db$vsr==1),]$count,na.rm = T)/sum(db[db$ndetect!=1 & db$covid_pcr!=1 & db$covid_rap!=1,]$count,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "VSR",color = "yellow",icon = icon("virus"))
  })
  
  output$top2.5 <- renderValueBox({
    db=database()
    
    valueBox(
      paste0(prettyNum(sum(db[which(db$hpvi==1),]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[which(db$hpvi==1),]$count,na.rm = T)/sum(db[db$ndetect!=1 & db$covid_pcr!=1 & db$covid_rap!=1,]$count,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "hPVI",color = "green",icon = icon("virus"))
  })
  
  output$top2.6 <- renderValueBox({
    db=database()
    
    valueBox(
      paste0(prettyNum(sum(db[which(db$adv==1),]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[which(db$adv==1),]$count,na.rm = T)/sum(db[db$ndetect!=1 & db$covid_pcr!=1 & db$covid_rap!=1,]$count,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "AdV",color = "navy",icon = icon("virus"))
  })
  
  output$top2.7 <- renderValueBox({
    db=database()
    
    valueBox(
      paste0(prettyNum(sum(db[which(db$hev==1),]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[which(db$hev==1),]$count,na.rm = T)/sum(db[db$ndetect!=1 & db$covid_pcr!=1 & db$covid_rap!=1,]$count,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "HEV",color = "teal",icon = icon("virus"))
  })
  
  output$top2.8 <- renderValueBox({
    db=database()
    
    valueBox(
      paste0(prettyNum(sum(db[which(db$mpv==1),]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[which(db$mpv==1),]$count,na.rm = T)/sum(db[db$ndetect!=1 & db$covid_pcr!=1 & db$covid_rap!=1,]$count,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "hMPV",color = "olive",icon = icon("virus"))
  })
  
  output$top2.9 <- renderValueBox({
    db=database()
    
    valueBox(
      paste0(prettyNum(sum(db[which(db$oc43==1),]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[which(db$oc43==1),]$count,na.rm = T)/sum(db[db$ndetect!=1 & db$covid_pcr!=1 & db$covid_rap!=1,]$count,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "β-coronavírus",color = "orange",icon = icon("virus"))
  })
  
  output$top2.10 <- renderValueBox({
    db=database()
    
    valueBox(
      paste0(prettyNum(sum(db[which(db$acorona==1),]$count,na.rm = T),big.mark = ".",decimal.mark = ",",scientific = FALSE),
             " (",prettyNum(round(sum(db[which(db$acorona==1),]$count,na.rm = T)/sum(db[db$ndetect!=1 & db$covid_pcr!=1 & db$covid_rap!=1,]$count,na.rm = T)*100,1),decimal.mark = ","),"%)"),
      "α-coronavírus",color = "purple",icon = icon("virus"))
  })
  
  output$plot1 <- renderGirafe({
    df=database() %>% 
      group_by(date=dateweek,result=resultado) %>% 
      summarise(cases=sum(count)) %>% 
      mutate(year=year(date)) %>% 
      mutate(result = forcats::fct_relevel(result, "Não detectado"))
    
    
    # df$result=factor(df$result,levels=c("Outros vírus respiratórios","COVID-19","COVID-19 + Outro vírus respiratórios",
    #                                        "Não detectado"))
    
    p=ggplot(df,aes(x = date,y=cases,fill=result))+
      geom_bar_interactive(stat='identity',aes(tooltip=cases),color="black")+
      theme_classic()+
      theme(legend.position = "bottom")+
      theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
            axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
            axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            legend.text=element_text(size = 20),
            text=element_text(size = 20),
            title = element_text(size = 20),
            strip.text = element_text(size=20))+
      guides(fill=guide_legend(nrow=1,byrow=TRUE))+
      scale_fill_manual("",values=c("COVID-19"="#DD4B39","COVID-19 + Outro vírus respiratórios"="#FF851B",
                                    "Outros vírus respiratórios"="#605CA8","Não detectado"="white"))+
      scale_x_date(date_breaks = "1 weeks",expand=c(0,0),date_labels = "%U")+
      facet_grid(~year, space="free_x", scales="free_x", switch="x") +
      theme(strip.placement = "outside",strip.background = element_rect(fill=NA,colour=NA),panel.spacing=unit(0.3,"cm"))+
      scale_y_continuous(expand=c(0,0))+
      labs(title="Distribuição dos resultados do painel de vírus respiratórios por semana epidemiológica.",x="Semana epidemiológica",y="Frequência",fill="")
    
    girafe(ggobj = p,width_svg = 20,height_svg = 9)
  })
  
  output$plot2 <- renderGirafe({
    df=database() %>% 
      dplyr::select(dateweek,covid_rap,covid_pcr,vsr,flu_a,flu_b,hpvi,acorona,adv,hev,mpv,oc43) %>% 
      #mutate(covid=ifelse(covid_rap==1 | covid_pcr==1,1,0)) %>% 
      dplyr::select(-covid_rap,-covid_pcr) %>% 
      mutate_at(.vars = c("vsr","flu_a","flu_b","hpvi","acorona","adv","hev","mpv","oc43"), funs(ifelse(. == 9, 0, .))) %>%
      group_by(dateweek) %>% 
      summarise(across(everything(), list(sum))) %>% 
      pivot_longer(!dateweek, names_to = "virus", values_to = "count",values_drop_na = TRUE) %>% 
      mutate(virus=stringr::str_replace_all(virus,"_1","")) %>% 
      mutate(year=lubridate::year(dateweek)) %>% 
      filter(virus %in% c(input$virus)) 
      # mutate(virus = recode(virus,"flu_a"="fluA","flu_b"="fluB","vsr"="VSR","hpvi"="hPVI","adv"="AdV",
      #                       "hev"="HEV","mpv"="hMPV","oc43"="β-coronavírus","acorona"="α-coronavírus"))
    
    p=ggplot(df,aes(x = dateweek,y=count,fill=virus))+
      geom_bar_interactive(stat='identity',aes(tooltip=paste0(virus,": ",count)),color="black")+
      #geom_bar(stat="identity",position = position_stack(reverse = T))+
      theme_classic()+
      theme(legend.position = "bottom")+
      theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
            axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
            axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            legend.text=element_text(size = 20),
            text=element_text(size = 20),
            title = element_text(size = 20),
            strip.text = element_text(size=20))+
      guides(fill=guide_legend(nrow=1,byrow=TRUE))+
      scale_fill_manual("",values=c("flu_a"="#0073B7","flu_b"="#DD4B39","vsr"="#F39C12","hpvi"="#00A65A","adv"="#001A36",
                                    "hev"="#39CCCC","mpv"="#3D9970","oc43"="#FF851B","acorona"="#605CA8",
                                    "HBoV"="#3C8DBC","HRV"="#D81B60"),
                        labels=c("flu_a"="fluA","flu_b"="fluB","vsr"="VSR","hpvi"="hPVI","adv"="AdV",
                                 "hev"="HEV","mpv"="hMPV","oc43"="β-coronavírus","acorona"="α-coronavírus"))+#,"COVID-19"="#528BB8"))+
      scale_x_date(date_breaks = "1 weeks",expand=c(0,0),date_labels = "%U")+
      facet_grid(~year, space="free_x", scales="free_x", switch="x") +
      theme(strip.placement = "outside",strip.background = element_rect(fill=NA,colour=NA),panel.spacing=unit(0.3,"cm"))+
      scale_y_continuous(expand=c(0,0))+
      labs(title="Distribuição dos resultados do painel de vírus respiratórios por semana epidemiológica.",x="Semana epidemiológica",y="Frequência",fill="")
    
    girafe(ggobj = p,width_svg = 20,height_svg = 9)
  })
  
  output$table <- DT::renderDataTable({
    
    db=database()
    
    tab=reshape2::dcast(data.frame(table(db[db$covid_pcr!=9 & db$covid_rap!=9,]$covid_pcr,db[db$covid_pcr!=9 & db$covid_rap!=9,]$covid_rap)),Var1~Var2,value.var="Freq") %>% 
      arrange(desc(Var1)) %>% 
      mutate(Var1=factor(Var1,levels=c(1,0),labels=c("Positivo","Negativo"))) %>% 
      relocate("1", .before = "0") %>% 
      rename("Teste Rápido"="Var1","Negativo"="0","Positivo"="1") %>% 
      mutate(Total=Positivo+Negativo) %>% 
      bind_rows(summarise(., across(where(is.numeric), sum),
                          across(where(is.factor), ~ 'Total')))
      
      my.options <- list(autoWidth = FALSE,
                         searching = FALSE,
                         ordering = FALSE,
                         lengthChange = FALSE,
                         lengthMenu = FALSE,
                         pageLength = FALSE,
                         paging = FALSE,
                         info = FALSE)

      header.style <- "th { font-family: 'Arial'; font-size:22px ;font-weight: bold; color: white; background-color: #3C8DBC;}"
      header.names <- names(tab)

      my.container <- withTags(table(
        style(type = "text/css", header.style),
        class = 'display',
        thead(
          tr(
            th(class = 'dt-left',rowspan = 2, 'Teste rápido',style = "border-right: solid 1px;"),
            th(class = 'dt-center',colspan = 10, 'PCR',style = "border-right: solid 1px;")
          ),
          tr(
           lapply(header.names[-1], th, style = "text-align: center; border-right-width: 1px; border-right-style: solid; border-right-color: white; border-bottom-width: 1px; border-bottom-style: solid; border-bottom-color: white")
          )
        )
      ))
      
      
      caption = htmltools::tags$caption(
        style = 'caption-side: bottom; text-align: left;',
        paste0("Teste de concordância kappa: ",
               prettyNum(sprintf("%.2f",vcd::Kappa(table(db[db$covid_pcr!=9 & db$covid_rap!=9,]$covid_pcr,db[db$covid_pcr!=9 & db$covid_rap!=9,]$covid_rap))$Unweighted[1]),decimal.mark = "."))
      )

      my.table <- datatable(tab,caption=caption,options = my.options,container = my.container,rownames = F,width = '100%') %>%
        formatStyle(columns = c(2:4),
                    width='100px',
                    fontFamily = "Arial",
                    fontSize = "20px",
                    borderRightWidth = "1px",
                    borderRightStyle = "solid",
                    borderRightColor = "white",
                    borderBottomColor = "#ffffff",
                    borderBottomStyle = "solid",
                    borderBottomWidth = "1px",
                    borderCollapse = "collapse",
                    verticalAlign = "middle",
                    textAlign = "center",
                    wordWrap = "break-word") %>%
        formatStyle(columns = c(1),
                    width='100px',
                    fontFamily = "Arial",
                    fontSize = "20px",
                    borderRightWidth = "1px",
                    borderRightStyle = "solid",
                    borderRightColor = "white",
                    borderBottomColor = "#ffffff",
                    borderBottomStyle = "solid",
                    borderBottomWidth = "1px",
                    borderCollapse = "collapse",
                    verticalAlign = "middle",
                    textAlign = "left",
                    wordWrap = "break-word")
      
      print(my.table)
  })
  
}


shinyApp(ui = ui, server = server)
