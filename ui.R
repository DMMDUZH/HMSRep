library(ggplot2)
library(shiny)
library(shinydashboard)
#library(XLConnect)
library("DT")
library(calibrate)
library(shinyjs)
library(plotly)
library(Hmisc)
library(stringr)

headers <- dashboardHeader(
        title = 'HMSRep v 0.4'
        )

body <- dashboardBody(
	fluidRow(
		column(width=12,style = 'overflow-x: scroll',
                        tabBox(width = 0,
                                id='tabvals',
                                tabPanel('EThcD 1',p(HTML("Distribution <a href='#freqdistr1'>ADPr</a>")),
                                        dataTableOutput("mytable1"), value=1,
					br(),
					hr(),
                                	fluidRow(
                                        column(3,
                                                dataTableOutput("mytable11")
                                                ),
					column(4,offset = 3,
                                                dataTableOutput("mytable12")
                                                )
					),	
					br(),
					fluidRow(
					column(3,
                                                dataTableOutput("mytable13")
                                                ),
					column(4,p(HTML("Main <a href='#mytable1'>table</a>")),offset = 3,
                                                plotOutput("freqdistr12", width = "500px", height = "600px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                plotlyOutput("freqdistr1", width = "500px", height = "600px")
                                                ),
					column(4,offset = 3,
						plotOutput("psmPepScore1", width = "500px", height = "600px")
						)
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                plotOutput("psmLocConf1", width = "500px", height = "600px")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(6,p(HTML("Main <a href='#mytable1'>table</a>")),
						tags$head(
                                        	tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
                                        	),
                                                DT::dataTableOutput("uniquePepTable1")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(6,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                DT::dataTableOutput("uniquePepTable11")
                                                )
                                        )
					#verbatimTextOutput("locresevent")

				),
				tabPanel('EThcD 2',p(HTML("Distribution <a href='#freqdistr2'>ADPr</a>")),
                                        dataTableOutput("mytable2"), value=1,
                                        br(),
                                        hr(),
                                        fluidRow(
                                        column(4,
                                                dataTableOutput("mytable21")
                                                )
					),
                                        br(),
                                        fluidRow(
                                        column(3,
                                                dataTableOutput("mytable22")
                                                ),
                                        column(4,p(HTML("Main <a href='#mytable2'>table</a>")),offset = 3,
                                                plotOutput("freqdistr22", width = "500px", height = "600px")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable2'>table</a>")),
                                                plotlyOutput("freqdistr2", width = "500px", height = "600px")
                                                ),
					column(4,offset = 3,
                                                plotOutput("psmPepScore2", width = "500px", height = "600px")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                plotOutput("psmLocConf2", width = "500px", height = "600px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                dataTableOutput("uniquePepTable2")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                dataTableOutput("uniquePepTable21")
                                                )
                                        )
                            	),
                                tabPanel('EThcD 3',p(HTML("Distribution <a href='#freqdistr3'>ADPr</a>")),
                                        dataTableOutput("mytable3"), value=1,
					br(),
                                        hr(),
                                        fluidRow(
                                        column(4,
                                                dataTableOutput("mytable31")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(3,
                                                dataTableOutput("mytable32")
                                                ),
                                        column(4,p(HTML("Main <a href='#mytable3'>table</a>")),offset = 3,
                                                plotOutput("freqdistr32", width = "500px", height = "600px")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable3'>table</a>")),
                                                plotlyOutput("freqdistr3", width = "500px", height = "600px")
                                                ),
					column(4,offset = 3,
                                                plotOutput("psmPepScore3", width = "500px", height = "600px")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                plotOutput("psmLocConf3", width = "500px", height = "600px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                dataTableOutput("uniquePepTable3")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                dataTableOutput("uniquePepTable31")
                                                )
                                        )
				),
                                tabPanel('EThcD 4',p(HTML("Distribution <a href='#freqdistr4'>ADPr</a>")),
                                        dataTableOutput("mytable4"), value=1,	
					br(),
                                        hr(),
                                        fluidRow(
                                        column(4,
                                                dataTableOutput("mytable41")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(3,
                                                dataTableOutput("mytable42")
                                                ),
                                        column(4,p(HTML("Main <a href='#mytable4'>table</a>")),offset = 3,
                                                plotOutput("freqdistr42", width = "500px", height = "600px")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable4'>table</a>")),
                                                plotlyOutput("freqdistr4", width = "500px", height = "600px")
                                                ),
					column(4,offset = 3,
                                                plotOutput("psmPepScore4", width = "500px", height = "600px")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                plotOutput("psmLocConf4", width = "500px", height = "600px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                dataTableOutput("uniquePepTable4")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                dataTableOutput("uniquePepTable41")
                                                )
                                        )
				),
				tabPanel('HCD 1',p(HTML("Distribution <a href='#freqdistrH1'>ADPr</a>")),
                                        dataTableOutput("mytableH1"), value=1,
					br(),
					hr(),
                                	fluidRow(
                                        column(3,
                                                dataTableOutput("mytableH11")
                                                ),
					column(4,offset = 3,
                                                dataTableOutput("mytableH12")
                                                )
					),	
					br(),
					fluidRow(
					column(3,
                                                dataTableOutput("mytableH13")
                                                ),
					column(4,p(HTML("Main <a href='#mytableH1'>table</a>")),offset = 3,
                                                plotOutput("freqdistrH12", width = "500px", height = "600px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH1'>table</a>")),
                                                plotlyOutput("freqdistrH1", width = "500px", height = "600px")
                                                ),
					column(4,offset = 3,
						plotOutput("psmPepScoreH1", width = "500px", height = "600px")
						)
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH1'>table</a>")),
                                                plotOutput("psmLocConfH1", width = "500px", height = "600px")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH1'>table</a>")),
                                                DT::dataTableOutput("uniquePepTableH1")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH1'>table</a>")),
                                                DT::dataTableOutput("uniquePepTableH11")
                                                )
                                        )
					#verbatimTextOutput("locresevent")

				),
				tabPanel('HCD 2',p(HTML("Distribution <a href='#freqdistrH2'>ADPr</a>")),
                                        dataTableOutput("mytableH2"), value=1,
					br(),
					hr(),
                                	fluidRow(
                                        column(3,
                                                dataTableOutput("mytableH21")
                                                ),
					column(4,offset = 3,
                                                dataTableOutput("mytableH22")
                                                )
					),	
					br(),
					fluidRow(
					column(3,
                                                dataTableOutput("mytableH23")
                                                ),
					column(4,p(HTML("Main <a href='#mytableH2'>table</a>")),offset = 3,
                                                plotOutput("freqdistrH22", width = "500px", height = "600px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH2'>table</a>")),
                                                plotlyOutput("freqdistrH2", width = "500px", height = "600px")
                                                ),
					column(4,offset = 3,
						plotOutput("psmPepScoreH2", width = "500px", height = "600px")
						)
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH2'>table</a>")),
                                                plotOutput("psmLocConfH2", width = "500px", height = "600px")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH2'>table</a>")),
                                                DT::dataTableOutput("uniquePepTableH2")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH2'>table</a>")),
                                                DT::dataTableOutput("uniquePepTableH21")
                                                )
                                        )
					#verbatimTextOutput("locresevent")

				),
				tabPanel('HCD 3',p(HTML("Distribution <a href='#freqdistrH3'>ADPr</a>")),
                                        dataTableOutput("mytableH3"), value=1,
					br(),
					hr(),
                                	fluidRow(
                                        column(3,
                                                dataTableOutput("mytableH31")
                                                ),
					column(4,offset = 3,
                                                dataTableOutput("mytableH32")
                                                )
					),	
					br(),
					fluidRow(
					column(3,
                                                dataTableOutput("mytableH33")
                                                ),
					column(4,p(HTML("Main <a href='#mytableH3'>table</a>")),offset = 3,
                                                plotOutput("freqdistrH32", width = "500px", height = "600px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH3'>table</a>")),
                                                plotlyOutput("freqdistrH3", width = "500px", height = "600px")
                                                ),
					column(4,offset = 3,
						plotOutput("psmPepScoreH3", width = "500px", height = "600px")
						)
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH3'>table</a>")),
                                                plotOutput("psmLocConfH3", width = "500px", height = "600px")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH3'>table</a>")),
                                                DT::dataTableOutput("uniquePepTableH3")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH3'>table</a>")),
                                                DT::dataTableOutput("uniquePepTableH31")
                                                )
                                        )
					#verbatimTextOutput("locresevent")

				),
				tabPanel('HCD 4',p(HTML("Distribution <a href='#freqdistrH4'>ADPr</a>")),
                                        dataTableOutput("mytableH4"), value=1,
					br(),
					hr(),
                                	fluidRow(
                                        column(3,
                                                dataTableOutput("mytableH41")
                                                ),
					column(4,offset = 3,
                                                dataTableOutput("mytableH42")
                                                )
					),	
					br(),
					fluidRow(
					column(3,
                                                dataTableOutput("mytableH43")
                                                ),
					column(4,p(HTML("Main <a href='#mytableH4'>table</a>")),offset = 3,
                                                plotOutput("freqdistrH42", width = "500px", height = "600px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH4'>table</a>")),
                                                plotlyOutput("freqdistrH4", width = "500px", height = "600px")
                                                ),
					column(4,offset = 3,
						plotOutput("psmPepScoreH4", width = "500px", height = "600px")
						)
                                        ),
					br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH4'>table</a>")),
                                                plotOutput("psmLocConfH4", width = "500px", height = "600px")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH4'>table</a>")),
                                                DT::dataTableOutput("uniquePepTableH4")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytableH4'>table</a>")),
                                                DT::dataTableOutput("uniquePepTableH41")
                                                )
                                        )
					#verbatimTextOutput("locresevent")

				),
				tabPanel('EThcD plots',p(HTML("Distribution <a href='#heatmapProt'>ADPr</a>")),
					value=1,
					br(),
                                        fluidRow(
                                        column(3,
                                                plotOutput("vennDiaProt", width = "400px", height = "450px")
                                                ),
					column(4,offset = 3,
                                                plotOutput("vennDiaPep", width = "400px", height = "450px")
                                                )
                                        ),
					br(),
					fluidRow(
                                        column(3,
                                                plotOutput("heatmapProt", width = "800px", height = "800px")
                                                )
					),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("heatmapProtr", width = "800px", height = "800px")
                                                )
                                        ),
					
					br(),
                                        fluidRow(
					column(3,
                                                plotOutput("psmBar1", width = "400px", height = "400px")
                                                ),
					column(4,offset = 3,
                                                plotOutput("psmBar2", width = "400px", height = "400px")
                                                )
					),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("psmBar1r", width = "400px", height = "400px")
                                                ),
                                        column(4,offset = 3,
                                                plotlyOutput("psmBar2r", width = "400px", height = "400px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotOutput("pepBar1", width = "400px", height = "400px")
                                                ),
                                        column(4,offset = 3,
                                                plotOutput("pepBar2", width = "400px", height = "400px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("pepBar1r", width = "400px", height = "400px")
                                                ),
                                        column(4,offset = 3,
                                                plotlyOutput("pepBar2r", width = "400px", height = "400px")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(3,
                                                plotOutput("protBar", width = "400px", height = "400px")
                                                ),
					column(4,offset = 3,
                                                plotlyOutput("protBarr", width = "400px", height = "400px")
                                                )
                                        )
                                ),
				tabPanel('HCD plots',p(HTML("Distribution <a href='#heatmapProtH'>ADPr</a>")),
					value=1,
					br(),
                                        fluidRow(
                                        column(3,
                                                plotOutput("vennDiaProtH", width = "400px", height = "450px")
                                                ),
					column(4,offset = 3,
                                                plotOutput("vennDiaPepH", width = "400px", height = "450px")
                                                )
                                        ),
					br(),
					fluidRow(
                                        column(3,
                                                plotOutput("heatmapProtH", width = "800px", height = "800px")
                                                )
					),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("heatmapProtHr", width = "800px", height = "800px")
                                                )
                                        ),
					
					br(),
                                        fluidRow(
					column(3,
                                                plotOutput("psmBar1H", width = "400px", height = "400px")
                                                ),
					column(4,offset = 3,
                                                plotOutput("psmBar2H", width = "400px", height = "400px")
                                                )
					),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("psmBar1Hr", width = "400px", height = "400px")
                                                ),
                                        column(4,offset = 3,
                                                plotlyOutput("psmBar2Hr", width = "400px", height = "400px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotOutput("pepBar1H", width = "400px", height = "400px")
                                                ),
                                        column(4,offset = 3,
                                                plotOutput("pepBar2H", width = "400px", height = "400px")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("pepBar1Hr", width = "400px", height = "400px")
                                                ),
                                        column(4,offset = 3,
                                                plotlyOutput("pepBar2Hr", width = "400px", height = "400px")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(3,
                                                plotOutput("protBarH", width = "400px", height = "400px")
                                                ),
					column(4,offset = 3,
                                                plotlyOutput("protBarHr", width = "400px", height = "400px")
                                                )
                                        )
                                ),
				tabPanel('HCD 1 & EThcD 1',p(HTML("Distribution <a href='#HCDEThcDTab14'>ADPr</a>")),
                                        dataTableOutput("HCDEThcDTab1"), value=1,
                                        br(),
                                        hr(),
                                        fluidRow(
                                        column(3,
                                                dataTableOutput("HCDEThcDTab11")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("HCDEThcDTab14")
                                                )
                                        )
				),
				tabPanel('HCD 2 & EThcD 2',p(HTML("Distribution <a href='#HCDEThcDTab24'>ADPr</a>")),
                                        dataTableOutput("HCDEThcDTab2"), value=1,
                                        br(),
                                        hr(),
                                        fluidRow(
                                        column(3,
                                                dataTableOutput("HCDEThcDTab21")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("HCDEThcDTab24")
                                                )
                                        ) 
                                ),
				tabPanel('HCD 3 & EThcD 3',p(HTML("Distribution <a href='#HCDEThcDTab34'>ADPr</a>")),
                                        dataTableOutput("HCDEThcDTab3"), value=1,
                                        br(),
                                        hr(),
                                        fluidRow(
                                        column(3,
                                                dataTableOutput("HCDEThcDTab31")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("HCDEThcDTab34")
                                                )
                                        ) 
                                ),
				tabPanel('HCD 4 & EThcD 4',p(HTML("Distribution <a href='#HCDEThcDTab44'>ADPr</a>")),
                                        dataTableOutput("HCDEThcDTab4"), value=1,
                                        br(),
                                        hr(),
                                        fluidRow(
                                        column(3,
                                                dataTableOutput("HCDEThcDTab41")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("HCDEThcDTab44")
                                                )
                                        ) 
                                ),				
				tabPanel('GO Enrichment BP EThcD',p(HTML("Distribution <a href='#goPlot1BP'>ADPr</a>")),
                                        value=1,
                                        br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot1BP", width = "600px", height = "750px")
                                                ),
					column(4,offset = 3,
                                                dataTableOutput("goBP1")
                                                )
					),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot2BP", width = "600px", height = "750px")
                                                )
					),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot3BP", width = "600px", height = "750px")
                                                )
					),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot4BP", width = "600px", height = "750px")
                                                )
					)
				),
				tabPanel('GO Enrichment MF EThcD',p(HTML("Distribution <a href='#goPlot1BP'>ADPr</a>")),
                                        value=1,	
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot1MF", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goMF1")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot2MF", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goMF2")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot3MF", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goMF3")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot4MF", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goMF4")
                                                )
                                        )
				),
				tabPanel('GO Enrichment CC EThcD',p(HTML("Distribution <a href='#goPlot1BP'>ADPr</a>")),
                                        value=1,
                                        br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot1CC", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goCC1")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot2CC", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goCC2")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot3CC", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goCC3")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot4CC", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goCC4")
                                                )
                                        )
				),
				tabPanel('GO Enrichment BP HCD',p(HTML("Distribution <a href='#goPlot1HBP'>ADPr</a>")),
                                        value=1,
                                        br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot1HBP", width = "600px", height = "750px")
                                                ),
					column(4,offset = 3,
                                                dataTableOutput("goBP1H")
                                                )
					),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot2HBP", width = "600px", height = "750px")
                                                )
					),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot3HBP", width = "600px", height = "750px")
                                                )
					),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot4HBP", width = "600px", height = "750px")
                                                )
					)
				),
				tabPanel('GO Enrichment MF HCD',p(HTML("Distribution <a href='#goPlot1HMF'>ADPr</a>")),
                                        value=1,	
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot1HMF", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goMF1H")
                                                )
                                        ),
                                        br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot2HMF", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goMF2H")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot3HMF", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goMF3H")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot4HMF", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goMF4H")
                                                )
                                        )
				),
				tabPanel('GO Enrichment CC HCD',p(HTML("Distribution <a href='#goPlot1HCC'>ADPr</a>")),
                                        value=1,
                                        br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot1HCC", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goCC1H")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot2HCC", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goCC2H")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot3HCC", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goCC3H")
                                                )
                                        ),
					br(),
                                        fluidRow(
                                        column(3,
                                                plotlyOutput("goPlot4HCC", width = "600px", height = "750px")
                                                ),
                                        column(4,offset = 3,
                                                dataTableOutput("goCC4H")
                                                )
                                        )
				)		
			)
		)
	)
)
dashboardPage(
        headers,
        dashboardSidebar(
                conditionalPanel(
                        condition = "input.tabvals == 1",
  			h2('Upload csv files'),
  			fileInput('csvfile1', 'Select a single or up to 4 EThcD CSV file(s) having "pep_var_mod_conf" column',
			multiple = TRUE,
            			accept = c(
              			'text/csv',
              			'text/comma-separated-values',
              			'text/tab-separated-values',
              			'text/plain',
              			'.csv',
              			'.tsv'
            			)
  			),
			fileInput('csvfile2', 'Select a single or up to 4 HCD CSV file(s) having "pep_var_mod_conf" column',
                        multiple = TRUE,
                                accept = c(
                                'text/csv',
                                'text/comma-separated-values',
                                'text/tab-separated-values',
                                'text/plain',
                                '.csv',
                                '.tsv'
                                )
                        ),
  			tags$hr(),
  			checkboxInput('header', 'Header', TRUE),
  			radioButtons('sep', 'Separator',
               			c(Comma=',',
                 		Semicolon=';',
                 		Tab='\t'),
               		','),
  			radioButtons('quote', 'Quote',
               			c(None='',
                 		'Double Quote'='"',
                 		'Single Quote'="'"),
               		''),
		#selectInput("searchtype1", "Select search type (EThcD)", choices = c('None','EThcD RSY search'='ADP-Ribosyl-EThcD (RSY) (RSY)','EThcD DEKRSY search'='ADP-Ribosyl-EThcD (SRDEKY) (DEKRSY)','EThcD DEKRS search'='ADP-Ribosyl-EThcD (KRDES) (DEKRS)')),
		#selectInput("searchorder1", "Select search order of the ADPr modification (EThcD) and then press Create Table button", choices = c('first'='1', 'second'='2', 'third'='3', 'fourth'='4', 'fifth'='5')),
		#selectInput("searchtype2", "Select search type (HCD)", choices = c('None','HCD RSY search'='ADP-Ribosyl-HCD (RSY) (RSY)', 'HCD DEKRSY search'='ADP-Ribosyl-HCD (SRDEKY) (DEKRSY)','HCD DEKRS search'='ADP-Ribosyl-HCD (KRDES) (DEKRS)', 'HCD DEKRSY search R spec'='ADP-Ribosyl-HCD_249_347_541_583 (DEKRSY)')),
		#selectInput("searchorder2", "Select search order of the ADPr modification (HCD) and then press Create Table button", choices = c('first'='1', 'second'='2', 'third'='3', 'fourth'='4', 'fifth'='5')),
		textAreaInput("exptype1", rows=1, cols=40, "Enter title for exp1"),
		textAreaInput("exptype2", rows=1, cols=40, "Enter title for exp2"),
		textAreaInput("exptype3", rows=1, cols=40, "Enter title for exp3"),
		textAreaInput("exptype4", rows=1, cols=40, "Enter title for exp4"),
		selectInput("organism", "Select organism for ontology analysis", choices = c('Mouse'='Mouse', 'Human'='Human')),
		numericInput("pval1", "Enter p-value for BP",max=1, min=0.00005, value = "0.01"),
		numericInput("pval2", "Enter p-value for MF",max=1, min=0.00005, value = "0.01"),
		numericInput("pval3", "Enter p-value for CC",max=1, min=0.00005, value = "0.01"),
		actionButton("choice", "Create Tables"),
		checkboxGroupInput('select1',
                        'Columns in msms_file1 to show:',
                        "test",
                        selected = NULL),
                        helpText('Here we can select the columns to show in the
                        table.')
		)
 	),
body
)

