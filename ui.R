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
        title = 'HMSRep v 0.1'
        )

body <- dashboardBody(
        fluidRow(
		column(width=12,style = 'overflow-x: scroll',
                        tabBox(width = 0,
                                id='tabvals',
                                tabPanel('msms_1',p(HTML("Distribution <a href='#freqdistr1'>ADPr</a>")),
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
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                DT::dataTableOutput("uniquePepTable1")
                                                )
                                        ),
					 br(),
                                        fluidRow(
                                        column(2,p(HTML("Main <a href='#mytable1'>table</a>")),
                                                DT::dataTableOutput("uniquePepTable11")
                                                )
                                        )
					#verbatimTextOutput("locresevent")

				),
				tabPanel('msms_2',p(HTML("Distribution <a href='#freqdistr2'>ADPr</a>")),
                                        dataTableOutput("mytable2"), value=2,
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
                                tabPanel('msms_3',p(HTML("Distribution <a href='#freqdistr3'>ADPr</a>")),
                                        dataTableOutput("mytable3"), value=3,
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
                                tabPanel('msms_4',p(HTML("Distribution <a href='#freqdistr4'>ADPr</a>")),
                                        dataTableOutput("mytable4"), value=4,	
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
				tabPanel('plots',p(HTML("Distribution <a href='#heatmapProt'>ADPr</a>")),
					value=5,
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
				tabPanel('GO Enrichment BP',p(HTML("Distribution <a href='#goPlot1BP'>ADPr</a>")),
                                        value=6,
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
				tabPanel('GO Enrichment MF',p(HTML("Distribution <a href='#goPlot1BP'>ADPr</a>")),
                                        value=6,	
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
				tabPanel('GO Enrichment CC',p(HTML("Distribution <a href='#goPlot1BP'>ADPr</a>")),
                                        value=6,
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
  			h2('Upload csv folder'),
  			fileInput('csvfile1', 'Select a single or up to 4 CSV file(s) having "pep_var_mod_conf" column',
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
		selectInput("searchtype1", "Select search type (EThcD)", choices = c('None','EThcD RSY search'='ADP-Ribosyl-EThcD (RSY) (RSY)','EThcD DEKRSY search'='ADP-Ribosyl-EThcD (SRDEKY) (DEKRSY)','EThcD DEKRS search'='ADP-Ribosyl-EThcD (KRDES) (DEKRS)')),
		selectInput("searchtype2", "Select search type (HCD)", choices = c('None','HCD RSY search'='ADP-Ribosyl-HCD (RSY) (RSY)', 'HCD DEKRSY search'='ADP-Ribosyl-HCD (SRDEKY) (DEKRSY)','HCD DEKRS search'='ADP-Ribosyl-HCD (KRDES) (DEKRS)')),
		selectInput("searchorder", "Select search order of the ADPr modification and then press Create Table button", choices = c('first'='1', 'second'='2', 'third'='3', 'fourth'='4', 'fifth'='5')),
		textAreaInput("exptype1", rows=1, cols=40, "Enter title for exp1"),
		textAreaInput("exptype2", rows=1, cols=40, "Enter title for exp2"),
		textAreaInput("exptype3", rows=1, cols=40, "Enter title for exp3"),
		textAreaInput("exptype4", rows=1, cols=40, "Enter title for exp4"),
		selectInput("organism", "Select organism fo ontology analysis", choices = c('Mouse'='Mouse', 'Human'='Human')),
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

