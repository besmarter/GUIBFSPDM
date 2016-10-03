####################################################
#           Posterior Model Probability            #
# Bayesian General Nesting Spatial Panel Data Model#
#                Contiguity Matrix                 #
#                September 12, 2016                #
#         Professor: Andres Ramirez Hassan         #
#             aramir21@eafit.edu.co                #
#            besmarter.team@gmail.com              #
#             www.besmarter-team.org               #
####################################################

####################################################
# Define UI for Bayesian Econometrics: simulations, models and applications to research, teaching and encoding with responsibility
image<- img(src='BEsmarterLogo.png', height = 200, width = 450) #Local variable

##### 1. Posterior Matrix Probability: First NavBar#####
file1am<- fileInput('file1a', 'Choose File: Contiguity Matrices',
                   accept=c('text/csv',
                            'text/comma-separated-values,text/plain',
                            '.csv'))
filech1am<- checkboxInput('header1a', 'Header', TRUE)
rb1am<- radioButtons('sep1a', 'Separator',
                    c(Comma=',',
                      Semicolon=';',
                      Tab='\t'),
                    selected=',')

file1bm<- fileInput('file1b', 'Choose File: Data set',
                    accept=c('text/csv',
                             'text/comma-separated-values,text/plain',
                             '.csv'))
filech1bm<- checkboxInput('header1b', 'Header', TRUE)
rb1bm<- radioButtons('sep1b', 'Separator',
                     c(Comma=',',
                       Semicolon=';',
                       Tab='\t'),
                     selected=',')

Formula1<- textInput('Formula1a', 'Main Equation', value = '')
HTForm1<- helpText('Introduce Formula. Example: y~x1+x2 where y is the dependent variable, and x1 and x2 are independent variables.')
Units<- textInput('Units1', 'Units', value = '')
HTMU<- helpText('Introduce number of cross sectional units')
Time<- textInput('Time1', 'Periods', value = '')
HTMT<- helpText('Introduce number of periods')

Formula2<- textInput('Formula2a', 'Main Equation', value = '')
HTForm1<- helpText('Introduce Formula. Example: y~x1+x2 where y is the dependent variable, and x1 and x2 are independent variables.')
Units2<- textInput('Units2a', 'Units', value = '')
HTMU<- helpText('Introduce number of cross sectional units')
Time2<- textInput('Time2a', 'Periods', value = '')
HTMT<- helpText('Introduce number of periods')
###################################
##### 2. General Nesting Spatial Panel Data Model: Second NavBar#####

file2am<- fileInput('file2a', 'Choose File: Contiguity Matrix',
                   accept=c('text/csv',
                            'text/comma-separated-values,text/plain',
                            '.csv'))
filech2am<- checkboxInput('header2a', 'Header', TRUE)
rb2am<- radioButtons('sep2a', 'Separator',
                    c(Comma=',',
                      Semicolon=';',
                      Tab='\t'),
                    selected=',')

file2bm<- fileInput('file2b', 'Choose File: Data set',
                   accept=c('text/csv',
                            'text/comma-separated-values,text/plain',
                            '.csv'))
filech2bm<- checkboxInput('header2b', 'Header', TRUE)
rb2bm<- radioButtons('sep2b', 'Separator',
                    c(Comma=',',
                      Semicolon=';',
                      Tab='\t'),
                    selected=',')

#######################################################

it1<- sliderInput('it',
                  'MCMC Iterations:',
                  value = 1000,
                  min = 1000,
                  max = 50000,
                  step = 1000)
it2<- sliderInput('burnin',
                  'Burn-in Sample:',
                  value = 100,
                  min = 100,
                  max = 10000,
                  step = 100)

it3<- selectInput('keep', 'Thinning parameter:',
                  choices = c('1', '10', '20', '50', '100'), selected = '1')

TPlambda<- textInput('TPlambda1', 'Tuning parameter lambda', value = '')
HTM<- helpText('Introduce tuning parameter for lambda in the Metropolis-Hastings Algorithm')
TPrho<- textInput('TPrho1', 'Tuning parameter rho', value = '')
HTM<- helpText('Introduce tuning parameter for rho in the Metropolis-Hastings Algorithm')

HT<- helpText('Click the button (Go!) after importing contiguity matrices and data set.')
BE<- helpText('Warning: Be patient this may take several minutes!!!')
PMPtext<- helpText('POSTERIOR MATRIX PROBABILITIES')

######Hyperparameters####
PMean<- textInput('PMeanL', 'Prior Mean Vector: Location Parameters', value = '')
HTM<- helpText('Introduce prior mean vector location parameters. Example: c(0,0).')
PVar<- textInput('PVarL', 'Prior Covariance Matrix: Location Parameters', value = '')
HTV<- helpText('Introduce prior covariances location parameters by row. It has to be symmetric. Example: c(100,0,0,100)')
Psh<-textInput('PShL', 'Prior Shape Parameter: Variance Parameter', value = '')
HTsh<- helpText('Introduce Prior shape Parameter. Example: 0.001')
Psc<- textInput('PScL', 'Prior Scale Parameter: Variance Parameter', value = '')
HTsc<- helpText('Introduce prior scale parameter. Example: 0.001')
TLambda<-textInput('TLambda1', 'Tuning parameter Lambda', value = '')
HTTLam<- helpText('Introduce tuning parameter for lambda in the Metropolis-Hastings algorimth. Example: 0.25')
TRho<- textInput('TRho1', 'Tuning parameter Rho', value = '')
HTTRho<- helpText('Introduce tuning parameter for rho in the Metropolis-Hastings algorimth. Example: 0.25')
##### 1.1 #########
go11<- actionButton('goButton11', 'Go!')
DL11<- downloadButton('download11', 'Download Posterior Matrix Probabilities')
DLP11<- downloadButton('multiDownload11', 'Download Posterior Graphs')

##### 2.1 #########
go21<- actionButton('goButton21', 'Go!')
DL21<- downloadButton('download21', 'Download Posterior Chains')
DLP21<- downloadButton('multiDownload21', 'Download Posterior Graphs')
pplot21<- plotOutput('plot21', height = 1)

shinyUI(
  navbarPage(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx'),

             sidebarLayout(
               sidebarPanel(h1(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx')),
                            h2('Bayesian Econometrics: simulations, models and applications to research, teaching and encoding with responsibility'),
                            image,
                            h4('Professor Andres Ramirez Hassan'),
                            h4(span('besmarter.team@gmail.com', style = 'color:blue')),
                            h4('This is free graphical user interface and comes with ABSOLUTELY NO WARRANTY.')),
               mainPanel(h3(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx'),' is a team of friends from ', a('Universidad EAFIT', href ='http://www.eafit.edu.co/Paginas/index.aspx'), ' (Medellin, Colombia) that promotes research, teaching and encoding of Bayesian Econometrics with social responsibility.'
               ),
               h3('Bayesian Econometrics allows establishing a framework that simultaneously unifies decision theory, statistical inference, and probability theory under a single philosophically and mathematically consistent structure.'),
               br(),
               h3(strong('VISION')),
               h4(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx'), 'envisions worldwide econometric research, teaching and applications based on a Bayesian framework that:'),
               h4(em('inspire'), ' new econometric ideas,'),
               h4(em('create'), ' a user friendly environment for applications of Bayesian econometrics,'),
               h4(em('transform'), ' classic econometric research, teaching and applications,'),
               h4('and where one of the main concerns of science is to solve social problems.'),
               br(),
               h3(strong('MISSION')),
               h4(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx'), 'leads and excels in the generation and dissemination of Bayesian Econometric knowledge through research, teaching and software.')
               )),

             tabPanel('Presentation',
                      sidebarLayout(
                        sidebarPanel(h1(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx')),
                                     h2('Bayesian Econometrics: simulations, models and applications to research, teaching and encoding with responsibility'),
                                     image,
                                     h4('Professor Andres Ramirez Hassan'),
                                     h4(span('besmarter.team@gmail.com', style = 'color:blue')),
                                     h4('This is free graphical user interface and comes with ABSOLUTELY NO WARRANTY.')
                      ),
                        mainPanel(h3(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx'),' is a team of friends from ', a('Universidad EAFIT', href ='http://www.eafit.edu.co/Paginas/index.aspx'), ' (Medellin, Colombia) that promotes research, teaching and encoding of Bayesian Econometrics with social responsibility.'
                        ),
                        h3('Bayesian Econometrics allows establishing a framework that simultaneously unifies decision theory, statistical inference, and probability theory under a single philosophically and mathematically consistent structure.'),
                        br(),
                        h3(strong('VISION')),
                        h4(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx'), 'envisions worldwide econometric research, teaching and applications based on a Bayesian framework that:'),
                        h4(em('inspire'), ' new econometric ideas.'),
                        h4(em('create'), ' a user friendly environment for applications of Bayesian econometrics.'),
                        h4(em('transform'), ' classic econometric research, teaching and applications.'),
                        h4('and where one of the main concerns of science is to solve social problems.'),
                        br(),
                        h3(strong('MISSION')),
                        h4(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx'), 'leads and excels in the generation and dissemination of Bayesian Econometric knowledge through research, teaching and software.')
                        ))),

             tabPanel('Posterior Matrix Probability',
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons('M11', 'Models',
                                       c('No Selection'='m110','Spatial Autoregressive with Autoregressive Disturbances'='m111', 'Spatial Autoregressive Model'='m112','Spatial Error Model'='m113')
                          ),
                          h1(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx')),
                          h2('Bayesian Econometrics: simulations, models and applications to research, teaching and encoding with responsibility'),
                          image,
                          h4('Professor Andres Ramirez Hassan'),
                          h4(span('besmarter.team@gmail.com', style = 'color:blue'))
                        ),
                          mainPanel(fluidRow(column(5,file1am),column(2,filech1am),column(2,rb1am)),fluidRow(column(5,file1bm),column(2,filech1bm),column(2,rb1bm)),fluidRow(column(6,Formula1),column(3,Units),column(3,Time)),fluidRow(column(6,HTForm1),column(3,HTMU),column(3,HTMT)),
                                    HT,go11,BE,DL11,PMPtext,verbatimTextOutput('summary11')))),

             tabPanel('Bayesian Spatial Fixed Effects Panel Data Models',
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons('M21', 'Models',
                                       c('No Selection'='m210','General Nesting Model'='m211', 'Spatial Autoregressive with Autoregressive Disturbances'='m212', 'Spatial Durbin Model'='m213','Spatial Durbin Error Model'='m214','Spatial Autoregressive Model'='m215','Spatial Error Model'='m216','Spatial Lag Model'='m217','No Spatial Effects'='m218')
                          ),
                          h1(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx')),
                          h2('Bayesian Econometrics: simulations, models and applications to research, teaching and encoding with responsibility'),
                          image,
                          h4('Professor Andres Ramirez Hassan'),
                          h4(span('besmarter.team@gmail.com', style = 'color:blue'))
                        ),
                        mainPanel(fluidRow(column(5,file2am),column(2,filech2am),column(2,rb2am)),fluidRow(column(5,file2bm),column(2,filech2bm),column(2,rb2bm)),fluidRow(column(4,it1),column(4,it2),column(4,it3)),fluidRow(column(3,TLambda),column(3,TRho)),fluidRow(column(3,HTTLam),column(3,HTTRho)),
                                  fluidRow(column(6,Formula2),column(3,Units2),column(3,Time2)),fluidRow(column(6,HTForm1),column(3,HTMU),column(3,HTMT)),
                                  fluidRow(column(3,PMean),column(3,PVar),column(3,Psh),column(3,Psc)),fluidRow(column(3,HTM),column(3,HTV),column(3,HTsh),column(3,HTsc)),HT,go21,BE,DL21,DLP21,verbatimTextOutput('summary21'),pplot21))),

             tabPanel('Help',
                      sidebarLayout(
                        sidebarPanel(
                          h1(a(em(strong('BEsmarter',style = 'color:light blue')),href = 'http://www.eafit.edu.co/docentes-investigadores/Paginas/andres-ramirez.aspx')),
                          h2('Bayesian Econometrics: simulations, models and applications to research, teaching and encoding with responsibility'),
                          image,
                          h4('Professor Andres Ramirez Hassan'),
                          h4(span('besmarter.team@gmail.com', style = 'color:blue'))),
                        mainPanel(h4(a(em(strong('Video tutorial: graphical user interface for Bayesian estimation of fixed effects spatial panel data models.',style = 'color:light blue')),href = 'https://www.youtube.com/watch?v=9vOyZQutjn8')
                                     ),
                                     h4(a(em(strong('Ramirez Hassan, A. (2016). "The Interplay Between the Bayesian and Frequentist Approaches: A General Nesting Spatial Panel Data Model" Spatial Economic Analysis. Accepted manuscript.',style = 'color:light blue')),href = 'http://www.tandfonline.com/toc/rsea20/current')
                        ))
                      ))
             ))

