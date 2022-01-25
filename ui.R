## ui.R ##


# PopHumanVar is licensed under the GNU General Public License (GPL) v3.0 (https://github.com/ainacolovila/PopHumanVar/blob/master/LICENSE)

ui <- dashboardPage(

	skin="blue",
	title= "PopHumanVar", # Show in tab
  
	# HEADER ------------------------------------------------------------------
  
	dashboardHeader(
		# Like them
		titleWidth='100%',
		title = span(
			column(3, tags$a(href="javascript:history.go(0)", tags$img(src="image.png", width = '50%'))),
			column(1),
			column(8, 
				class="title-box",
				tags$h1(class="primary-title",  
				   'PopHumanVar'), 
			 	tags$h3(class="primary-subtitle",  
					'Functional Characterization & Prioritization of Genomic Variants')
			)
		),
		dropdownMenu(
			type = "notifications", 
			headerText = strong("INSTRUCTIONS"), 
			icon = icon("question"), 
			badgeStatus = NULL,
			notificationItem(
				text = help$text[1],
				icon = icon(help$icon[1])
		  	),
			notificationItem(
				text = help$text[2],
				icon = icon(help$icon[2])
		  	),
			notificationItem(
				text = help$text[3],
				icon = icon(help$icon[3])
		  	),
			notificationItem(
				text = help$text[4],
				icon = icon(help$icon[4])
		  	),
			notificationItem(
				text = help$text[5],
				icon = icon(help$icon[5])
			),
			notificationItem(
				text = help$text[6],
				icon = icon(help$icon[6])
			),
			notificationItem(
				text = strong(help$text[7]),
				icon = icon(help$icon[7])
		  	),
			notificationItem(
				text = help$text[8],
				icon = icon(help$icon[8])
		  	),
		  	notificationItem(
				text = help$text[9],
				icon = icon(help$icon[9])
			)
		),
		tags$li(
		  a(tags$img(src='PHcara.png',height='20'),
			href = "https://pophuman.uab.cat/",
			title = "",
			target = "_blank"
		  ),
		  class = "dropdown"
		),
		tags$li(
		  a(tags$img(src='PHSprint.png',height='20'),
			href = "https://pophumanscan.uab.cat/",
			title = "",
			target = "_blank"
		  ),
		  class = "dropdown"
		),
		tags$li(
		  a(tags$img(src='GHw.png',height='20'),
			href = "https://github.com/ainacolovila",
			title = "",
			target = "_blank"
		  ),
		  class = "dropdown"
		)
		
	),
  
  # SIDEBAR -----------------------------------------------------------------
  
	dashboardSidebar(
		width = 325,
		useShinyjs(),
		# Adjust the sidebar
		sidebarMenu(
			id="navMenu",
			h5(" "),
			menuItem("Stats Visualization", 
				icon = icon("chart-bar"),
				tabName = 'principal',
				menuItem("Selection (iHS & nSL)",
					icon = icon("stream"), 
					tabName = "selection"),
				menuItem("Favored Mutation (iSAFE)",
					icon = icon("asterisk"), 
					tabName = "isafe"),
				menuItem("Functional Description",
					icon = icon("layer-group"), 
					tabName = "functional"),
				menuItem("Age Information",
					icon = icon("hourglass"), 
					tabName = "age"),
				menuItem("Summary Report",
					icon = icon("file-alt"), 
					tabName = "report")
	  		),
			menuItem("Download",
				icon = icon("file-download"),
				tabName = "download",
				menuItem("Current Region",
					icon = icon("thumbtack"),
					tabName = "region_download"),
				menuItem("Batch Download",
					icon = icon("file-archive"),
					tabName = "batch_download")
			),
			menuItem("Tutorial", 
					icon = icon("question-circle"), 
					tabName = "tutorial"),
			menuItem("About Us", 
					icon = icon("info"), 
					tabName = "Information"),
			br()
		),
		shiny::conditionalPanel(
		  condition = "(input.navMenu !== 'tutorial' && input.navMenu !== 'Information')",
		  sidebarMenu(
			id="filtersMenu",
			h5("FILTERS MENU"),
			menuItem("Coordinates",
				icon = icon("search"),
				tabName = "dashboard",
				br(),
				"  GRCh37/hg19",
				br(),
					# Button for quick search
					useShinyalert(),  # Set up shinyalert
					tags$div(align = 'right',
						actionBttn(inputId = "SearchID",
							label = " Quick search",
							icon = icon("searchengin"),
							style = "simple",
							color="primary",
							size="sm"
						)
					),
					# useSweetAlert("borderless"),
					useSweetAlert(),    
					selectInput('chrInput', 'Chromosome:',
									choices = c(1:22),
									selected = 2),
								textInput('minPosInput', 'Start Position', 
									value = 109500927),
								textInput('maxPosInput', 'End Position', 
									value = 109615828),
								h5(""),
								#Input per fer zoom out
								tags$div(align = 'right', 
									  radioGroupButtons(
										inputId = 'zoomCoordenates', 
										choices =c(`<i class="fa fa-search-minus fa-2x"  aria-hidden="true"></i>`= "-1", 
											 `<i class="fa fa-search-plus fa-2x"  aria-hidden="true"></i>`= "1"), 
										selected = '', 
										size='xs',
										direction = "horizontal",
										individual=T, disabled = F))
			),
			menuItem("Populations",
				  icon=icon("globe-africa"),
				  tabName="population",
				  #Input per triar les metapoblacions on mirem el valor de IHS
				  checkboxGroupButtons('metapopsInput', 'Metapopulations:', 
						   choices = as.vector(metapopsData$meta), selected = ''),
				  tags$script("$(\"input:checkbox[name='metapopsInput'][value='AFR']\").parent().css('background-color', '#F7F14A');"),
				  tags$script("$(\"input:checkbox[name='metapopsInput'][value='EAS']\").parent().css('background-color', '#33B033');"),
				  tags$script("$(\"input:checkbox[name='metapopsInput'][value='EUR']\").parent().css('background-color', '#5691C4');"),
				  tags$script("$(\"input:checkbox[name='metapopsInput'][value='SAS']\").parent().css('background-color', '#A965BA');"),
				  #Input per triar les poblacions on mirem el valor de IHS
				  checkboxGroupButtons('popsInput', 'Populations:', 
						   choices = sort(as.vector(popsData$pops)), 
						   selected = c('CEU','CHS','YRI','BEB'), 
						   direction = "horizontal",
						   status = "primary", size="sm"),
				  actionButton('selectAllPopsInput', 'Deselect All')
			),
			menuItem("Selection",
				  icon = icon("stream"),
				  tabName = "dashboard",
				  #Input per filtrar per extreme values
				  awesomeCheckbox(inputId = 'extremeValuesInput_ihs',
						  label="Only extreme iHS Values",
						  value = FALSE,
						  status="primary"),
				  awesomeCheckbox(inputId = 'extremeValuesInput_nsl',
						  label="Only extreme nSL Values",
						  value = FALSE,
						  status="primary")
				  ),
					menuItem("Favored mutation",
				  icon = icon("asterisk"),
				  tabName = "dashboard",
				  #Input per filtrar per extreme values
				  awesomeCheckbox(inputId = 'extremeISAFE',
						  label="Only extreme iSAFE Values",
						  value = FALSE,
						  status="primary")
				  ),
			menuItem("Functional Description",
				  icon = icon("layer-group"), 
				  tabName = "dashboard",
				  # SNPEFF
				  menuItem("Variant effect (SnpEff)",
						awesomeCheckboxGroup(inputId = 'VEPInput', 
								 'SnpEff Consequences', 
								 choices = as.vector(VEPinfo$outLabels),
								 selected = as.vector(VEPinfo[rank < 20]$outLabels),
								 status='primary'),
						  actionButton('selectAllVEPInput', 'Select All')
				  ),
				  # REGULOMEDB
				  menuItem("RegulomeDB filter",
						sliderTextInput('RegulomeInput', 
							   h5('Regulome Score:'),
							   choices = regulomeLabels, 
							   selected = c('1a','7'),
							   grid = TRUE)),
				  # CLINVAR
				  menuItem("ClinVar filter",
						awesomeCheckboxGroup(inputId = 'ClinVarInput', 
								 'Clinical Significances:', 
								 choices = as.vector(ClinVarLabels$ClinicalSign),
								 selected = as.vector(ClinVarLabels$ClinicalSign),
								 status='primary'),
						actionButton('selectAllClinVarInput', 'Select All')),
				  # GWAS CAT
				  menuItem("GWAS Catalog",
						sliderInput(
						  'GWASInput',
						  h5('Filter by Association Counts:'),
						  min = 0,
						  max = 100,
						  value = c(0,10),
						  step = 1)
				  ),
				  # DISGENET
				  menuItem("DisGeNET filters",
						sliderInput('DSIInput',
							  h5('Disease Specificity Index:'),
							  min = 0, max = 1, step=0.05, value = c(0, 1)),
						sliderInput('EInput', 
							  h5('Evidence Index:'),
							  min = 0, max = 1, step=0.05, value = c(0, 1)),
						sliderInput('NofPmidsInput',
							  h5('Number of Pubmed IDs:'),
							  min = 0, max = 20,
							  value = c(0,20))
				  )
			),
			menuItem("Age Information",
				  icon = icon("hourglass"),
				  tabName = "dashboard",
				  pickerInput(
					'AgeModelInput',
					'Clock Model:',
					choices = c("Mutation clock"='Mut', 
						  "Recombination clock"='Rec', 
						  "Joint clock"='Jnt'),
					selected = 'Jnt'
				  ),
				  numericInput("GenYears",
						 "Years per generation", 
						 value =0, min = 0, max = 40, step=1),
				  sliderInput(
					'AgeModeInput',
					'Age Mode (generations ago):',
					min = 0,
					max = 50000,
					value = c(0,30000)),
				  sliderInput(
					'qualScoreInput',
					h5('Quality Score:'),
					min = 0,
					max = 1,
					value = c(0.5, 1),
					step=0.05
				  ),
				  awesomeCheckbox(inputId = 'errorBars',
						  label="Show error bars",
						  value = FALSE,
						  status="primary")
			),
			br(),
			useShinyjs(),  # Set up shinyjs
			bsButton(inputId = "UpdateButt",
				class= "UpdateButt",
				  label = "Update", 
				  icon = icon("redo-alt"), 
				  style = "primary",
				  size='large',
				  width='85%',
				ignoreInit = FALSE)
		  )
		)
	),
	  
  
  
  # MAIN BODY ---------------------------------------------------------------
  
	dashboardBody(
		tags$head(
			tags$link(rel="icon",
				type = "image/gif",
				href="logoBW.png")
		),
		tags$head(
			tags$link(
				rel = "stylesheet", 
				type = "text/css", 
				href = "appStyle.css")
		),	 
		tabItems(
			tabItem(
				tabName = "selection",
				fluidRow(align="center",HTML('<h1><i class="fas fa-stream"></i> <b>Selection (iHS & nSL)</b></h1>')),
				fluidRow(
					div(
						id = "Selection_panel", 
						column(
							width = 12,
							h2("General Distribution"),
							uiOutput("box_general")
						),
						column(
							width = 12,
							h2("Distribution by population"),
							uiOutput("box_bypop")
						)
					)
				)
	  		),
			# Functional Characht
			tabItem(
				tabName = "functional",
				fluidRow(align="center",HTML('<h1><i class="fas fa-layer-group"></i> <b>Functional Description</b></h1>')),
				# Table functional
				fluidRow(
					div(
						id = "Functional_panel", 
						column(
							width = 12,
							h2("Variant Effect (SnpEff)"),
							uiOutput("box_SNPEFF")
						),
						column(
							width = 6,
							h2("Regulome DB"),
							uiOutput("box_ReguDB") 
						),
						column(
							 width = 6,
							h2("ClinVar"),
							uiOutput("box_ClinVar")
						),
						column(
							width = 6,
							h2("GWAS Catalog"),
							uiOutput("box_GWAScat")
						),
						column(
							width = 6,
							h2("DisGeNET"),
							uiOutput("box_DisGeNET")
						)
					)
				)
			),
			# Age information , 
			tabItem(
				tabName = "age",
				fluidRow(align="center",HTML('<h1><i class="fas fa-hourglass"></i> <b>Age Information</b></h1>')),
				fluidRow(
					div(
						id = "GEVA_panel", 
						column(
							width = 12,
							h2("Atlas of Variant Age"),
						uiOutput("box_geva")
						)
					)
		  		)
		  	),
			# iSAFE information
			tabItem(
				tabName = "isafe",
				fluidRow(align="center",HTML('<h1><i class="fas fa-asterisk"></i> <b>Favored Mutation (iSAFE)</b></h1>')),
				h2('General Distribution'),
				fluidRow(
					div(
						id = "iSAFE_panel", 
						column(
							width = 12,
							uiOutput("box_isafe")
						),
						column(
							width = 12,
							h2("Distribution by population"),
							"To reduce complexity, we only show those positions reporting an iSAFE score higher than 0.05. Drag the mouse arround to see if a position belongs to the top 0.01% of values.",
							uiOutput("scatter_isafe")
						)
					)
				)
			),
			# General Table 
			tabItem(
				tabName = "report",
				fluidRow(align="center",HTML('<h1><i class="fas fa-file-alt"></i> <b>Summary Report</b></h1>')),br(),
				fluidRow(
					align = "center",
					column(6,
						htmlOutput("DescriptionRegionTitle")
						),
					column(6,
						# add some buttons according to the region
						div(style="display:inline-block;width:49%;text-align: right;",uiOutput("pophumanscanLink")),
						div(style="display:inline-block;width:49%;text-align: right;",uiOutput("pophumanLink"))						
						)
					),
				br(),
				fluidRow(
					div(
						column(12,
							uiOutput("box_jbrows")
						)
					)
				),
				br(),br(),
				fluidRow(
					div(
						id = "report_panel", 
						column(
							width = 9,
							uiOutput("box_report")
						), 
						column(
							width =3,
							box(
								title = "DISPLAY OPTIONS",
								width = NULL,
								height = NULL,
								status = "primary",
								br(),
								"Options in this menu apply to this page only.\nIn the plot, iSAFE scores (y-axis) are represented for all SNVs across the genomic region of itnerest (x-axis). Color represents the strongest SNpEFF functional effect of each variant, and size represents its combined iHS + nSL value. The rest of information is displayed when dragging the mouse arround the plot.", 
								br(),br(), 
								h3("Filters:"),br(),
								pickerInput(
									inputId  = "Overview_Pop",
									label    = "Population:",
									choices  = list(
									AFR      = AFRpops,
									EUR      = EURpops,
									EAS      = EASpops,
									SAS      = SASpops),
									multiple = FALSE,
									selected = "CHS"
						 		),
								br(),
								pickerInput(
									'Overview_Clock',
									'Clock Model:',
									choices = c("Mutation clock"='Mut', 
										"Recombination clock"='Rec', 
										"Joint clock"='Jnt'),
									selected = 'Jnt'
								),
								br(),
								sliderInput(
									'Overview_AgeMode',
									'Age Mode (generations ago):',
									min = 55.444,
									max = 109013,
									value = c(55.444,2000)),
								br(),
								sliderInput(
									'Overview_qualScore',
									h5('Quality Score:'),
									min = 0,
									max = 1,
									value = c(0.75, 1),
									step=0.05
								),
								br(),br(),
								tags$div(align = 'center',
									bsButton(inputId = "ButtonNum2",
									label = "Refresh", 
									icon = icon("angle-double-left"),
									width= '100%',
									size='large',
									style = "primary")
								),
								br(),
								br()
							)
						)
					)
				),
				br(),
				br(),
				fluidRow(
					div(
						tags$head(tags$style(HTML('.info-box {min-height: 85px;} .info-box-icon {height: 85px; line-height: 85px;} .info-box-content {padding-top: 0px; padding-bottom: 0px;}'))),
						column(
							width=12,
							infoBoxOutput("info_isafe"),
							infoBoxOutput("info_iHS"),
							infoBoxOutput("info_nsl")
						),
						column(
							width=12,
							infoBoxOutput("info_snpeff"),
							infoBoxOutput("info_ReguDB"),
							infoBoxOutput("info_ClinVar")
						),
						column(
							width=12,
							infoBoxOutput("info_GWAScat"),
							infoBoxOutput("info_DGN")
						)
					)
				),
				br(),br(),
				fluidRow(
					div(
						id = "report_panel", 
						column(
							width = 12,
							uiOutput("box_prior")
						)
					)
				)
			),
			# Download
			tabItem(
				tabName="region_download",
				fluidRow(align="center",HTML('<h1><i class="fas fa-thumbtack"></i> <b>Current Region</b></h1>')),br(),
				fluidRow(
					column(1),
					column(8,
						h2("Region table"),
						fluidPage(
							withSpinner(
								DTOutput('dtDownload'),
								type  = 5,
								color = "#d33724",
								size  = 0.7
							)
						)
					),
					column(3,
						box(title = "DISPLAY OPTIONS",
								width = NULL,
								height = NULL,
								status = "primary",
								br(),
							# Filter AND/OR
							radioGroupButtons(
								inputId = "FilterType",
								label = "Combine activated filters in the FILTERS MENU with (boolean):",
								choices = c("AND","OR"),
								selected="OR",
								justified = FALSE
							),
							br(),
							pickerInput(
								inputId = "headersDownload",
								label = "Include the following columns (click the Update button in the FILTERS MENU first)", 
								choices = c('Position', 'rsid'),
								selected = NULL,
								multiple = TRUE,
								options = pickerOptions(
									dropupAuto=FALSE,
									virtualScroll=TRUE
								)
							),
							actionButton('selectAllheaders', 'Select All'),
							br(),br(),
							downloadButton(
								outputId = "bt_download",
								label    = "Download",
								icon     = icon("download"),
								style="width:100%;
								border: auto;
								color: white;
								background-color: #219086;
								padding: 15px 32px;
								text-align: center;
								display: inline-block;
								font-size: 18px;
								margin: 10px 5px;
								cursor: pointer;
								border-radius: 5px;"
							),
							br(),
							br()
							)
					)
				)
			),

			tabItem(
				tabName="batch_download",
				fluidRow(align="center",HTML('<h1><i class="fas fa-file-archive"></i> <b>Batch Download</b></h1>')),br(),
				fluidRow(
					column(1),
					column(10,
						box(
							width = NULL,
							height = NULL,
							align="center",
							HTML("<br>
								<p style='text-align:left;'>&emsp;Input a <b>CSV file</b> with one or more sets of coordinates to download PopHumanVar raw data.</p>
								<p style='text-align:left;'>&emsp;Format of the input file: <span style='background-color: #fcf3cf '>&#60;chr&#62;,&#60;start&#62;,&#60;end&#62;</span></p>
								<p style='text-align:left;'>&emsp;&emsp;Example: 2,109510927,109705828</p>
								<p style='text-align:left;'>&emsp;Data for a total maximum of <b>50 Mbp</b> can be downloaded at a time.</p>
								<p style='text-align:left;'>&emsp;Raw files for all the chromosomes are available from <a href='https://pophumanscan.uab.cat/data/phv/' target='_blank'>here</a></p>"),
							br(),
							fileInput("batch_regions", 
								"Choose CSV File",
								multiple = FALSE,
								accept = c("text/csv",".csv"),
								width = "60%"
							),
							uiOutput("hidden_download"),
							br()
						),
					column(1)
					),
				),
				fluidRow(
					column(1),
					column(11,
						withSpinner(
							DTOutput('batch_download'),
							type  = 5,
							color = "#d33724",
							size  = 0.7
							)
					)
				)
			),
		 			
			
			# Tutorial
			tabItem(tabName="tutorial",
				fluidRow(column(1),
					column(11,
						includeHTML("tutorial_text.html")
					)
				)
			),          
	  
			# About Us
			tabItem(
				tabName = "Information",
				fluidRow(
					align="center",
					HTML('<h1><i class="fas fa-info"></i> <b>About Us</b></h1>')),
				br(),
				fluidRow(column(1),
					column(10,
						box(
							width = NULL, 
							height = NULL,
							includeHTML("aboutUs_text.html"))),
					column(1)
				)
			)
		)
	)
)

