# SYSTEM PARAMETERS
# HOSTNAME	  = '137.227.224.220'
# PORT		  = 16002
HOSTNAME	  = 'pubavo1.wr.usgs.gov'
PORT		  = 16022
OUT_WEB_DIR   = '/www/avosouth.wr.usgs.gov/htdocs/infrasound'	# html & image ouptut directory
OUT_ASCII_DIR = '/www/avosouth.wr.usgs.gov/htdocs/infrasound/ascii_output'	# ascii output directory (delete if undesired)
LOGS_DIR	  = ''

# DEFAUlT PROCESSING PARAMETERS
DURATION      = 600  # DON'T CHANGE THIS!
LATENCY       = 60	 # how long to pause and wait for data latency to catch up before processing (seconds)
EXTRA_PAUSE   = 0 	 # extra pause in case some arrays are extra latent (Dillingham)
FREQMIN		  = 0.8  # minimum frequency
FREQMAX		  = 5.0  # maximum frequency
TAPER	      = 5.0  # seconds to taper beginning and end of trace before filtering
WINDOW_LENGTH = 30   # window length for each calculation [seconds]
OVERLAP       = 15   # amount of overlap between windows [seconds]
MIN_CHAN      = 3	 # minimum number of channels needed to perform inversion

# DEFAULT PLOTTING PARAMETERS
MCTHRESH      = 0.6  # where to draw the threshold line in MCCM plot
AZ_MIN		  = 0	 # minimum azimuth to plot
AZ_MAX		  = 360	 # maximum azimuth to plot
VEL_MIN		  = 0.25 # minimum sound speed range
VEL_MAX		  = 0.45 # maximum sound speed range
ARRAY_LABEL	  = 'Infrasound'

# Infrasound channels list
NETWORKS=[
	{
		'Name' : 'AVO',
		'ARRAYS' : [

			dict({'Name':'Kenai',
				  'SCNL':[
							{'scnl':'KENI.HDF.AV.01'	, 'sta_lat': 60.6413700	, 'sta_lon': -151.070200},
							{'scnl':'KENI.HDF.AV.02'	, 'sta_lat': 60.6404567 , 'sta_lon': -151.070330},
							{'scnl':'KENI.HDF.AV.03'	, 'sta_lat': 60.6406033	, 'sta_lon': -151.072020},
							{'scnl':'KENI.HDF.AV.04'	, 'sta_lat': 60.6412000	, 'sta_lon': -151.073000},
							{'scnl':'KENI.HDF.AV.05'	, 'sta_lat': 60.6415300	, 'sta_lon': -151.072000},
							{'scnl':'KENI.HDF.AV.06'	, 'sta_lat': 60.6409167 , 'sta_lon': -151.071170},
						],
				'digouti': (1/419430.0)/(0.0275),
				'volcano':[
							{'name': 'Wrangell',  'v_lat': 62.00572,  'v_lon': -144.01935},
							{'name': 'Spurr',     'v_lat': 61.29897,  'v_lon': -152.25122},
							{'name': 'Redoubt',   'v_lat': 60.48576,  'v_lon': -152.74282},
							{'name': 'Iliamna',   'v_lat': 60.03220,  'v_lon': -153.09002},
							{'name': 'Augustine', 'v_lat': 59.36107,  'v_lon': -153.42938},
						  ],
				'AZ_MIN': 200,
				'AZ_MAX': 80
			}),

			dict({'Name':'Sand Point',
				  'SCNL':[
							{'scnl':'SDPI.HDF.AV.01'	, 'sta_lat': 55.34900	, 'sta_lon': -160.47640},
							{'scnl':'SDPI.HDF.AV.02'	, 'sta_lat': 55.34870	, 'sta_lon': -160.47683},
							# {'scnl':'SDPI.HDF.AV.03'	, 'sta_lat': 55.34934	, 'sta_lon': -160.47732},
							{'scnl':'SDPI.HDF.AV.04'	, 'sta_lat': 55.34952	, 'sta_lon': -160.47661},
							{'scnl':'SDPI.HDF.AV.05'	, 'sta_lat': 55.34922	, 'sta_lon': -160.47650},
							{'scnl':'SDPI.HDF.AV.06'	, 'sta_lat': 55.34919	, 'sta_lon': -160.47710},
						 ],
				'digouti': (1/419430.0)/0.0275,
				'volcano':[
			   				{'name': 'Pavlof',    'v_lat': 55.417622,'v_lon': -161.893669},
			   				{'name': 'Veniaminof','v_lat': 56.195825,'v_lon': -159.389536},
			   				{'name': 'Shishaldin','v_lat': 54.755856,'v_lon': -163.969961},
			   				{'name': 'Trident',	 'v_lat': 58.234389,	'v_lon': -155.103988},
		   				  ],
		   		'AZ_MIN': 240,
				'AZ_MAX': 60,
				'FREQMIN': 1,
			}),

			dict({'Name':'Akutan',
					 'SCNL':[
								{'scnl':'AKS.HDF.AV.01'	, 'sta_lat': 54.11048	, 'sta_lon': -165.69774},
								{'scnl':'AKS.HDF.AV.02'	, 'sta_lat': 54.11105	, 'sta_lon': -165.69705},
								{'scnl':'AKS.HDF.AV.03'	, 'sta_lat': 54.11028	, 'sta_lon': -165.69616},
								{'scnl':'AKS.HDF.AV.04'	, 'sta_lat': 54.11051	, 'sta_lon': -165.69681},
								# {'scnl':'AKS.HDF.AV.04'	, 'sta_lat': 54.11051	, 'sta_lon': -165.69683},
								{'scnl':'AKS.HDF.AV.06'	, 'sta_lat': 54.11005	, 'sta_lon': -165.69720},
							],
				'digouti': (1/400000)/0.0275,	# convert counts to Pressure in Pa (Centaur + Chaparral Vx2 mics)
				'volcano':[
							{'name': 'Makushin', 'v_lat': 53.8900,   'v_lon': -166.9200},
							{'name': 'Akutan',   'v_lat': 54.1300,   'v_lon': -165.9900},
							# {'name': 'Westdahl', 'v_lat': 54.5200,   'v_lon': -164.6500}
			   				{'name': 'Shishaldin','v_lat': 54.755856,'v_lon': -163.969961},
						  ],
				'FREQMIN': 1,
			}),

			dict({'Name':'Okmok',
				  'SCNL':[
							{'scnl':'OKIF.HDF.AV.01'	, 'sta_lat': 53.41083004	, 'sta_lon': -167.91426701},
							{'scnl':'OKIF.HDF.AV.02'	, 'sta_lat': 53.41001901	, 'sta_lon': -167.91366301},
							{'scnl':'OKIF.HDF.AV.03'	, 'sta_lat': 53.40998297	, 'sta_lon': -167.91499598},
							{'scnl':'OKIF.HDF.AV.04'	, 'sta_lat': 53.41029796	, 'sta_lon': -167.91431696},
							{'scnl':'OKIF.HDF.AV.05'	, 'sta_lat': 53.41038496	, 'sta_lon': -167.91331901},
							{'scnl':'OKIF.HDF.AV.06'	, 'sta_lat': 53.41045604	, 'sta_lon': -167.91544802},
						],
				'digouti': 1/(419430*.03),	# convert counts to Pressure in Pa (Q330 + Chaparral 64Vx original)
				'volcano':[
							{'name': 'Bogoslof', 'v_lat': 53.9310,   'v_lon': -168.0360},
							{'name': 'Cleveland','v_lat': 52.8222,   'v_lon': -169.9464},
							{'name': 'Okmok',	 'v_lat': 53.428865, 'v_lon': -168.131632},
							{'name': 'Makushin', 'v_lat': 53.8900,   'v_lon': -166.9200},
			   				{'name': 'Shishaldin','v_lat': 54.755856,'v_lon': -163.969961},

						  ],
			}),

			dict({'Name':'Cleveland',
				  'SCNL':[
							{'scnl':'CLCO.HDF.AV.01'	, 'sta_lat': 52.7864125000	, 'sta_lon': -169.7229250000},
							{'scnl':'CLCO.HDF.AV.02'	, 'sta_lat': 52.7871866667	, 'sta_lon': -169.7244333330},
							{'scnl':'CLCO.HDF.AV.03'	, 'sta_lat': 52.7874600000	, 'sta_lon': -169.7210166667},
							{'scnl':'CLCO.HDF.AV.04'	, 'sta_lat': 52.7861266667	, 'sta_lon': -169.7203966667},
							{'scnl':'CLCO.HDF.AV.05'	, 'sta_lat': 52.7851866667	, 'sta_lon': -169.7250066667},
							{'scnl':'CLCO.HDF.AV.06'	, 'sta_lat': 52.785861	, 'sta_lon': -169.723179},
						 ],
				'digouti': (1/419430.0)/(0.0275),
				'volcano':[
							{'name': 'Cleveland',   'v_lat': 52.8222,   'v_lon': -169.9464},
							{'name': 'Bogoslof', 'v_lat': 53.9310,   'v_lon': -168.0360},
			   				{'name': 'Shishaldin','v_lat': 54.755856,'v_lon': -163.969961},
						  ],
			}),

			dict({'Name':'Adak',
				  'SCNL':[
							{'scnl':'ADKI.HDF.AV.01'	, 'sta_lat': 51.86190727	, 'sta_lon': -176.6438598},
							{'scnl':'ADKI.HDF.AV.02'	, 'sta_lat': 51.86324162	, 'sta_lon': -176.6435998},
							{'scnl':'ADKI.HDF.AV.03'	, 'sta_lat': 51.86226962	, 'sta_lon': -176.6446503},
							{'scnl':'ADKI.HDF.AV.04'	, 'sta_lat': 51.86246609	, 'sta_lon': -176.6457851},
							{'scnl':'ADKI.HDF.AV.05'	, 'sta_lat': 51.86326916	, 'sta_lon': -176.6461231},
							# {'scnl':'ADKI.HDF.AV.06'	, 'sta_lat': 51.86157572	, 'sta_lon': -176.6469340},
						 ],
				'digouti': (1/419430.0)/(0.05),
				'volcano':[
							{'name':	'Cleveland',   'v_lat': 52.8222,   'v_lon': -169.9464},
							{'name':	'Great Sitkin','v_lat': 52.077282, 'v_lon': -176.131317},
							{'name':	'Moffett',     'v_lat': 51.931876, 'v_lon': -176.740191},
							{'name': 'Semisopochnoi',	 'v_lat': 51.947,	'v_lon': 179.623},
			   				{'name': 'Shishaldin','v_lat': 54.755856,'v_lon': -163.969961},
						  ],
			}),

			dict({'Name':'Amchitka',
				  'SCNL':[
							{'scnl':'AMKA.HDF.AV.01'	, 'sta_lat': 51.379130	, 'sta_lon': 179.30130},
							{'scnl':'AMKA.HDF.AV.02'	, 'sta_lat': 51.378489	, 'sta_lon': 179.30183},
							{'scnl':'AMKA.HDF.AV.03'	, 'sta_lat': 51.378105	, 'sta_lon': 179.301225},
							{'scnl':'AMKA.HDF.AV.04'	, 'sta_lat': 51.37831	, 'sta_lon': 179.30028},
							{'scnl':'AMKA.HDF.AV.05'	, 'sta_lat': 51.379055	, 'sta_lon': 179.30026},
							# {'scnl':'AMKA.HDF.AV.06'	, 'sta_lat': 51.37871	, 'sta_lon': 179.30093},
						],
				'digouti': (1/400000)/0.0275,
				'volcano':[
			   			   {'name':	'Little Sitkin','v_lat': 51.95,	'v_lon': 178.543},
						   {'name': 'Semisopochnoi','v_lat': 51.947,'v_lon': 179.623},
						   {'name': 'Gareloi',	    'v_lat': 51.79,	'v_lon':-178.794},
						 ],
				'AZ_MIN': 270,
				'AZ_MAX': 90,
			}),

			dict({'Name':'Dillingham',
				  'SCNL':[
							{'scnl':'DLL.HDF.AV.01'	, 'sta_lat': 59.13988781	, 'sta_lon': -158.6209290},
							{'scnl':'DLL.HDF.AV.02'	, 'sta_lat': 59.13620003	, 'sta_lon': -158.6053376},
							{'scnl':'DLL.HDF.AV.03'	, 'sta_lat': 59.12904044	, 'sta_lon': -158.6146964},
							{'scnl':'DLL.HDF.AV.04'	, 'sta_lat': 59.13602776	, 'sta_lon': -158.6142354},
							{'scnl':'DLL.HDF.AV.05'	, 'sta_lat': 59.13509488	, 'sta_lon': -158.6136803},
							{'scnl':'DLL.HDF.AV.06'	, 'sta_lat': 59.13532733	, 'sta_lon': -158.6155527},
						 ],
				'digouti': 1/33505.5968,   # convert counts to Pressure in Pa 01-Jun-2021
				'volcano':[
							{'name': 'Bogoslof', 'v_lat': 53.9310,   'v_lon': -168.0360},
							{'name': 'Veniaminof','v_lat': 56.195825,'v_lon': -159.389536},
							{'name': 'Semisopochnoi',	 'v_lat': 51.947,	'v_lon': 179.623},
							{'name': 'Trident',	 'v_lat': 58.234389,	'v_lon': -155.103988},
			   				{'name': 'Shishaldin','v_lat': 54.755856,'v_lon': -163.969961},
						  ],
				'EXTRA_PAUSE': 90,
			}),			
		]
	},

	{
		'Name' : 'CNMI',
		'ARRAYS' : [

			dict({'Name':'Saipan',
				  'SCNL':[
							{'scnl':'FLX.HDF.MI.01'	, 'sta_lat': 15.23388	, 'sta_lon': 145.79172},
							{'scnl':'FLX.HDF.MI.02'	, 'sta_lat': 15.23320	, 'sta_lon': 145.79234},
							{'scnl':'FLX.HDF.MI.03'	, 'sta_lat': 15.23475	, 'sta_lon': 145.79216},
							# {'scnl':'FLX.HDF.MI.04'	, 'sta_lat': 15.23206	, 'sta_lon': 145.79389},
							{'scnl':'FLX.HDF.MI.05'	, 'sta_lat': 15.23215	, 'sta_lon': 145.79112},
							{'scnl':'FLX.HDF.MI.06'	, 'sta_lat': 15.23647	, 'sta_lon': 145.79045},
						 ],
				'digouti': (1/400000)/0.0275, #new sensors and locations installed Feb 20, 2024
				'volcano':[
							{'name': 'Anatahan', 'v_lat': 18.141555, 'v_lon': 145.786260},
							{'name':	'Pagan', 'v_lat': 18.141555, 'v_lon': 145.786260},
						  ],
				'AZ_MIN' : 205,
				'AZ_MAX' :   15,
			}),

			# dict({'Name':'Sarigan',
			# 	  'SCNL':[
			# 				{'scnl':'SRN1.BDF.MI.--'	, 'sta_lat': 16.70035	, 'sta_lon': 145.77809},
			# 				{'scnl':'SRN2.BDF.MI.--'	, 'sta_lat': 16.69977	, 'sta_lon': 145.77798},
			# 				{'scnl':'SRN3.BDF.MI.--'	, 'sta_lat': 16.69960	, 'sta_lon': 145.77841},
			# 				{'scnl':'SRN4.BDF.MI.--'	, 'sta_lat': 16.69990	, 'sta_lon': 145.77876},
			# 				# {'scnl':'SRN5.BDF.MI.--'	, 'sta_lat': 16.70036	, 'sta_lon': 145.77872},
			# 				{'scnl':'SRN6.BDF.MI.--'	, 'sta_lat': 16.70003	, 'sta_lon': 145.77838},
			# 			 ],
			# 	'digouti': (1/419430.0)/(1e-2),
			# 	'volcano':[
			# 				{'name':	'Pagan?',   'v_lat': 18.141555,	'v_lon': 145.786260},
			# 				{'name':	'Anatahan?','v_lat': 16.351323,	'v_lon': 145.687602},
			# 			  ],
			# }),

			dict({'Name':'Wake Island North',
				  'SCNL':[
							{'scnl':'H11N1.EDH.IM.--'	, 'sta_lat': 19.71356	, 'sta_lon':  166.891083},
							{'scnl':'H11N2.EDH.IM.--'	, 'sta_lat': 19.730801	, 'sta_lon':  166.897675},
							{'scnl':'H11N3.EDH.IM.--'	, 'sta_lat': 19.71718	, 'sta_lon':  166.909988},
						 ],

				'HOSTNAME' : 'IRIS',
			 	'digouti': 1/1855.04,
				'volcano':[
							{'name':	'Saipan',	'v_lat': 15.19,  'v_lon': 145.74},
							{'name':	'Anatahan', 'v_lat': 16.35,  'v_lon': 145.67},
							# {'name':	'Sarigan', 	'v_lat': 16.708, 'v_lon': 145.78},
							{'name':	'Pagan', 	'v_lat': 18.13,  'v_lon': 145.80},
							{'name':	'Ahyi', 	'v_lat': 20.42,  'v_lon': 145.03},
							# {'name':	'FdP', 		'v_lat': 20.546, 'v_lon': 144.893}
						  ],
				'FREQMIN': 1,
				'FREQMAX': 10,
				'VEL_MIN': 1.4,
				'VEL_MAX': 1.6,
				'AZ_MIN' : 250,
				'AZ_MAX' : 290,
				'ARRAY_LABEL' : 'Hydroacoustic'
			}),

			dict({'Name':'Wake Island South',
				  'SCNL':[
							{'scnl':'H11S1.EDH.IM.--'	, 'sta_lat': 18.50827	, 'sta_lon':  166.700272},
							{'scnl':'H11S2.EDH.IM.--'	, 'sta_lat': 18.49082	, 'sta_lon':  166.705002},
							{'scnl':'H11S3.EDH.IM.--'	, 'sta_lat': 18.49568	, 'sta_lon':  166.686462},
						 ],
 
			 	'HOSTNAME' : 'IRIS',
				'digouti': 1/1860.86,
				'volcano':[
							{'name':	'Saipan',	'v_lat': 15.19,  'v_lon': 145.74},
							{'name':	'Anatahan', 'v_lat': 16.35,  'v_lon': 145.67},
							# {'name':	'Sarigan', 	'v_lat': 16.708, 'v_lon': 145.78},
							{'name':	'Pagan', 	'v_lat': 18.13,  'v_lon': 145.80},
							{'name':	'Ahyi', 	'v_lat': 20.42,  'v_lon': 145.03},
							# {'name':	'FdP', 		'v_lat': 20.546, 'v_lon': 144.893}
						  ],
				'FREQMIN': 1,
				'FREQMAX': 10,
				'VEL_MIN': 1.4,
				'VEL_MAX': 1.6,
				'AZ_MIN' : 250,
				'AZ_MAX' : 290,	
				'ARRAY_LABEL' : 'Hydroacoustic'			
			}),
		]
	}
]
