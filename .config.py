out_dir='/path/to/output/directory'
winston_address='winston address'
winston_port=11111
winston_address_cnmi='CNMI winston address'
winston_port_cnmi=11111
working_dir='/path/to/working/directory'


duration      = 600  # DON'T CHANGE THIS!
latency       = 30.0 # seconds between timestamps and end of data window
taper_val     = 5.0  # seconds to taper beginning and end of trace before filtering
f1		      = 0.5  # minimum frequency for bandpass filter
f2		      = 10.0 # maximum frequency for bandpass filter
window_length = 30   # window length for each calculation [seconds]
overlap       = 15   # seconds
min_chan      = 3
mcthresh      = 0.6  # where to draw the threshold line in MCCM plot
network       = 'AVO'

# Infrasound channels list
arrays=[
	dict({'Name':'Okmok',
		  'SCNL':[
					{'scnl':'OK01.BDF.AV.--'	, 'sta_lat': 53.41084	, 'sta_lon': -167.91431},
					{'scnl':'OK02.BDF.AV.--'	, 'sta_lat': 53.41002	, 'sta_lon': -167.91367},
					{'scnl':'OK03.BDF.AV.--'	, 'sta_lat': 53.40995	, 'sta_lon': -167.91505},
					{'scnl':'OK04.BDF.AV.--'	, 'sta_lat': 53.41031	, 'sta_lon': -167.91440}],
		'digouti':  (1/419430.0)/0.05,
		'volcano':[
					{'name': 'Bogoslof', 'v_lat': 53.9310,   'v_lon': -168.0360},
					{'name': 'Cleveland','v_lat': 52.8222,   'v_lon': -169.9464},
					{'name': 'Okmok',	 'v_lat': 53.428865, 'v_lon': -168.131632},
					{'name': 'Makushin', 'v_lat': 53.8900,   'v_lon': -166.9200}]}),



	dict({'Name':'Sand Point',
		  'SCNL':[
					{'scnl':'SDPI.BDF.AV.01'	, 'sta_lat': 55.34900	, 'sta_lon': -160.47640},
					{'scnl':'SDPI.BDF.AV.02'	, 'sta_lat': 55.34870	, 'sta_lon': -160.47683},
					{'scnl':'SDPI.BDF.AV.03'	, 'sta_lat': 55.34934	, 'sta_lon': -160.47732},
					{'scnl':'SDPI.BDF.AV.04'	, 'sta_lat': 55.34952	, 'sta_lon': -160.47661},
					{'scnl':'SDPI.BDF.AV.05'	, 'sta_lat': 55.34922	, 'sta_lon': -160.47650},
					{'scnl':'SDPI.BDF.AV.06'	, 'sta_lat': 55.34919	, 'sta_lon': -160.47710},
				 ],
	   'digouti': (1/419430.0)/(1e-2),
	   'volcano':[
	   				{'name': 'Pavlof',    'v_lat': 55.417622,'v_lon': -161.893669},
	   				{'name': 'Veniaminof','v_lat': 56.195825,'v_lon': -159.389536},
	   				{'name': 'Shishaldin','v_lat': 54.755856,'v_lon': -163.969961}]}),



	dict({'Name':'Adak',
		  'SCNL':[
					{'scnl':'ADKI.BDF.AV.01'	, 'sta_lat': 51.86190727	, 'sta_lon': -176.6438598},
					{'scnl':'ADKI.BDF.AV.02'	, 'sta_lat': 51.86324162	, 'sta_lon': -176.6435998},
					{'scnl':'ADKI.BDF.AV.03'	, 'sta_lat': 51.86226962	, 'sta_lon': -176.6446503},
					{'scnl':'ADKI.BDF.AV.04'	, 'sta_lat': 51.86246609	, 'sta_lon': -176.6457851},
					{'scnl':'ADKI.BDF.AV.05'	, 'sta_lat': 51.86326916	, 'sta_lon': -176.6461231},
					{'scnl':'ADKI.BDF.AV.06'	, 'sta_lat': 51.86157572	, 'sta_lon': -176.6469340}],
	   'digouti': (1/419430.0)/(0.05),
	   'volcano':[
				   {'name':	'Cleveland',   'v_lat': 52.8222,   'v_lon': -169.9464},
				   {'name':	'Great Sitkin','v_lat': 52.077282, 'v_lon': -176.131317},
				   {'name':	'Moffett',     'v_lat': 51.931876, 'v_lon': -176.740191}]}),




	dict({'Name':'Cleveland',
		  'SCNL':[
					{'scnl':'CLCO1.BDF.AV.--'	, 'sta_lat': 52.7864125000	, 'sta_lon': -169.7229250000},
					{'scnl':'CLCO2.BDF.AV.--'	, 'sta_lat': 52.7871866667	, 'sta_lon': -169.7244333330},
					{'scnl':'CLCO3.BDF.AV.--'	, 'sta_lat': 52.7874600000	, 'sta_lon': -169.7210166667},
					{'scnl':'CLCO4.BDF.AV.--'	, 'sta_lat': 52.7861266667	, 'sta_lon': -169.7203966667},
					{'scnl':'CLCO5.BDF.AV.--'	, 'sta_lat': 52.7851866667	, 'sta_lon': -169.7250066667},
				 ],
	   'digouti': (1/419430.0)/(1.62e-2),
	   'volcano':[
				   {'name':	'Cleveland',   'v_lat': 52.8222,   'v_lon': -169.9464},
				   {'name': 'Bogoslof', 'v_lat': 53.9310,   'v_lon': -168.0360}]}),


	dict({'Name':'Akutan',
  		  'SCNL':[
					{'scnl':'AKS.BDF.AV.--'	, 'sta_lat': 54.11050	, 'sta_lon': -165.69773},
					{'scnl':'AKS.BDG.AV.--'	, 'sta_lat': 54.11028	, 'sta_lon': -165.69618},
					{'scnl':'AKS.BDH.AV.--'	, 'sta_lat': 54.11105	, 'sta_lon': -165.69700},
					{'scnl':'AKS.BDI.AV.--'	, 'sta_lat': 54.11053	, 'sta_lon': -165.69683},
			    ],
		  'digouti': (1/419430.0)/0.05,
		  'volcano':[
				   {'name': 'Makushin', 'v_lat': 53.8900,   'v_lon': -166.9200},
				   {'name': 'Akutan',   'v_lat': 54.1300,   'v_lon': -165.9900},
				   {'name': 'Westdahl', 'v_lat': 54.5200,   'v_lon': -164.6500}]}),


	dict({'Name':'Saipan',
		  'SCNL':[
					{'scnl':'FLX2.BDF.MI.01'	, 'sta_lat': 15.23481	, 'sta_lon': 145.79219},
					{'scnl':'FLX2.BDF.MI.02'	, 'sta_lat': 15.23325	, 'sta_lon': 145.79242},
					{'scnl':'FLX2.BDF.MI.03'	, 'sta_lat': 15.23461	, 'sta_lon': 145.79097},
					{'scnl':'FLX2.BDF.MI.04'	, 'sta_lat': 15.23206	, 'sta_lon': 145.79389},
					{'scnl':'FLX2.BDF.MI.05'	, 'sta_lat': 15.23217	, 'sta_lon': 145.79119},
					{'scnl':'FLX2.BDF.MI.06'	, 'sta_lat': 15.23650	, 'sta_lon': 145.79053},
					# {'scnl':'DPS.BDF.MI.--'		, 'sta_lat': 15.24082	, 'sta_lon': 145.78909},
					# {'scnl':'HTSP.BDF.MI.--'	, 'sta_lat': 15.23166	, 'sta_lon': 145.79851},
					# {'scnl':'GOLF.BDF.MI.--'	, 'sta_lat': 15.22511	, 'sta_lon': 145.78625}, # still has IML sensor
					# {'scnl':'FLX.BDF.MI.--'		, 'sta_lat': 15.23389	, 'sta_lon': 145.79172},
				 ],
	   'digouti': (1/419430.0)/(1e-2),
	   'volcano':[
				   {'name':	'Anatahan?', 'v_lat': 18.141555, 'v_lon': 145.786260}]}),



	dict({'Name':'Sarigan',
		  'SCNL':[
					{'scnl':'SRN1.BDF.MI.--'	, 'sta_lat': 16.70035	, 'sta_lon': 145.77809},
					{'scnl':'SRN2.BDF.MI.--'	, 'sta_lat': 16.69977	, 'sta_lon': 145.77798},
					{'scnl':'SRN3.BDF.MI.--'	, 'sta_lat': 16.69960	, 'sta_lon': 145.77841},
					{'scnl':'SRN4.BDF.MI.--'	, 'sta_lat': 16.69990	, 'sta_lon': 145.77876},
					# {'scnl':'SRN5.BDF.MI.--'	, 'sta_lat': 16.70036	, 'sta_lon': 145.77872},
					{'scnl':'SRN6.BDF.MI.--'	, 'sta_lat': 16.70003	, 'sta_lon': 145.77838},
				 ],
	   'digouti': (1/419430.0)/(1e-2),
	   'volcano':[
				   {'name':	'Pagan?',   'v_lat': 18.141555,	'v_lon': 145.786260},
				   {'name':	'Anatahan?','v_lat': 16.351323,	'v_lon': 145.687602}]}),

]
