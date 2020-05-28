out_web_dir = '/www/avosouth.wr.usgs.gov/htdocs/infrasound'

#######
# TODO: Determine if we need these
out_valve_dir = '/lamp/valve3/seismic/infrasound/raw'
working_dir = '/usr/local/ipensive'
winston_address = 'hvo-wws.wr.usgs.gov'
winston_port = 16022
##########

duration = 600  # DON'T CHANGE THIS!
latency = 30.0  # seconds between timestamps and end of data window
taper_val = 5.0  # seconds to taper beginning and end of trace before filtering
f1 = 0.5  # minimum frequency for bandpass filter
f2 = 5.0  # maximum frequency for bandpass filter
window_length = 30   # window length for each calculation [seconds]
overlap = .5   # percent
min_chan = 3
mcthresh = 0.5  # where to draw the threshold line in MCCM plot

# Infrasound channels list
arrays = [
    # AVO Arrays
    {
        'network': 'AV',
        'display name': 'AVO',
        'arrays': [
            {'Name': 'Dillingham',
             'id': 'DLL',
             'Alpha': 0.75,
             'channel': 'HDF',
             },
            {'Name': 'Sand Point',
             'id': 'SDPI',
             'Alpha': 0.75,
             'channel': 'HDF',
             },
            {'Name': 'Akutan',
             'id': 'AKS',
             'Alpha': 0.75,
             'channel': 'BD*',
             },
            {'Name': 'Okmok',
             'id': 'OKIF',
             'Alpha': 0.75,
             'channel': 'HDF',
             },
            {'Name': 'Cleveland',
             'id': 'CLC*',
             'Alpha': 0.75,
             'channel': 'BDF',

             },
            {'Name': 'Adak',
             'id': 'ADKI',
             'Alpha': 0.75,
             'channel': 'HDF',
             },
        ]
    },
    # CNMI arrays
    #     {
    #         'network': 'CNMI',
    #         'arrays': [
    #             {'Name': 'Saipan',
    #              'network': 'CNMI',
    #              },
    #             {'Name': 'Sarigan',
    #              'network': 'CNMI',
    #              },
    #         ],
    #     }
]

# Infrasound channels list
# arrays=[
# 	dict({'Name':'Ahua',
# 		  'SCNL':[
# 					{'scnl':'AHUD.BDF.HV.01'	, 'sta_lat': 19.371500	, 'sta_lon': -155.263400},
# 					{'scnl':'AHUD.BDF.HV.02'	, 'sta_lat': 19.371935  , 'sta_lon': -155.263505},
# 					{'scnl':'AHUD.BDF.HV.03'	, 'sta_lat': 19.371370  , 'sta_lon': -155.262850},
# 					{'scnl':'AHUD.BDF.HV.04'	, 'sta_lat': 19.371195  , 'sta_lon': -155.263880}],
# 		'digouti':  (1/6291.3),
# 		'volcano':[
# 					{'name': 'Halemaumau', 'v_lat': 19.405360,   'v_lon': -155.280972},
# 					{'name': 'PuuOo','v_lat': 19.388485,   'v_lon': -155.106042},
# 					{'name': 'MLsummit',	 'v_lat': 19.472753, 'v_lon': -155.590970},
# 					{'name': 'MLswrz', 'v_lat': 19.346941,   'v_lon': -155.667472}]}),
#
# 	dict({'Name':'Ainapo',
# 		  'SCNL':[
# 					{'scnl':'AIND.BDF.HV.01'	, 'sta_lat': 19.372062	, 'sta_lon': -155.457170},
# 					{'scnl':'AIND.BDF.HV.02'	, 'sta_lat': 19.372244  , 'sta_lon': -155.457710},
# 					{'scnl':'AIND.BDF.HV.03'	, 'sta_lat': 19.372158  , 'sta_lon': -155.456800},
# 					{'scnl':'AIND.BDF.HV.04'	, 'sta_lat': 19.371567  , 'sta_lon': -155.457410}],
# 		'digouti':  (1/6291.3),
# 		'volcano':[
# 					{'name': 'Halemaumau', 'v_lat': 19.405360,   'v_lon': -155.280972},
# 					{'name': 'PuuOo','v_lat': 19.388485,   'v_lon': -155.106042},
# 					{'name': 'MLsummit',	 'v_lat': 19.472753, 'v_lon': -155.590970},
# 					{'name': 'MLnerz',	 'v_lat': 19.525019, 'v_lon': -155.486837},
# 					{'name': 'MLswrz', 'v_lat': 19.346941,   'v_lon': -155.667472}]}),
#
# 	dict({'Name':'Nanawale',
# 		'SCNL':[
# 					{'scnl':'WALE.HDF.HV.01'        , 'sta_lat': 19.4942275  , 'sta_lon': -154.91405},
# 					{'scnl':'WALE.HDF.HV.02'        , 'sta_lat': 19.4943302  , 'sta_lon': -154.91453},
# 					{'scnl':'WALE.HDF.HV.03'        , 'sta_lat': 19.4946066  , 'sta_lon': -154.91475},
# 					{'scnl':'WALE.HDF.HV.04'        , 'sta_lat': 19.4948083  , 'sta_lon': -154.91434},
# 					{'scnl':'WALE.HDF.HV.05'        , 'sta_lat': 19.4949789  , 'sta_lon': -154.91403},
# 					{'scnl':'WALE.HDF.HV.06'        , 'sta_lat': 19.4946475  , 'sta_lon': -154.91388}],
# 		'digouti':  (2.6686e-4),
#         'volcano': [
# 					{'name': 'Fissure 8', 'v_lat': 19.460222,   'v_lon': -154.908751},
# 					{'name': 'Fissure 10', 'v_lat': 19.455894,   'v_lon': -154.919373},
# 					{'name': 'Fissure 22', 'v_lat': 19.474969,   'v_lon': -154.884739},
# 					{'name': 'Fissure 16/18', 'v_lat': 19.478863,   'v_lon': -154.876357},
# 					{'name': 'Kapoho Cone',  'v_lat': 19.500132, 'v_lon': -154.839707},
# 					{'name': 'Cracks on 130', 'v_lat': 19.446489,   'v_lon': -154.941694},
# 					{'name': 'PuuOo','v_lat': 19.388485,   'v_lon': -155.106042},
# 					{'name': 'Ocean Entry', 'v_lat': 19.498129, 'v_lon': -154.814505}]}),
#
# #	dict({'Name':'North Pit',
# #		  'SCNL':[
# #					{'scnl':'NPT.HDF.HV.01'	, 'sta_lat': 19.411640	, 'sta_lon': -155.28107},
# #					#{'scnl':'NPT.HDF.HV.02'	, 'sta_lat': 19.412110  , 'sta_lon': -155.28110},
# #					{'scnl':'NPT.HDF.HV.03'	, 'sta_lat': 19.411972  , 'sta_lon': -155.280642},
# #					{'scnl':'NPT.HDF.HV.04'	, 'sta_lat': 19.412015  , 'sta_lon': -155.280997}],
# #		'digouti':  (1/906.0),
# #		'volcano':[
# #					{'name': 'Halemaumau', 'v_lat': 19.405360,   'v_lon': -155.280972},
# #					{'name': 'PuuOo','v_lat': 19.388485,   'v_lon': -155.106042},
# #					{'name': 'MLsummit',	 'v_lat': 19.472753, 'v_lon': -155.590970},
# #					{'name': 'MLnerz',	 'v_lat': 19.525019, 'v_lon': -155.486837},
# #					{'name': 'MLswrz', 'v_lat': 19.346941,   'v_lon': -155.667472}]}),
#
# 	dict({'Name':'ERZ2',
# 	  'SCNL':[
# 					{'scnl':'ERZ2.HDF.HV.01'	, 'sta_lat': 19.474760	, 'sta_lon': -154.833630},
# 					{'scnl':'ERZ2.HDF.HV.02'	, 'sta_lat': 19.474380  , 'sta_lon': -154.834030},
# 					{'scnl':'ERZ2.HDF.HV.03'	, 'sta_lat': 19.475210  , 'sta_lon': -154.833980}],
# 		'digouti':  (1/6291.3),
# 		'volcano':[
# 					{'name': 'Fissure 17', 'v_lat': 19.483665,   'v_lon': -154.871967},
# 					{'name': 'Kapoho Cone',	 'v_lat': 19.500132, 'v_lon': -154.839707},
#                                         {'name': 'Fissure 8', 'v_lat': 19.460222,   'v_lon': -154.908751},
# 					{'name': 'Kapoho Ocean Entry', 'v_lat': 19.480619, 'v_lon': -154.820470}]}),
#
# ]
