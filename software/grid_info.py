valid_regions = ['polstere-wide-nh', 'polstere-nh', 'polstere-arc',
                 'polstere-can', 'polstere-ala', 'polstere-eur',
                 'polstere-sib', 'polstere-bfg', 'polstere-fra',
                 'polstere-wide-sh', 'polstere-sh', 'polstere-wed',
                 'polstere-dml', 'polstere-ross', 'polstere-lross',
                 'polstere-bell', 'polstere-wil',
                 'ease-nh', 'ease-nh-wide', 'ease-nh-very-wide',
                 'ease-eur', 'ease-eur-trim', 'ease-arc', 'ease-can',
                 'ease-sh', 'ease-sh-wide', 'ease-sh-very-wide',
                 'ease-wed', 'ease-ross', 'ease-bell', 'ease-wil']

def region_params(region):
    '''Setting up parameters according to the plotting region'''

    rp = {}
    # Defaults
    rp['scale'] = 1
    rp['skip'] = 1

    if region == 'polstere-wide-nh':
        rp['lllat'] = 22.0
        rp['lllon'] = -80.0
        rp['urlat'] = 22.0
        rp['urlon'] = 100.0
        rp['long_name'] = 'Wide Northern Hemisphere'
        rp['scale'] = 5
        rp['skip']  = 3
    elif region == 'polstere-nh':
        rp['lllat'] = 34.9407
        rp['lllon'] = -80.5377
        rp['urlat'] = 32.4423
        rp['urlon'] = 102.771
        rp['long_name'] = 'Northern Hemisphere'
        rp['scale'] = 5
        rp['skip']  = 3
    elif region == 'polstere-arc':
        rp['lllat'] = 63.98820
        rp['lllon'] = -105.75117
        rp['urlat'] = 60.80976
        rp['urlon'] = 88.39744
        rp['long_name'] = 'Arctic Ocean'
        rp['scale'] = 8
        rp['labelpos'] = 0.85
    elif region == 'polstere-can':
        rp['lllat'] = 51.
        rp['lllon'] = -84.
        rp['urlat'] = 88.
        rp['urlon'] = -30
        rp['long_name'] = 'West Greenland and Canada'
        rp['scale'] = 8
        rp['logo_left'] = 0.25
        rp['logo_bottom'] = 0.44
        rp['labelpos'] = 0.76
    elif region == 'polstere-ala':
        rp['lllat'] = 65.
        rp['lllon'] = -130.
        rp['urlat'] = 75.
        rp['urlon'] = 150
        rp['long_name'] = 'Alaska'
        rp['scale'] = 1
    elif region == 'polstere-eur':
        rp['lllat'] = 78.
        rp['lllon'] = -55.
        rp['urlat'] = 69.5
        rp['urlon'] = 63.
        rp['long_name'] = 'European Arctic'
        rp['scale'] = 2
    elif region == 'polstere-sib':
        rp['lllat'] = 75.
        rp['lllon'] = -145.
        rp['urlat'] = 68.
        rp['urlon'] = 110.
        rp['long_name'] = 'Siberian Arctic'
        rp['scale'] = 2
    elif region == 'polstere-bfg':
        rp['lllat'] = 67.
        rp['lllon'] = -120.
        rp['urlat'] = 80.
        rp['urlon'] = 170.
        rp['long_name'] = 'Beaufort Sea'
        rp['scale'] = 1
    elif region == 'polstere-fra':
        rp['lllat'] = 75.
        rp['lllon'] = -30
        rp['urlat'] = 81.
        rp['urlon'] = 32.
        rp['long_name'] = 'Fram Strait'
        rp['scale'] = 1
    elif region == 'polstere-wide-sh':
        rp['lllat'] = -15.000
        rp['lllon'] = -135.000
        rp['urlat'] = -15.000
        rp['urlon'] = 45.000
        rp['long_name'] = 'Wide Southern Hemisphere'
        rp['scale'] = 5
        rp['skip']  = 3
    elif region == 'polstere-sh':
        rp['lllat'] = -42.2681
        rp['lllon'] = -135.000
        rp['urlat'] = -40.1782
        rp['urlon'] = 42.357
        rp['long_name'] = 'Southern Hemisphere'
        rp['scale'] = 5
        rp['skip']  = 3
    elif region == 'polstere-wed':
        rp['lllat'] = -66.0
        rp['lllon'] = -75.0
        rp['urlat'] = -52.5
        rp['urlon'] = 20.0
        rp['long_name'] = 'Weddell and Lazarev Sea'
        rp['scale'] = 6
        rp['logo_left'] = 0.63
        rp['logo_bottom'] = 0.08
    elif region == 'polstere-dml':
        rp['lllat'] = -75.0
        rp['lllon'] = -20.0
        rp['urlat'] = -50.0
        rp['urlon'] = 20.0
        rp['long_name'] = 'Dronning Maud Land'
        rp['scale'] = 6
    elif region == 'polstere-ross' or region == 'polstere-lross':
        if region == 'polstere-ross':
            rp['lllat'] = -55.0
            rp['lllon'] =-150.0
        else:
            rp['lllat'] = -50.0
            rp['lllon'] =-135.0
        rp['urlat'] = -80.0
        rp['urlon'] = 150.0
        rp['long_name'] = 'Ross Sea'
        rp['scale'] = 6
        rp['logo_left'] = 0.65
        rp['logo_bottom'] = 0.75
    elif region == 'polstere-bell':
        rp['lllat'] = -60.0
        rp['lllon'] = -100.0
        rp['urlat'] = -68.0
        rp['urlon'] = -45.0
        rp['long_name'] = 'Bellinghausen Sea'
        rp['scale'] = 6
        rp['labelpos'] = 0.72
        rp['logo_left'] = 0.32
        rp['logo_bottom'] = 0.08
    elif region == 'polstere-wil':
        rp['lllat'] = -55.0
        rp['lllon'] = 180.0
        rp['urlat'] = -45.0
        rp['urlon'] = 70.0
        rp['long_name'] = 'Wilkes Land Coast'
        rp['scale'] = 6
        rp['labelpos'] = 0.72
        rp['logo_left'] = 0.32
        rp['logo_bottom'] = 0.75
    elif region == 'ease-nh':
        rp['lllat'] = 60.
        rp['lllon'] = -140.
        rp['urlat'] = 65.
        rp['urlon'] = 55.
        rp['long_name'] = 'Northern Hemisphere EASE2'
        rp['scale'] = 5
        rp['skip']  = 3
    elif region == 'ease-nh-wide':
        rp['lllat'] = 40.
        rp['lllon'] = -50.
        rp['urlat'] = 50.
        rp['urlon'] = 145.
        rp['long_name'] = 'Northern Hemisphere EASE2'
        rp['scale'] = 5
        rp['skip']  = 3
        rp['labelpos'] = 0.75
        rp['logo_left'] = 0.27
        rp['logo_bottom'] = 0.65
    elif region == 'ease-nh-very-wide':
        rp['lllat'] = 10.
        rp['lllon'] = -45.
        rp['urlat'] = 10.
        rp['urlon'] = 135.
        rp['long_name'] = 'Northern Hemisphere EASE2'
        rp['scale'] = 5
        rp['skip']  = 3
    elif region == 'ease-eur':
        rp['lllat'] = 65.
        rp['lllon'] = -30
        rp['urlat'] = 70.
        rp['urlon'] = 95.
        rp['long_name'] = 'European Arctic'
        rp['scale'] = 5
        rp['skip']  = 3
        rp['labelpos'] = 0.78
        rp['logo_left'] = 0.15
        rp['logo_bottom'] = 0.5
    elif region == 'ease-eur-trim':
        rp['lllat'] = 67.
        rp['lllon'] = -25.
        rp['urlat'] = 72.
        rp['urlon'] = 95.
        rp['long_name'] = 'European Arctic'
        rp['scale'] = 5
        rp['skip']  = 3
        rp['labelpos'] = 0.82
        rp['logo_left'] = 0.18
        rp['logo_bottom'] = 0.5
    elif region == 'ease-arc':
        rp['lllat'] = 58.
        rp['lllon'] = -50.
        rp['urlat'] = 58.
        rp['urlon'] = 140.
        rp['long_name'] = 'Arctic Ocean'
        rp['scale'] = 8
        rp['logo_left'] = 0.64
        rp['logo_bottom'] = 0.75
        rp['labelpos'] = 0.8
    elif region == 'ease-can':
        rp['lllat'] = 45.
        rp['lllon'] = -55.
        rp['urlat'] = 78.
        rp['urlon'] = -110
        rp['long_name'] = 'West Greenland and Canada'
        rp['scale'] = 8
        rp['logo_left'] = 0.27
        rp['logo_bottom'] = 0.44
        rp['labelpos'] = 0.76
    elif region == 'ease-sh':
        rp['lllat'] = -45.
        rp['lllon'] = -130.
        rp['urlat'] = -40.
        rp['urlon'] = 50.
        rp['long_name'] = 'Southern Hemisphere EASE2'
        rp['scale'] = 5
        rp['skip']  = 3
    elif region == 'ease-sh-wide':
        rp['lllat'] = -40.
        rp['lllon'] = -130.
        rp['urlat'] = -35.
        rp['urlon'] = 40.
        rp['long_name'] = 'Southern Hemisphere EASE2'
        rp['scale'] = 5
        rp['skip']  = 3
        rp['labelpos'] = 0.75
        rp['logo_left'] = 0.62
        rp['logo_bottom'] = 0.75
    elif region == 'ease-sh-very-wide':
        rp['lllat'] = -10.
        rp['lllon'] = -135.
        rp['urlat'] = -10.
        rp['urlon'] = 45.
        rp['long_name'] = 'Southern Hemisphere EASE2'
        rp['scale'] = 5
        rp['skip']  = 3
    elif region == 'ease-wed':
        rp['lllat'] = -65.
        rp['lllon'] = -82.
        rp['urlat'] = -48.
        rp['urlon'] = 25.
        rp['long_name'] = 'Weddell and Lazarev Sea'
        rp['scale'] = 5
        rp['skip']  = 3
        rp['logo_left'] = 0.52
        rp['logo_bottom'] = 0.23
        rp['labelpos'] = 0.83
    elif region == 'ease-ross':
        rp['lllat'] = -55.0
        rp['lllon'] =-150.0
        rp['urlat'] = -80.0
        rp['urlon'] = 150.0
        rp['long_name'] = 'Ross Sea'
        rp['scale'] = 6
        rp['logo_left'] = 0.63
        rp['logo_bottom'] = 0.75
    elif region == 'ease-bell':
        rp['lllat'] = -62.0
        rp['lllon'] = -102.0
        rp['urlat'] = -68.0
        rp['urlon'] = -45.0
        rp['long_name'] = 'Bellinghausen Sea'
        rp['scale'] = 6
        rp['labelpos'] = 0.65
        rp['logo_left'] = 0.36
        rp['logo_bottom'] = 0.08
    elif region == 'ease-wil':
        rp['lllat'] = -55.0
        rp['lllon'] = 180.0
        rp['urlat'] = -45.0
        rp['urlon'] = 70.0
        rp['long_name'] = 'Wilkes Land Coast'
        rp['scale'] = 6
        rp['labelpos'] = 0.75
        rp['logo_left'] = 0.29
        rp['logo_bottom'] = 0.75
    else:
        print("Region {} is not valid.".format(region))
        return None

    return rp
