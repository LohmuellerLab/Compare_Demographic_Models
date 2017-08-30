# November 28, 2016

# Author: Beichman, Annabel
# Contact: annabel.beichman@gmail.com

# Purpose: Use dadi (Gutenkunst et al. 2009) functions to calculate
# the expected proportional site frequency spectrum
# under a stepwise demographic model inferred by Schiffels & Durbin (2014)
# using MSMC on two, four or eight human genome haplotypes.


from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot
# set number of haploid individuals
ns = np.array([20])
# set extrapolation points
pts = [40,50,60]
# set sequence length
Length = 4040000
# Set integration timescale factor
Integration.set_timescale_factor = 0.0001

#### dadi model for 2-Haplotype and 4-Haplotype MSMC models: (38 steps)
# nuA is ancestral size
# T0 is oldest time interval
# Time intervals and Population sizes are in terms of 2*nuA generations

def msmc_model_38steps((nu38, nu37, nu36, nu35, nu34, nu33, nu32, nu31, nu30, nu29,nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19,nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9,nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nuA, T38, T37, T36, T35, T34, T33, T32, T31, T30, T29, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12,T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
    xx = Numerics.default_grid(pts)
    # initial phi with ancestral pop: (nuA = 1)
    phi = PhiManip.phi_1D(xx,nu=nuA)
    # stays at nuA for T0 duration of time:
    phi = Integration.one_pop(phi,xx,T0,nuA)
    # followed by a number of time steps, with associated pop changes:
    phi = Integration.one_pop(phi, xx, T1, nu1)
    phi = Integration.one_pop(phi, xx, T2, nu2)
    phi = Integration.one_pop(phi, xx, T3, nu3)
    phi = Integration.one_pop(phi, xx, T4, nu4)
    phi = Integration.one_pop(phi, xx, T5, nu5)
    phi = Integration.one_pop(phi, xx, T6, nu6)
    phi = Integration.one_pop(phi, xx, T7, nu7)
    phi = Integration.one_pop(phi, xx, T8, nu8)
    phi = Integration.one_pop(phi, xx, T9, nu9)
    phi = Integration.one_pop(phi, xx, T10, nu10)
    phi = Integration.one_pop(phi, xx, T11, nu11)
    phi = Integration.one_pop(phi, xx, T12, nu12)
    phi = Integration.one_pop(phi, xx, T13, nu13)
    phi = Integration.one_pop(phi, xx, T14, nu14)
    phi = Integration.one_pop(phi, xx, T15, nu15)
    phi = Integration.one_pop(phi, xx, T16, nu16)
    phi = Integration.one_pop(phi, xx, T17, nu17)
    phi = Integration.one_pop(phi, xx, T18, nu18)
    phi = Integration.one_pop(phi, xx, T19, nu19)
    phi = Integration.one_pop(phi, xx, T20, nu20)
    phi = Integration.one_pop(phi, xx, T21, nu21)
    phi = Integration.one_pop(phi, xx, T22, nu22)
    phi = Integration.one_pop(phi, xx, T23, nu23)
    phi = Integration.one_pop(phi, xx, T24, nu24)
    phi = Integration.one_pop(phi, xx, T25, nu25)
    phi = Integration.one_pop(phi, xx, T26, nu26)
    phi = Integration.one_pop(phi, xx, T27, nu27)
    phi = Integration.one_pop(phi, xx, T28, nu28)
    phi = Integration.one_pop(phi, xx, T29, nu29)
    phi = Integration.one_pop(phi, xx, T30, nu30)
    phi = Integration.one_pop(phi, xx, T31, nu31)
    phi = Integration.one_pop(phi, xx, T32, nu32)
    phi = Integration.one_pop(phi, xx, T33, nu33)
    phi = Integration.one_pop(phi, xx, T34, nu34)
    phi = Integration.one_pop(phi, xx, T35, nu35)
    phi = Integration.one_pop(phi, xx, T36, nu36)
    phi = Integration.one_pop(phi, xx, T37, nu37)
    phi = Integration.one_pop(phi, xx, T38, nu38)
    # get sfs:
    fs = Spectrum.from_phi(phi,ns,(xx,))
    return fs

# wrap your function in Numerics.make_extrap_log_func to extrapolate over pts
msmc_extrap_function_38 = Numerics.make_extrap_log_func(msmc_model_38steps)

# usage: msmc_extrap_function_38(Parameters,ns,pts)
# where Parameters is a vector of MSMC output parameters
# converted to population sizes across time intervals, in units of 2*nuA (see below)

###### dadi model for 8 Haplotype MSMC model (28 steps) #########################
# need 28 steps for 8 haplotype model
def msmc_model_28steps((nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19,nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9,nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nuA, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12,T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
    xx = Numerics.default_grid(pts)
    # initial phi with ancestral pop: (nuA = 1)
    phi = PhiManip.phi_1D(xx,nu=nuA)
    # stays at nuA for T0 duration of time:
    phi = Integration.one_pop(phi,xx,T0,nuA)
    # followed by a number of time steps, with associated pop changes:
    phi = Integration.one_pop(phi, xx, T1, nu1)
    phi = Integration.one_pop(phi, xx, T2, nu2)
    phi = Integration.one_pop(phi, xx, T3, nu3)
    phi = Integration.one_pop(phi, xx, T4, nu4)
    phi = Integration.one_pop(phi, xx, T5, nu5)
    phi = Integration.one_pop(phi, xx, T6, nu6)
    phi = Integration.one_pop(phi, xx, T7, nu7)
    phi = Integration.one_pop(phi, xx, T8, nu8)
    phi = Integration.one_pop(phi, xx, T9, nu9)
    phi = Integration.one_pop(phi, xx, T10, nu10)
    phi = Integration.one_pop(phi, xx, T11, nu11)
    phi = Integration.one_pop(phi, xx, T12, nu12)
    phi = Integration.one_pop(phi, xx, T13, nu13)
    phi = Integration.one_pop(phi, xx, T14, nu14)
    phi = Integration.one_pop(phi, xx, T15, nu15)
    phi = Integration.one_pop(phi, xx, T16, nu16)
    phi = Integration.one_pop(phi, xx, T17, nu17)
    phi = Integration.one_pop(phi, xx, T18, nu18)
    phi = Integration.one_pop(phi, xx, T19, nu19)
    phi = Integration.one_pop(phi, xx, T20, nu20)
    phi = Integration.one_pop(phi, xx, T21, nu21)
    phi = Integration.one_pop(phi, xx, T22, nu22)
    phi = Integration.one_pop(phi, xx, T23, nu23)
    phi = Integration.one_pop(phi, xx, T24, nu24)
    phi = Integration.one_pop(phi, xx, T25, nu25)
    phi = Integration.one_pop(phi, xx, T26, nu26)
    phi = Integration.one_pop(phi, xx, T27, nu27)
    phi = Integration.one_pop(phi, xx, T28, nu28)
    # get sfs:
    fs = Spectrum.from_phi(phi,ns,(xx,))
    return fs

# wrap your function in Numerics.make_extrap_log_func to extrapolate over pts

msmc_extrap_function_28 = Numerics.make_extrap_log_func(msmc_model_28steps)

# usage: msmc_extrap_function_28(Parameters,ns,pts)
# where Parameters is a vector of MSMC output parameters
# converted to population sizes across time intervals, in units of 2*nuA

############################# transform MSMC parameters #############################
# "LY" is the "left years" column of MSMC output
# "Dips" is the corresponding population size in diploids
# YRI, CHB, CEU are the three populations, and 2, 4 and 8 refer to the
# number of haplotypes used in the inference
#YRI : 2
yri_2_LY = (0.0000,22173.1680,44922.2400,68278.3200,92274.0000,116946.0000,142332.9600,168478.0800,195427.6800,223232.8800,251949.6000,281640.0000,312374.4000,344224.8000,377277.6000,411626.4000,447376.8000,484651.2000,523581.6000,564324.0000,607053.6000,651976.8000,699328.8000,749385.6000,802480.8000,859003.2000,919428.0000,984331.2000,1054430.4000,1130635.2000,1214107.2000,1306382.4000,1409534.4000,1526481.6000,1661484.0000,1821160.8000,2016588.0000,2268540.0000,2623632.0000,3230688.0000)
yri_2_Dips = (17385.0306,16432.9085,14094.7804,11228.8582,13596.1931,18897.6973,24303.5514,28720.9828,31233.4073,32311.2217,32899.6068,33255.7366,33011.7439,33011.7439,31332.0017,31332.0017,28537.5911,28537.5911,25268.7968,25268.7968,22146.0644,22146.0644,19524.3860,19524.3860,17473.4295,17473.4295,15927.2443,15927.2443,14790.1645,14790.1645,13997.6134,13997.6134,13578.9770,13578.9770,13957.3185,13957.3185,18543.1569,18543.1569,41100.4227,41100.4227)
#CHB : 2
chb_2_LY = (0.0000,21185.1360,42920.6400,65235.8400,88162.3200,111734.8800,135990.7200,160970.6400,186719.5200,213285.6000,240722.4000,269090.4000,298454.4000,328886.4000,360465.6000,393283.2000,427442.4000,463056.0000,500251.2000,539176.8000,580003.2000,622924.8000,668164.8000,715994.4000,766723.2000,820725.6000,878457.6000,940468.8000,1007445.6000,1080254.4000,1160006.4000,1248168.0000,1346726.4000,1458460.8000,1587448.8000,1740009.6000,1926729.6000,2167452.0000,2506728.0000,3086736.0000)
chb_2_Dips = (61086.7019,5121.4815,2279.8908,3293.1568,4924.8953,7712.7762,11940.5480,17004.2723,21922.6132,25886.6166,28725.9331,30666.1453,32109.9444,32109.9444,31517.1572,31517.1572,28813.6692,28813.6692,25523.0633,25523.0633,22497.9471,22497.9471,20023.2269,20023.2269,18113.4810,18113.4810,16662.9175,16662.9175,15534.1012,15534.1012,14643.8613,14643.8613,14054.2706,14054.2706,14289.3887,14289.3887,18881.9970,18881.9970,41126.2844,41126.2844)
#CEU : 2
ceu_2_LY = (0.0000,22396.2240,45374.1600,68964.9600,93202.3200,118122.2400,143764.8000,170172.7200,197393.5200,225478.5600,254484.0000,284474.4000,315516.0000,347685.6000,381072.0000,415766.4000,451879.2000,489525.6000,528849.6000,570000.0000,613159.2000,658533.6000,706363.2000,756924.0000,810554.4000,867645.6000,928675.2000,994231.2000,1065038.4000,1142008.8000,1226320.8000,1319522.4000,1423713.6000,1541836.8000,1678197.6000,1839480.0000,2036875.2000,2291359.2000,2650032.0000,3263184.0000)
ceu_2_Dips = (20589.0528,5770.6984,2996.7485,3714.2963,5541.2330,9115.2925,14000.6510,18990.0160,23260.9530,26686.0585,29403.5490,31284.2171,32060.2733,32060.2733,30416.4005,30416.4005,27366.8943,27366.8943,24236.2550,24236.2550,21533.5088,21533.5088,19346.6632,19346.6632,17618.1961,17618.1961,16256.2637,16256.2637,15207.7377,15207.7377,14463.9306,14463.9306,14068.2592,14068.2592,14447.7875,14447.7875,19054.3334,19054.3334,41261.2746,41261.2746)
#YRI : 4
yri_4_LY = (0.0000,5016.6720,10163.6640,15447.9120,20876.9760,26459.0400,32202.7200,38118.0000,44215.4400,50506.3200,57003.6000,63721.2000,70674.4800,77880.4800,85358.6400,93130.3200,101219.0400,109652.1600,118460.1600,127678.0800,137345.7600,147509.2800,158222.6400,169548.4800,181561.2000,194349.3600,208020.0000,222704.4000,238564.8000,255806.4000,274692.0000,295567.2000,318907.2000,345364.8000,375909.6000,412036.8000,456252.0000,513256.8000,593599.2000,730944.0000)
yri_4_Dips = (36822.2406,11856.0322,10732.6078,13292.3044,14356.4712,14278.8808,13980.6857,13547.3361,13030.7591,12500.5469,12639.1491,12863.5149,13476.0446,13476.0446,15151.6299,15151.6299,17828.0132,17828.0132,20871.4890,20871.4890,24560.2186,24560.2186,28248.3881,28248.3881,31695.7211,31695.7211,33906.0633,33906.0633,34981.5472,34981.5472,35069.5693,35069.5693,33498.0320,33498.0320,28924.5142,28924.5142,22568.5237,22568.5237,205844.9679,205844.9679)
#CHB : 4
chb_4_LY = (0.0000,3503.0160,7097.0400,10786.9200,14577.8880,18475.6560,22486.4400,26616.9600,30874.5600,35267.2800,39804.2400,44495.0400,49350.2400,54382.0800,59604.0000,65030.6400,70678.8000,76567.4400,82717.9200,89154.4800,95905.2000,103002.2400,110483.0400,118391.7600,126779.7600,135709.4400,145255.4400,155509.2000,166584.0000,178623.1200,191810.4000,206388.2400,222684.9600,241161.6000,262490.4000,287716.8000,318590.4000,358394.4000,414494.4000,510400.8000)
chb_4_Dips = (1543352.7796,48196.0809,63813.5623,25405.6972,11547.8775,6332.6009,4383.6914,3523.0806,3076.3315,2755.8458,2706.1585,2745.3672,2918.0679,2918.0679,3331.9450,3331.9450,3922.1454,3922.1454,4620.4155,4620.4155,5870.0862,5870.0862,7892.1615,7892.1615,10973.8466,10973.8466,14786.1734,14786.1734,19331.2359,19331.2359,24615.5361,24615.5361,28892.5486,28892.5486,30091.4781,30091.4781,27128.7599,27128.7599,187513.4775,187513.4775)
#CEU : 4
ceu_4_LY = (0.0000,3730.3680,7557.6480,11487.0000,15524.0160,19674.7680,23945.8560,28344.4800,32878.3200,37556.4000,42387.6000,47382.7200,52553.2800,57911.7600,63472.3200,69251.2800,75265.9200,81536.8800,88086.4800,94940.8800,102129.6000,109687.2000,117653.7600,126075.3600,135007.9200,144517.2000,154682.8800,165601.9200,177395.7600,190216.0800,204259.2000,219783.3600,237137.7600,256812.0000,279525.6000,306388.8000,339266.4000,381655.2000,441396.0000,543525.6000)
ceu_4_Dips = (613285.2927,14521.5863,13292.3044,14309.5294,11374.9467,7539.5449,5660.5589,4713.3572,4095.1466,3570.0261,3359.3402,3286.5817,3371.6863,3371.6863,3840.9464,3840.9464,4670.0221,4670.0221,5791.0167,5791.0167,7476.6495,7476.6495,9826.6095,9826.6095,12957.6899,12957.6899,16490.2213,16490.2213,20389.1265,20389.1265,24601.7590,24601.7590,28383.8921,28383.8921,29324.0083,29324.0083,25618.1992,25618.1992,191238.4121,191238.4121)
#YRI : 8
yri_8_LY = (0.0000,1483.4112,3018.8880,4610.2080,6261.5760,7977.7440,9763.9680,11626.2240,13571.2800,15606.8160,17741.7120,19986.1200,22351.9200,24852.9600,27505.6800,30329.7600,33348.4800,36591.1200,40093.6800,43901.0400,48071.2800,52681.6800,57835.4400,63678.2400,70423.2000,78401.0400,88164.9600,100752.9600,118494.7200,148824.2400)
yri_8_Dips = (1062487.5490,53906.6122,20514.9246,13146.2188,10684.5312,9804.9309,9696.7111,10075.3637,10685.1020,11294.2027,12007.3966,12007.3966,12613.4024,12613.4024,12898.4409,12898.4409,12884.6885,12884.6885,12705.0274,12705.0274,12670.5375,12670.5375,12946.5341,12946.5341,13928.6921,13928.6921,16715.6994,16715.6994,132016.6869,132016.6869)
#CHB : 8
chb_8_LY = (0.0000,981.2184,1996.8696,3049.4640,4141.7760,5276.9520,6458.4720,7690.2720,8976.8640,10323.2880,11735.4240,13220.0160,14784.8880,16439.2320,18193.8960,20061.8640,22058.7120,24203.7600,26520.2400,29038.8000,31797.3600,34846.8000,38255.7600,42120.4800,46582.0800,51859.2000,58317.6000,66643.9200,78379.4400,98441.2800)
chb_8_Dips = (16559924.1555,55514151.2511,1234179.3633,288817.0055,154615.6641,91088.3694,54805.2564,33788.3498,21339.1376,14102.5325,8474.0018,8474.0018,5251.6526,5251.6526,3856.3509,3856.3509,3303.9284,3303.9284,2940.2684,2940.2684,2884.2133,2884.2133,2937.0732,2937.0732,3171.0295,3171.0295,3727.1016,3727.1016,5666.1718,5666.1718)
#CEU : 8
ceu_8_LY = (0.0000,1044.3936,2125.4400,3245.8080,4408.4640,5616.7200,6874.3200,8185.4400,9554.8320,10987.9680,12491.0400,14071.2000,15736.8240,17497.7040,19365.3360,21353.5440,23478.9840,25762.0800,28227.8400,30908.4000,33844.5600,37090.3200,40718.8800,44832.4800,49581.3600,55198.0800,62072.4000,70934.8800,83425.9200,104779.4400)
ceu_8_Dips = (4084053.9136,25633.3028,11027.4252,11224.0374,12481.5118,13371.5760,13630.2455,13042.1882,11702.0654,10052.3476,7749.5350,7749.5350,5639.0288,5639.0288,4372.8990,4372.8990,3759.9639,3759.9639,3371.2032,3371.2032,3297.5820,3297.5820,3369.6980,3369.6980,3616.1788,3616.1788,3929.2731,3929.2731,4036.5225,4036.5225)

########################### conversion to dadi parameters: ###########################
# need to divide by 30yrs/gen (value msmc used) and scale by 2*Na (ancestral size)
# then convert times into intervals

## YRI 2 ###
yri_2_LY_dadi = tuple(x/(30*2*yri_2_Dips[-1]) for x in yri_2_LY)
# then need to turn into intervals using numpy.diff()
yri_2_LY_dadi_intervals = tuple(np.diff(yri_2_LY_dadi))
# scale Dips relative to Na (oldest ancestral size)
yri_2_Dips_dadi = tuple(x/yri_2_Dips[-1] for x in yri_2_Dips)

### YRI 4 ###
yri_4_LY_dadi = tuple(x/(30*2*yri_4_Dips[-1]) for x in yri_4_LY)
# then need to turn into intervals using numpy.diff()
yri_4_LY_dadi_intervals = tuple(np.diff(yri_4_LY_dadi))
# scale Dips relative to Na (oldest ancestral size)
yri_4_Dips_dadi = tuple(x/yri_4_Dips[-1] for x in yri_4_Dips)

### YRI 8 ###
yri_8_LY_dadi = tuple(x/(30*2*yri_8_Dips[-1]) for x in yri_8_LY)
# then need to turn into intervals using numpy.diff()
yri_8_LY_dadi_intervals = tuple(np.diff(yri_8_LY_dadi))
# scale Dips relative to Na (oldest ancestral size)
yri_8_Dips_dadi = tuple(x/yri_8_Dips[-1] for x in yri_8_Dips)

### CHB dadi parameters ###
### CHB 2 ###
chb_2_LY_dadi = tuple(x/(30*2*chb_2_Dips[-1]) for x in chb_2_LY)
# then need to turn into intervals using numpy.diff()
chb_2_LY_dadi_intervals = tuple(np.diff(chb_2_LY_dadi))
chb_2_Dips_dadi = tuple(x/chb_2_Dips[-1] for x in chb_2_Dips)

### CHB 4 ###
chb_4_LY_dadi = tuple(x/(30*2*chb_4_Dips[-1]) for x in chb_4_LY)
# then need to turn into intervals using numpy.diff()
chb_4_LY_dadi_intervals = tuple(np.diff(chb_4_LY_dadi))
chb_4_Dips_dadi = tuple(x/chb_4_Dips[-1] for x in chb_4_Dips)

### CHB 8 ###
chb_8_LY_dadi = tuple(x/(30*2*chb_8_Dips[-1]) for x in chb_8_LY)
# then need to turn into intervals using numpy.diff()
chb_8_LY_dadi_intervals = tuple(np.diff(chb_8_LY_dadi))
chb_8_Dips_dadi = tuple(x/chb_8_Dips[-1] for x in chb_8_Dips)

### CEU dadi parameters ###
## CEU 2 ###
ceu_2_LY_dadi = tuple(x/(30*2*ceu_2_Dips[-1]) for x in ceu_2_LY)
# then need to turn into intervals using numpy.diff()
ceu_2_LY_dadi_intervals = tuple(np.diff(ceu_2_LY_dadi))
# scale Dips relative to Na (oldest ancestral size)
ceu_2_Dips_dadi = tuple(x/ceu_2_Dips[-1] for x in ceu_2_Dips)

### CEU 4 ###
ceu_4_LY_dadi = tuple(x/(30*2*ceu_4_Dips[-1]) for x in ceu_4_LY)
# then need to turn into intervals using numpy.diff()
ceu_4_LY_dadi_intervals = tuple(np.diff(ceu_4_LY_dadi))
ceu_4_Dips_dadi = tuple(x/ceu_4_Dips[-1] for x in ceu_4_Dips)

### CEU 8 ###
ceu_8_LY_dadi = tuple(x/(30*2*ceu_8_Dips[-1]) for x in ceu_8_LY)
# then need to turn into intervals using numpy.diff()
ceu_8_LY_dadi_intervals = tuple(np.diff(ceu_8_LY_dadi))
ceu_8_Dips_dadi = tuple(x/ceu_8_Dips[-1] for x in ceu_8_Dips)


### concatenate the sizes and times to make your dadi paramteers (Nus & Thetas)
# (cut off last Nu)
yri_2_dadi_params_30yr_125mu = yri_2_Dips_dadi[:-1]+yri_2_LY_dadi_intervals
yri_4_dadi_params_30yr_125mu = yri_4_Dips_dadi[:-1]+yri_4_LY_dadi_intervals
yri_8_dadi_params_30yr_125mu = yri_8_Dips_dadi[:-1]+yri_8_LY_dadi_intervals


ceu_2_dadi_params_30yr_125mu = ceu_2_Dips_dadi[:-1]+ceu_2_LY_dadi_intervals
ceu_4_dadi_params_30yr_125mu = ceu_4_Dips_dadi[:-1]+ceu_4_LY_dadi_intervals
ceu_8_dadi_params_30yr_125mu = ceu_8_Dips_dadi[:-1]+ceu_8_LY_dadi_intervals

chb_2_dadi_params_30yr_125mu = chb_2_Dips_dadi[:-1]+chb_2_LY_dadi_intervals
chb_4_dadi_params_30yr_125mu = chb_4_Dips_dadi[:-1]+chb_4_LY_dadi_intervals
chb_8_dadi_params_30yr_125mu = chb_8_Dips_dadi[:-1]+chb_8_LY_dadi_intervals

################# Get your expected SFSs #####################################

# these are relative to theta = 1
msmc_fs_yri_relT_2_new = msmc_extrap_function_38(yri_2_dadi_params_30yr_125mu,ns,pts)
msmc_fs_yri_relT_4_new = msmc_extrap_function_38(yri_4_dadi_params_30yr_125mu,ns,pts)
msmc_fs_yri_relT_8_new = msmc_extrap_function_28(yri_8_dadi_params_30yr_125mu,ns,pts)

# to get proportional SFS, divide by the sum of the SFS
msmc_fs_yri_freq_2_new = msmc_fs_yri_relT_2_new/np.sum(msmc_fs_yri_relT_2_new)
msmc_fs_yri_freq_4_new = msmc_fs_yri_relT_4_new/np.sum(msmc_fs_yri_relT_4_new)
msmc_fs_yri_freq_8_new = msmc_fs_yri_relT_8_new/np.sum(msmc_fs_yri_relT_8_new)

msmc_fs_ceu_relT_2_new = msmc_extrap_function_38(ceu_2_dadi_params_30yr_125mu,ns,pts)
msmc_fs_ceu_relT_4_new = msmc_extrap_function_38(ceu_4_dadi_params_30yr_125mu,ns,pts)
msmc_fs_ceu_relT_8_new = msmc_extrap_function_28(ceu_8_dadi_params_30yr_125mu,ns,pts)

msmc_fs_ceu_freq_2_new = msmc_fs_ceu_relT_2_new/np.sum(msmc_fs_ceu_relT_2_new)
msmc_fs_ceu_freq_4_new = msmc_fs_ceu_relT_4_new/np.sum(msmc_fs_ceu_relT_4_new)
msmc_fs_ceu_freq_8_new = msmc_fs_ceu_relT_8_new/np.sum(msmc_fs_ceu_relT_8_new)

msmc_fs_chb_relT_2_new = msmc_extrap_function_38(chb_2_dadi_params_30yr_125mu,ns,pts)
msmc_fs_chb_relT_4_new = msmc_extrap_function_38(chb_4_dadi_params_30yr_125mu,ns,pts)
msmc_fs_chb_relT_8_new = msmc_extrap_function_28(chb_8_dadi_params_30yr_125mu,ns,pts)

msmc_fs_chb_freq_2_new = msmc_fs_chb_relT_2_new/np.sum(msmc_fs_chb_relT_2_new)
msmc_fs_chb_freq_4_new = msmc_fs_chb_relT_4_new/np.sum(msmc_fs_chb_relT_4_new)
msmc_fs_chb_freq_8_new = msmc_fs_chb_relT_8_new/np.sum(msmc_fs_chb_relT_8_new)

# To get the absolute SFS based on SNP Counts
# you multiply the SFS relative to theta = 1 ("relT")
# by your desired theta (4*nuA*mutation rate*sequence length)

# you can calculate liklihoods using the two functions in:
# multinomial_poisson_likelihood_Functions.py
########################## write out proportional SFS: ###################################
msmc_fs_yri_freq_2_new.tofile("msmc_dadi_YRI_30yr_freq_2.txt")
msmc_fs_ceu_freq_2_new.tofile("msmc_dadi_CEU_30yr_freq_2.txt")
msmc_fs_chb_freq_2_new.tofile("msmc_dadi_CHB_30yr_freq_2.txt")

msmc_fs_yri_freq_4_new.tofile("msmc_dadi_YRI_30yr_freq_4.txt")
msmc_fs_ceu_freq_4_new.tofile("msmc_dadi_CEU_30yr_freq_4.txt")
msmc_fs_chb_freq_4_new.tofile("msmc_dadi_CHB_30yr_freq_4.txt")

msmc_fs_yri_freq_8_new.tofile("msmc_dadi_YRI_30yr_freq_8.txt")
msmc_fs_ceu_freq_8_new.tofile("msmc_dadi_CEU_30yr_freq_8.txt")
msmc_fs_chb_freq_8_new.tofile("msmc_dadi_CHB_30yr_freq_8.txt")
