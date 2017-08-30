# January 3, 2017

# Author: Beichman, Annabel
# Contact: annabel.beichman@gmail.com

# Purpose: Use dadi (Gutenkunst et al. 2010) functions to calculate
# the expected proportional site frequency spectrum
# under a stepwise demographic model inferred by Terhorst, Kamm & Song (2017) in SMC++
# for CEU, CHB, YRI populations

from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot
import pandas as pd
ns = np.array([20])
pts = [40,50,60]
Length = 4040000
Integration.set_timescale_factor = 0.0001


############  SMC++ RESULTS ##################################################
## Got data from Yun Song, results of SMC++ on CEU/CHB/YRI many genomes ##########
# gen time: 29 yr/gen
song_gen = 29
song_data = pd.read_csv("/Users/annabelbeichman/Documents/UCLA/PSMC_Project/SMCPLUS_DADI_20170103/data.csv")
song_data.head(10)
song_data.shape


####### CEU Song data ############
song_ceu_yrs = song_data.loc[song_data['label']=="ceu",'x']
song_ceu_yrs.head()
song_ceu_yrs.tail()
song_ceu_yrs.shape

# need to scale by 2 * Na for dadi (times as well as sizes!)
song_ceu_dips = tuple(song_data.loc[song_data['label']=="ceu",'y'])
song_ceu_dips_2Na = tuple(x/song_ceu_dips[-1] for x in song_ceu_dips)

# need to go from years --> 2Na generations
song_ceu_gen = song_ceu_yrs / (song_gen * 2 * song_ceu_dips[-1])
song_ceu_gen.head

# need intervals: (this is intervals from past to present -- goes to 1000ya, not to present day)
song_ceu_gen_int = tuple(np.diff(song_ceu_gen))

len(song_ceu_dips_2Na)
len(song_ceu_gen_int)

### dadi params (Nus & Thetas) (cut off first Nu -- the most recent pop size)
# you cut off the first Nu because you don't have an interval for that one; it's the end point
song_ceu_dadi = song_ceu_dips_2Na[1:]+song_ceu_gen_int
len(song_ceu_dadi) # should be 199+199= 398


# note: oldest time and size should be at end of each tuple (T0 and nuA)
# and most recent times should be at the beginning (nu198 T198)

####### CHB Song data ############
song_chb_yrs = song_data.loc[song_data['label']=="chb",'x']
song_chb_yrs.head()
song_chb_yrs.tail()
song_chb_yrs.shape

# dips: effective size:
# need to scale by 2 * Na for dadi (times as well as sizes!)
song_chb_dips = tuple(song_data.loc[song_data['label']=="chb",'y'])
song_chb_dips_2Na = tuple(x/song_chb_dips[-1] for x in song_chb_dips)

# need to go from years --> 2Na generations
song_chb_gen = song_chb_yrs / (song_gen * 2 * song_chb_dips[-1])
song_chb_gen.head
# need intervals: (this is intervals from past to present -- goes to 1000ya, not to present day)
song_chb_gen_int = tuple(np.diff(song_chb_gen))

len(song_chb_dips_2Na)
len(song_chb_gen_int)

### dadi params (Nus & Thetas) (cut off first Nu -- the most recent pop size)
# you cut off the first Nu because you don't have an interval for that one; it's the end point
song_chb_dadi = song_chb_dips_2Na[1:]+song_chb_gen_int
len(song_chb_dadi) # should be 199+199= 398


# note: oldest time and size should be at end of each tuple (T0 and nuA)
# and most recent times should be at the beginning (nu198 T198)


####### YRI Song data ############
song_yri_yrs = song_data.loc[song_data['label']=="yri",'x']
song_yri_yrs.head()
song_yri_yrs.tail()
song_yri_yrs.shape

# dips: effective size:
# need to scale by 2 * Na for dadi (times as well as sizes!)
song_yri_dips = tuple(song_data.loc[song_data['label']=="yri",'y'])
song_yri_dips_2Na = tuple(x/song_yri_dips[-1] for x in song_yri_dips)

# need to go from years --> 2Na generations
song_yri_gen = song_yri_yrs / (song_gen * 2 * song_yri_dips[-1])
song_yri_gen.head
# need intervals: (this is intervals from past to present -- goes to 1000ya, not to present day)
song_yri_gen_int = tuple(np.diff(song_yri_gen))

len(song_yri_dips_2Na)
len(song_yri_gen_int)

### dadi params (Nus & Thetas) (cut off first Nu -- the most recent pop size)
# you cut off the first Nu because you don't have an interval for that one; it's the end point
song_yri_dadi = song_yri_dips_2Na[1:]+song_yri_gen_int
len(song_yri_dadi) # should be 199+199= 398


# note: oldest time and size should be at end of each tuple (T0 and nuA)
# and most recent times should be at the beginning (nu198 T198)
######### dadi model for SMC++ (199 steps) ###################################
def SMCPLUS_model_198steps((nu198, nu197, nu196, nu195, nu194, nu193, nu192, nu191, nu190, nu189, nu188, nu187, nu186, nu185, nu184, nu183, nu182, nu181, nu180, nu179, nu178, nu177, nu176, nu175, nu174, nu173, nu172, nu171, nu170, nu169, nu168, nu167, nu166, nu165, nu164, nu163, nu162, nu161, nu160, nu159, nu158, nu157, nu156, nu155, nu154, nu153, nu152, nu151, nu150, nu149, nu148, nu147, nu146, nu145, nu144, nu143, nu142, nu141, nu140, nu139, nu138, nu137, nu136, nu135, nu134, nu133, nu132, nu131, nu130, nu129, nu128, nu127, nu126, nu125, nu124, nu123, nu122, nu121, nu120, nu119, nu118, nu117, nu116, nu115, nu114, nu113, nu112, nu111, nu110, nu109, nu108, nu107, nu106, nu105, nu104, nu103, nu102, nu101, nu100, nu99, nu98, nu97, nu96, nu95, nu94, nu93, nu92, nu91, nu90, nu89, nu88, nu87, nu86, nu85, nu84, nu83, nu82, nu81, nu80, nu79, nu78, nu77, nu76, nu75, nu74, nu73, nu72, nu71, nu70, nu69, nu68, nu67, nu66, nu65, nu64, nu63, nu62, nu61, nu60, nu59, nu58, nu57, nu56, nu55, nu54, nu53, nu52, nu51, nu50, nu49, nu48, nu47, nu46, nu45, nu44, nu43, nu42, nu41, nu40, nu39, nu38, nu37, nu36, nu35, nu34, nu33, nu32, nu31, nu30, nu29, nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nuA, T198, T197, T196, T195, T194, T193, T192, T191, T190, T189, T188, T187, T186, T185, T184, T183, T182, T181, T180, T179, T178, T177, T176, T175, T174, T173, T172, T171, T170, T169, T168, T167, T166, T165, T164, T163, T162, T161, T160, T159, T158, T157, T156, T155, T154, T153, T152, T151, T150, T149, T148, T147, T146, T145, T144, T143, T142, T141, T140, T139, T138, T137, T136, T135, T134, T133, T132, T131, T130, T129, T128, T127, T126, T125, T124, T123, T122, T121, T120, T119, T118, T117, T116, T115, T114, T113, T112, T111, T110, T109, T108, T107, T106, T105, T104, T103, T102, T101, T100, T99, T98, T97, T96, T95, T94, T93, T92, T91, T90, T89, T88, T87, T86, T85, T84, T83, T82, T81, T80, T79, T78, T77, T76, T75, T74, T73, T72, T71, T70, T69, T68, T67, T66, T65, T64, T63, T62, T61, T60, T59, T58, T57, T56, T55, T54, T53, T52, T51, T50, T49, T48, T47, T46, T45, T44, T43, T42, T41, T40, T39, T38, T37, T36, T35, T34, T33, T32, T31, T30, T29, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
    xx = Numerics.default_grid(pts)
    # intial phi with ancestral pop: (nuA = 1)
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
    phi = Integration.one_pop(phi, xx, T39, nu39)
    phi = Integration.one_pop(phi, xx, T40, nu40)
    phi = Integration.one_pop(phi, xx, T41, nu41)
    phi = Integration.one_pop(phi, xx, T42, nu42)
    phi = Integration.one_pop(phi, xx, T43, nu43)
    phi = Integration.one_pop(phi, xx, T44, nu44)
    phi = Integration.one_pop(phi, xx, T45, nu45)
    phi = Integration.one_pop(phi, xx, T46, nu46)
    phi = Integration.one_pop(phi, xx, T47, nu47)
    phi = Integration.one_pop(phi, xx, T48, nu48)
    phi = Integration.one_pop(phi, xx, T49, nu49)
    phi = Integration.one_pop(phi, xx, T50, nu50)
    phi = Integration.one_pop(phi, xx, T51, nu51)
    phi = Integration.one_pop(phi, xx, T52, nu52)
    phi = Integration.one_pop(phi, xx, T53, nu53)
    phi = Integration.one_pop(phi, xx, T54, nu54)
    phi = Integration.one_pop(phi, xx, T55, nu55)
    phi = Integration.one_pop(phi, xx, T56, nu56)
    phi = Integration.one_pop(phi, xx, T57, nu57)
    phi = Integration.one_pop(phi, xx, T58, nu58)
    phi = Integration.one_pop(phi, xx, T59, nu59)
    phi = Integration.one_pop(phi, xx, T60, nu60)
    phi = Integration.one_pop(phi, xx, T61, nu61)
    phi = Integration.one_pop(phi, xx, T62, nu62)
    phi = Integration.one_pop(phi, xx, T63, nu63)
    phi = Integration.one_pop(phi, xx, T64, nu64)
    phi = Integration.one_pop(phi, xx, T65, nu65)
    phi = Integration.one_pop(phi, xx, T66, nu66)
    phi = Integration.one_pop(phi, xx, T67, nu67)
    phi = Integration.one_pop(phi, xx, T68, nu68)
    phi = Integration.one_pop(phi, xx, T69, nu69)
    phi = Integration.one_pop(phi, xx, T70, nu70)
    phi = Integration.one_pop(phi, xx, T71, nu71)
    phi = Integration.one_pop(phi, xx, T72, nu72)
    phi = Integration.one_pop(phi, xx, T73, nu73)
    phi = Integration.one_pop(phi, xx, T74, nu74)
    phi = Integration.one_pop(phi, xx, T75, nu75)
    phi = Integration.one_pop(phi, xx, T76, nu76)
    phi = Integration.one_pop(phi, xx, T77, nu77)
    phi = Integration.one_pop(phi, xx, T78, nu78)
    phi = Integration.one_pop(phi, xx, T79, nu79)
    phi = Integration.one_pop(phi, xx, T80, nu80)
    phi = Integration.one_pop(phi, xx, T81, nu81)
    phi = Integration.one_pop(phi, xx, T82, nu82)
    phi = Integration.one_pop(phi, xx, T83, nu83)
    phi = Integration.one_pop(phi, xx, T84, nu84)
    phi = Integration.one_pop(phi, xx, T85, nu85)
    phi = Integration.one_pop(phi, xx, T86, nu86)
    phi = Integration.one_pop(phi, xx, T87, nu87)
    phi = Integration.one_pop(phi, xx, T88, nu88)
    phi = Integration.one_pop(phi, xx, T89, nu89)
    phi = Integration.one_pop(phi, xx, T90, nu90)
    phi = Integration.one_pop(phi, xx, T91, nu91)
    phi = Integration.one_pop(phi, xx, T92, nu92)
    phi = Integration.one_pop(phi, xx, T93, nu93)
    phi = Integration.one_pop(phi, xx, T94, nu94)
    phi = Integration.one_pop(phi, xx, T95, nu95)
    phi = Integration.one_pop(phi, xx, T96, nu96)
    phi = Integration.one_pop(phi, xx, T97, nu97)
    phi = Integration.one_pop(phi, xx, T98, nu98)
    phi = Integration.one_pop(phi, xx, T99, nu99)
    phi = Integration.one_pop(phi, xx, T100, nu100)
    phi = Integration.one_pop(phi, xx, T101, nu101)
    phi = Integration.one_pop(phi, xx, T102, nu102)
    phi = Integration.one_pop(phi, xx, T103, nu103)
    phi = Integration.one_pop(phi, xx, T104, nu104)
    phi = Integration.one_pop(phi, xx, T105, nu105)
    phi = Integration.one_pop(phi, xx, T106, nu106)
    phi = Integration.one_pop(phi, xx, T107, nu107)
    phi = Integration.one_pop(phi, xx, T108, nu108)
    phi = Integration.one_pop(phi, xx, T109, nu109)
    phi = Integration.one_pop(phi, xx, T110, nu110)
    phi = Integration.one_pop(phi, xx, T111, nu111)
    phi = Integration.one_pop(phi, xx, T112, nu112)
    phi = Integration.one_pop(phi, xx, T113, nu113)
    phi = Integration.one_pop(phi, xx, T114, nu114)
    phi = Integration.one_pop(phi, xx, T115, nu115)
    phi = Integration.one_pop(phi, xx, T116, nu116)
    phi = Integration.one_pop(phi, xx, T117, nu117)
    phi = Integration.one_pop(phi, xx, T118, nu118)
    phi = Integration.one_pop(phi, xx, T119, nu119)
    phi = Integration.one_pop(phi, xx, T120, nu120)
    phi = Integration.one_pop(phi, xx, T121, nu121)
    phi = Integration.one_pop(phi, xx, T122, nu122)
    phi = Integration.one_pop(phi, xx, T123, nu123)
    phi = Integration.one_pop(phi, xx, T124, nu124)
    phi = Integration.one_pop(phi, xx, T125, nu125)
    phi = Integration.one_pop(phi, xx, T126, nu126)
    phi = Integration.one_pop(phi, xx, T127, nu127)
    phi = Integration.one_pop(phi, xx, T128, nu128)
    phi = Integration.one_pop(phi, xx, T129, nu129)
    phi = Integration.one_pop(phi, xx, T130, nu130)
    phi = Integration.one_pop(phi, xx, T131, nu131)
    phi = Integration.one_pop(phi, xx, T132, nu132)
    phi = Integration.one_pop(phi, xx, T133, nu133)
    phi = Integration.one_pop(phi, xx, T134, nu134)
    phi = Integration.one_pop(phi, xx, T135, nu135)
    phi = Integration.one_pop(phi, xx, T136, nu136)
    phi = Integration.one_pop(phi, xx, T137, nu137)
    phi = Integration.one_pop(phi, xx, T138, nu138)
    phi = Integration.one_pop(phi, xx, T139, nu139)
    phi = Integration.one_pop(phi, xx, T140, nu140)
    phi = Integration.one_pop(phi, xx, T141, nu141)
    phi = Integration.one_pop(phi, xx, T142, nu142)
    phi = Integration.one_pop(phi, xx, T143, nu143)
    phi = Integration.one_pop(phi, xx, T144, nu144)
    phi = Integration.one_pop(phi, xx, T145, nu145)
    phi = Integration.one_pop(phi, xx, T146, nu146)
    phi = Integration.one_pop(phi, xx, T147, nu147)
    phi = Integration.one_pop(phi, xx, T148, nu148)
    phi = Integration.one_pop(phi, xx, T149, nu149)
    phi = Integration.one_pop(phi, xx, T150, nu150)
    phi = Integration.one_pop(phi, xx, T151, nu151)
    phi = Integration.one_pop(phi, xx, T152, nu152)
    phi = Integration.one_pop(phi, xx, T153, nu153)
    phi = Integration.one_pop(phi, xx, T154, nu154)
    phi = Integration.one_pop(phi, xx, T155, nu155)
    phi = Integration.one_pop(phi, xx, T156, nu156)
    phi = Integration.one_pop(phi, xx, T157, nu157)
    phi = Integration.one_pop(phi, xx, T158, nu158)
    phi = Integration.one_pop(phi, xx, T159, nu159)
    phi = Integration.one_pop(phi, xx, T160, nu160)
    phi = Integration.one_pop(phi, xx, T161, nu161)
    phi = Integration.one_pop(phi, xx, T162, nu162)
    phi = Integration.one_pop(phi, xx, T163, nu163)
    phi = Integration.one_pop(phi, xx, T164, nu164)
    phi = Integration.one_pop(phi, xx, T165, nu165)
    phi = Integration.one_pop(phi, xx, T166, nu166)
    phi = Integration.one_pop(phi, xx, T167, nu167)
    phi = Integration.one_pop(phi, xx, T168, nu168)
    phi = Integration.one_pop(phi, xx, T169, nu169)
    phi = Integration.one_pop(phi, xx, T170, nu170)
    phi = Integration.one_pop(phi, xx, T171, nu171)
    phi = Integration.one_pop(phi, xx, T172, nu172)
    phi = Integration.one_pop(phi, xx, T173, nu173)
    phi = Integration.one_pop(phi, xx, T174, nu174)
    phi = Integration.one_pop(phi, xx, T175, nu175)
    phi = Integration.one_pop(phi, xx, T176, nu176)
    phi = Integration.one_pop(phi, xx, T177, nu177)
    phi = Integration.one_pop(phi, xx, T178, nu178)
    phi = Integration.one_pop(phi, xx, T179, nu179)
    phi = Integration.one_pop(phi, xx, T180, nu180)
    phi = Integration.one_pop(phi, xx, T181, nu181)
    phi = Integration.one_pop(phi, xx, T182, nu182)
    phi = Integration.one_pop(phi, xx, T183, nu183)
    phi = Integration.one_pop(phi, xx, T184, nu184)
    phi = Integration.one_pop(phi, xx, T185, nu185)
    phi = Integration.one_pop(phi, xx, T186, nu186)
    phi = Integration.one_pop(phi, xx, T187, nu187)
    phi = Integration.one_pop(phi, xx, T188, nu188)
    phi = Integration.one_pop(phi, xx, T189, nu189)
    phi = Integration.one_pop(phi, xx, T190, nu190)
    phi = Integration.one_pop(phi, xx, T191, nu191)
    phi = Integration.one_pop(phi, xx, T192, nu192)
    phi = Integration.one_pop(phi, xx, T193, nu193)
    phi = Integration.one_pop(phi, xx, T194, nu194)
    phi = Integration.one_pop(phi, xx, T195, nu195)
    phi = Integration.one_pop(phi, xx, T196, nu196)
    phi = Integration.one_pop(phi, xx, T197, nu197)
    phi = Integration.one_pop(phi, xx, T198, nu198)
    # simulate sfs:
    fs = Spectrum.from_phi(phi,ns,(xx,))
    return fs

SMCPLUS_extrap_function_198 = Numerics.make_extrap_log_func(SMCPLUS_model_198steps)


#################################### CEU ##################################
# get SFS relative to theta = 1
smc_song_fs_ceu_relT = SMCPLUS_extrap_function_198(song_ceu_dadi,ns,pts)
# get the proportional sfs:
smc_song_fs_ceu_freq = smc_song_fs_ceu_relT / np.sum(smc_song_fs_ceu_relT)
# Plotting.plot_1d_comp_Poisson(obs_fs_CEU/np.sum(obs_fs_CEU),smc_song_fs_ceu_freq)
# get num snp, need a theta:
# Length is from Gutenkunst so that SFSes are comparable
Length = 4040000
# wsong mu: 1.25 × 10−8
song_mu= 1.25e-08
# ancestral size (oldest nu)
song_ceu_nuA = song_ceu_dips[-1]
# can also try it with the gutenkunst mu
song_theta_ceu_125 = 4 * song_ceu_nuA * song_mu

smc_song_fs_ceu_NUMSNP_125 = smc_song_fs_ceu_relT*song_theta_ceu_125*Length


# write out: Jan 3 2017

#################################### CHB ##################################
# get SFS relative to theta = 1
smc_song_fs_chb_relT = SMCPLUS_extrap_function_198(song_chb_dadi,ns,pts)
# get the proportional sfs:
smc_song_fs_chb_freq = smc_song_fs_chb_relT / np.sum(smc_song_fs_chb_relT)
# Plotting.plot_1d_comp_Poisson(obs_fs_CHB/np.sum(obs_fs_CHB),smc_song_fs_chb_freq)
# get num snp, need a theta:
# Length is from Gutenkunst so that SFSes are comparable
Length = 4040000
# what mu to use? song mu: 1.25 × 10−8 # cool this is the same as MSMC ; also use gutenkunst 2.35 mu
song_mu= 1.25e-08
# ancestral size (oldest nu)
song_chb_nuA = song_chb_dips[-1]
# can also try it with the gutenkunst mu
song_theta_chb_125 = 4 * song_chb_nuA * song_mu

smc_song_fs_chb_NUMSNP_125 = smc_song_fs_chb_relT*song_theta_chb_125*Length


#################################### YRI ##################################
# get SFS relative to theta = 1
smc_song_fs_yri_relT = SMCPLUS_extrap_function_198(song_yri_dadi,ns,pts)
# get the proportional sfs:
smc_song_fs_yri_freq = smc_song_fs_yri_relT / np.sum(smc_song_fs_yri_relT)
# Plotting.plot_1d_comp_Poisson(obs_fs_YRI/np.sum(obs_fs_YRI),smc_song_fs_yri_freq)
# get num snp, need a theta:
# Length is from Gutenkunst so that SFSes are comparable
Length = 4040000
# what mu to use? song mu: 1.25 × 10−8
song_mu= 1.25e-08
# ancestral size (oldest nu)
song_yri_nuA = song_yri_dips[-1]
# can also try it with the gutenkunst mu
song_theta_yri_125 = 4 * song_yri_nuA * song_mu

smc_song_fs_yri_NUMSNP_125 = smc_song_fs_yri_relT*song_theta_yri_125*Length

############ WRITE OUT SFSes #######################################
## CEU:
smc_song_fs_ceu_freq.to_file("smcPlus_dadi_CEU_SFS_freq.txt")
smc_song_fs_ceu_NUMSNP_125.to_file("smcPlus_dadi_CEU_SFS_NUMSNP_125mu.txt")
#
### CHB:
smc_song_fs_chb_freq.to_file("smcPlus_dadi_CHB_SFS_freq.txt")
smc_song_fs_chb_NUMSNP_125.to_file("smcPlus_dadi_CHB_SFS_NUMSNP_125mu.txt")
#
### YRI:
smc_song_fs_yri_freq.to_file("smcPlus_dadi_YRI_SFS_freq.txt")
smc_song_fs_yri_NUMSNP_125.to_file("smcPlus_dadi_YRI_SFS_NUMSNP_125mu.txt")
