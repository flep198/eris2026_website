input_idifits=["rsm07_1_1.IDI1","rsm07_1_1.IDI2"]
msfile="rsm07.ms" #name of the measurement set to be created
sources=["3C345","J1848+3219","3C395"] #sources to export

#load the data
importfitsidi(input_idifits,vis=msfile,constobsid=True,scanreindexgap_s=15)

#generate tsys calibration table
gencal(vis=msfile,caltable="tsys.cal",caltype="tsys",uniform=False)

#generate gain calibration table
gencal(vis=msfile,caltable="gain.cal",caltype="gc")

#apply flag data from stations
flagdata(msfile, mode="list", inpfile="rsm07.casa.flag")

#apply flags from PI letter information
flagdata(msfile, mode='manual', antenna='O8,TR,CM,DA,KN,PI')
flagdata(msfile, mode='manual', antenna='DE', timerange='13:45:00~15:03:00')

#flag edge channels
flagdata(msfile,mode="manual",spw="*:0~3;60~63")

#WB requires some additional flagging at the beginning of each scan ("quacking")
flagdata(msfile, mode='quack', antenna='WB', quackinterval=40, quackmode='beg')


#run single-band delay fringefit
fringefit(vis=msfile,
          caltable="sbd.cal",
          timerange="15:48:00~15:50:00",
          solint='inf',
          zerorates=True,
          refant='EF',
          minsnr=10,
          gaintable=['gain.cal','tsys.cal'],
          interp=['nearest','nearest,nearest'],
          parang=True)

#run multi-band delay fringefit
fringefit(vis=msfile,
          caltable="mbd.cal",
          field="J1848+3219,3C345",
          solint="inf",
          zerorates=False,
          refant="EF,MC",
          combine="spw",
          minsnr=7,
          gaintable=["gain.cal","tsys.cal","sbd.cal"],
          interp=["nearest","nearest,nearest","nearest"],
          parang=True)

#bandpass correction
bandpass(vis=msfile,
         caltable="bpass.cal",
         field="J1848+3219,3C345",
         gaintable=["gain.cal","tsys.cal","sbd.cal","mbd.cal"],
         interp=["nearest","nearest,nearest","nearest","linear"],
         solnorm=True,
         solint="inf",
         refant="EF",
         bandtype="B",
         spwmap=[[],[],[],[0,0,0,0]],
         parang=True)

#apply calibration
applycal(vis=msfile,
         field="",
         gaintable=["gain.cal","tsys.cal","sbd.cal","mbd.cal","bpass.cal"],
         interp=["nearest","nearest,nearest","nearest","linear","nearest,nearest"],
         spwmap=[[],[],[],[0,0,0,0],[]],
         parang=True)

#average and export calibrated data
for source in sources:
    try:
        mstransform(vis=msfile,
                field=source,
                outputvis=source+"_calibrated.ms",
                datacolumn="corrected",
                keepflags=True,
                chanaverage=True,
                chanbin=64)

        exportuvfits(vis=source+"_calibrated.ms",
                     field=source,
                     fitsfile=source+".uvf",
                     datacolumn="corrected",
                     combinespw=True,
                     multisource=False)
    except:
        print(f"Error with exporting source {source}")
