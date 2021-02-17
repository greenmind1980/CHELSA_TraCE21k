### =========================================================================
### CHELSA-TraCE21k v1.0. Downscaled transient temperature and precipitation 
### data since the last glacial maximum
### Dirk Nikolaus Karger1, Michael Nobis1, Signe Normand2, Catherine H. Graham1, Niklaus E. Zimmermann1
### 1Swiss Federal Research Institute WSL, Zürcherstrasse 111, 8903 Birmensdorf, Switzerland
### 2Aarhus University, Ny Munkegade 116, 8000 Aarhus, Denmark
### Correspondence to: Dirk Nikolaus Karger (dirk.karger@wsl.ch)
### THIS MODULE CALCULATES AND DOWNSCALES GLACIER HEIGHT, 
### AN OROGRAPHY WITH SEA LEVEL RISE, AND MEAN ANNUAL TEMPERATURE FROM TRACE AND ICE6G
### dependencies SAGA-GIS 2.7.1, R 3.4.4
### code version 1.0
### =========================================================================

### =========================================================================
### load required libraries
### =========================================================================

library(ncdf4)
library(raster)

### =========================================================================
### create vectors of focal timesteps (can be adapted to purpose)
### =========================================================================

yearvector<-seq(1,21000,100)            # vector of years 1,endyear starting from LGM, resolution (not adaptable)
tvect<-seq(0,length(yearvector)/10,0.5) # vector converting the timesteps of the ICE6G data
gvect<-rev(rep(0:42,by=0.5,each=5)/2)   # vector containing the timesteps of the ICE6G data
decades<-rep(1:2200,by=1)               # vector of decades
centuries<-rep(1:220,by=10,each=10)     # vector of centuries

### =========================================================================
### set the path to the input, temporary, and output directories
### =========================================================================

TEMP  <-"/mnt/md0/data/TESTICEMELT/TEMP/"   # path to temporary directory
INPUT <-"/mnt/md0/data/TESTICEMELT/INPUT/"  # path to input directory
OUTPUT<-"/mnt/md0/data/TESTICEMELT/OUTPUT/" # path to output directory

### =========================================================================
### loop over timesteps
### =========================================================================
for (t in 1:lenght(yearvector))    
{ 
  ### =========================================================================
  ### calculations for t
  ### =========================================================================
  
  # set the names of the input files
  tas_glacieroutline_t0_interpol <- paste(INPUT,"glacier_boarder_heights_t0_interpol_biascor.sgrd",sep="") # file containing the annual mean 2m air temperature at the glacier outline at t0 (LGM) at 0.5° resolution.
  nc1                            <- nc_open(paste0(INPUT,"b30.00_4kaDVTd.cam2.ncrcat.ann.nc")) # ncdf4 containing the TraCE21k data
  #tas_name_t1_highres_temp       <- paste0(TEMP,"tas_t1_highres_temp")
  glacier_lgm                    <- paste0(INPUT,"glacier_heights.sgrd")  # file containing the high resolution glacier at LGM
  biasgrid                       <- paste0(INPUT,"t0_t210_bias.sgrd") # file containing the bias between current and LGM glacier temperatures
  gmted2010                      <- paste0(INPUT,"GMTED_lgm_sea_zero.sgrd") # file containing digital elevation model
  altgridname                    <- paste0(INPUT,"GMTED_bathymetry.sgrd") # file containing the bathymetry
  lookup1                        <- paste0(INPUT,"lookup1.txt") # lookup table for classification
  randompoints                   <- paste0(INPUT,"randompoints.shp") # random points for sample of glacier heights
  ice6gname                      <- paste0(INPUT,"ICE6G/I6_C.VM5a_1deg.",format(round(gvect[t], 1), nsmall = 1),".nc")
  orogbias                       <- paste0(INPUT,"orog_bias.sdat",sep="")
  sealevel                       <- read.csv(paste(INPUT,"sealevelchange.csv"))
  lookup2                        <- paste0(INPUT,"lookup2.txt")
  traceinput                     <- paste0(INPUT,"trace_template.sdat")
  
  # set the names of the output files (names can be arbritray, but have to include 't')
  tas_name_t0                    <- paste0(OUTPUT,"CHELSA_TraCE_tas_",t,".sgrd") # ANNUAL TEMPERATURE AT T0 <t of timestep AT TRACE RESOLUTION
  tas_name_t1                    <- paste0(OUTPUT,"CHELSA_TraCE_tas_",t,".sgrd") # ANNUAL TEMPERATURE AT T1 <t of timestep AT TRACE RESOLUTION
  glacier_t                      <- paste0(OUTPUT,"glacierheigth_",t,".sgrd")
  gmted_tx<-paste(OUTPUT,"gmted2010_",t,".sgrd",sep="")
  glacierheights<-paste(OUTPUT,"glacierheigth_",t,".sgrd",sep="")
  glacier_dem<-paste(OUTPUT,"glacier_plus_dem_",t,".sgrd",sep="")
  
  # set the names of the temporary files (names are arbritrary)
  glacier_outline_line           <- paste0(TEMP,"glacier_outline_temp_line.shp")
  glacier_outline_grid           <- paste0(TEMP,"glacier_outline_temp.sgrd")
  glacier_outline_shape          <- paste0(TEMP,"glacier_outline_temp.shp")
  ice6grandom_points             <- paste0(TEMP,"ice6g_randompoints.shp")
  glacier_outline_points         <- paste0(TEMP,"glacier_outline_temp_points.shp")
  combpoints                     <- paste0(TEMP,"combindedpoints.shp")
  orogname                       <- paste0(TEMP,"orog.sgrd")
  
  glacier_t <- glacier_lgm
  
  #### ---- get glacier outline
  systemcommand<-paste("saga_cmd grid_tools 12 -INPUT=",glacier_t," -OUTPUT=",glacier_outline_grid," -METHOD=1 -RANGE=",lookup1,sep="")
  system(systemcommand)
  
  #### ---- vectorizing grid classes
  systemcommand<-paste("saga_cmd shapes_grid 6 -GRID=",glacier_outline_grid," -POLYGONS=",glacier_outline_shape," -CLASS_ALL=1 -CLASS_ID=1.000000 -SPLIT=0 -ALLVERTICES=0",sep="")
  system(systemcommand)
  
  #### ---- polygon to lines
  systemcommand<-paste("saga_cmd shapes_lines 0 -POLYGONS=",glacier_outline_shape," -LINES=",glacier_outline_line,sep="")
  system(systemcommand)
  
  #### ---- clip random points to glacier extent
  systemcommand<-paste("saga_cmd shapes_points 8 -POINTS=",randompoints," -POLYGONS=",glacier_outline_shape," -CLIPS=",ice6grandom_points," -FIELD=-1 -METHOD=0",sep="")
  system(systemcommand)
  
  #### ---- convert lines to points
  systemcommand<-paste("saga_cmd shapes_points 5 -LINES=",glacier_outline_line," -POINTS=",glacier_outline_points," -ADD=0 -METHOD_INSERT=0 -DIST=1.000000 -ADD_POINT_ORDER=0",sep="")
  system(systemcommand)
  
  #### ---- add grid values to points
  systemcommand<-paste("saga_cmd shapes_grid 0 -SHAPES=",glacier_outline_points," -GRIDS=",gmted2010," -RESULT=",glacier_outline_points," -RESAMPLING=3",sep="")
  system(systemcommand)
  
  #### ---- get ICE6G data
  r1<-raster(paste(ice6gname),varname="orog")
  r2<-raster(orogbias)
  r1<-r1-r2
  writeRaster(r1,file=orogname,format="SAGA",overwrite=TRUE)
  
  # CHANGE LONGITUDINAL RANGE FOR GRIDS
  systemcommand<-paste("saga_cmd pj_proj4 13 -INPUT=",orogname," -OUTPUT=",orogname," -DIRECTION=0 -PATCH=1",sep = "")
  system(systemcommand)
  
  #### ---- add grid values to points
  s1<-shapefile(ice6grandom_points)
  s1<-s1[,1:2]
  shapefile(s1,filename=ice6grandom_points,overwrite=TRUE)
  systemcommand<-paste("saga_cmd shapes_grid 0 -SHAPES=", ice6grandom_points," -GRIDS=",orogname," -RESULT=",ice6grandom_points," -RESAMPLING=3",sep="")
  system(systemcommand)
  
  #### --- combine random points and outline
  s1<-shapefile(ice6grandom_points)
  s2<-shapefile(glacier_outline_points)
  names(s1)<-c("ID","VALUE","elev")
  s2<-s2[,-3]
  names(s2)<-c("ID","VALUE","elev")
  ss<-rbind(as.data.frame(s1),as.data.frame(s2))
  ss[ss>=99999]<--99999
  ss[ss<=-9999]<--99999
  coordinates(ss)<-~coords.x1+coords.x2
  combined_points <-combpoints 
  shapefile(ss,filename=combined_points,overwrite=T)
  
  #### ---- Multilevel b spline
  systemcommand<-paste("saga_cmd grid_spline 4 -SHAPES=",combined_points," -FIELD=3 -TARGET_DEFINITION=1 -TARGET_TEMPLATE=",gmted2010," -TARGET_OUT_GRID=",glacierheights," -TARGET_USER_FITS=0 -METHOD=1 -EPSILON=0.000100 -LEVEL_MAX=14",sep="")
  system(systemcommand)
  grids<-paste("\"",glacier_outline_grid,";",glacierheights,"\"", sep="")
  glacier_t_name <- paste("glacierheigth_",t,sep="")
  equ<-paste("\"a*b\"", sep="")
  systemcommand<-paste("saga_cmd grid_calculus 1 -GRIDS=",grids," -RESULT=",glacierheights," -FORMULA=",equ," -NAME=",glacier_t_name," -TYPE=7",sep="")
  system(systemcommand)
  
  #### ---- add glaciers to dem    
  grids<-paste("\"",gmted2010,";",glacierheights,"\"", sep="")
  equ<-paste("\"ifelse(b>a,b*1,a*1)\"", sep="")
  systemcommand<-paste("saga_cmd grid_calculus 1 -GRIDS=",grids," -RESULT=",gmted_tx," -FORMULA=",equ," -NAME=Calculation -TYPE=7 -USE_NODATA=1",sep="")
  system(systemcommand)
  
  #### ---- adjust DEM by sealevel and glaciers (The Phanerozoic Record of Global Sea-Level Change, Miller et al Science)
  sealevel_t<-sealevel$sealevel[t+1]    
  grids<-paste("\"",gmted_tx,"\"", sep="")
  equ<-paste("\"a-(",sealevel_t,")\"", sep="")    
  glacier_dem_name<-paste("glacier_plus_dem_",t,sep="")
  systemcommand<-paste("saga_cmd grid_calculus 1 -GRIDS=",grids," -RESULT=",glacier_dem," -FORMULA=",equ," -NAME=",glacier_dem_name," -TYPE=7",sep="")
  system(systemcommand)
  systemcommand<-paste("saga_cmd grid_tools 12 -INPUT=",glacier_dem," -OUTPUT=",glacier_dem," -METHOD=1 -RANGE=",lookup2,sep="")
  system(systemcommand)
  
  #### ---- add glaciers to dem
  glacier_dem<-paste(OUTPUT,"glacier_plus_dem_",t,".sgrd",sep="")
  grids<-paste("\"",gmted_tx,";",glacierheights,"\"", sep="")
  equ<-paste("\"ifelse(b>a,b*1,a*1)\"", sep="")
  systemcommand<-paste("saga_cmd grid_calculus 1 -GRIDS=",grids," -RESULT=",glacier_dem," -FORMULA=",equ," -NAME=Calculation -TYPE=7 -USE_NODATA=1",sep="")
  system(systemcommand)
  
  ### =========================================================================
  ### calculations for t+1
  ### =========================================================================
  
  tt=t+1
  nc_in<-nc1
  
  ### =========================================================================
  ### set filenames
  ### =========================================================================
  
  biasgrid                       <- paste0(INPUT,"t0_t210_bias.sgrd")
  bio12                          <- paste0(INPUT,"chelsa_bio1_V1.2.sdat")
  tas_glacieroutline_t0_interpol <- paste0(INPUT,"tas_t0_interpol_biascor.sgrd")      
  
  bias_t          <- paste0(TEMP,"biastemp.sgrd")
  lowresdem       <- paste0(TEMP,"loresdem.sgrd")
  tasname_tt      <- paste0(TEMP,"tas",tt,".tif")
  lowresdem       <- paste0(TEMP,"loresdem.sgrd")
  tasname_t       <- paste0(TEMP,"tas",t,".tif")
  
  tmean_t1_name   <- paste0(OUTPUT,"CHELSA_TraCE_tas_",tt,".sgrd")
  glacier_null    <- paste0(OUTPUT,"glacierheigth_0.sgrd")
  glacier_tt      <- paste0(OUTPUT,"glacierheigth_",tt,".sgrd")
  resultname      <- paste0(OUTPUT,"CHELSA_TraCE_tas_",tt,".sgrd")
  lapsername      <- paste0(OUTPUT,"CHELSA_TraCE_t_lapserate_",t,".tif")
  lapsetname      <- paste0(OUTPUT,"CHELSA_TraCE_t_lapserate_",t,".sdat")
  coeff1          <- paste0(OUTPUT,"CHELSA_TraCE_t_intercept_",t,".sgrd")
  coeff2          <- paste0(OUTPUT,"CHELSA_TraCE_t_lapserate_",t,".sgrd")
  resultname      <- paste0(OUTPUT,"CHELSA_TraCE_tas_",t,".sgrd")
  coeff2tif       <- paste0(OUTPUT,"CHELSA_TraCE_t_lapserate_",t,".tif")
  
  #### ---- downscale temp for t+1 on that dem
  nc1  <-nc_open(nc_in)
  times<- nc1$dim$time$vals
  rt   <-raster(tracetemplate)
  pos  <-decades[centuries==tt+10]
  
  names      <-c()
  names_z    <-c()
  names_180  <-c()
  names_z_180<-c()
   
  for (l in 16:26)
  {
    r1<-stack()
    for (tt_dec in 1:10)
    {
      r1_x<-raster(nc_in  ,var="T",band=pos[tt_dec],lev=l)
      r1_x<-rotate(r1_x)
      r1<-stack(r1,r1_x)
    }
    r1<-calc(r1,mean)
    r2<-stack()
    for (tt_dec in 1:10)
    {
      r2_x<-raster(nc_in  ,var="Z3",band=pos[tt_dec],lev=l)
      r2_x<-rotate(r2_x)
      r2<-stack(r2,r2_x)
    }
    r2<-calc(r2,mean)
    fname   <-paste(TEMP,"t_",l,".tif",sep="")
    fname180<-paste(TEMP,"t_",l,"_180.tif",sep="")
    zname   <-paste(TEMP,"z_",l,".tif",sep="")
    zname180<-paste(TEMP,"z_",l,"_180.tif",sep="")
    
    NAvalue(r1)<--99999
    NAvalue(r2)<--99999
    
    writeRaster(r1,filename=fname,format="GTiff",overwrite=TRUE)
    writeRaster(r2,filename=zname,format="GTiff",overwrite=TRUE)
    names  <-paste(names,";",fname,sep="")
    names_z<-paste(names_z,";",zname,sep="")
    names_180  <-paste(names_180,";",fname180,sep="")
    names_z_180<-paste(names_z_180,";",zname180,sep="")
  }
  names<-substring(names,2)
  names<-paste("\"",names,"\"",sep="")
  names_z<-substring(names_z,2)
  names_z<-paste("\"",names_z,"\"",sep="")
  
  #### ---- CALCULATE TEMPERATURE LAPSERATE
  coeffs<-paste("\"",coeff1,";",coeff2,"\"", sep="")
  systemcommand<-paste("saga_cmd statistics_regression 9 -Y_GRIDS=",names," -X_GRIDS=",names_z," -COEFF=",coeffs," -ORDER=1 -XSOURCE=2",sep="")
  system(systemcommand)
  equ2<-paste("\"","a*100*(-1)","\"")
  systemcommand<-paste("saga_cmd grid_calculus 1 -GRIDS=",coeff2," -RESULT=",coeff2," -FORMULA=",equ2," -NAME=Calculation -TYPE=7",sep="")
  system(systemcommand)
  
  #### ---- DOWNSCALE TAS USING LAPSERATES --- FOR TIME t + 1 (needed for next iteration)
  #### ---- Calculate temp in last timestep
  posx<-decades[centuries==210+10]
  tas_ref<-stack()
  for (tt_dec in 1:10)
  {
    tas_x  <-raster(nc_in  ,var="TREFHT",band=posx[tt_dec])
    tas_x  <-tas_x-273.15
    tas_x  <-rotate(tas_x)
    tas_ref<-stack(tas_ref,tas_x)
  }
  tas_ref<-calc(tas_ref,mean)
  rx<-raster(nc_in,var="TREFHT",band=2204)
  rx<-rx-273.15
  rx<-rotate(rx)
  tas_ref<-rx
  
  tas<-stack()
  for (tt_dec in 1:10)
  {
    tas_x <-raster(nc_in  ,var="TREFHT",band=pos[tt_dec])
    tas_x <-tas_x-273.15
    tas_x <-rotate(tas_x)
    tas   <-stack(tas,tas_x)
  } 
  tas<-calc(tas,mean)
  tasxx<-tas_ref-tas
  tasxx<-resample(tasxx,tas_biasgrid,method="bilinear")
  tas<-tas_biasgrid-tasxx
  
  lapse           <-raster(lapsetname)
  lapse<-resample(lapse,tas_biasgrid,method="bilinear")
  writeRaster(lapse,file=lapsername,format="GTiff",overwrite=TRUE)

  ##### 
  writeRaster(tas,filename=tasname_tt,format="GTiff",overwrite=TRUE)
  systemcommand<-paste("saga_cmd grid_tools 0 -INPUT=",glacier_dem," -OUTPUT=",lowresdem," -KEEP_TYPE=0 -SCALE_UP=4 -TARGET_DEFINITION=1 -TARGET_TEMPLATE=",tasname_tt,sep="")
  system(systemcommand)
  systemcommand<-paste("saga_cmd climate t_downscale -LORES_DEM=",lowresdem," -LORES_T=",tasname_tt," -LAPSE_RATES=",coeff2tif," -HIRES_DEM=",glacier_dem," -HIRES_T=",resultname,sep="")
  system(systemcommand)
  
  #### ---- DOWNSCALE TAS USING LAPSERATES--- FOR TIME t
  pos<-decades[centuries==t+10]
  tas<-stack()
  for (tt_dec in 1:10)
  {
    tas_x<-raster(nc_in  ,var="TREFHT",band=pos[tt_dec])
    tas_x <-tas_x-273.15
    tas_x <-rotate(tas_x)
    tas<-stack(tas,tas_x)
  } 
  tas<-calc(tas,mean)
  tasxx<-tas_ref-tas
  tas_biasgrid<-raster(bio12)
  tasxx<-resample(tasxx,tas_biasgrid,method="bilinear")
  tas<-tas_biasgrid-tasxx
  writeRaster(tas,filename=tasname_t,format="GTiff",overwrite=TRUE)

  #lowresdem<-paste(INPUT,"gmted2010_lowres.sgrd",sep="")
  systemcommand<-paste("saga_cmd grid_tools 0 -INPUT=",glacier_dem," -OUTPUT=",lowresdem," -KEEP_TYPE=0 -SCALE_UP=4 -TARGET_DEFINITION=1 -TARGET_TEMPLATE=",tasname_t,sep="")
  system(systemcommand)
  systemcommand<-paste("saga_cmd climate t_downscale -LORES_DEM=",lowresdem," -LORES_T=",tasname_t," -LAPSE_RATES=",coeff2," -HIRES_DEM=",glacier_dem," -HIRES_T=",resultname,sep="")
  system(systemcommand)
  
  #### ---- ADJUST BIAS TO LENGTH OF TIMESTEP
  grids <-paste("\"",biasgrid,"\"", sep="")
  equx<-paste("\"","ifelse(a<0,a*100/21000*",tt,",a*1)\"", sep="")
  systemcommand<-paste("saga_cmd grid_calculus 1 -GRIDS=",grids,"  -RESULT=",bias_t," -FORMULA=",equx," -NAME=bias -TYPE=7",sep="")
  system(systemcommand)
  
  #### ---- CALCULATE GLACIERS AT T+1

  grids                    <- paste("\"",glacier_null,";",tmean_t1_name,"\"", sep="")
  xgrids                   <- paste("\"",tas_glacieroutline_t0_interpol,";",bias_t,"\"", sep="")
  glacier_tt_name          <- paste("glacierheigth_",tt,sep="")
  equx                     <- paste("\"","ifelse(b>(c+d),a-a-99999,a*1)","\"", sep="")
  systemcommand            <- paste("saga_cmd grid_calculus 1 -GRIDS=",grids," -XGRIDS=",xgrids," -RESULT=",glacier_tt," -FORMULA=",equ," -NAME=",glacier_tt_name," -TYPE=7",sep="")
  system(systemcommand)
  
  #### ---- END OF LOOP
}    

